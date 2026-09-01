#pragma once
#include "lenscore/color/spectable.hpp"
#include "lenscore/color/transfer.hpp"
#include "lenscore/optics/dispersion.hpp"
#include "lenscore/optics/distortion.hpp"
#include "lenscore/optics/pupil.hpp"
#include "lenscore/optics/vignette.hpp"
#include <utility>
#include <vector>

namespace lens {

struct Params {
    int schema = 1;

    // Stage switches. A user who wants only chromatic aberration pays for nothing else.
    bool doLateralCa = true;
    bool doPsf       = true;
    bool doVignette  = true;

    // Colour
    bool         highlightRecovery = true;
    color::Knee  knee{};
    int          bands = 11;                 // 3 is the preview tier
    std::vector<std::pair<float, float>> equalizer{};   // (lambda, weight); empty means flat

    // Lens
    optics::Dispersion     dispersion{};
    optics::Distortion     distortion{};
    optics::PupilParams    pupil{};
    optics::VignetteParams vignette{};
    float lateralK   = 0.0f;    // k_l
    float lambdaHat  = 650.0f;  // reference wavelength, keeps magnification positive

    // Wavefront coefficients in waves, at the corner of the frame
    float petzval   = 0.0f;
    float astig     = 0.0f;
    float coma      = 0.0f;
    float spherical = 0.0f;

    // Sampling
    float focal_mm     = 32.0f;
    float fNumberWide  = 2.0f;
    float pixelPitchUm = 5.0f;
    float pupilFill    = 0.25f;  // fraction of the FFT grid the pupil fills at working stop;
                                  // 1.0 would put the first Airy zero at 1.22 samples --
                                  // critically undersampled. Also scales the sample spacing.
    // psfGrid and psfKernel must satisfy, at every sampled wavelength:
    //   psfKernel * pixelPitchUm <= psfGrid * psfSampleSpacingUm(lambda, fNumberWide) * pupilFill
    // (see lenscore/pipeline.hpp, "spp" and the psfKernel-vs-grid validation it does per
    // band) -- the kernel's physical footprint (psfKernel pixels wide) cannot exceed the
    // PSF ring's own FFT grid footprint (psfGrid samples wide, each psfSampleSpacingUm *
    // pupilFill micrometres -- pupilFill MULTIPLIES here: psfSampleSpacingUm is
    // calibrated for a pupil filling the whole grid, and pupilFill is how much smaller
    // than that the actual working aperture is; see pipeline.hpp for the reasoning), or
    // psfAtField's box average (see psfrings.hpp) draws from clamped edge samples and
    // energy conservation breaks. render() throws std::invalid_argument if a caller's
    // psfKernel violates this for the wavelength actually being rendered, rather than
    // silently truncating.
    //
    // At these defaults an unaberrated PSF is genuinely sub-pixel (samplesPerPixel is
    // large, ~18 at 550nm: gridSampleUm = 550e-3 * 2.0 * 0.25 = 0.275um, vs a 5um
    // pixel) -- that is real, not a bug, and is exactly why psfAtField box-averages
    // rather than point-samples once samplesPerPixel > 1 (see psfrings.hpp).
    //
    // gridSampleUm is proportional to lambda, so the binding case for the invariant
    // above is the SHORTEST wavelength actually sampled, not the 550nm reference --
    // sampleBands (pipeline.hpp) can draw as low as kLambdaMin = 380nm. At 380nm,
    // gridSampleUm = 380e-3 * 2.0 * 0.25 = 0.19um, so maxKernel = floor(256 * 0.19 /
    // 5.0) = 9 (odd already); psfKernel = 7 clears that with margin (128 does not --
    // floor(128 * 0.19 / 5.0) = 4, which floors to an even 4 and rounds down to 3,
    // well under 7 -- an earlier default of 128 threw std::invalid_argument on a
    // plain `Params{}` render as soon as sampleBands touched the blue end of the
    // spectrum; 256 is the smallest power of two that stays clear of that all the way
    // down to 380nm at the default pixelPitchUm/pupilFill/fNumberWide). 7 pixels
    // covers 7 * 18.18 ~= 127 ring samples at 550nm -- ample for the aberrated PSFs
    // this project actually renders (a several-pixel-wide defocused spot is a much
    // easier target than the unaberrated case the bound is computed from). Measured
    // directly (see psfrings.hpp's psfAtField and test_psfrings.cpp): with the
    // box-average fix, this configuration holds the on-axis kernel's energy within a
    // couple of percent of 1 from bare Airy through 20 waves of defocus.
    //
    // This is still tied to fNumberWide: a lens file with a faster wide-open stop
    // (smaller t_stop) shrinks gridSampleUm further and can still violate the
    // invariant at these defaults -- render() throws a clear, actionable message
    // rather than silently truncating (see pipeline.hpp), so that is a loud failure
    // for the caller to raise psfGrid or shrink psfKernel, not a silent one.
    int   psfGrid      = 256;   // pupil FFT grid, power of two
    int   psfRings     = 12;
    int   psfKernel    = 7;     // odd; see the derivation above
    int   effPatch     = 64;    // power of two
};

}  // namespace lens
