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
    //   psfKernel * pixelPitchUm <= psfGrid * psfSampleSpacingUm(lambda, fNumberWide) / pupilFill
    // (see lenscore/pipeline.hpp, "spp" and the psfKernel-vs-grid validation it does per
    // band) -- the kernel's physical footprint (psfKernel pixels wide) cannot exceed the
    // PSF ring's own FFT grid footprint (psfGrid samples wide, each psfSampleSpacingUm /
    // pupilFill micrometres -- pupilFill DIVIDES here because it shrinks the illuminated
    // pupil into a fraction of the grid, i.e. represents the same physical aperture at
    // higher sample density, not lower; see pipeline.hpp for the verification), or the
    // convolution reads clamped edge samples of the ring and energy conservation breaks.
    // render() throws std::invalid_argument if a caller's psfKernel violates this for the
    // wavelength actually being rendered, rather than silently truncating.
    //
    // Getting the coverage condition right is necessary but not sufficient: the kernel
    // also has to be wide enough, in OUTPUT PIXELS, to actually capture the diffraction
    // pattern's energy, or it just conserves energy for a kernel that's throwing most of
    // that energy away. At a representative mid-visible wavelength (550nm) with the
    // sampling parameters above, one grid sample is psfSampleSpacingUm(550,2.0)/0.25 =
    // 4.4um, so the pattern spans a few output pixels (4.4/5.0 grid samples per pixel) --
    // not sub-pixel, not dozens of pixels wide. Measured against the energy-conservation
    // test's own aberration-free configuration across its 3 sampled bands (447-713nm),
    // psfKernel=65 recovers >99% of the on-axis kernel's energy at every one of them,
    // comfortably inside the coverage bound above (grid=128 clears up to 111 pixels at
    // 550nm, 91 at the bluest sampled band) -- so 65 was chosen as a default with real
    // margin on both sides, not the tightest value that merely passes.
    int   psfGrid      = 128;   // pupil FFT grid, power of two
    int   psfRings     = 12;
    int   psfKernel    = 65;    // odd; see the derivation above
    int   effPatch     = 128;   // power of two; must be >= psfKernel for effConvolve's
                                  // per-patch window to comfortably contain the kernel
};

}  // namespace lens
