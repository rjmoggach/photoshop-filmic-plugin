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
    int   psfGrid      = 128;   // pupil FFT grid, power of two
    int   psfRings     = 12;
    int   psfKernel    = 33;    // odd
    int   effPatch     = 64;    // power of two
};

}  // namespace lens
