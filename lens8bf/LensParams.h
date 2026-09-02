#pragma once
#include <cmath>

// The controls the dialog exposes.
//
// Split into the ones a photographer already has words for and the ones that are
// lens-design vocabulary. Everything runs -100..+100 with 0 meaning off, except
// aperture (a real T-stop), squeeze (a real ratio), grain size (pixels) and
// quality (a band count).
//
// Kept free of both Photoshop and lenscore headers so the Objective-C++ dialog
// and the C++ filter entry point can share it without dragging either in.

struct LensControls {
    // Plain language
    double distortion          = 0.0;   // minus barrel, plus pincushion
    double chromaticAberration = 0.0;
    double vignette            = 0.0;
    double edgeBlur            = 0.0;

    // Lens character
    double astigmatism         = 0.0;
    double coma                = 0.0;
    double tStop               = 2.0;   // 1.4 .. 16
    double squeeze             = 1.0;   // 1.0 spherical, 2.0 classic scope

    // Film
    double grain               = 0.0;   // 0 .. 100
    double grainSize           = 1.5;   // pixels
    double grainColour         = 35.0;  // 0 monochrome, 100 independent per channel

    // 7 is the lowest count that reaches the colour model's own accuracy floor,
    // so it is Draft rather than a wrong answer -- below it a render is not lower
    // quality, it is incorrect: three bands put a saturated red out by 0.43 of
    // full scale. 11 is the default because it is comfortably past the floor.
    int    bands               = 11;    // 7 draft, 11 normal, 15 best

    bool anyOptics() const {
        return distortion != 0.0 || chromaticAberration != 0.0 || vignette != 0.0 ||
               edgeBlur != 0.0 || astigmatism != 0.0 || coma != 0.0 ||
               std::abs(squeeze - 1.0) > 1e-6;
    }
    bool any() const { return anyOptics() || grain > 0.0; }
};
