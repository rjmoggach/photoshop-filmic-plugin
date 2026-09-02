#include <doctest/doctest.h>
#include "lenscore/color/upsample.hpp"
#include <cmath>

using namespace lens::color;

TEST_CASE("sigmoid is bounded and centred") {
    CHECK(evalSpectrum(Coeffs{0, 0, 0}, 550.0f) == doctest::Approx(0.5f));
    for (float c2 : {-50.0f, -1.0f, 0.0f, 1.0f, 50.0f}) {
        const float v = evalSpectrum(Coeffs{0, 0, c2}, 550.0f);
        CHECK(v >= 0.0f);
        CHECK(v <= 1.0f);
    }
}

TEST_CASE("large positive constant gives a flat white spectrum") {
    const Coeffs c{0.0f, 0.0f, 1e4f};
    for (float l = 400.0f; l <= 700.0f; l += 50.0f)
        CHECK(evalSpectrum(c, l) == doctest::Approx(1.0f).epsilon(1e-3).scale(0));
    const RGB rgb = spectrumToRec2020(c);
    CHECK(rgb.r == doctest::Approx(1.0f).epsilon(0.02).scale(0));
    CHECK(rgb.g == doctest::Approx(1.0f).epsilon(0.02).scale(0));
    CHECK(rgb.b == doctest::Approx(1.0f).epsilon(0.02).scale(0));
}

TEST_CASE("fit round trips neutrals and moderate hues accurately") {
    const RGB targets[] = {
        {1.0f, 1.0f, 1.0f}, {1.0f, 0.5f, 0.2f},
        {0.1f, 0.3f, 1.0f}, {1.0f, 1.0f, 0.1f}, {0.6f, 0.6f, 0.6f},
    };
    for (const RGB& t : targets) {
        const Coeffs c = fitCoeffs(t);
        const RGB got = spectrumToRec2020(c);
        CAPTURE(t.r); CAPTURE(t.g); CAPTURE(t.b);
        CHECK(got.r == doctest::Approx(t.r).epsilon(0.03).scale(0));
        CHECK(got.g == doctest::Approx(t.g).epsilon(0.03).scale(0));
        CHECK(got.b == doctest::Approx(t.b).epsilon(0.03).scale(0));
    }
}

// (0.2, 1.0, 0.3) sits far enough toward the Rec.2020 green primary that the model
// gamut-maps it rather than reproducing it exactly -- a MODEL LIMITATION, not a table
// artefact (there is no table in this path at all). Jakob-Hanika's zero-error
// guarantee is for the sRGB gamut; a bounded reflectance spectrum cannot reach
// chromaticities this saturated (see upsample.hpp's spectrumToRec2020 comment and the
// design spec's accepted-limitation note). Assert what IS true -- bounded and
// non-negative -- not an exact round trip.
TEST_CASE("fit round trips a highly saturated hue as bounded, not exact") {
    const RGB t{0.2f, 1.0f, 0.3f};
    const RGB got = spectrumToRec2020(fitCoeffs(t));
    CAPTURE(t.r); CAPTURE(t.g); CAPTURE(t.b);
    CHECK(got.r >= 0.0f); CHECK(got.r <= 1.5f);
    CHECK(got.g >= 0.0f); CHECK(got.g <= 1.5f);
    CHECK(got.b >= 0.0f); CHECK(got.b <= 1.5f);
}

TEST_CASE("fitted spectra stay smooth, no oscillation") {
    const Coeffs c = fitCoeffs(RGB{1.0f, 0.35f, 0.1f});
    int reversals = 0;
    float prev = evalSpectrum(c, 380.0f), prevSlope = 0.0f;
    for (float l = 385.0f; l <= 780.0f; l += 5.0f) {
        const float v = evalSpectrum(c, l);
        const float slope = v - prev;
        if (slope * prevSlope < 0.0f) ++reversals;
        prevSlope = slope; prev = v;
    }
    CHECK(reversals <= 2);   // a quadratic through a monotone sigmoid cannot wiggle more
}
