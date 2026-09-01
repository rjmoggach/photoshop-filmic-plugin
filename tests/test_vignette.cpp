#include <doctest/doctest.h>
#include "lenscore/optics/vignette.hpp"
#include <cmath>

using namespace lens::optics;

static VignetteParams defaults() {
    VignetteParams p;
    p.focal_mm = 32.0f; p.sensorHalfDiag_mm = 14.0f; p.naturalExp = 4.0f;
    p.tStopWide = 2.0f; p.tStop = 2.0f; p.rEntrance = 1.0f; p.sepNorm = 0.5f;
    return p;
}

TEST_CASE("overlap area handles containment and disjointness") {
    const float pi = 3.14159265f;
    CHECK(circleOverlapArea(5.0f, 1.0f, 1.0f) == doctest::Approx(0.0f));            // disjoint
    CHECK(circleOverlapArea(0.0f, 1.0f, 2.0f) == doctest::Approx(pi));              // contained
    CHECK(circleOverlapArea(0.0f, 1.0f, 1.0f) == doctest::Approx(pi));              // coincident
}

TEST_CASE("overlap area of two unit circles at unit separation is the known value") {
    // 2*acos(1/2) - sqrt(3)/2 = 1.22837
    CHECK(circleOverlapArea(1.0f, 1.0f, 1.0f) == doctest::Approx(1.22837f).epsilon(1e-4).scale(0));
}

TEST_CASE("everything is exactly unity at the optical centre") {
    const VignetteParams p = defaults();
    CHECK(naturalFalloff(p, 0.0f) == doctest::Approx(1.0f));
    CHECK(mechanicalFraction(p, 0.0f) == doctest::Approx(1.0f));
    CHECK(vignette(p, 0.0f) == doctest::Approx(1.0f));
}

TEST_CASE("natural falloff matches the cosine law analytically") {
    VignetteParams p = defaults();
    const float r = 1.0f * p.sensorHalfDiag_mm;
    const float expected = std::pow(std::cos(std::atan(r / p.focal_mm)), 4.0f);
    CHECK(naturalFalloff(p, 1.0f) == doctest::Approx(expected).epsilon(1e-5).scale(0));
}

TEST_CASE("the exponent is free, not pinned at four") {
    VignetteParams a = defaults(); a.naturalExp = 4.0f;
    VignetteParams b = defaults(); b.naturalExp = 3.2f;
    CHECK(naturalFalloff(b, 1.0f) > naturalFalloff(a, 1.0f));
}

TEST_CASE("natural falloff decreases monotonically toward the corner") {
    const VignetteParams p = defaults();
    float prev = 1.0f;
    for (float t = 0.1f; t <= 1.0f; t += 0.1f) {
        const float v = naturalFalloff(p, t);
        CHECK(v < prev);
        prev = v;
    }
}

TEST_CASE("mechanical vignetting bites wide open") {
    const VignetteParams p = defaults();
    CHECK(mechanicalFraction(p, 1.0f) < 0.99f);
    CHECK(mechanicalFraction(p, 1.0f) > 0.0f);
}

TEST_CASE("mechanical vignetting vanishes above f/4, as measured") {
    VignetteParams p = defaults();
    CHECK(mechanicalVanishStop(p) == doctest::Approx(4.0f));
    p.tStop = 4.0f;  CHECK(mechanicalFraction(p, 1.0f) == doctest::Approx(1.0f).epsilon(1e-4).scale(0));
    p.tStop = 8.0f;  CHECK(mechanicalFraction(p, 1.0f) == doctest::Approx(1.0f).epsilon(1e-4).scale(0));
    p.tStop = 2.8f;  CHECK(mechanicalFraction(p, 1.0f) < 1.0f);
}

TEST_CASE("stopping down monotonically reduces mechanical clipping") {
    VignetteParams p = defaults();
    float prev = 0.0f;
    bool sawStrictIncrease = false;
    for (float s : {2.0f, 2.4f, 2.8f, 3.4f, 4.0f}) {
        p.tStop = s;
        const float m = mechanicalFraction(p, 1.0f);
        CHECK(m >= prev - 1e-6f);
        if (m > prev + 1e-6f) sawStrictIncrease = true;
        prev = m;
    }
    // A constant function would satisfy the non-decreasing checks above too;
    // require a genuine increase somewhere in the sweep to rule that out.
    CHECK(sawStrictIncrease);
}
