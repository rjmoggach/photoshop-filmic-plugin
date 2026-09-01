#include <doctest/doctest.h>
#include "lenscore/color/cie.hpp"

using namespace lens::color;

TEST_CASE("ybar peaks near 555nm") {
    float best = 0.0f, arg = 0.0f;
    for (float l = 380.0f; l <= 780.0f; l += 0.5f)
        if (cmf(l).y > best) { best = cmf(l).y; arg = l; }
    CHECK(arg == doctest::Approx(555.0f).epsilon(0.02));
    CHECK(best == doctest::Approx(1.0f).epsilon(0.05));
}

TEST_CASE("equal energy spectrum is illuminant E, chromaticity one third") {
    XYZ sum{0, 0, 0};
    for (float l = 380.0f; l <= 780.0f; l += 1.0f) {
        const XYZ c = cmf(l);
        sum.x += c.x; sum.y += c.y; sum.z += c.z;
    }
    const float t = sum.x + sum.y + sum.z;
    CHECK(sum.x / t == doctest::Approx(1.0f / 3.0f).epsilon(0.02));
    CHECK(sum.y / t == doctest::Approx(1.0f / 3.0f).epsilon(0.02));
}

TEST_CASE("Rec2020 matrices are mutual inverses") {
    const RGB in{0.3f, 0.7f, 0.2f};
    const RGB out = xyzToRec2020(rec2020ToXyz(in));
    CHECK(out.r == doctest::Approx(in.r).epsilon(1e-4));
    CHECK(out.g == doctest::Approx(in.g).epsilon(1e-4));
    CHECK(out.b == doctest::Approx(in.b).epsilon(1e-4));
}

TEST_CASE("Rec2020 white maps to a neutral chromaticity") {
    const XYZ w = rec2020ToXyz(RGB{1.0f, 1.0f, 1.0f});
    const float t = w.x + w.y + w.z;
    CHECK(w.x / t == doctest::Approx(0.3127f).epsilon(0.01));  // D65
    CHECK(w.y / t == doctest::Approx(0.3290f).epsilon(0.01));
}

TEST_CASE("cieYNormalisation returns neutral divisor on empty input") {
    CHECK(cieYNormalisation(nullptr, nullptr, 0) == 1.0f);
}
