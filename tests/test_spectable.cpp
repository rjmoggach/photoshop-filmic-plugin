#include <doctest/doctest.h>
#include "lenscore/color/spectable.hpp"
#include <cstdio>

using namespace lens::color;

static const SpecTable& smallTable() {
    static SpecTable t = buildTable(12);   // small, so the suite stays fast
    return t;
}

TEST_CASE("table has the documented size") {
    const SpecTable& t = smallTable();
    CHECK(t.res == 12);
    CHECK(t.data.size() == size_t(3) * 12 * 12 * 12 * 3);
}

TEST_CASE("lookup round trips grey at several brightnesses") {
    for (float v : {0.15f, 0.4f, 0.75f, 1.0f}) {
        const RGB in{v, v, v};
        const RGB got = spectrumToRec2020(lookup(smallTable(), in));
        CAPTURE(v);
        CHECK(got.r / v == doctest::Approx(1.0f).epsilon(0.06));
        CHECK(got.g / v == doctest::Approx(1.0f).epsilon(0.06));
        CHECK(got.b / v == doctest::Approx(1.0f).epsilon(0.06));
    }
}

TEST_CASE("lookup is scale equivariant, so HDR values pass through") {
    // The table only spans [0,1]^3; anything at or beyond the gamut boundary along a
    // ray clamps to the same cell, so the coefficients must agree exactly there. Below
    // the boundary the model is NOT scale-invariant -- a muted, less-bright colour
    // legitimately needs a different spectral shape than the same hue pushed to full
    // saturation, so this deliberately compares two points already at/past the boundary
    // (1.0x and 10x), not an in-gamut point against an out-of-gamut one.
    const Coeffs a = lookup(smallTable(), RGB{1.0f, 0.5f, 0.25f});
    const Coeffs b = lookup(smallTable(), RGB{10.0f, 5.0f, 2.5f});
    CHECK(a.c0 == doctest::Approx(b.c0).epsilon(1e-4));
    CHECK(a.c1 == doctest::Approx(b.c1).epsilon(1e-4));
    CHECK(a.c2 == doctest::Approx(b.c2).epsilon(1e-4));
}

TEST_CASE("lookup round trips saturated hues within tolerance") {
    const RGB targets[] = {{1.0f, 0.2f, 0.1f}, {0.1f, 1.0f, 0.2f}, {0.15f, 0.2f, 1.0f}};
    for (const RGB& t : targets) {
        const RGB got = spectrumToRec2020(lookup(smallTable(), t));
        CAPTURE(t.r); CAPTURE(t.g); CAPTURE(t.b);
        CHECK(got.r == doctest::Approx(t.r).epsilon(0.10));
        CHECK(got.g == doctest::Approx(t.g).epsilon(0.10));
        CHECK(got.b == doctest::Approx(t.b).epsilon(0.10));
    }
}

TEST_CASE("table serialises and reloads bit exactly") {
    const std::string p = "tbl.bin";
    REQUIRE(writeTable(p, smallTable()));
    auto back = readTable(p);
    REQUIRE(back.has_value());
    CHECK(back->res == smallTable().res);
    CHECK(back->data == smallTable().data);
    std::remove(p.c_str());
}
