#include <doctest/doctest.h>
#include "lenscore/color/spectable.hpp"
#include <cstdio>
#include <unistd.h>

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

TEST_CASE("lookup round trips near-black brightnesses without a warp/bracket mismatch") {
    // Regression for a bug where buildTable floored the brightness axis at 1e-4 (so
    // several low-z slices all fit the SAME target) while lookup bracketed against the
    // unfloored axis (so it believed those slices sat at distinct brightnesses). Queries
    // in the affected range reconstructed at up to ~100x the correct brightness. These
    // values sit below the old 1e-4 floor, so they land squarely in what used to be the
    // collapsed region.
    for (float v : {0.001f, 0.0002f}) {
        const RGB in{v, v, v};
        const RGB got = spectrumToRec2020(lookup(smallTable(), in));
        CAPTURE(v);
        // A 100x error would put the ratio at ~100 or ~0.01; res=12 is coarse this close
        // to black, so the tolerance is loose in absolute terms but nowhere near that.
        CHECK(got.r / v == doctest::Approx(1.0f).epsilon(0.6));
        CHECK(got.g / v == doctest::Approx(1.0f).epsilon(0.6));
        CHECK(got.b / v == doctest::Approx(1.0f).epsilon(0.6));
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

TEST_CASE("readTable rejects a truncated file without throwing") {
    const std::string p = "truncated.bin";
    REQUIRE(writeTable(p, smallTable()));
    // Chop the file down to just past the header: the res the header claims no longer
    // matches the data actually present.
    const long full = [&] {
        std::FILE* f = std::fopen(p.c_str(), "rb");
        REQUIRE(f != nullptr);
        std::fseek(f, 0, SEEK_END);
        const long n = std::ftell(f);
        std::fclose(f);
        return n;
    }();
    REQUIRE(truncate(p.c_str(), full / 2) == 0);
    CHECK_NOTHROW(auto back = readTable(p));
    CHECK_FALSE(readTable(p).has_value());
    std::remove(p.c_str());
}

TEST_CASE("readTable rejects an absurd res header without throwing or allocating huge memory") {
    const std::string p = "absurd_res.bin";
    {
        std::FILE* f = std::fopen(p.c_str(), "wb");
        REQUIRE(f != nullptr);
        std::fwrite("LSPT", 1, 4, f);
        const int hostileRes = 2000000000;
        std::fwrite(&hostileRes, sizeof(int), 1, f);
        // No data follows -- a genuine file of that claimed size would be enormous.
        std::fclose(f);
    }
    CHECK_NOTHROW(auto back = readTable(p));
    CHECK_FALSE(readTable(p).has_value());
    std::remove(p.c_str());
}
