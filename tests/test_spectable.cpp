#include <doctest/doctest.h>
#include "lenscore/color/spectable.hpp"
#include <cstdio>
#include <stdexcept>
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

// buildTable(res) sizes a 3*res^3*3-float allocation directly from res; res <= 1 also
// divides by (res - 1) in axisScale/axisBrightness. Both the shipped generator
// (rgb2spec/main.cpp) and readTable's own bound feed from untrusted-ish external input
// (a CLI argument, a file header respectively), so buildTable must reject the same
// range readTable already does rather than sizing an allocation from it unchecked.
TEST_CASE("buildTable rejects a resolution outside [2, kMaxTableRes]") {
    CHECK_THROWS_AS(buildTable(0), std::invalid_argument);
    CHECK_THROWS_AS(buildTable(1), std::invalid_argument);
    CHECK_THROWS_AS(buildTable(-5), std::invalid_argument);
    CHECK_THROWS_AS(buildTable(kMaxTableRes + 1), std::invalid_argument);
    CHECK_NOTHROW(buildTable(2));
}

TEST_CASE("lookup round trips grey at several brightnesses") {
    // res=12 is a deliberate speed compromise: measured grey round-trip error is
    // ~11.8% on neutrals at res=12 (worst case here, v=0.15), vs. ~2.5% at res=24 and
    // ~0.3% at res=48/64 -- the shipping table is res=64 and is accurate. This states
    // the measured res=12 behaviour rather than raising the test table's resolution
    // (which would slow the suite for no benefit).
    for (float v : {0.15f, 0.4f, 0.75f, 1.0f}) {
        const RGB in{v, v, v};
        const RGB got = spectrumToRec2020(lookup(smallTable(), in));
        CAPTURE(v);
        CHECK(got.r / v == doctest::Approx(1.0f).epsilon(0.13).scale(0));
        CHECK(got.g / v == doctest::Approx(1.0f).epsilon(0.13).scale(0));
        CHECK(got.b / v == doctest::Approx(1.0f).epsilon(0.13).scale(0));
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
        CHECK(got.r / v == doctest::Approx(1.0f).epsilon(0.6).scale(0));
        CHECK(got.g / v == doctest::Approx(1.0f).epsilon(0.6).scale(0));
        CHECK(got.b / v == doctest::Approx(1.0f).epsilon(0.6).scale(0));
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
    CHECK(a.c0 == doctest::Approx(b.c0).epsilon(1e-4).scale(0));
    CHECK(a.c1 == doctest::Approx(b.c1).epsilon(1e-4).scale(0));
    CHECK(a.c2 == doctest::Approx(b.c2).epsilon(1e-4).scale(0));
}

TEST_CASE("lookup round trips saturated hues within tolerance") {
    // (0.15, 0.2, 1.0) is saturated but sits close enough to the model's reachable
    // gamut that the round trip is accurate -- keep a tight check on it, and this holds
    // at every table resolution, not just here.
    const RGB accurate{0.15f, 0.2f, 1.0f};
    const RGB got = spectrumToRec2020(lookup(smallTable(), accurate));
    CHECK(got.r == doctest::Approx(accurate.r).epsilon(0.10).scale(0));
    CHECK(got.g == doctest::Approx(accurate.g).epsilon(0.10).scale(0));
    CHECK(got.b == doctest::Approx(accurate.b).epsilon(0.10).scale(0));
}

// (1.0, 0.2, 0.1) and (0.1, 1.0, 0.2) sit far out toward the Rec.2020 red and green
// primaries. This is a MODEL LIMITATION, not a table-resolution artefact: it fails
// identically at every table resolution and in the direct fitCoeffs path with no table
// at all. Jakob-Hanika report zero error only on the much narrower sRGB gamut; a
// bounded reflectance spectrum cannot reach chromaticities this saturated, so the
// model gamut-maps them (measured: up to +73% on a channel, and a channel driven
// slightly negative before spectrumToRec2020's clamp). This is a recorded, accepted
// limitation (see upsample.hpp and the design spec) -- assert what IS true, bounded
// and non-negative, not an exact round trip.
TEST_CASE("lookup round trips near-primary saturated hues as bounded, not exact") {
    const RGB targets[] = {{1.0f, 0.2f, 0.1f}, {0.1f, 1.0f, 0.2f}};
    for (const RGB& t : targets) {
        const RGB got = spectrumToRec2020(lookup(smallTable(), t));
        CAPTURE(t.r); CAPTURE(t.g); CAPTURE(t.b);
        CHECK(got.r >= 0.0f); CHECK(got.r <= 1.5f);
        CHECK(got.g >= 0.0f); CHECK(got.g <= 1.5f);
        CHECK(got.b >= 0.0f); CHECK(got.b <= 1.5f);
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
