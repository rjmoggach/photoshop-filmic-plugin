#include <doctest/doctest.h>
#include "lenscore/color/transfer.hpp"

using namespace lens::color;

TEST_CASE("sRGB transfer round trips") {
    for (float v : {0.0f, 0.001f, 0.04f, 0.5f, 1.0f})
        CHECK(linearToSrgb(srgbToLinear(v)) == doctest::Approx(v).epsilon(1e-5));
}

TEST_CASE("sRGB transfer hits the known anchor") {
    CHECK(srgbToLinear(0.5f) == doctest::Approx(0.2140f).epsilon(1e-3).scale(0));
}

TEST_CASE("knee is identity below the threshold") {
    Knee k;
    CHECK(expandHighlights(0.5f, k) == doctest::Approx(0.5f));
    CHECK(expandHighlights(k.threshold, k) == doctest::Approx(k.threshold));
}

TEST_CASE("knee lifts white to the peak") {
    Knee k;
    CHECK(expandHighlights(1.0f, k) == doctest::Approx(k.peak).epsilon(1e-4).scale(0));
}

TEST_CASE("knee has unit slope at the threshold so there is no crease") {
    Knee k;
    const float h = 1e-4f;
    const float slope = (expandHighlights(k.threshold + h, k) - k.threshold) / h;
    CHECK(slope == doctest::Approx(1.0f).epsilon(1e-2).scale(0));
}

TEST_CASE("knee is strictly increasing") {
    Knee k;
    float prev = -1.0f;
    for (int i = 0; i <= 200; ++i) {
        const float y = expandHighlights(float(i) / 100.0f, k);  // spans 0..2
        CHECK(y > prev);
        prev = y;
    }
}

TEST_CASE("knee inverts exactly") {
    Knee k;
    for (float v : {0.1f, 0.85f, 0.9f, 0.99f, 1.0f, 1.5f})
        CHECK(compressHighlights(expandHighlights(v, k), k) == doctest::Approx(v).epsilon(1e-4).scale(0));
}

TEST_CASE("degenerate knee: threshold at 1.0 returns identity") {
    Knee k{1.0f, 8.0f};
    for (float v : {0.0f, 0.5f, 1.0f, 1.5f}) {
        float expanded = expandHighlights(v, k);
        float compressed = compressHighlights(v, k);
        CHECK(std::isfinite(expanded));
        CHECK(std::isfinite(compressed));
        CHECK(expanded == v);
        CHECK(compressed == v);
    }
}

TEST_CASE("degenerate knee: threshold at 0.0 returns identity") {
    Knee k{0.0f, 8.0f};
    for (float v : {0.0f, 0.5f, 1.0f, 1.5f}) {
        float expanded = expandHighlights(v, k);
        float compressed = compressHighlights(v, k);
        CHECK(std::isfinite(expanded));
        CHECK(std::isfinite(compressed));
        CHECK(expanded == v);
        CHECK(compressed == v);
    }
}

TEST_CASE("degenerate knee: peak equals threshold returns identity") {
    Knee k{0.85f, 0.85f};
    for (float v : {0.0f, 0.5f, 1.0f, 1.5f}) {
        float expanded = expandHighlights(v, k);
        float compressed = compressHighlights(v, k);
        CHECK(std::isfinite(expanded));
        CHECK(std::isfinite(compressed));
        CHECK(expanded == v);
        CHECK(compressed == v);
    }
}

TEST_CASE("degenerate knee: peak below 1.0 returns identity") {
    Knee k{0.85f, 0.5f};
    for (float v : {0.0f, 0.5f, 1.0f, 1.5f}) {
        float expanded = expandHighlights(v, k);
        float compressed = compressHighlights(v, k);
        CHECK(std::isfinite(expanded));
        CHECK(std::isfinite(compressed));
        CHECK(expanded == v);
        CHECK(compressed == v);
    }
}
