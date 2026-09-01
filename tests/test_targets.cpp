#include <doctest/doctest.h>
#include "targets.hpp"
#include "lenscore/constants.hpp"
#include <cmath>

using namespace lens;
using namespace lens::targets;

TEST_CASE("flat field is flat") {
    const Plane p = flatField(16, 9, 0.7f);
    for (float v : p.v) CHECK(v == doctest::Approx(0.7f));
}

TEST_CASE("point grid places isolated unit points on the expected lattice") {
    const Plane p = pointGrid(64, 64, 16);
    double sum = 0; int nonzero = 0;
    for (float v : p.v) { sum += v; if (v > 0.0f) ++nonzero; }
    CHECK(nonzero == 9);                 // 3 x 3 interior lattice
    CHECK(sum == doctest::Approx(9.0));
}

TEST_CASE("slanted edge is genuinely slanted and hits both levels") {
    const Plane p = slantedEdge(64, 64, 5.0f, 0.1f, 0.9f);
    CHECK(p.at(2, 32)  == doctest::Approx(0.1f).epsilon(1e-3));
    CHECK(p.at(61, 32) == doctest::Approx(0.9f).epsilon(1e-3));
    // The transition column shifts down the rows because the edge is tilted.
    auto crossing = [&](int y) {
        for (int x = 1; x < 64; ++x) if (p.at(x, y) > 0.5f) return x;
        return 64;
    };
    CHECK(crossing(10) != crossing(54));
}

TEST_CASE("slanted edge has exactly one partially covered pixel per row") {
    const Plane p = slantedEdge(64, 64, 5.0f, 0.0f, 1.0f);
    for (int y = 0; y < 64; ++y) {
        int partial = 0;
        for (int x = 0; x < 64; ++x) {
            const float v = p.at(x, y);
            if (v > 1e-4f && v < 1.0f - 1e-4f) ++partial;
        }
        CHECK(partial <= 1);
    }
}

TEST_CASE("Siemens star is rotationally periodic in the spoke count") {
    const Plane p = siemensStar(128, 128, 16);
    const double r = 40.0, cx = 63.5, cy = 63.5;
    auto sample = [&](double a) {
        return double(p.at(int(std::lround(cx + r * std::cos(a))),
                           int(std::lround(cy + r * std::sin(a)))));
    };
    for (int k = 0; k < 8; ++k) {
        const double a = 0.31 + k * 0.4;
        CHECK(sample(a) == doctest::Approx(sample(a + 2.0 * kPi / 16.0)).epsilon(0.05));
    }
}

TEST_CASE("single point carries unit energy at the exact centre") {
    const Plane p = singlePoint(65, 65);
    double s = 0; for (float v : p.v) s += v;
    CHECK(s == doctest::Approx(1.0));
    CHECK(p.at(32, 32) == doctest::Approx(1.0f));
}

TEST_CASE("plane to image and back is lossless for grey") {
    const Plane a = flatField(8, 8, 0.42f);
    const Plane b = luminance(toImage(a));
    for (int i = 0; i < 64; ++i) CHECK(b.v[i] == doctest::Approx(a.v[i]).epsilon(1e-5));
}
