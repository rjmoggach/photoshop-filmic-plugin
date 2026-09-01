#include <doctest/doctest.h>
#include "lenscore/conv/eff.hpp"
#include <cmath>

using namespace lens;
using namespace lens::conv;

static Plane noise(int w, int h) {
    Plane p(w, h);
    unsigned s = 12345u;
    for (auto& v : p.v) { s = s * 1664525u + 1013904223u; v = float(s >> 8 & 0xFFFF) / 65535.0f; }
    return p;
}

static Plane deltaKernel(int k) { Plane p(k, k); p.at(k / 2, k / 2) = 1.0f; return p; }

static Plane boxKernel(int k) {
    Plane p(k, k);
    for (auto& v : p.v) v = 1.0f / float(k * k);
    return p;
}

TEST_CASE("Hann windows at fifty percent overlap sum to one") {
    const int P = 32;
    const std::vector<float> w = hannWindow(P);
    for (int i = 0; i < P / 2; ++i)
        CHECK(w[i] + w[i + P / 2] == doctest::Approx(1.0f).epsilon(1e-5));
}

TEST_CASE("a delta kernel everywhere is the identity") {
    const Plane src = noise(64, 64);
    const Plane out = effConvolve(src, 32, [](float, float) { return deltaKernel(9); });
    for (int y = 8; y < 56; ++y)
        for (int x = 8; x < 56; ++x)
            CHECK(out.at(x, y) == doctest::Approx(src.at(x, y)).epsilon(1e-3));
}

TEST_CASE("a uniform kernel reduces exactly to plain convolution") {
    const Plane src = noise(64, 64);
    const Plane k = boxKernel(7);
    const Plane want = convolveDirect(src, k);
    const Plane got  = effConvolve(src, 32, [&](float, float) { return k; });
    for (int y = 12; y < 52; ++y)
        for (int x = 12; x < 52; ++x)
            CHECK(got.at(x, y) == doctest::Approx(want.at(x, y)).epsilon(2e-3));
}

TEST_CASE("a normalised kernel preserves energy in the interior") {
    Plane src(64, 64);
    for (int y = 16; y < 48; ++y) for (int x = 16; x < 48; ++x) src.at(x, y) = 1.0f;
    const Plane out = effConvolve(src, 32, [](float, float) { return boxKernel(5); });
    double a = 0, b = 0;
    for (float v : src.v) a += v;
    for (float v : out.v) b += v;
    CHECK(b / a == doctest::Approx(1.0).epsilon(0.01));
}

TEST_CASE("a varying kernel really does vary across the frame") {
    Plane src(128, 128);
    for (int y = 0; y < 128; ++y)
        for (int x = 0; x < 128; ++x) src.at(x, y) = ((x / 4) % 2) ? 1.0f : 0.0f;   // vertical bars

    // Sharp on the left, blurry on the right.
    const Plane out = effConvolve(src, 32, [](float cx, float) {
        return (cx < 64.0f) ? deltaKernel(11) : boxKernel(11);
    });

    auto contrast = [&](const Plane& p, int x0, int x1) {
        double lo = 1e9, hi = -1e9;
        for (int y = 40; y < 88; ++y) for (int x = x0; x < x1; ++x) {
            lo = std::min(lo, double(p.at(x, y))); hi = std::max(hi, double(p.at(x, y)));
        }
        return hi - lo;
    };
    CHECK(contrast(out, 10, 40) > 0.9);
    CHECK(contrast(out, 90, 120) < 0.5);
}
