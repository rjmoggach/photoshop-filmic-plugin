#include <doctest/doctest.h>
#include "lenscore/optics/lateralca.hpp"
#include <cmath>

using namespace lens;
using namespace lens::optics;

static Plane dot(int w, int h, int px, int py) {
    Plane p(w, h);
    p.at(px, py) = 1.0f;
    return p;
}

TEST_CASE("frame puts t = 1 at the corner and 0 at the centre") {
    const Frame f = frameOf(101, 51);
    CHECK(f.cx == doctest::Approx(50.0f));
    CHECK(f.cy == doctest::Approx(25.0f));
    const float t = std::sqrt(f.cx * f.cx + f.cy * f.cy) / f.halfDiag;
    CHECK(t == doctest::Approx(1.0f));
}

TEST_CASE("inverse radius exactly inverts the forward magnification") {
    for (float K : {-0.05f, -0.001f, 0.0f, 0.002f, 0.08f}) {
        for (float tOut : {0.0f, 0.25f, 0.5f, 1.0f}) {
            const float ts = inverseLateralRadius(K, tOut);
            CHECK(ts * (1.0f + K * ts) == doctest::Approx(tOut).epsilon(1e-5));
        }
    }
}

TEST_CASE("magnification is unity at the reference wavelength") {
    CHECK(lateralMagnification(0.02f, 650.0f, 650.0f, 1.0f) == doctest::Approx(1.0f));
}

TEST_CASE("blue is magnified differently from red, and only off axis") {
    const float mBlue = lateralMagnification(1e-4f, 450.0f, 650.0f, 1.0f);
    CHECK(mBlue != doctest::Approx(1.0f));
    CHECK(lateralMagnification(1e-4f, 450.0f, 650.0f, 0.0f) == doctest::Approx(1.0f));
}

TEST_CASE("identity warp reproduces the source") {
    Plane src(33, 21);
    for (int y = 0; y < 21; ++y)
        for (int x = 0; x < 33; ++x) src.at(x, y) = float((x * 7 + y * 3) % 11) / 11.0f;
    const Plane out = warpPlane(src, Distortion{}, 0.0f);
    for (int y = 0; y < 21; ++y)
        for (int x = 0; x < 33; ++x) CHECK(out.at(x, y) == doctest::Approx(src.at(x, y)).epsilon(1e-5));
}

TEST_CASE("positive K pushes a feature outward from the centre") {
    const int w = 129, h = 129;
    const Plane src = dot(w, h, 100, 64);          // right of centre
    // K = 0.01 (the original brief value) cannot move the argmax pixel at all here:
    // the largest possible shift anywhere on this row is K * 0.5 * halfDiag =~ 0.45px
    // (t maxes out at cx/halfDiag =~ 0.707 at the row's own edges, short of the corner's
    // t = 1), which never crosses the 0.5px needed to flip the discrete argmax -- true
    // for every point on this row, not just x = 100. K = 0.08 (already exercised by the
    // inverseLateralRadius test above) clears that bound with margin.
    const Plane out = warpPlane(src, Distortion{}, 0.08f);
    float best = -1.0f; int bx = 0;
    for (int x = 0; x < w; ++x) if (out.at(x, 64) > best) { best = out.at(x, 64); bx = x; }
    CHECK(bx > 100);
}

TEST_CASE("barrel distortion pulls the corners inward") {
    Distortion d; d.k1 = -0.10f;
    float dx = 0, dy = 0;
    applyDistortion(d, 0.7071f, 0.7071f, dx, dy);
    CHECK(std::sqrt(dx * dx + dy * dy) < 1.0f);
}

TEST_CASE("clamped edge sampling stays bounded, and inflates energy in a documented direction") {
    // The previous version of this test used a 25x25 block spanning [20,45)
    // in a 65x65 plane (halfDiag =~ 45.25) with K = 0.005f -- a region and a
    // K small enough that sampleBilinear's edge clamp never actually
    // triggers, so a passing result proved nothing about clamping. Positive K
    // can never trigger it either way: the forward map t_out = t_src(1+Kt_src)
    // is increasing for K >= 0, so the source radius pulled for any in-frame
    // output pixel (t_out <= 1) is always <= t_out <= 1, i.e. still in frame.
    //
    // To exercise the clamp for real, this block instead reaches the plane's
    // own edge (touching both the right and bottom boundary), and K is
    // negative and large enough that, near that corner, the inverted source
    // radius exceeds 1 -- outside the frame -- forcing sampleBilinear to read
    // the boundary pixel instead of the (nonexistent) pixel beyond it.
    Plane src(65, 65);
    for (int y = 50; y < 65; ++y) for (int x = 50; x < 65; ++x) src.at(x, y) = 1.0f;
    const Plane out = warpPlane(src, Distortion{}, -0.1f);

    // This holds unconditionally, clamp or no clamp: sampleBilinear returns a
    // convex combination of up to four source samples, so no output value can
    // exceed the brightest source value. That is the property that is
    // actually true here, and it is the one this test leans on.
    for (float v : out.v) CHECK(v <= 1.0f + 1e-5f);

    double a = 0, b = 0;
    for (float v : src.v) a += v;
    for (float v : out.v) b += v;
    // Exact conservation is NOT true once clamping engages: the block sits
    // right on the boundary, so every clamped query near that corner reads
    // back into the block (1.0) rather than off the edge into nothing. That
    // can only add energy, never remove it, so total energy is expected to
    // measurably exceed the source's -- and, since this is duplication of a
    // bounded region rather than a runaway effect, to stay within a modest
    // multiple of it.
    CHECK(b > a);
    CHECK(b < 2.0 * a);
}
