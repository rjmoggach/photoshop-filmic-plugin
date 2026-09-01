#include <doctest/doctest.h>
#include "lenscore/optics/pupil.hpp"
#include "lenscore/optics/vignette.hpp"
#include <cmath>

using namespace lens;
using namespace lens::optics;

static PupilParams circular() {
    PupilParams p;
    p.blades = 0; p.curvature = 1.0f; p.rotationRad = 0.0f;
    p.apertureRadius = 1.0f; p.rEntrance = 1e6f; p.rExit = 1e6f;
    p.sepNorm = 0.0f; p.apodizationSlope = 0.0f;
    return p;
}

TEST_CASE("a circular aperture has unit edge radius in every direction") {
    for (float th = 0.0f; th < 6.28f; th += 0.3f)
        CHECK(apertureEdgeRadius(0, 1.0f, 0.0f, th) == doctest::Approx(1.0f));
}

TEST_CASE("a straight-sided hexagon touches 1 at vertices and cos(pi/6) at flats") {
    float lo = 2.0f, hi = 0.0f;
    for (float th = 0.0f; th < 6.283f; th += 0.001f) {
        const float r = apertureEdgeRadius(6, 0.0f, 0.0f, th);
        lo = std::min(lo, r); hi = std::max(hi, r);
    }
    CHECK(hi == doctest::Approx(1.0f).epsilon(1e-3));
    CHECK(lo == doctest::Approx(std::cos(3.14159265f / 6.0f)).epsilon(1e-3));
}

TEST_CASE("curvature interpolates the polygon back to a circle") {
    const float straight = apertureEdgeRadius(6, 0.0f, 0.0f, 3.14159265f / 6.0f);
    const float round    = apertureEdgeRadius(6, 1.0f, 0.0f, 3.14159265f / 6.0f);
    const float half     = apertureEdgeRadius(6, 0.5f, 0.0f, 3.14159265f / 6.0f);
    CHECK(round == doctest::Approx(1.0f));
    CHECK(half > straight);
    CHECK(half < round);
}

TEST_CASE("an unclipped circular pupil has full energy on axis") {
    CHECK(pupilEnergyFraction(circular(), 0.0f, 256) == doctest::Approx(1.0f).epsilon(0.01));
}

TEST_CASE("a hexagonal aperture passes the known polygon area ratio") {
    PupilParams p = circular();
    p.blades = 6; p.curvature = 0.0f;
    // regular hexagon in a unit circle: area = 3*sin(60) = 2.598, over pi = 0.827
    CHECK(pupilEnergyFraction(p, 0.0f, 512) == doctest::Approx(0.827f).epsilon(0.02));
}

TEST_CASE("the cat's-eye clip removes energy off axis and not on axis") {
    PupilParams p = circular();
    p.rEntrance = 1.0f; p.rExit = 1.0f; p.sepNorm = 0.5f;
    CHECK(pupilEnergyFraction(p, 0.0f, 256) == doctest::Approx(1.0f).epsilon(0.01));
    CHECK(pupilEnergyFraction(p, 1.0f, 256) < 0.95f);
}

TEST_CASE("stopping down removes the cat's-eye clip entirely") {
    PupilParams p = circular();
    p.rEntrance = 1.0f; p.rExit = 1.0f; p.sepNorm = 0.5f;
    p.apertureRadius = 0.5f;                       // f/4 with a wide-open f/2
    CHECK(pupilEnergyFraction(p, 1.0f, 256) == doctest::Approx(1.0f).epsilon(0.01));
}

TEST_CASE("rasterised pupil energy agrees with the independent closed form") {
    // One clipping circle only, so the two models describe the same geometry.
    VignetteParams v;
    v.tStopWide = 2.0f; v.rEntrance = 1.0f; v.sepNorm = 0.5f;
    for (float stop : {2.0f, 2.4f, 2.8f, 3.2f, 4.0f}) {
        v.tStop = stop;
        PupilParams p = circular();
        p.apertureRadius = v.tStopWide / v.tStop;
        p.rEntrance = v.rEntrance; p.rExit = 1e6f; p.sepNorm = v.sepNorm;
        CAPTURE(stop);
        // Boundary-pixel quantisation at N=512 is on the order of 1/N (~0.2%),
        // so 0.5% is already loose relative to rasterisation error; it is not
        // a relaxation for either model's own imprecision.
        CHECK(pupilEnergyFraction(p, 1.0f, 512) ==
              doctest::Approx(mechanicalFraction(v, 1.0f)).epsilon(0.005));
    }
}

TEST_CASE("apodization shifts the pupil centroid outward without adding energy") {
    PupilParams flat = circular();
    PupilParams ramp = circular(); ramp.apodizationSlope = 0.31f;
    const Plane a = rasterPupil(flat, 1.0f, 256);
    const Plane b = rasterPupil(ramp, 1.0f, 256);
    double sa = 0, sb = 0, cx = 0;
    for (int y = 0; y < 256; ++y)
        for (int x = 0; x < 256; ++x) {
            const double u = (2.0 * x) / 255.0 - 1.0;
            sa += a.at(x, y); sb += b.at(x, y); cx += b.at(x, y) * u;
        }
    CHECK(sb / sa == doctest::Approx(1.0).epsilon(0.01));   // energy preserved
    CHECK(cx / sb > 0.02);                                   // centroid moved outward
}

TEST_CASE("apodization stays energy-neutral at the maximum allowed slope") {
    // At kMaxApodizationSlope (0.9) the ramp's floor, 1 - 0.9 = 0.1, is still
    // strictly positive at t=1, u=-1, so the clamp on out.at() never fires and
    // the ramp's mean over the disc is still exactly 1 -- energy neutrality
    // must hold right up to the boundary, not just at the small slope (0.31)
    // the other apodization test uses.
    PupilParams flat = circular();
    PupilParams ramp = circular(); ramp.apodizationSlope = kMaxApodizationSlope;
    const Plane a = rasterPupil(flat, 1.0f, 256);
    const Plane b = rasterPupil(ramp, 1.0f, 256);
    double sa = 0, sb = 0;
    for (int y = 0; y < 256; ++y)
        for (int x = 0; x < 256; ++x) {
            sa += a.at(x, y); sb += b.at(x, y);
        }
    CHECK(sb / sa == doctest::Approx(1.0).epsilon(0.01));   // energy preserved at the boundary
}
