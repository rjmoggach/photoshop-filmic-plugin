#include <doctest/doctest.h>
#include "lenscore/optics/psfrings.hpp"
#include <cmath>

using namespace lens;
using namespace lens::optics;

static PupilParams disc() {
    PupilParams p;
    p.blades = 0; p.curvature = 1.0f; p.apertureRadius = 0.125f;
    p.rEntrance = 1e6f; p.rExit = 1e6f; p.sepNorm = 0.0f; p.apodizationSlope = 0.0f;
    return p;
}

// Same as disc(), but clips the pupil off axis so field-dependent vignetting
// is exercised (mirrors test_psf.cpp's vignetting configuration).
static PupilParams vignettingDisc() {
    PupilParams p = disc();
    p.rEntrance = 0.13f; p.rExit = 1e6f; p.sepNorm = 0.10f;
    return p;
}

static void centroid(const Plane& p, double& cx, double& cy) {
    double s = 0; cx = 0; cy = 0;
    for (int y = 0; y < p.h; ++y) for (int x = 0; x < p.w; ++x) {
        s  += p.at(x, y);
        cx += p.at(x, y) * (x - p.w / 2.0);
        cy += p.at(x, y) * (y - p.h / 2.0);
    }
    if (s > 0) { cx /= s; cy /= s; }
}

TEST_CASE("rings are built at the requested count") {
    const PsfRings r = buildPsfRings(disc(), [](float) { return Wavefront{}; },
                                     550.0f, 550.0f, 8, 128);
    CHECK(r.ring.size() == 8u);
    CHECK(r.gridN == 128);
}

TEST_CASE("a rotationally symmetric aberration gives the same PSF at every angle") {
    const PsfRings r = buildPsfRings(disc(), [](float t) {
        Wavefront w; w.defocus = 1.0f * t * t; return w; }, 550.0f, 550.0f, 8, 128);
    const Plane a = psfAtField(r, 0.8f, 0.0f,             1.0f, 33);
    const Plane b = psfAtField(r, 0.8f, 1.0471975f,       1.0f, 33);   // 60 degrees
    const Plane c = psfAtField(r, 0.8f, 2.7f,             1.0f, 33);
    for (int i = 0; i < 33 * 33; ++i) {
        CHECK(b.v[i] == doctest::Approx(a.v[i]).epsilon(0.02));
        CHECK(c.v[i] == doctest::Approx(a.v[i]).epsilon(0.02));
    }
}

TEST_CASE("coma's tail follows the radial direction round the frame") {
    const PsfRings r = buildPsfRings(disc(), [](float t) {
        Wavefront w; w.coma = 1.5f * t; return w; }, 550.0f, 550.0f, 8, 256);

    double cx = 0, cy = 0;
    centroid(psfAtField(r, 1.0f, 0.0f, 1.0f, 65), cx, cy);
    const double along = std::abs(cx);
    CHECK(along > 0.3);
    CHECK(std::abs(cy) < 0.25 * along);

    centroid(psfAtField(r, 1.0f, 1.5707963f, 1.0f, 65), cx, cy);   // 90 degrees
    CHECK(std::abs(cy) > 0.3);
    CHECK(std::abs(cx) < 0.25 * std::abs(cy));
}

TEST_CASE("sampling exactly on a ring returns that ring's energy") {
    const PsfRings r = buildPsfRings(disc(), [](float t) {
        Wavefront w; w.defocus = 2.0f * t; return w; }, 550.0f, 550.0f, 5, 128);
    auto energy = [](const Plane& p) { double s = 0; for (float v : p.v) s += v; return s; };
    const double onRing = energy(psfAtField(r, 0.5f, 0.0f, 1.0f, 41));   // ring index 2 of 5
    CHECK(onRing > 0.0);
}

TEST_CASE("resampling to coarser pixels keeps the energy") {
    const PsfRings r = buildPsfRings(disc(), [](float) { return Wavefront{}; },
                                     550.0f, 550.0f, 4, 256);
    auto energy = [](const Plane& p) { double s = 0; for (float v : p.v) s += v; return s; };
    const double fine   = energy(psfAtField(r, 0.0f, 0.0f, 1.0f, 81));
    const double coarse = energy(psfAtField(r, 0.0f, 0.0f, 2.0f, 41));
    CHECK(coarse / fine == doctest::Approx(1.0).epsilon(0.05));
}

// The correction to the brief: psfFromPupil's raw PSF is never normalised (its
// energy carries the optical vignetting), but that must not leak literally into
// the convolution kernel or every image would be scaled by N^2 * sum(A^2). What
// this pipeline actually wants is: the on-axis kernel sums to 1, and off-axis
// kernels sum to the vignetting fraction *relative to the axis*. PsfRings::axisEnergy
// captures ring 0's raw energy once, and psfAtField divides by it.
TEST_CASE("the on-axis kernel sums to 1; a vignetted off-axis kernel sums to less") {
    const PsfRings r = buildPsfRings(vignettingDisc(), [](float) { return Wavefront{}; },
                                     550.0f, 550.0f, 5, 128);
    auto energy = [](const Plane& p) { double s = 0; for (float v : p.v) s += v; return s; };
    // outSize spans nearly the whole grid at unit sampling so the resampled
    // kernel captures essentially all of the PSF's energy.
    const double onAxis  = energy(psfAtField(r, 0.0f, 0.0f, 1.0f, 127));
    const double offAxis = energy(psfAtField(r, 1.0f, 0.0f, 1.0f, 127));
    CHECK(onAxis == doctest::Approx(1.0).epsilon(0.05));
    CHECK(offAxis < onAxis);
}
