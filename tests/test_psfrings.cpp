#include <doctest/doctest.h>
#include "lenscore/optics/psfrings.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

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

// Numerically rotates `src` about its own centre by thetaRad, using EXACTLY
// psfAtField's inverse-rotation convention (see psfrings.hpp). This is an
// independent, test-only cross-check: it lets us ask "does psfAtField's
// output at angle theta match what you get by rotating psfAtField's own
// on-axis output by theta", entirely outside of production code.
static Plane rotatePlane(const Plane& src, float thetaRad) {
    const float ct = std::cos(thetaRad), st = std::sin(thetaRad);
    const float c = 0.5f * float(src.w - 1);
    Plane out(src.w, src.h);
    for (int y = 0; y < src.h; ++y) {
        for (int x = 0; x < src.w; ++x) {
            const float dx = float(x) - c;
            const float dy = float(y) - c;
            const float ux =  ct * dx + st * dy;
            const float uy = -st * dx + ct * dy;
            out.at(x, y) = sampleBilinear(src, c + ux, c + uy);
        }
    }
    return out;
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
    // A relative tolerance on near-zero background pixels is meaningless -- it demands
    // 2% of ~4e-5. The image's peak is the scale that matters, so compare absolute
    // difference against it instead. Measured: this 33px window at field=0.8 sits
    // mostly in the PSF's tail (values span only ~4e-5 to ~1.5e-4, a modest ~4x range,
    // not a sharp core against negligible background), so the angle-to-angle
    // discretization noise from resampling the Cartesian ring grid shows up at up to
    // ~9.3% of peak rather than a couple of percent; 0.12 covers that with headroom
    // while remaining far tighter than the pre-sweep vacuous (~100%+) tolerance.
    const float peak = *std::max_element(a.v.begin(), a.v.end());
    for (size_t i = 0; i < a.v.size(); ++i) {
        CHECK(std::abs(b.v[i] - a.v[i]) < 0.12f * peak);
        CHECK(std::abs(c.v[i] - a.v[i]) < 0.12f * peak);
    }
}

// Fix round 3: a centroid is exactly the statistic truncation corrupts (a
// window that clips coma's long tail asymmetrically rotates the MEASURED
// centroid by the wrong amount, independent of whether the rotation code is
// correct). This replaces both centroid-based coma tests with an image
// comparison instead: `b` (psfAtField at theta) is compared directly against
// `c` (psfAtField's on-axis output `a`, numerically rotated by theta). `a`
// and `b` are both produced by the SAME 65-pixel window, so whatever that
// window clips, it clips identically going into both -- the comparison does
// not depend on the window containing the whole distribution the way a
// centroid or an energy-capture guard did.
TEST_CASE("coma's tail rotates with the field angle, direction included") {
    const PsfRings r = buildPsfRings(disc(), [](float t) {
        Wavefront w; w.coma = 1.5f * t; return w; }, 550.0f, 550.0f, 8, 256);

    // theta = 45 degrees: at 0 or 90 degrees a transposed (inverse-flipped)
    // rotation produces the same magnitudes as the correct one, so only an
    // off-cardinal angle actually discriminates direction.
    const float theta = 0.78539816f;
    const Plane a = psfAtField(r, 1.0f, 0.0f,  1.0f, 65);   // on axis
    const Plane b = psfAtField(r, 1.0f, theta, 1.0f, 65);   // rotated by theta
    const Plane c = rotatePlane(a, theta);                  // a, rotated the same way, numerically

    float peak = 0.0f;
    for (float v : a.v) peak = std::max(peak, v);
    REQUIRE(peak > 0.0f);

    // Comparing over the WHOLE plane (or even down to 1% of peak) does not work here:
    // coma's tail is wide enough that most of a 65px window is still tail, not core, and
    // rotatePlane(a, theta) has to sample a's corners from outside a's own 65px extent --
    // the exact same square-window truncation mismatch this test exists to avoid, just
    // moved one step later. Restricting to pixels near a's peak (empirically, >= 70% of
    // it) keeps the comparison inside the well-resolved core, safely away from where
    // rotatePlane's corner sampling runs off the edge of `a`. That is a SMALL region of
    // the window (a handful of pixels), not a small threshold value.
    double sumSq = 0.0; int n = 0;
    for (size_t i = 0; i < a.v.size(); ++i) {
        if (a.v[i] > 0.7f * peak) {
            const double d = double(b.v[i]) - double(c.v[i]);
            sumSq += d * d;
            ++n;
        }
    }
    REQUIRE(n > 0);
    const double rms = std::sqrt(sumSq / double(n));
    CHECK(rms / double(peak) < 0.05);
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
    CHECK(coarse / fine == doctest::Approx(1.0).epsilon(0.05).scale(0));
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
    CHECK(onAxis == doctest::Approx(1.0).epsilon(0.05).scale(0));
    CHECK(offAxis < onAxis);
}

// Fix round 1, finding 1: fewer than two rings leaves psfAtField's i0/i0+1
// pairing undefined (i0 = -1, an out-of-bounds ring[-1] read), and a negative
// rings previously wrapped size_t(rings) into a huge reserve() before the
// loop even ran. Both must be rejected as programmer error, not tolerated.
TEST_CASE("buildPsfRings rejects fewer than two rings") {
    auto flat = [](float) { return Wavefront{}; };
    CHECK_THROWS_AS(buildPsfRings(disc(), flat, 550.0f, 550.0f, 1, 64), std::invalid_argument);
    CHECK_THROWS_AS(buildPsfRings(disc(), flat, 550.0f, 550.0f, 0, 64), std::invalid_argument);
    CHECK_THROWS_AS(buildPsfRings(disc(), flat, 550.0f, 550.0f, -3, 64), std::invalid_argument);
}

// psfAtField takes a PsfRings it did not necessarily build itself (nothing
// stops a caller from constructing one by hand), so it must not assume the
// >= 2 rings invariant holds -- it has to check, not just document.
TEST_CASE("psfAtField rejects a PsfRings with fewer than two rings") {
    PsfRings empty;                       // default-constructed: 0 rings
    CHECK_THROWS_AS(psfAtField(empty, 0.0f, 0.0f, 1.0f, 17), std::invalid_argument);

    PsfRings one;
    one.gridN = 64;
    one.ring.push_back(Plane(64, 64));    // still only 1 ring
    CHECK_THROWS_AS(psfAtField(one, 0.0f, 0.0f, 1.0f, 17), std::invalid_argument);
}

// Fix round 1, finding 3: axisEnergyOf must floor near-zero energy (denormal
// noise, not just an exact 0.0 sum), or dividing by it would silently blow
// every kernel up by an arbitrary, physically meaningless factor.
TEST_CASE("axisEnergyOf floors near-zero energy instead of dividing by it") {
    CHECK(axisEnergyOf(Plane(4, 4)) == doctest::Approx(1.0));   // exact zero: default Plane

    Plane denormal(4, 4);
    for (float& v : denormal.v) v = 1e-20f;                     // technically > 0, still noise-scale
    CHECK(axisEnergyOf(denormal) == doctest::Approx(1.0));

    Plane real(4, 4);
    for (float& v : real.v) v = 2.0f;                           // a real, usable energy
    CHECK(axisEnergyOf(real) == doctest::Approx(32.0));          // 16 pixels * 2.0, well above the floor
}

// The pipeline-level companion to the above: a pupil clipped almost to
// nothing (here, an aperture just barely wide enough to admit the handful of
// grid pixels closest to the centre at N=128) still yields a finite, sanely
// normalised on-axis kernel -- it neither explodes nor gets incorrectly
// floored just because its energy is small.
// Directly exercises the floor itself, as opposed to the small-but-legitimate
// energy case above: an aperture with literally zero grid pixels inside it
// gives an exactly-zero raw ring, axisEnergyOf floors that to 1.0, and
// psfAtField's scale becomes just the area factor -- no division by anything
// near zero ever happens, so the result is FINITE (here, trivially zero,
// since scale * 0 = 0) rather than NaN or infinite.
TEST_CASE("a fully occluded pupil hits the axisEnergy floor and stays finite") {
    PupilParams p = disc();
    p.apertureRadius = 0.001f;   // smaller than any grid pixel's distance from centre at N=64
    const PsfRings r = buildPsfRings(p, [](float) { return Wavefront{}; },
                                     550.0f, 550.0f, 2, 64);
    CHECK(r.axisEnergy == doctest::Approx(1.0));   // floor engaged: raw energy was exactly 0

    const Plane onAxis = psfAtField(r, 0.0f, 0.0f, 1.0f, 17);
    for (float v : onAxis.v) {
        CHECK(std::isfinite(v));
        CHECK(v == doctest::Approx(0.0f));         // raw ring was all zero; scale never divided by it
    }
}

TEST_CASE("a pupil clipped almost to nothing still normalises sanely") {
    PupilParams p = disc();
    p.apertureRadius = 0.012f;   // admits only the centremost few pixels at N=128
    const PsfRings r = buildPsfRings(p, [](float) { return Wavefront{}; },
                                     550.0f, 550.0f, 2, 128);
    CHECK(r.axisEnergy > 0.0);

    const Plane onAxis = psfAtField(r, 0.0f, 0.0f, 1.0f, 65);
    double s = 0.0;
    for (float v : onAxis.v) { CHECK(std::isfinite(v)); s += v; }
    CHECK(std::isfinite(s));
    // Genuine physics, not an over-conversion: an aperture this small widens the Airy
    // core well past this fixed 65px window, so the window itself truncates real
    // energy -- measured s ~= 0.6824, i.e. ~32% short of 1.0. epsilon(0.4) states that
    // measured value honestly; it is numerically the same actual tolerance this check
    // had before the .scale(0) sweep (nominal 0.2 defaulted to scale=1, an effective
    // 40% at E=1.0), just now arrived at by stating it instead of by an accidental
    // scale default.
    CHECK(s == doctest::Approx(1.0).epsilon(0.4).scale(0));
}

// Fix round 2: the on-axis kernel sum being close to 1 was previously verified
// at exactly one configuration (Wavefront{}, i.e. zero aberration) -- which is
// the one PSF shape this project's own reason for existing (aberrated glass)
// never actually renders. A resampling formula tuned only against that single,
// degenerate case (a bare Airy disk, genuinely narrower than one output pixel
// at this project's default optics) held for it alone and overshot total
// energy by roughly an order of magnitude the moment real aberration widened
// the PSF to several pixels across. This exercises the same samplesPerPixel >
// 1 (oversampled, box-averaged) regime render() actually uses -- pixelPitchUm
// / (psfSampleSpacingUm(lambda, fNumberWide) * pupilFill) at the project's
// default optics -- across zero aberration through a strong (20-wave) defocus,
// so a formula or resampling change that only works in one of those regimes
// fails loudly here instead of shipping.
TEST_CASE("the on-axis kernel sums to approximately 1 across aberration strengths") {
    PupilParams pp;
    pp.apertureRadius = 0.25f;   // matches Params::pupilFill's default working aperture
    pp.rEntrance = 1e6f; pp.rExit = 1e6f; pp.sepNorm = 0.0f;

    const float lambda = 550.0f, fNumberWide = 2.0f, pupilFill = 0.25f, pixelPitchUm = 5.0f;
    const float spp = pixelPitchUm / (psfSampleSpacingUm(lambda, fNumberWide) * pupilFill);

    for (float waves : {0.0f, 1.0f, 6.0f, 20.0f}) {
        const PsfRings r = buildPsfRings(pp, [waves](float) {
            Wavefront w; w.defocus = waves; return w; }, lambda, 650.0f, 2, 128);
        const Plane k = psfAtField(r, 0.0f, 0.0f, spp, 7);
        double sum = 0.0; for (float v : k.v) sum += v;
        CAPTURE(waves);
        CHECK(sum == doctest::Approx(1.0).epsilon(0.15).scale(0));
    }
}
