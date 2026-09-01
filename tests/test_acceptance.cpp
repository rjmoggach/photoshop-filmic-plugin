#include <doctest/doctest.h>
#include "lenscore/constants.hpp"
#include "lenscore/pipeline.hpp"
#include "metrics.hpp"
#include "targets.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace lens;
using namespace lens::metrics;
using namespace lens::targets;

static const color::SpecTable& table() {
    static color::SpecTable t = color::buildTable(12);
    return t;
}

// psfGrid=256/psfKernel=7 is the combination already proven, across
// test_pipeline.cpp, to satisfy the psfKernel-vs-grid consistency condition
// documented in params.hpp (psfKernel*pixelPitchUm <= psfGrid *
// psfSampleSpacingUm(lambda, fNumberWide) * pupilFill) at every sampled
// wavelength in the visible range, at the default pixelPitchUm/pupilFill.
// The brief this file implements predates that runtime-enforced condition
// and its literal psfGrid=64/psfKernel=21 throws std::invalid_argument.
static Params base() {
    Params p;
    p.bands = 7;
    p.highlightRecovery = false;
    p.doLateralCa = false; p.doPsf = false; p.doVignette = false;
    p.psfGrid = 256; p.psfRings = 8; p.psfKernel = 7; p.effPatch = 64;
    p.dispersion.residual = 0.0f;
    return p;
}

// Anchor-clamp: keeps the whole roi x roi box inside the frame, rather than
// letting per-pixel writes/reads silently drop or repeat whatever falls
// outside it. frameOf's corner convention puts t=1.0 exactly on the image's
// last pixel (see geometry.hpp), so a box of any nonzero size CENTRED there
// always runs off-canvas on two sides; clamping per pixel (rather than the
// box's placement) would then write only a fraction of the edge pattern and
// read back a plateau of repeated boundary pixels -- corrupting exactly the
// corner measurement the acceptance test cares about most. Clamping the
// anchor instead only changes anything once the box would otherwise clip
// (verified: at t in {0, 0.35, 0.7} for these sizes it is a no-op).
static int clampedAnchor(int centre, int box, int full) {
    return std::clamp(centre - box / 2, 0, full - box);
}

// An edge patch placed at a chosen normalised field radius. azimuthRad
// defaults to the diagonal (45 degrees) that every other acceptance test
// uses -- t=1.0 there is the true image corner. The astigmatism test below
// passes 0.0 (the +x axis) instead: see its comment for why.
static Image edgeAt(const Params& p, int full, float t, bool vertical,
                     float azimuthRad = kPi * 0.25f) {
    Plane big(full, full);
    const Frame f = frameOf(full, full);
    const int roi = 64;
    const int cx = int(f.cx + t * f.halfDiag * std::cos(azimuthRad));
    const int cy = int(f.cy + t * f.halfDiag * std::sin(azimuthRad));
    const int ax = clampedAnchor(cx, roi, full), ay = clampedAnchor(cy, roi, full);
    const Plane e = vertical ? slantedEdge(roi, roi, 5.0f, 0.05f, 0.95f)
                             : slantedEdge(roi, roi, 85.0f, 0.05f, 0.95f);
    for (int y = 0; y < roi; ++y)
        for (int x = 0; x < roi; ++x)
            big.at(ax + x, ay + y) = e.at(x, y);
    return render(toImage(big), p, table());
}

// Must use the same azimuth and clamped anchor as edgeAt's write, or the read
// window disagrees with where the edge actually landed (see clampedAnchor's
// comment) and samples a mismatched region instead of the rendered edge.
static Plane roiAt(const Image& im, float t, int roi, float azimuthRad = kPi * 0.25f) {
    const Frame f = frameOf(im.w, im.h);
    const int cx = int(f.cx + t * f.halfDiag * std::cos(azimuthRad));
    const int cy = int(f.cy + t * f.halfDiag * std::sin(azimuthRad));
    const int ax = clampedAnchor(cx, roi, im.w), ay = clampedAnchor(cy, roi, im.h);
    return crop(luminance(im), ax, ay, roi, roi);
}

TEST_CASE("ACCEPTANCE: lateral CA fringe width grows linearly with field radius") {
    Params p = base();
    p.doLateralCa = true;
    p.lateralK = 6e-5f;

    float prev = 1.0f;   // sentinel above any real value; growth here means decreasing (more negative)
    std::vector<float> widths;
    for (float t : {0.0f, 0.35f, 0.7f, 1.0f}) {
        const Image im = edgeAt(p, 320, t, true);
        const Frame f = frameOf(320, 320);
        const int cx = int(f.cx + t * f.halfDiag * 0.70710678f);
        const int cy = int(f.cy + t * f.halfDiag * 0.70710678f);
        const int ax = clampedAnchor(cx, 48, 320), ay = clampedAnchor(cy, 48, 320);
        Image roi(48, 48);
        for (int y = 0; y < 48; ++y)
            for (int x = 0; x < 48; ++x)
                for (int c = 0; c < 3; ++c)
                    roi.at(x, y, c) = im.at(ax + x, ay + y, c);
        // No abs(): fringeWidthPx's sign is the direction of the fringe (R
        // edge position minus B edge position), and a sign-flipped lateralK
        // would still grow in MAGNITUDE with radius and pass a magnitude-
        // only check. At this lateralK/lambdaHat the true signed value is
        // consistently NEGATIVE (measured: 0.0006, -0.16, -0.64, -0.88 at
        // t=0/0.35/0.7/1.0 -- growing in magnitude, same sign throughout,
        // not a mid-series flip), so growth means more negative, not more
        // positive; asserting a positive-going magnitude here would pass
        // just as readily for a sign-reversed k_l as for the real one.
        const float wpx = fringeWidthPx(roi);
        CAPTURE(t); CAPTURE(wpx);
        CHECK(wpx < prev + 0.02f);      // grows (more negative) with radius
        CHECK(wpx <= 0.02f);            // and does not flip direction
        widths.push_back(wpx);
        prev = wpx;
    }
    CHECK(std::abs(widths.front()) < 0.10f);   // no fringe on axis
    CHECK(widths.back() < -0.40f);             // real fringe at the corner, same direction
}

TEST_CASE("ACCEPTANCE: field curvature softens the corners relative to the centre") {
    Params p = base();
    p.doPsf = true;
    p.petzval = 1.2f;

    const float centre = mtf50(roiAt(edgeAt(p, 256, 0.0f, true), 0.0f, 48));
    const float corner = mtf50(roiAt(edgeAt(p, 256, 1.0f, true), 1.0f, 48));
    CAPTURE(centre); CAPTURE(corner);
    CHECK(corner < centre);
    CHECK(corner < 0.85f * centre);
}

TEST_CASE("ACCEPTANCE: astigmatism separates sagittal from tangential resolution") {
    // fitEdge/mtf50 scan rows and locate the edge along x, so they only apply
    // to a near-vertical edge; the 85-degree ("horizontal") ROI must go
    // through transpose first (verified standalone in test_metrics.cpp
    // against the analytic 0.1874/sigma result) or the number it returns is
    // the metric failing on an input it cannot read, not a blur measurement.
    //
    // Pure astigmatism at the medial focus (petzval == 0) is four-fold
    // symmetric -- the two meridians only separate once defocus displaces
    // the image plane toward one line focus. p.petzval supplies that offset.
    Params p = base();
    p.doPsf = true;
    p.astig = 1.2f;
    p.petzval = 0.8f;

    // Measuring at the shared 45-degree diagonal (like every other test
    // here) does NOT work, even with a valid metric and defocus present:
    // astigmatism's principal axes are RADIAL/TANGENTIAL, which at 45
    // degrees fall exactly on the diagonal, not on absolute x/y. Defocus
    // (Z(2,0), no theta dependence) cannot break that alignment, so a
    // vertical-edge (x-resolution) and a transposed horizontal-edge
    // (y-resolution) measurement there are, by that symmetry, nearly equal
    // regardless of how much defocus is added. Measured directly: sag =
    // 0.1357, tan = 0.1360 -- a 0.23% gap, nowhere near a real separation
    // (the tiny residual is the pupil's 9-blade aperture, which is not
    // itself 45-degree-symmetric). That is a genuine finding about this
    // measurement geometry, not a defect in the astigmatism physics --
    // recorded in the report rather than papered over here.
    //
    // Placing the field point on the +x axis instead (azimuthRad = 0) makes
    // radial = x and tangential = y directly, so vertical/horizontal edges
    // measure the real meridians without needing an interpolated image
    // rotation (which mtf50 was never validated against). t=0.7 is as far
    // out along +x as this frame allows without edgeAt's anchor-clamp
    // kicking in and shifting the measured position inward.
    constexpr float azimuth = 0.0f;
    const float sag = mtf50(roiAt(edgeAt(p, 256, 0.7f, true, azimuth), 0.7f, 48, azimuth));
    const float tan = mtf50(transpose(roiAt(edgeAt(p, 256, 0.7f, false, azimuth), 0.7f, 48, azimuth)));
    CAPTURE(sag); CAPTURE(tan);
    CHECK(std::abs(sag - tan) / std::max(sag, tan) > 0.08f);
}

TEST_CASE("ACCEPTANCE: an achromat fringes green and magenta, not red and blue") {
    // Corrected at the F and C lines: red and blue focus together, green does not.
    //
    // eF == eG == 0 to within float noise is a mathematical guarantee of
    // Lagrange interpolation through its own nodes (focusError's correction
    // term is built to match rawFocusError exactly AT 486.1nm and 656.3nm,
    // for any B/C coefficients at all) -- it is not, by itself, evidence
    // this is a real achromat. The physical content is eG being genuinely
    // nonzero (green is NOT a node, so nothing forces it to zero) and the
    // singlet-vs-achromat comparison below; both are checked and both are
    // real measurements, not restatements of the interpolation identity.
    Params p = base();
    p.doPsf = true;
    p.dispersion.correction_nm = {486.1f, 656.3f};
    p.dispersion.residual = 1.0f;
    p.bands = 11;

    const float eF = optics::focusError(p.dispersion, 486.1f, p.lambdaHat);
    const float eC = optics::focusError(p.dispersion, 656.3f, p.lambdaHat);
    const float eG = optics::focusError(p.dispersion, 546.1f, p.lambdaHat);
    CAPTURE(eF); CAPTURE(eC); CAPTURE(eG);

    CHECK(std::abs(eF) < 1e-6f);                       // guaranteed by construction, see above
    CHECK(std::abs(eC) < 1e-6f);                       // guaranteed by construction, see above
    CHECK(std::abs(eG) > 1e-5f);                       // green is the odd one out -- the real test

    // The uncorrected singlet does the opposite: a monotonic red-to-blue ramp.
    Params s = p;
    s.dispersion.correction_nm.clear();
    const float sF = optics::focusError(s.dispersion, 486.1f, s.lambdaHat);
    const float sC = optics::focusError(s.dispersion, 656.3f, s.lambdaHat);
    CHECK(std::abs(sF - sC) > std::abs(eF - eC));
}

TEST_CASE("ACCEPTANCE: stopping down removes mechanical vignetting from the render") {
    Params p = base();
    p.doVignette = true;
    p.vignette.tStop = 2.0f;
    const float wide = radialMean(luminance(render(toImage(flatField(128, 128, 1.0f)), p, table())), 8).back();
    p.vignette.tStop = 8.0f;
    const float stopped = radialMean(luminance(render(toImage(flatField(128, 128, 1.0f)), p, table())), 8).back();
    CHECK(stopped > wide);
    // naturalFalloff(...) at t=0.94 is well below 1.0; the Approx epsilon must
    // be scaled to the value itself, not doctest's default scale=1.0, or a
    // nominal 5% tolerance is really far looser (see the global .scale(0) rule).
    CHECK(stopped == doctest::Approx(optics::naturalFalloff(p.vignette, 0.94f)).epsilon(0.05).scale(0));
}
