#include <doctest/doctest.h>
#include "metrics.hpp"
#include "targets.hpp"
#include <cmath>
#include <stdexcept>

using namespace lens;
using namespace lens::metrics;
using namespace lens::targets;

TEST_CASE("gaussian blur preserves energy and is symmetric") {
    Plane p(64, 64); p.at(32, 32) = 1.0f;
    const Plane b = gaussianBlur(p, 2.0f);
    CHECK(totalEnergy(b) == doctest::Approx(1.0).epsilon(0.01));
    CHECK(b.at(30, 32) == doctest::Approx(b.at(34, 32)).epsilon(1e-4));
    CHECK(b.at(32, 30) == doctest::Approx(b.at(32, 34)).epsilon(1e-4));
}

TEST_CASE("mtf50 of a sharp edge is near the Nyquist limit") {
    const Plane e = slantedEdge(96, 96, 5.0f, 0.0f, 1.0f);
    CHECK(mtf50(e) > 0.35f);
}

TEST_CASE("mtf50 recovers the analytic value for a known gaussian") {
    for (float sigma : {1.0f, 1.5f, 2.0f, 3.0f}) {
        const Plane e = gaussianBlur(slantedEdge(128, 128, 5.0f, 0.0f, 1.0f), sigma);
        const float got  = mtf50(crop(e, 24, 24, 80, 80));
        const float want = 0.1874f / sigma;
        CAPTURE(sigma); CAPTURE(got); CAPTURE(want);
        CHECK(got == doctest::Approx(want).epsilon(0.08));
    }
}

TEST_CASE("mtf50 recovers the analytic value with a nonzero dark level") {
    // Contrast is hi - lo, and MTF50 is contrast-normalised, so a nonzero
    // dark level should not change the answer. This is the regression test
    // for the leading-bin fill: with lo == 0.0 (every other test target in
    // this file) a wrong fill that pins empty bins to zero is invisible,
    // because zero happens to be the right answer too.
    for (float sigma : {1.0f, 1.5f, 2.0f, 3.0f}) {
        const Plane e = gaussianBlur(slantedEdge(128, 128, 5.0f, 0.2f, 0.9f), sigma);
        const float got  = mtf50(crop(e, 24, 24, 80, 80));
        const float want = 0.1874f / sigma;
        CAPTURE(sigma); CAPTURE(got); CAPTURE(want);
        CHECK(got == doctest::Approx(want).epsilon(0.08));
    }
}

TEST_CASE("mtf50 recovers the analytic value with a narrow ROI that leaves ESF bins empty") {
    // The 80x80 nonzero-dark-level test above never exercises the leading-
    // bin backfill: the ESF spans +/-32px at 4x oversampling over 256 bins,
    // and pixels projected from an 80px-wide ROI span roughly +/-40px --
    // wider than the ESF range -- so every bin gets a real sample there and
    // the fix and the old last=0.0 bug are bit-identical. A 40px ROI
    // projects only about +/-20px, leaving roughly 48 leading (and 48
    // trailing) ESF bins empty, so this is the test that actually walks the
    // backfill path. Measured directly against a standalone build of both
    // versions (crop=40x40, lo=0.4): the fix stays within ~4.4% of the
    // analytic value at every sigma below, matching the base 80x80 test's
    // accuracy almost exactly -- the narrower ROI is not meaningfully
    // noisier once the dark level is filled correctly. The reverted
    // last=0.0 code was off by 17-34% on the same inputs. 10% tolerance
    // (vs the base test's 8%) leaves headroom for the fix while sitting
    // nowhere near the broken result, so this test passes with the fix and
    // fails without it -- see task-17-report.md for the recorded numbers.
    for (float sigma : {1.0f, 1.2f, 1.5f}) {
        const Plane e = gaussianBlur(slantedEdge(128, 128, 5.0f, 0.4f, 0.9f), sigma);
        const float got  = mtf50(crop(e, 44, 44, 40, 40));
        const float want = 0.1874f / sigma;
        CAPTURE(sigma); CAPTURE(got); CAPTURE(want);
        // .scale(0) makes epsilon a true relative tolerance: doctest::Approx
        // defaults scale to 1.0, and threshold = epsilon*(scale + max(|a|,|b|));
        // with values in the 0.06-0.2 range that default floor is 50-140% of
        // the target, which would pass almost anything. Verified below.
        CHECK(got == doctest::Approx(want).epsilon(0.10).scale(0));
    }
}

TEST_CASE("mtf50 falls monotonically as blur grows") {
    float prev = 1.0f;
    for (float sigma : {0.8f, 1.2f, 1.8f, 2.6f}) {
        const float m = mtf50(crop(gaussianBlur(slantedEdge(128, 128, 5.0f, 0, 1), sigma), 24, 24, 80, 80));
        CHECK(m < prev);
        prev = m;
    }
}

TEST_CASE("edge position finds the true sub-pixel intercept") {
    const Plane e = slantedEdge(64, 64, 5.0f, 0.0f, 1.0f);
    CHECK(edgePosition(e) == doctest::Approx(32.0f).epsilon(0.02));
}

TEST_CASE("radial mean of a flat field is flat") {
    const std::vector<float> r = radialMean(flatField(64, 64, 0.6f), 8);
    CHECK(r.size() == 8u);
    for (float v : r) CHECK(v == doctest::Approx(0.6f).epsilon(1e-4));
}

TEST_CASE("radial mean detects a synthetic vignette") {
    Plane p(128, 128);
    const Frame f = frameOf(128, 128);
    for (int y = 0; y < 128; ++y)
        for (int x = 0; x < 128; ++x) {
            const float dx = (x - f.cx) / f.halfDiag, dy = (y - f.cy) / f.halfDiag;
            p.at(x, y) = 1.0f - 0.5f * std::sqrt(dx * dx + dy * dy);
        }
    const std::vector<float> r = radialMean(p, 10);
    CHECK(r.front() > 0.95f);
    CHECK(r.back() < 0.60f);
    for (size_t i = 1; i < r.size(); ++i) CHECK(r[i] < r[i - 1]);
}

TEST_CASE("fringe width is zero with no chromatic offset and nonzero with one") {
    const Plane e = slantedEdge(96, 96, 5.0f, 0.0f, 1.0f);
    CHECK(fringeWidthPx(toImage(e)) == doctest::Approx(0.0f).epsilon(0.05));

    Image shifted = toImage(e);                              // move blue by one pixel
    for (int y = 0; y < 96; ++y)
        for (int x = 95; x > 0; --x) shifted.at(x, y, 2) = shifted.at(x - 1, y, 2);
    CHECK(fringeWidthPx(shifted) == doctest::Approx(-1.0f).epsilon(0.10));
}

TEST_CASE("rotational asymmetry is zero for a symmetric image and large otherwise") {
    Plane sym(65, 65);
    const Frame f = frameOf(65, 65);
    for (int y = 0; y < 65; ++y)
        for (int x = 0; x < 65; ++x) {
            const float dx = x - f.cx, dy = y - f.cy;
            sym.at(x, y) = std::exp(-0.01f * (dx * dx + dy * dy));
        }
    CHECK(rot90Asymmetry(sym) < 1e-4f);

    Plane skew = sym;
    for (int y = 0; y < 65; ++y) skew.at(10, y) += 1.0f;
    CHECK(rot90Asymmetry(skew) > 0.1f);
}

TEST_CASE("rot90Asymmetry throws on a non-square image") {
    const Plane rect(64, 32);
    CHECK_THROWS_AS(rot90Asymmetry(rect), std::invalid_argument);
}
