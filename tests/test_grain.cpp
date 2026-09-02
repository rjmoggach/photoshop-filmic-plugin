#include "lenscore/film/grain.hpp"
#include <doctest/doctest.h>
#include <cmath>

namespace {
lens::Image flat(int w, int h, float v) {
    lens::Image img(w, h);
    for (int i = 0; i < w * h * 3; ++i) img.px[i] = v;
    return img;
}
float stddev(const lens::Image& a, const lens::Image& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.px.size(); ++i) { const double d = a.px[i] - b.px[i]; s += d * d; }
    return float(std::sqrt(s / double(a.px.size())));
}
}  // namespace

TEST_CASE("grain is strongest in the mid tones and absent at both ends") {
    lens::film::GrainParams g;
    g.amount = 1.0f;
    g.size = 1.0f;

    const lens::Image black = flat(64, 64, 0.0f);
    const lens::Image mid   = flat(64, 64, 0.5f);
    const lens::Image white = flat(64, 64, 1.0f);

    lens::Image gb = black, gm = mid, gw = white;
    lens::film::applyGrain(gb, g);
    lens::film::applyGrain(gm, g);
    lens::film::applyGrain(gw, g);

    // Clear film has nothing to develop and solid black nothing left undeveloped;
    // the fluctuation belongs in between. Uniform noise would fail this.
    CHECK(stddev(gb, black) < 1e-6f);
    CHECK(stddev(gw, white) < 1e-6f);
    // Full travel now moves a mid grey by about 3% rather than 10%: the old
    // constant read as sand rather than as grain.
    CHECK(stddev(gm, mid) > 0.01f);
    CHECK(stddev(gm, mid) < 0.06f);
}

TEST_CASE("grain is deterministic, so a preview matches the render") {
    lens::film::GrainParams g;
    g.amount = 1.0f;
    lens::Image a = flat(32, 32, 0.5f), b = flat(32, 32, 0.5f);
    lens::film::applyGrain(a, g);
    lens::film::applyGrain(b, g);
    CHECK(stddev(a, b) == 0.0f);   // bit-for-bit, not merely close
}

TEST_CASE("grain amount zero leaves the image untouched") {
    lens::film::GrainParams g;
    g.amount = 0.0f;
    const lens::Image src = flat(32, 32, 0.5f);
    lens::Image out = src;
    lens::film::applyGrain(out, g);
    CHECK(stddev(out, src) == 0.0f);
}

TEST_CASE("bigger grain makes coarser clumps, not just more noise") {
    // Measure how much neighbouring pixels agree. Fine grain decorrelates from
    // one pixel to the next; coarse grain does not. Comparing only the overall
    // deviation would pass on noise of any size at all.
    auto neighbourDiff = [](float size) {
        lens::film::GrainParams g;
        g.amount = 1.0f;
        g.size = size;
        lens::Image img = flat(96, 96, 0.5f);
        lens::film::applyGrain(img, g);
        double s = 0.0;
        int n = 0;
        for (int y = 0; y < 96; ++y)
            for (int x = 0; x + 1 < 96; ++x, ++n) {
                const double d = img.at(x, y, 0) - img.at(x + 1, y, 0);
                s += d * d;
            }
        return std::sqrt(s / double(n));
    };
    CHECK(neighbourDiff(6.0f) < neighbourDiff(1.0f) * 0.5);
}

// ---------------------------------------------------------------------------

#include "lenscore/color/spectable.hpp"
#include "lenscore/pipeline.hpp"

TEST_CASE("an anamorphic squeeze widens the blur sideways, not evenly") {
    // The point of a squeeze is that the blur stops being round. Measuring only
    // that the picture got softer would pass on a plain spherical blur.
    static const lens::color::SpecTable tbl = lens::color::buildTable(8);

    auto spread = [&](float squeeze, float& sx, float& sy) {
        const int w = 65, h = 65;
        lens::Image dot(w, h);
        for (int i = 0; i < w * h * 3; ++i) dot.px[i] = 0.0f;
        for (int c = 0; c < 3; ++c) dot.at(w / 2, h / 2, c) = 1.0f;

        lens::Params p{};
        p.bands = 7;
        p.doPsf = true;
        p.doLateralCa = false;
        p.doVignette = false;
        p.highlightRecovery = false;
        p.petzval = 3.0f;                 // something to spread
        p.psfSqueeze = squeeze;
        p.dispersion.correction_nm = {486.1f, 587.6f, 656.3f};
        p.dispersion.residual = 0.0f;

        const lens::Image out = lens::render(dot, p, tbl);
        double s = 0, mx = 0, my = 0;
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) { const double v = out.at(x, y, 0); s += v; mx += v * x; my += v * y; }
        mx /= s; my /= s;
        double vx = 0, vy = 0;
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) {
                const double v = out.at(x, y, 0) / s;
                vx += v * (x - mx) * (x - mx);
                vy += v * (y - my) * (y - my);
            }
        sx = float(std::sqrt(vx));
        sy = float(std::sqrt(vy));
    };

    float rx = 0, ry = 0, ax = 0, ay = 0;
    spread(1.0f, rx, ry);
    spread(2.0f, ax, ay);

    CAPTURE(rx); CAPTURE(ry); CAPTURE(ax); CAPTURE(ay);
    CHECK(rx == doctest::Approx(ry).epsilon(0.15).scale(0));   // round without it
    CHECK(ax > rx * 1.2f);                                     // wider across
    CHECK(ax > ay * 1.2f);                                     // and oval, not just bigger
}
