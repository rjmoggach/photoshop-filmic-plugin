// A flat field must survive the blur stage flat.
//
// This is the regression test for the second defect that shipped: psfAtField's
// kernel energy depends on the ORIENTATION it is asked for, not only on the
// field angle. Measured at the frame corner on a perfect lens, the kernel sums
// 0.607 along the frame axes against 0.538 on the diagonals -- a 12% swing with
// the four-fold symmetry of the square resampling window. Orientation varies
// with position, so that painted a fixed cross of shading across the picture,
// which showed up as a wedge of shading through the centre of the frame.
//
// The radial falloff was NOT the problem and these tests must not claim it was:
// kernel sum over pupil throughput held flat at 0.95 across the field, so the
// mechanical vignetting was being reproduced correctly all along.
//
// Nothing in the suite caught this, because every existing vignetting test
// measures a RADIAL mean -- averaging around each ring is precisely the
// operation that destroys an angular artefact before it can be asserted on.
// These tests therefore hold radius fixed and vary the angle.
#include "lenscore/color/spectable.hpp"
#include "lenscore/pipeline.hpp"
#include <doctest/doctest.h>
#include <cmath>

namespace {

lens::Image grey(int w, int h, float v) {
    lens::Image img(w, h);
    for (int i = 0; i < w * h; ++i) {
        img.px[i * 3 + 0] = v; img.px[i * 3 + 1] = v; img.px[i * 3 + 2] = v;
    }
    return img;
}

const lens::color::SpecTable& table() {
    static const lens::color::SpecTable t = lens::color::buildTable(8);
    return t;
}

}  // namespace

TEST_CASE("the blur stage conserves energy on a flat field") {
    struct Case { const char* name; float petzval, astig, coma; };
    const Case cases[] = {
        {"field curvature", 2.0f, 0.0f, 0.0f},
        {"astigmatism",     0.0f, 2.0f, 0.0f},
        {"coma",            0.0f, 0.0f, 2.0f},
        {"all three",       1.5f, 1.5f, 1.5f},
    };

    const int w = 96, h = 54;      // deliberately not square: t contours must not
    for (const auto& c : cases) {  // line up with the patch grid and hide a seam
        CAPTURE(c.name);
        lens::Params p{};
        p.bands = 7;
        p.doPsf = true;
        p.doLateralCa = false;
        p.doVignette = false;      // vignetting is the ONLY thing allowed to darken
        p.highlightRecovery = false;
        p.petzval = c.petzval; p.astig = c.astig; p.coma = c.coma;

        const lens::Image out = lens::render(grey(w, h, 0.5f), p, table());

        float lo = 1e9f, hi = -1e9f;
        for (size_t i = 0; i < out.px.size(); i += 3) {
            lo = std::min(lo, out.px[i]);
            hi = std::max(hi, out.px[i]);
        }
        // Before the fix this ran 0.33 to 0.55 against an input of 0.50.
        //
        // The band left is EFF's own approximation, measured at 3.4% in the
        // interior at two waves of coma. Overlapping patches carry DIFFERENT
        // kernels, and while their Hann windows sum to exactly one, those windows
        // convolved with different kernels do not -- so a flat field picks up a
        // gentle ripple wherever the kernel changes quickly across a patch. It is
        // inherent to the method, it scales with how hard the aberration is
        // driven, and a smaller patch would reduce it at proportionate cost.
        // 4% covers the measured worst case; anything beyond that is a defect.
        CHECK(lo == doctest::Approx(0.5f).epsilon(0.04).scale(0));
        CHECK(hi == doctest::Approx(0.5f).epsilon(0.04).scale(0));
    }
}

TEST_CASE("the blur stage is isotropic: equal radius means equal shading") {
    // Vignetting ON, so the only thing allowed to vary the brightness is
    // DISTANCE from the centre. Two points at the same distance but different
    // angles must match. A radial mean cannot see this; that is the point.
    lens::Params p{};
    p.bands = 7;
    p.doPsf = true;
    p.doLateralCa = false;
    p.doVignette = true;
    p.highlightRecovery = false;
    p.coma = 2.0f;

    const int w = 129, h = 129;      // square and odd, so the centre is a pixel
    const lens::Image out = lens::render(grey(w, h, 0.5f), p, table());

    const float cx = 0.5f * float(w - 1), cy = 0.5f * float(h - 1);
    for (float radius : {20.0f, 40.0f, 60.0f}) {
        CAPTURE(radius);
        float lo = 1e9f, hi = -1e9f;
        for (int deg = 0; deg < 360; deg += 5) {
            const float a = float(deg) * 3.14159265f / 180.0f;
            const int x = int(std::lround(cx + radius * std::cos(a)));
            const int y = int(std::lround(cy + radius * std::sin(a)));
            const float v = out.px[(size_t(y) * w + x) * 3];
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
        // Before the fix the square window's four-fold ripple put this well
        // outside 2%; the residue left is EFF's patch blending.
        CAPTURE(lo); CAPTURE(hi);
        CHECK(hi - lo < 0.02f * 0.5f);
    }
}

TEST_CASE("a warp reaching past the border yields black, not a smeared stripe") {
    // Distortion deliberately samples outside the frame. Clamping returns the
    // same edge pixel for every such sample, smearing one pixel into a stripe
    // down the side of the picture. Black says plainly where the frame ended.
    const int w = 120, h = 80;
    lens::Image tex(w, h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            const float v = 0.5f + 0.35f * std::sin(float(x) * 0.7f) * std::cos(float(y) * 0.9f);
            for (int c = 0; c < 3; ++c) tex.at(x, y, c) = v;
        }

    lens::Params p{};
    p.bands = 7;
    p.highlightRecovery = false;
    p.doPsf = false;
    p.doVignette = false;
    p.doLateralCa = false;
    p.distortion.k1 = 0.25f;      // pincushion: pulls content in, exposing outside
    p.dispersion.correction_nm = {486.1f, 587.6f, 656.3f};
    p.dispersion.residual = 0.0f;

    const lens::Image out = lens::render(tex, p, table());

    // The corners are what a pincushion pulls furthest in, so they must be black.
    CHECK(out.at(0, 0, 0) == doctest::Approx(0.0f).epsilon(0.02).scale(1));
    CHECK(out.at(w - 1, 0, 0) == doctest::Approx(0.0f).epsilon(0.02).scale(1));
    CHECK(out.at(0, h - 1, 0) == doctest::Approx(0.0f).epsilon(0.02).scale(1));
    CHECK(out.at(w - 1, h - 1, 0) == doctest::Approx(0.0f).epsilon(0.02).scale(1));

    // And the middle must still be picture, or "all black" would pass this.
    float spread = 0.0f;
    for (int y = h / 3; y < 2 * h / 3; ++y)
        for (int x = w / 3; x < 2 * w / 3; ++x) spread = std::max(spread, out.at(x, y, 0));
    CHECK(spread > 0.4f);
}
