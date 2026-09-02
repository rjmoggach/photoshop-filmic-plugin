// Colour fidelity through the full pipeline.
//
// This is the regression test for the defect that shipped: reconstructing RGB
// from a few-band quadrature is not colour-accurate, because a handful of
// samples cannot integrate the sharply peaked CMFs. The white-point scale made
// a FLAT spectrum exact, so every neutral came out right and the whole suite
// stayed green -- while an orange sky rendered green in Photoshop.
//
// Neutrals are therefore not enough. These cases are saturated on purpose, and
// they run at every band count the panel can ask for, because the old code was
// wrong at all of them and merely less wrong as bands grew.
#include <doctest/doctest.h>
#include <cmath>
#include "lenscore/color/spectable.hpp"
#include "lenscore/pipeline.hpp"

namespace {

lens::Image flat(int w, int h, float r, float g, float b) {
    lens::Image img(w, h);
    for (int i = 0; i < w * h; ++i) {
        img.px[i * 3 + 0] = r; img.px[i * 3 + 1] = g; img.px[i * 3 + 2] = b;
    }
    return img;
}

const lens::color::SpecTable& table() {
    // res=8 keeps the test fast; the defect was in the reconstruction, not the
    // table's resolution, and it reproduces identically at any res.
    static const lens::color::SpecTable t = lens::color::buildTable(8);
    return t;
}

}  // namespace

TEST_CASE("a flat field keeps its colour through the pipeline") {
    // Measured worst-channel error against band count, with the table fitted to
    // the band set: 0.43 at 3 bands, 0.21 at 5, then 0.10 from 7 upward. The
    // floor is not the quadrature -- it is the sigmoid model's own gamut limit,
    // which the dense integral shares: the table alone maps a saturated red
    // (0.80, 0.10, 0.10) to (0.79, 0.13, 0.00). That limit is in the design spec.
    //
    // Three bands is therefore not a quality setting, it is a wrong answer, and
    // the plug-in no longer offers it. These are the counts it does offer.
    struct Patch { const char* name; float r, g, b; float tol; };
    const Patch patches[] = {
        // In-gamut: these must be close, and they are the ones that were wrong.
        {"sunset orange", 0.98f, 0.55f, 0.22f, 0.02f},
        {"dusk sky",      0.35f, 0.38f, 0.60f, 0.02f},
        {"deep shadow",   0.04f, 0.03f, 0.06f, 0.03f},
        {"mid grey",      0.50f, 0.50f, 0.50f, 0.01f},
        {"saturated blue",0.10f, 0.10f, 0.80f, 0.03f},
        // Outside what a bounded reflectance spectrum can reach. The tolerances
        // here record the model's measured gamut limit; they are not slack.
        {"saturated red",  0.80f, 0.10f, 0.10f, 0.11f},
        {"saturated green",0.10f, 0.80f, 0.10f, 0.07f},
    };

    const int w = 12, h = 12;
    for (int bands : {7, 11, 15}) {
        CAPTURE(bands);
        for (const auto& p : patches) {
            CAPTURE(p.name);
            lens::Params par{};
            par.bands = bands;
            par.lateralK = 1.0e-6f;      // the smallest the panel can ask for
            par.doPsf = false;           // isolate the colour path
            par.doVignette = false;
            par.highlightRecovery = false;   // the knee is a look, not a colour error

            const lens::Image out = lens::render(flat(w, h, p.r, p.g, p.b), par, table());
            const int c = (h / 2) * w + (w / 2);

            // A uniform field has nothing for the optics to displace, so the
            // output must be the input. Before the fix, "sunset orange" came back
            // as (0.97, 0.98, -0.05) -- green, with negative blue.
            CHECK(std::abs(out.px[c * 3 + 0] - p.r) < p.tol);
            CHECK(std::abs(out.px[c * 3 + 1] - p.g) < p.tol);
            CHECK(std::abs(out.px[c * 3 + 2] - p.b) < p.tol);
        }
    }
}

TEST_CASE("no channel goes negative on a saturated flat field") {
    lens::Params par{};
    par.bands = 7;
    par.lateralK = 1.0e-6f;
    par.doPsf = false;
    par.doVignette = false;
    par.highlightRecovery = false;

    const lens::Image out = lens::render(flat(12, 12, 0.98f, 0.55f, 0.22f), par, table());
    for (float v : out.px) CHECK(v >= -1e-4f);
}

TEST_CASE("geometric distortion moves every channel identically") {
    // Distortion is achromatic: it bends the picture, it does not colour it.
    //
    // Regression test for the third defect that shipped. Distortion was applied
    // per spectral band, inside the loop whose result is reconstructed as a
    // DIFFERENCE against an unaberrated reference. That difference cancels the
    // quadrature's colour bias only while the two images are close; a warp of
    // tens of pixels makes them entirely different, so the bias stopped
    // cancelling and became the output. A barrel setting produced violent green
    // and red fringing along every edge and barely any visible bending.
    //
    // It is applied to RGB ahead of the band loop now. Lateral chromatic
    // aberration stays inside it, being the part that really is per-wavelength.
    const int w = 96, h = 96;
    lens::Image grid(w, h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            const float v = (x % 16 == 0 || y % 16 == 0) ? 0.9f : 0.25f;
            const size_t o = (size_t(y) * w + x) * 3;
            grid.px[o] = grid.px[o + 1] = grid.px[o + 2] = v;
        }

    for (float amount : {-0.15f, 0.15f}) {
        CAPTURE(amount);
        lens::Params p{};
        p.bands = 3;
        p.highlightRecovery = false;
        p.doPsf = false;
        p.doVignette = false;
        p.doLateralCa = false;         // no lateral colour: nothing may split
        p.distortion.k1 = amount;

        const lens::Image out = lens::render(grid, p, table());

        float split = 0.0f, changed = 0.0f;
        for (size_t i = 0; i < out.px.size(); i += 3) {
            split = std::max(split, std::abs(out.px[i]     - out.px[i + 1]));
            split = std::max(split, std::abs(out.px[i + 1] - out.px[i + 2]));
            changed = std::max(changed, std::abs(out.px[i] - grid.px[i]));
        }
        // A grey input must stay grey. Before the fix this reached 0.5 and more.
        CHECK(split < 0.002f);
        // And it must actually bend something, or the test above passes on a
        // filter that quietly does nothing at all.
        CHECK(changed > 0.1f);
    }
}
