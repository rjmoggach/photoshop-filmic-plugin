#include <doctest/doctest.h>
#include "lenscore/pipeline.hpp"
#include "metrics.hpp"
#include "targets.hpp"
#include <cmath>

using namespace lens;
using namespace lens::metrics;
using namespace lens::targets;

static const color::SpecTable& table() {
    static color::SpecTable t = color::buildTable(12);
    return t;
}

static Params nullParams() {
    Params p;
    p.bands = 3;
    p.highlightRecovery = false;
    p.doLateralCa = false; p.doPsf = false; p.doVignette = false;
    return p;
}

TEST_CASE("bands are inside the visible range and their weights sum to one") {
    Params p; p.bands = 11;
    const auto b = sampleBands(p);
    CHECK(b.size() == 11u);
    double w = 0.0;
    for (auto& [l, wt] : b) { CHECK(l >= 380.0f); CHECK(l <= 780.0f); w += wt; }
    CHECK(w == doctest::Approx(1.0).epsilon(1e-4));
}

TEST_CASE("the preview tier uses three bands") {
    Params p; p.bands = 3;
    CHECK(sampleBands(p).size() == 3u);
}

TEST_CASE("every stage disabled is the identity") {
    const Image src = toImage(siemensStar(64, 64, 12));
    const Image out = render(src, nullParams(), table());
    for (size_t i = 0; i < src.px.size(); ++i)
        CHECK(out.px[i] == doctest::Approx(src.px[i]).epsilon(0.03));
}

TEST_CASE("the spectral round trip keeps grey grey") {
    Params p = nullParams(); p.bands = 11;
    for (float v : {0.2f, 0.5f, 0.9f}) {
        const Image out = render(toImage(flatField(32, 32, v)), p, table());
        CAPTURE(v);
        CHECK(out.at(16, 16, 0) == doctest::Approx(v).epsilon(0.05));
        CHECK(out.at(16, 16, 1) == doctest::Approx(v).epsilon(0.05));
        CHECK(out.at(16, 16, 2) == doctest::Approx(v).epsilon(0.05));
    }
}

TEST_CASE("chromatic defocus in waves falls quadratically with the f-number") {
    const float a = chromaticDefocusWaves(1e-3f, 32.0f, 550.0f, 2.0f);
    const float b = chromaticDefocusWaves(1e-3f, 32.0f, 550.0f, 4.0f);
    CHECK(a / b == doctest::Approx(4.0f).epsilon(1e-4));
}

TEST_CASE("wavefront field dependence follows the documented powers") {
    Params p;
    p.petzval = 1.0f; p.astig = 1.0f; p.coma = 1.0f; p.spherical = 1.0f;
    p.dispersion.residual = 0.0f;                    // isolate the geometric terms
    const optics::Wavefront a = wavefrontAt(p, 0.5f, 550.0f);
    const optics::Wavefront b = wavefrontAt(p, 1.0f, 550.0f);
    CHECK(b.defocus   / a.defocus   == doctest::Approx(4.0f).epsilon(1e-3));   // t^2
    CHECK(b.astig     / a.astig     == doctest::Approx(4.0f).epsilon(1e-3));   // t^2
    CHECK(b.coma      / a.coma      == doctest::Approx(8.0f).epsilon(1e-3));   // t^3
    CHECK(b.spherical == doctest::Approx(a.spherical));                        // constant
}

TEST_CASE("vignetting reproduces the model's radial curve") {
    Params p = nullParams();
    p.doVignette = true;
    p.vignette.tStop = 2.0f;
    const Image out = render(toImage(flatField(128, 128, 1.0f)), p, table());
    const std::vector<float> got = radialMean(luminance(out), 8);
    for (int i = 0; i < 8; ++i) {
        const float t = (i + 0.5f) / 8.0f;
        CAPTURE(i);
        CHECK(got[i] == doctest::Approx(optics::vignette(p.vignette, t)).epsilon(0.02));
    }
}

TEST_CASE("a rotationally symmetric input stays symmetric through the whole pipeline") {
    Params p;
    p.bands = 3;
    p.highlightRecovery = false;
    p.doLateralCa = true; p.doPsf = true; p.doVignette = true;
    p.lateralK = 2e-5f; p.petzval = 0.8f; p.spherical = 0.3f;
    p.psfGrid = 64; p.psfRings = 6; p.psfKernel = 17; p.effPatch = 32;

    Plane sym(65, 65);
    const Frame f = frameOf(65, 65);
    for (int y = 0; y < 65; ++y)
        for (int x = 0; x < 65; ++x) {
            const float dx = x - f.cx, dy = y - f.cy;
            sym.at(x, y) = std::exp(-0.004f * (dx * dx + dy * dy));
        }
    CHECK(rot90Asymmetry(luminance(render(toImage(sym), p, table()))) < 0.02f);
}

TEST_CASE("energy is conserved across the spectral and convolution stages") {
    Params p;
    p.bands = 3;
    p.highlightRecovery = false;
    p.doLateralCa = false; p.doVignette = false; p.doPsf = true;
    p.petzval = 0.0f; p.astig = 0.0f; p.coma = 0.0f; p.spherical = 0.0f;
    p.dispersion.residual = 0.0f;
    p.pupil.rEntrance = 1e6f; p.pupil.rExit = 1e6f; p.pupil.sepNorm = 0.0f;
    p.psfGrid = 64; p.psfRings = 4; p.psfKernel = 17; p.effPatch = 32;

    Plane blob(64, 64);
    for (int y = 20; y < 44; ++y) for (int x = 20; x < 44; ++x) blob.at(x, y) = 1.0f;
    const Image src = toImage(blob);
    const Image out = render(src, p, table());
    CHECK(totalEnergy(out) / totalEnergy(src) == doctest::Approx(1.0).epsilon(0.02));
}
