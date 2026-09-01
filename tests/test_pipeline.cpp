#include <doctest/doctest.h>
#include "lenscore/pipeline.hpp"
#include "metrics.hpp"
#include "targets.hpp"
#include <cmath>
#include <stdexcept>

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
    p.psfGrid = 256; p.psfRings = 6; p.psfKernel = 7; p.effPatch = 64;

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
    p.psfGrid = 256; p.psfRings = 4; p.psfKernel = 7; p.effPatch = 64;

    Plane blob(64, 64);
    for (int y = 20; y < 44; ++y) for (int x = 20; x < 44; ++x) blob.at(x, y) = 1.0f;
    const Image src = toImage(blob);
    const Image out = render(src, p, table());
    CHECK(totalEnergy(out) / totalEnergy(src) == doctest::Approx(1.0).epsilon(0.02));
}

// The self-normalisation bug this pipeline once had (each of X, Y, Z divided by its
// OWN weighted-CMF sum) forced a flat spectrum to reconstruct exactly, at the cost of
// applying three different quadrature corrections to a coloured spectrum -- a
// distortion no test caught, because every other case here is achromatic. This checks
// a clearly non-grey colour through the FULL pipeline (spectral round trip only; PSF,
// lateral CA and vignette are orthogonal to colour and are covered elsewhere) against
// the dense-integral ground truth spectrumToRec2020(lookup(...)) computes, at both the
// preview tier and the default band count. 3 bands is a genuinely coarse quadrature of
// a bimodal CMF -- some error there is real and expected, not a bug -- so its tolerance
// is deliberately looser than 11's; both must still land in the right neighbourhood of
// the true colour, not some other hue entirely.
TEST_CASE("a coloured input keeps its hue through the spectral stage at both band counts") {
    const color::RGB in{0.6f, 0.45f, 0.35f};
    const color::RGB truth = color::spectrumToRec2020(color::lookup(table(), in));

    Image src(8, 8);
    for (int y = 0; y < 8; ++y)
        for (int x = 0; x < 8; ++x) {
            src.at(x, y, 0) = in.r; src.at(x, y, 1) = in.g; src.at(x, y, 2) = in.b;
        }

    for (int n : {3, 11}) {
        Params p = nullParams(); p.bands = n;
        const Image out = render(src, p, table());
        const float eps = (n == 3) ? 0.30f : 0.08f;
        CAPTURE(n);
        CHECK(out.at(4, 4, 0) == doctest::Approx(truth.r).epsilon(eps).scale(0));
        CHECK(out.at(4, 4, 1) == doctest::Approx(truth.g).epsilon(eps).scale(0));
        CHECK(out.at(4, 4, 2) == doctest::Approx(truth.b).epsilon(eps).scale(0));
    }
}

TEST_CASE("render rejects a psfKernel the grid cannot cover") {
    // psfGrid=8 covers only a handful of micrometres at these defaults; psfKernel=7
    // (the Params default) demands far more than that grid can supply at any sampled
    // wavelength. Programmer error should fail loudly here, the same policy
    // psfrings.hpp already applies to "rings must be >= 2", rather than silently
    // truncating a request the caller's own parameters cannot support.
    Params p = nullParams();
    p.doPsf = true;
    p.psfGrid = 8;
    p.psfRings = 2;
    CHECK_THROWS_AS(render(toImage(flatField(8, 8, 0.5f)), p, table()), std::invalid_argument);
}
