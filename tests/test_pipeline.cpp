#include <doctest/doctest.h>
#include "lenscore/pipeline.hpp"
#include "metrics.hpp"
#include "targets.hpp"
#include <algorithm>
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
    // Not doctest::Approx(...).scale(0): siemensStar's pattern is exactly 0.0
    // or 1.0, and doctest's Approx compares with a strict "<", so a target of
    // exactly 0 with scale(0) makes the tolerance exactly 0 too -- an EXACT
    // 0-vs-0 match would then fail the "<" check despite being a perfect
    // match. This is the same class of case the project's own sweep already
    // carved out for near-zero per-pixel data (test_psfrings.cpp's coma-tail
    // check): a plain absolute bound states the real intent correctly.
    const Image src = toImage(siemensStar(64, 64, 12));
    const Image out = render(src, nullParams(), table());
    for (size_t i = 0; i < src.px.size(); ++i)
        CHECK(std::abs(out.px[i] - src.px[i]) < 0.03f);
}

TEST_CASE("the spectral round trip keeps grey grey") {
    Params p = nullParams(); p.bands = 11;
    for (float v : {0.2f, 0.5f, 0.9f}) {
        const Image out = render(toImage(flatField(32, 32, v)), p, table());
        CAPTURE(v);
        CHECK(out.at(16, 16, 0) == doctest::Approx(v).epsilon(0.05).scale(0));
        CHECK(out.at(16, 16, 1) == doctest::Approx(v).epsilon(0.05).scale(0));
        CHECK(out.at(16, 16, 2) == doctest::Approx(v).epsilon(0.05).scale(0));
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

// Critical 1b: stage 4 (the PSF's pupil clip) owns mechanical vignetting now
// that Critical 1 makes it actually fire; stage 6 applies natural falloff
// ONLY, or the two would double-count. So the rendered radial curve is no
// longer naturalFalloff * mechanicalFraction (optics::vignette) -- it is
// naturalFalloff * (the PSF's OWN energy ratio at that field radius),
// because the mechanical component now comes from the convolution kernel's
// energy, not from a second, independent closed-form multiply. Measuring
// that ratio via pupilEnergyFraction on the SAME pupil geometry render()
// derives internally (pupilFill-scaled, per Critical 1) is the most direct
// way to predict it without depending on mechanicalFraction's simplified,
// entrance-pupil-only, circular-aperture model, which is not what the PSF
// stage actually rasterises (this pupil has 9 blades and both an entrance
// AND an exit clip circle) and is not expected to match it exactly -- the
// two remain independently cross-checked instead in test_pupil.cpp.
TEST_CASE("vignetting reproduces the model's radial curve") {
    Params p = nullParams();
    p.bands = 11;   // isolate the vignetting-curve comparison from 3-band spectral
                     // quadrature noise (measured several percent on its own at the
                     // preview tier; not what this test is about).
    p.doVignette = true;
    p.doPsf = true;
    p.vignette.tStop = 2.0f;
    const Image out = render(toImage(flatField(128, 128, 1.0f)), p, table());
    const std::vector<float> got = radialMean(luminance(out), 8);

    // Mirror render()'s own Critical-1 pupil-clip scaling so the prediction
    // uses exactly the geometry the PSF convolution actually rasterises.
    optics::PupilParams predPupil = p.pupil;
    const float stopFactor = optics::apertureRadius(p.vignette);
    predPupil.apertureRadius = p.pupilFill * stopFactor;
    predPupil.rEntrance *= p.pupilFill;
    predPupil.rExit     *= p.pupilFill;
    predPupil.sepNorm   *= p.pupilFill;
    const float axisEnergyFraction = optics::pupilEnergyFraction(predPupil, 0.0f, p.psfGrid);

    for (int i = 0; i < 8; ++i) {
        const float t = (i + 0.5f) / 8.0f;
        const float psfRatio = optics::pupilEnergyFraction(predPupil, t, p.psfGrid) / axisEnergyFraction;
        const float expected = optics::naturalFalloff(p.vignette, t) * psfRatio;
        CAPTURE(i);
        // Measured, not assumed: the point-wise prediction above under-states
        // the rendered value by up to ~17% at the outer bins, always in the
        // same direction (rendered > predicted). That is EFF itself, not
        // slack in the prediction -- effPatch's 50%-overlap Hann windows
        // (conv/eff.hpp) blend each patch's kernel into its neighbours so
        // the partition of unity holds, and a patch nearer the centre
        // carries a LESS vignetted kernel than one nearer the corner; that
        // blend dilutes the sharp point-wise falloff into a shallower one,
        // more so where the falloff's gradient is steepest (the corner).
        // 0.2 covers the measured worst case with headroom.
        CHECK(got[i] == doctest::Approx(expected).epsilon(0.2).scale(0));
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
    CHECK(totalEnergy(out) / totalEnergy(src) == doctest::Approx(1.0).epsilon(0.02).scale(0));
}

// Important 5: the case above disables everything that COULD break energy
// conservation -- zero aberration, zero dispersion residual, rEntrance/rExit
// pushed out to 1e6 so the cat's-eye clip never fires, sepNorm zeroed so
// there is no field-dependent offset either. That proves the PSF stage
// conserves energy in the trivial, unclipped case, but says nothing about
// the case this project actually renders: real aberration, with clipping
// enabled. With clipping and natural falloff both enabled, energy is NOT
// supposed to be conserved -- vignetting is deliberate, modelled light loss
// -- so this asserts the measured loss matches the PREDICTED vignetting
// (naturalFalloff * the PSF's own pupil-clip energy ratio, the same
// construction as "vignetting reproduces the model's radial curve" above)
// rather than asserting conservation.
TEST_CASE("energy loss with realistic aberration and clipping matches predicted vignetting") {
    Params p;
    p.bands = 3;
    p.highlightRecovery = false;
    p.doLateralCa = false; p.doVignette = true; p.doPsf = true;
    p.petzval = 0.5f; p.astig = 0.2f; p.coma = 0.1f; p.spherical = 0.05f;
    p.psfGrid = 256; p.psfRings = 8; p.psfKernel = 7; p.effPatch = 64;
    // Sync vignette to the pupil's own clip geometry -- exactly what loading
    // a real .lens file does (lensdata/lensfile.hpp) -- so the prediction
    // below and the pupil clip render() actually rasterises describe the
    // same barrel, the way a real preset's parameters would.
    p.vignette.rEntrance = p.pupil.rEntrance;
    p.vignette.sepNorm   = p.pupil.sepNorm;

    Plane blob(64, 64);
    for (int y = 20; y < 44; ++y) for (int x = 20; x < 44; ++x) blob.at(x, y) = 1.0f;
    const Image src = toImage(blob);
    const Image out = render(src, p, table());
    const double measured = totalEnergy(out) / totalEnergy(src);

    // Predict the loss by averaging naturalFalloff(t) * psfEnergyRatio(t)
    // over the blob's own footprint, mirroring render()'s own Critical-1
    // pupil-clip scaling. A modest N (128, versus the render's own 256) and
    // a stride over the footprint keep this cheap; pupilEnergyFraction is a
    // smooth function of t so a coarser probe grid does not change the
    // predicted average meaningfully.
    optics::PupilParams predPupil = p.pupil;
    const float stopFactor = optics::apertureRadius(p.vignette);
    predPupil.apertureRadius = p.pupilFill * stopFactor;
    predPupil.rEntrance *= p.pupilFill;
    predPupil.rExit     *= p.pupilFill;
    predPupil.sepNorm   *= p.pupilFill;
    const int probeGrid = 128;
    const float axisEnergyFraction = optics::pupilEnergyFraction(predPupil, 0.0f, probeGrid);

    const Frame frame = frameOf(64, 64);
    double predictedSum = 0.0; int n = 0;
    for (int y = 20; y < 44; y += 2)
        for (int x = 20; x < 44; x += 2) {
            const float dx = (float(x) - frame.cx) / frame.halfDiag;
            const float dy = (float(y) - frame.cy) / frame.halfDiag;
            const float t = std::clamp(std::sqrt(dx * dx + dy * dy), 0.0f, 1.0f);
            const float psfRatio = optics::pupilEnergyFraction(predPupil, t, probeGrid) / axisEnergyFraction;
            predictedSum += optics::naturalFalloff(p.vignette, t) * psfRatio;
            ++n;
        }
    REQUIRE(n > 0);
    const double predicted = predictedSum / n;
    CHECK(measured == doctest::Approx(predicted).epsilon(0.1).scale(0));
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
