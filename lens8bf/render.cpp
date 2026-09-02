#include "render.h"

#include "lenscore/color/spectable.hpp"
#include "lenscore/film/grain.hpp"
#include "lenscore/pipeline.hpp"

#include <cmath>
#include <mutex>
#include <optional>

namespace {

// The 9.4MB spectral table is read once and reused for the life of the process.
// The dialog re-renders a proxy on every slider move, so reloading it per call
// would dominate the preview's cost entirely.
std::mutex gTableMutex;
std::optional<lens::color::SpecTable> gTable;

const lens::color::SpecTable& tableFrom(const std::string& path) {
    std::lock_guard<std::mutex> lock(gTableMutex);
    if (!gTable) {
        auto t = lens::color::readTable(path);
        if (!t) throw std::runtime_error("could not read the spectral table at " + path);
        gTable = std::move(*t);
    }
    return *gTable;
}

}  // namespace

std::string lensRender(float* rgb, int w, int h, const LensControls& c,
                       double pixelScale, const std::string& tablePath) {
    try {
        if (w <= 0 || h <= 0 || rgb == nullptr) return {};
        if (!c.any()) return {};

        lens::Params p{};
        p.bands = c.bands;
        p.highlightRecovery = true;

        p.lateralK = float(c.chromaticAberration * 1.0e-6);   // +-1e-4 at the extremes

        // Brown-Conrady k1, +-0.25 at the extremes: roughly a fisheye-ish bulge
        // through to a strong telephoto pinch. Beyond that the corners fold over
        // themselves and stop looking like a lens.
        //
        // Sign matters and was inverted here. applyDistortion maps an OUTPUT
        // position to the SOURCE position it came from, so a positive k1 samples
        // further out than it draws: the picture is squeezed inward at the edges
        // and straight lines bow in -- pincushion. Negative k1 is the bulge. The
        // control reads "minus bulges outward (barrel)", so minus must be
        // negative k1, which is what this now does. Measured before the fix:
        // slider -60 pulled content in and blacked the corners, i.e. pincushion
        // under a label saying barrel.
        p.distortion.k1 = float(c.distortion * 0.0025);
        p.petzval  = float(c.edgeBlur    * 0.02);             // +-2 waves
        p.astig    = float(c.astigmatism * 0.02);
        p.coma     = float(c.coma        * 0.02);

        // lenscore's default Dispersion is an UNCORRECTED BK7 singlet. That is
        // physically honest for a singlet and completely wrong for a photographic
        // lens: measured, it puts 70 waves of longitudinal colour into the
        // wavefront at 408nm, identical at the centre of the frame and at the
        // corner. Two waves of astigmatism or coma against a 70-wave field-flat
        // defocus is invisible -- which is exactly how it looked, astigmatism
        // appearing to soften the whole picture uniformly and coma appearing to
        // do nothing at all. It also smeared the point spread function across the
        // whole FFT grid, turning every blur into a featureless box.
        //
        // Model a corrected lens instead, and let the chromatic aberration
        // control set how much colour error is left: an apochromat's residual is
        // about 2 waves at full travel, and none at zero.
        p.dispersion.correction_nm = {486.1f, 587.6f, 656.3f};   // F, d, C lines
        p.dispersion.residual = float(std::abs(c.chromaticAberration) / 100.0);

        p.vignette.tStop     = float(c.tStop);
        // The shape is the lens's; the slider says how far to push it, and which
        // way. Previously the slider only moved the natural-falloff exponent
        // while the mechanical term applied at FULL strength regardless -- so
        // half travel already darkened the corners to 0.50.
        p.vignette.naturalExp = 4.0f;
        p.vignetteAmount = float(c.vignette / 100.0);

        // Distortion is applied to RGB ahead of the band loop (see pipeline.hpp),
        // so this gate is about lateral colour alone.
        p.doLateralCa = p.lateralK != 0.0f;
        p.doPsf = (p.petzval != 0.0f) || (p.astig != 0.0f) || (p.coma != 0.0f) ||
                  (std::abs(c.squeeze - 1.0) > 1e-6);
        p.doVignette  = c.vignette != 0.0;

        // Blur radii are in pixels. A proxy pixel stands for `pixelScale` real
        // pixels, so telling the model the sensor's pixels are that much bigger
        // makes the preview's blur cover the same FRACTION of the frame as the
        // final render will. Without this the preview systematically under-blurs.
        if (pixelScale > 0.0) p.pixelPitchUm = float(p.pixelPitchUm * pixelScale);

        p.psfSqueeze = float(c.squeeze);

        lens::Image img(w, h);
        std::copy(rgb, rgb + size_t(w) * h * 3, img.px.begin());

        if (c.anyOptics()) img = lens::render(img, p, tableFrom(tablePath));

        // Grain goes on last, as it does on film: the emulsion records whatever
        // the lens delivered, so it is not blurred by the lens's own aberrations.
        if (c.grain > 0.0) {
            lens::film::GrainParams g;
            g.amount = float(c.grain / 100.0);
            g.size   = float(c.grainSize);
            g.colour = float(c.grainColour / 100.0);
            lens::film::applyGrain(img, g, pixelScale);
        }

        std::copy(img.px.begin(), img.px.end(), rgb);
        return {};
    } catch (const std::exception& e) {
        return e.what();
    } catch (...) {
        return "the renderer failed";
    }
}
