#pragma once
#include "lenscore/color/cie.hpp"
#include "lenscore/color/spectable.hpp"
#include "lenscore/color/transfer.hpp"
#include "lenscore/color/upsample.hpp"
#include "lenscore/conv/eff.hpp"
#include "lenscore/geometry.hpp"
#include "lenscore/image.hpp"
#include "lenscore/optics/lateralca.hpp"
#include "lenscore/optics/psf.hpp"
#include "lenscore/optics/psfrings.hpp"
#include "lenscore/optics/vignette.hpp"
#include "lenscore/optics/zernike.hpp"
#include "lenscore/params.hpp"
#include "lenscore/plane.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace lens {

// Inverse-CDF importance sampling of the equalizer, so emphasising a band
// costs no extra noise. A flat equalizer reduces to uniform spacing.
inline std::vector<std::pair<float, float>> sampleBands(const Params& p) {
    const int N = std::max(1, p.bands);
    const int R = 512;
    const size_t Rz = size_t(R);
    std::vector<float> pdf(Rz), cdf(Rz + 1, 0.0f);
    for (int i = 0; i < R; ++i) {
        const float l = color::kLambdaMin + (color::kLambdaMax - color::kLambdaMin) * (i + 0.5f) / R;
        float w = 1.0f;
        if (!p.equalizer.empty()) {
            w = p.equalizer.back().second;
            for (size_t k = 1; k < p.equalizer.size(); ++k)
                if (l <= p.equalizer[k].first) {
                    const auto& a = p.equalizer[k - 1]; const auto& b = p.equalizer[k];
                    const float f = (l - a.first) / std::max(1e-6f, b.first - a.first);
                    w = a.second + f * (b.second - a.second);
                    break;
                }
        }
        pdf[i] = std::max(1e-6f, w);
    }
    for (int i = 0; i < R; ++i) cdf[i + 1] = cdf[i] + pdf[i];
    for (float& c : cdf) c /= cdf[R];

    std::vector<std::pair<float, float>> out;
    out.reserve(size_t(N));
    for (int k = 0; k < N; ++k) {
        const float u = (k + 0.5f) / float(N);
        int i = 0;
        while (i + 1 < R && cdf[i + 1] < u) ++i;
        const float l = color::kLambdaMin + (color::kLambdaMax - color::kLambdaMin) * (i + 0.5f) / R;
        out.emplace_back(l, 1.0f / float(N));   // importance sampling makes weights equal
    }
    return out;
}

// Relative focal-length error to waves of defocus: W = F e / (8 lambda N^2).
inline float chromaticDefocusWaves(float relFocalError, float focal_mm, float lambdaNm, float fNumber) {
    const float lambda_mm = lambdaNm * 1e-6f;
    return focal_mm * relFocalError / (8.0f * lambda_mm * fNumber * fNumber);
}

inline optics::Wavefront wavefrontAt(const Params& p, float t, float lambdaNm) {
    optics::Wavefront w;
    const float e = optics::focusError(p.dispersion, lambdaNm, p.lambdaHat);
    w.defocus   = p.petzval * t * t + chromaticDefocusWaves(e, p.focal_mm, lambdaNm, p.fNumberWide);
    w.astig     = p.astig * t * t;
    w.coma      = p.coma * t * t * t;
    w.spherical = p.spherical;
    return w;
}

// color::equalEnergyWhitePointScale() computes its correction from a DENSE 5nm
// integral, but render() integrates over only a handful of importance-sampled
// bands. Applying the dense-integral correction to a coarse-quadrature error is
// two different quadratures fighting each other -- that mismatch is exactly why
// a flat spectrum tinted at low band counts. This mirrors that function's own
// logic but replaces the dense loop with THIS render call's own band list, so
// the correction and the error it is correcting are computed the same way. A
// flat spectrum then maps to white exactly, at any band count, by construction
// (the numerator and denominator are the same finite sum).
inline std::array<float, 2> bandWhitePointScale(const float* lambdas, const float* weights, int n) {
    color::XYZ acc{0.0f, 0.0f, 0.0f};
    for (int i = 0; i < n; ++i) {
        const color::XYZ m = color::cmf(lambdas[i]);
        acc.x += weights[i] * m.x;
        acc.y += weights[i] * m.y;
        acc.z += weights[i] * m.z;
    }
    const float norm = color::cieYNormalisation(lambdas, weights, n);
    if (norm <= 0.0f || acc.x <= 0.0f || acc.z <= 0.0f) return {1.0f, 1.0f};
    acc.x /= norm;
    acc.z /= norm;
    const color::XYZ ref = color::rec2020ToXyz(color::RGB{1.0f, 1.0f, 1.0f});
    return {ref.x / acc.x, ref.z / acc.z};
}

inline Image render(const Image& src, const Params& p, const color::SpecTable& tbl) {
    const int w = src.w, h = src.h;
    const Frame frame = frameOf(w, h);

    // Local copy: the pupil aperture is derived from the working stop, so
    // stopping down (vignette.tStop) narrows the PSF too, not just the
    // corners. pupilFill keeps the pupil from filling the whole FFT grid
    // (which would put the first Airy zero at a critically undersampled
    // 1.22 samples) at the wide-open stop.
    Params lp = p;
    lp.pupil.apertureRadius = lp.pupilFill * (lp.vignette.tStopWide / lp.vignette.tStop);

    // Per-pixel spectral coefficients plus the scale that rides outside the model.
    std::vector<color::Coeffs> coeff(size_t(w) * h);
    std::vector<float> scale(size_t(w) * h, 0.0f);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            color::RGB c{src.at(x, y, 0), src.at(x, y, 1), src.at(x, y, 2)};
            if (lp.highlightRecovery) {
                c.r = color::expandHighlights(c.r, lp.knee);
                c.g = color::expandHighlights(c.g, lp.knee);
                c.b = color::expandHighlights(c.b, lp.knee);
            }
            const size_t i = size_t(y) * w + x;
            // The table reproduces the colour directly for in-gamut inputs, so only the
            // EXCESS above the gamut boundary rides outside the model. Scaling by the raw
            // maximum here would square the brightness of every in-gamut pixel.
            scale[i] = std::max(1.0f, std::max({c.r, c.g, c.b, 0.0f}));
            coeff[i] = color::lookup(tbl, c);
        }

    const auto bands = sampleBands(lp);
    std::vector<float> X(size_t(w) * h, 0.0f), Y(X.size(), 0.0f), Z(X.size(), 0.0f);

    // One authority for the CIE Y normalisation -- do not recompute it inline.
    std::vector<float> lambdas(bands.size()), weights(bands.size());
    for (size_t k = 0; k < bands.size(); ++k) { lambdas[k] = bands[k].first; weights[k] = bands[k].second; }
    const float normY = color::cieYNormalisation(lambdas.data(), weights.data(), int(bands.size()));

    for (const auto& [lambda, weight] : bands) {
        Plane band(w, h);
        for (size_t i = 0; i < band.v.size(); ++i)
            band.v[i] = scale[i] * color::evalSpectrum(coeff[i], lambda);

        if (lp.doLateralCa)
            band = optics::warpPlane(band, lp.distortion, lp.lateralK * (lp.lambdaHat - lambda));

        if (lp.doPsf) {
            // psfAtField's samplesPerPixel means "grid samples spanned by one output
            // pixel" (psfrings.hpp: dx = (x - c) * samplesPerPixel indexes INTO the
            // ring grid as x steps by one output pixel) -- i.e. pixel size divided by
            // grid-sample size (verified against test_psfrings.cpp's "resampling to
            // coarser pixels keeps the energy": a LARGER samplesPerPixel is explicitly
            // the coarser-pixel case).
            //
            // The grid-sample size is psfSampleSpacingUm(lambda, fNumberWide) DIVIDED
            // by pupilFill, not multiplied: psfSampleSpacingUm is calibrated for a
            // pupil that fills the whole psfGrid FFT (apertureRadius == 1). Shrinking
            // the pupil to pupilFill of the grid (see the apertureRadius line above)
            // packs the SAME physical aperture into a smaller fraction of the same
            // N-sample grid, i.e. represents it at higher sample density -- one grid
            // sample now covers LESS physical distance, not more. Verified two ways:
            // multiplying (matching an apertureRadius==1 reference literally) pushes
            // samplesPerPixel past 14 at these parameters and every kernel from 5 to
            // grid-filling overshoots total energy by 1-2 orders of magnitude,
            // independent of psfGrid -- i.e. it is not a grid-too-small/clamping
            // artifact, the ratio itself is wrong. Dividing lands samplesPerPixel
            // near 1-2 (the diffraction pattern really does span a few output pixels
            // at these test parameters) and, with a kernel sized generously enough
            // to cover it (see maxKernel below), reproduces energy conservation to
            // within a fraction of a percent, independent of psfGrid.
            const float gridSampleUm = optics::psfSampleSpacingUm(lambda, lp.fNumberWide) / lp.pupilFill;
            const float spp = lp.pixelPitchUm / gridSampleUm;

            // The kernel's physical footprint (psfKernel pixels) cannot exceed the
            // ring's own grid footprint (psfGrid samples), or psfAtField reads
            // clamped edge samples for its outer pixels and energy conservation
            // breaks silently. Fail loudly instead -- same policy as psfrings.hpp's
            // "rings must be >= 2" check -- rather than quietly truncating a
            // request the caller's own parameters cannot support.
            int maxKernel = int(std::floor(lp.psfGrid * gridSampleUm / lp.pixelPitchUm));
            if (maxKernel % 2 == 0) --maxKernel;
            if (lp.psfKernel > maxKernel) {
                throw std::invalid_argument(
                    "render: psfKernel (" + std::to_string(lp.psfKernel) +
                    ") exceeds what psfGrid=" + std::to_string(lp.psfGrid) +
                    " covers at " + std::to_string(lambda) + "nm (max " +
                    std::to_string(maxKernel) + "); shrink psfKernel or grow psfGrid/pupilFill");
            }

            // Build the PSF rings ONCE per band, outside the convolution: the
            // effConvolve callback below only interpolates and rotates.
            const optics::PsfRings rings = optics::buildPsfRings(
                lp.pupil, [&](float t) { return wavefrontAt(lp, t, lambda); },
                lambda, lp.lambdaHat, lp.psfRings, lp.psfGrid);
            band = conv::effConvolve(band, lp.effPatch, [&](float cx, float cy) {
                const float dx = (cx - frame.cx) / frame.halfDiag;
                const float dy = (cy - frame.cy) / frame.halfDiag;
                const float t  = std::clamp(std::sqrt(dx * dx + dy * dy), 0.0f, 1.0f);
                return optics::psfAtField(rings, t, std::atan2(dy, dx), spp, lp.psfKernel);
            });
        }

        const color::XYZ m = color::cmf(lambda);
        for (size_t i = 0; i < X.size(); ++i) {
            X[i] += weight * band.v[i] * m.x;
            Y[i] += weight * band.v[i] * m.y;
            Z[i] += weight * band.v[i] * m.z;
        }
    }

    // Diagonal white-point scale that reconciles the equal-energy spectral
    // integration above with the D65-referenced Rec.2020 matrices. Computed
    // from THIS render's own band set (see bandWhitePointScale above), not
    // the dense continuous integral color::equalEnergyWhitePointScale() uses
    // for spectrumToRec2020 -- the correction must be integrated the same way
    // as the error it corrects, or a flat spectrum comes out tinted instead
    // of white at low band counts.
    const auto wb = bandWhitePointScale(lambdas.data(), weights.data(), int(bands.size()));

    Image out(w, h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            const size_t i = size_t(y) * w + x;
            color::RGB c = color::xyzToRec2020(color::XYZ{
                X[i] / normY * wb[0], Y[i] / normY, Z[i] / normY * wb[1]});

            if (lp.doVignette) {
                const float dx = (float(x) - frame.cx) / frame.halfDiag;
                const float dy = (float(y) - frame.cy) / frame.halfDiag;
                const float v = optics::vignette(lp.vignette, std::sqrt(dx * dx + dy * dy));
                c.r *= v; c.g *= v; c.b *= v;
            }
            if (lp.highlightRecovery) {
                c.r = color::compressHighlights(c.r, lp.knee);
                c.g = color::compressHighlights(c.g, lp.knee);
                c.b = color::compressHighlights(c.b, lp.knee);
            }
            out.at(x, y, 0) = c.r; out.at(x, y, 1) = c.g; out.at(x, y, 2) = c.b;
        }
    return out;
}

}  // namespace lens
