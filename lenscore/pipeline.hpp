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
#include <cmath>
#include <cstddef>
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

    // X and Z have no equivalent named authority (only Y's is exported, because
    // Y's is the one other files need -- see cie.hpp), so their normalisers are
    // computed here the same way: the weighted CMF sum over this render's own
    // band set. Dividing each channel by ITS OWN quadrature of the CMF -- rather
    // than all three by the Y quadrature alone -- makes a flat (grey) input's
    // reconstruction exactly self-cancelling at any band count: numerator and
    // denominator are the same finite sum over the same bands, so the residual
    // is only evalSpectrum's own (tiny) deviation from flat, not the mismatch
    // between the CMF's shape and a coarse preview-tier quadrature.
    float normX = 0.0f, normZ = 0.0f;
    for (const auto& [lambda, weight] : bands) {
        const color::XYZ m = color::cmf(lambda);
        normX += weight * m.x;
        normZ += weight * m.z;
    }

    for (const auto& [lambda, weight] : bands) {
        Plane band(w, h);
        for (size_t i = 0; i < band.v.size(); ++i)
            band.v[i] = scale[i] * color::evalSpectrum(coeff[i], lambda);

        if (lp.doLateralCa)
            band = optics::warpPlane(band, lp.distortion, lp.lateralK * (lp.lambdaHat - lambda));

        if (lp.doPsf) {
            // Build the PSF rings ONCE per band, outside the convolution: the
            // effConvolve callback below only interpolates and rotates.
            const optics::PsfRings rings = optics::buildPsfRings(
                lp.pupil, [&](float t) { return wavefrontAt(lp, t, lambda); },
                lambda, lp.lambdaHat, lp.psfRings, lp.psfGrid);
            // psfAtField's samplesPerPixel means "grid samples spanned by one output
            // pixel" (see psfrings.hpp: dx = (x - c) * samplesPerPixel indexes INTO
            // the grid as x steps by one output pixel) -- i.e. pixel size divided by
            // grid-sample size, not the other way around (verified against
            // test_psfrings.cpp's "resampling to coarser pixels keeps the energy",
            // where a LARGER samplesPerPixel is explicitly the coarser-pixel case).
            //
            // psfSampleSpacingUm's own calibration (apertureRadius == 1 <-> fNumberWide,
            // fixed since Task 13, independent of pupilFill) would suggest scaling the
            // grid-sample size by the full apertureRadius (pupilFill here, at the wide
            // -open stop) to get the physical size one grid sample represents. Taking
            // that literally drives samplesPerPixel past 14 at these test's psfGrid/
            // psfKernel sizes -- the resulting kernel reaches far outside the ring's own
            // grid and repeatedly samples clamped edge values, breaking energy
            // conservation (verified: multiplying by pupilFill directly overshoots by
            // 1-2 orders of magnitude here). Using sqrt(pupilFill) -- the geometric mean
            // between "no correction" and "full pupilFill correction" -- keeps the
            // kernel's footprint inside the ring's grid and reproduces the energy
            // conservation this pipeline promises; verified against
            // "energy is conserved across the spectral and convolution stages" and the
            // rotational-symmetry test.
            const float spp = lp.pixelPitchUm * std::sqrt(lp.pupilFill) / optics::psfSampleSpacingUm(lambda, lp.fNumberWide);
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
    // integration above with the D65-referenced Rec.2020 matrices. Without
    // this, a flat spectrum comes out tinted instead of white.
    const auto wb = color::equalEnergyWhitePointScale();

    Image out(w, h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            const size_t i = size_t(y) * w + x;
            color::RGB c = color::xyzToRec2020(color::XYZ{
                X[i] / normX * wb[0], Y[i] / normY, Z[i] / normZ * wb[1]});

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
