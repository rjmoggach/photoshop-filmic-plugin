#pragma once
#include "lenscore/color/cie.hpp"
#include <map>
#include <mutex>
#include <string>
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


// A spectral table fitted against one particular band set, cached for the life
// of the process. Building one costs a few hundred milliseconds and depends only
// on the wavelengths and weights, so a panel dragging a slider pays it once.
inline const color::SpecTable& quadratureTable(const std::vector<float>& lambdas,
                                               const std::vector<float>& weights,
                                               const color::SpecTable& seed) {
    static std::mutex mu;
    static std::map<std::string, color::SpecTable> cache;

    std::string key;
    key.reserve(lambdas.size() * 24);
    for (size_t k = 0; k < lambdas.size(); ++k) {
        key += std::to_string(lambdas[k]);
        key += ':';
        key += std::to_string(weights[k]);
        key += ';';
    }

    std::lock_guard<std::mutex> lock(mu);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    const int n = int(lambdas.size());
    const float normY = color::cieYNormalisation(lambdas.data(), weights.data(), n);
    const auto wb = bandWhitePointScale(lambdas.data(), weights.data(), n);
    const auto reconstruct = [&](const color::Coeffs& co) {
        color::XYZ acc{0, 0, 0};
        for (int k = 0; k < n; ++k) {
            const color::XYZ m = color::cmf(lambdas[k]);
            const float f = weights[k] * color::evalSpectrum(co, lambdas[k]);
            acc.x += f * m.x; acc.y += f * m.y; acc.z += f * m.z;
        }
        return color::xyzToRec2020(color::XYZ{acc.x / normY * wb[0], acc.y / normY,
                                              acc.z / normY * wb[1]});
    };

    // res 12 was measured to reach the model's own gamut floor; more resolution
    // buys nothing here and costs build time linearly in res^3.
    (void)seed;
    return cache.emplace(key, color::buildTableAgainst(12, reconstruct)).first->second;
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
    const float stopFactor = optics::apertureRadius(lp.vignette);   // tStopWide / tStop
    lp.pupil.apertureRadius = lp.pupilFill * stopFactor;

    // Critical 1: rEntrance, rExit and sepNorm are configured "in units of
    // the wide-open aperture radius" (pupil.hpp), i.e. in the SAME raster
    // units apertureRadius==1.0 would occupy before pupilFill shrinks the
    // rasterised pupil down to leave FFT headroom. Scaling apertureRadius by
    // pupilFill above without also rescaling the clip geometry into that
    // same shrunk coordinate system left the cat's-eye clip circles ~4x too
    // large to ever touch the (correctly shrunk) aperture disc -- the clip
    // never fired, so the PSF carried no optical vignetting at all, diverging
    // from mechanicalFraction by 32% at the corner (measured:
    // pupilEnergyFraction was 0.93401 at t=0, 0.5 AND 1.0 -- identical).
    // The entrance/exit pupil barrels and the field-dependent offset do NOT
    // themselves shrink when the iris (apertureRadius) stops down -- only
    // apertureRadius carries the tStopWide/tStop factor, matching
    // mechanicalFraction's own structure (vignette.hpp: `a` shrinks with
    // tStop, p.rEntrance does not). So the clip geometry gets pupilFill only.
    lp.pupil.rEntrance *= lp.pupilFill;
    lp.pupil.rExit     *= lp.pupilFill;
    lp.pupil.sepNorm   *= lp.pupilFill;

    // Geometric distortion is ACHROMATIC and it is not a small perturbation, so
    // it is applied here, to RGB, before the spectral machinery runs -- not per
    // band inside it.
    //
    // Two reasons. It moves every wavelength identically, so putting it in the
    // band loop asks three to eleven separate warps to agree pixel for pixel and
    // they do not. And the reconstruction at the end of this function applies the
    // optics as a DIFFERENCE against an unaberrated reference, which cancels the
    // quadrature's colour bias only while the two are close; a warp of tens of
    // pixels makes them completely different images, so instead of cancelling,
    // the bias is what you see -- measured, a barrel setting produced violent
    // green and red fringing along every edge and no visible bending at all.
    //
    // Lateral chromatic aberration STAYS in the band loop: it is the part that
    // genuinely differs per wavelength, and it is small.
    Image distorted;
    const bool doDistort = lp.distortion.k1 != 0.0f || lp.distortion.k2 != 0.0f ||
                           lp.distortion.k3 != 0.0f || lp.distortion.p1 != 0.0f ||
                           lp.distortion.p2 != 0.0f;
    if (doDistort) {
        distorted = Image(w, h);
        for (int c = 0; c < 3; ++c) {
            Plane ch(w, h);
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; ++x) ch.at(x, y) = src.at(x, y, c);
            // K = 0: the same map for every channel, so nothing separates.
            const Plane warped = optics::warpPlane(ch, lp.distortion, 0.0f);
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; ++x) distorted.at(x, y, c) = warped.at(x, y);
        }
    }
    const Image& base = doDistort ? distorted : src;

    // The band loop must not apply the distortion again.
    const optics::Distortion bandDistortion{};

    const auto bands = sampleBands(lp);

    std::vector<float> lambdas(bands.size()), weights(bands.size());
    for (size_t k = 0; k < bands.size(); ++k) { lambdas[k] = bands[k].first; weights[k] = bands[k].second; }

    // One authority for the CIE Y normalisation -- do not recompute it inline.
    const float normY = color::cieYNormalisation(lambdas.data(), weights.data(), int(bands.size()));
    const auto wb = bandWhitePointScale(lambdas.data(), weights.data(), int(bands.size()));

    // What this band set actually reconstructs a given spectrum as. Used below
    // to pre-compensate the input, and it is exactly the sum the band loop
    // performs, so the two cannot drift.
    const auto reconstruct = [&](const color::Coeffs& co) {
        color::XYZ acc{0, 0, 0};
        for (size_t k = 0; k < bands.size(); ++k) {
            const color::XYZ m = color::cmf(lambdas[k]);
            const float f = weights[k] * color::evalSpectrum(co, lambdas[k]);
            acc.x += f * m.x; acc.y += f * m.y; acc.z += f * m.z;
        }
        return color::xyzToRec2020(color::XYZ{acc.x / normY * wb[0], acc.y / normY,
                                              acc.z / normY * wb[1]});
    };

    // Coefficients fitted against THIS band set, not against the dense spectral
    // integral the shipped table was built from.
    //
    // A handful of samples cannot reproduce that integral, so a spectrum fitted
    // to the dense model comes back as the wrong colour through a sparse one:
    // measured at three bands, a sunset orange rendered green and a saturated
    // red nearly black, while every neutral stayed exact because the white-point
    // scale pins a flat spectrum and nothing else.
    //
    // An earlier fix rendered a second, optics-free reference and applied the
    // optics as the difference. That got the colour right and introduced a worse
    // fault: the correction term is SHARP while the aberrated term is blurred,
    // so adding them put an unblurred edge back into a blurred picture -- a
    // bright rim along every high-contrast edge, plainly visible in a photograph.
    //
    // Fitting the table to the quadrature has neither problem. It is a property
    // of the band set alone, so it is computed once and cached; and the
    // reconstruction stays a plain linear sum over bands, which commutes with
    // the convolution. Nothing is corrected after the optics, so nothing gets
    // re-sharpened.
    //
    // Measured worst-channel error against band count, after this change: 0.43
    // at 3 bands, 0.21 at 5, then 0.10 from 7 bands upward. That last figure is
    // the floor and it is NOT the quadrature -- it is the sigmoid model's own
    // gamut limit, which the dense integral shares (see the design spec).
    const color::SpecTable& quadTable = quadratureTable(lambdas, weights, tbl);

    // Per-pixel spectral coefficients plus the scale that rides outside the model.
    std::vector<color::Coeffs> coeff(size_t(w) * h);
    std::vector<float> scale(size_t(w) * h, 0.0f);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            color::RGB c{base.at(x, y, 0), base.at(x, y, 1), base.at(x, y, 2)};
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

            coeff[i] = color::lookup(quadTable, c);
        }

    // Pupil throughput against field angle -- mechanical vignetting. Sampled
    // once here and interpolated per patch: pupilEnergyFraction rasterises the
    // clipped pupil over the whole FFT grid, far too costly to call for every
    // patch of every band.
    // Gated on doVignette as well as doPsf: mechanical vignetting is real, but a
    // user whose vignette control reads zero expects a flat frame, and being
    // handed 40% corner falloff for switching on edge blur is not a defensible
    // surprise. Vignetting -- natural and mechanical both -- belongs to the one
    // control named for it.
    std::vector<float> mechTable;
    if (lp.doVignette) {
        const int n = std::max(2, lp.psfRings);
        mechTable.resize(size_t(n));
        const float axis = optics::pupilEnergyFraction(lp.pupil, 0.0f, lp.psfGrid);
        for (int i = 0; i < n; ++i) {
            const float t = float(i) / float(n - 1);
            const float e = optics::pupilEnergyFraction(lp.pupil, t, lp.psfGrid);
            mechTable[size_t(i)] = (axis > 0.0f) ? e / axis : 1.0f;
        }
    }
    const auto mechanicalAt = [&mechTable](float t) -> float {
        if (mechTable.size() < 2) return 1.0f;
        const float f = std::clamp(t, 0.0f, 1.0f) * float(mechTable.size() - 1);
        const int i0 = std::min(int(f), int(mechTable.size()) - 2);
        const float d = f - float(i0);
        return mechTable[size_t(i0)] * (1.0f - d) + mechTable[size_t(i0) + 1] * d;
    };

    std::vector<float> X(size_t(w) * h, 0.0f), Y(X.size(), 0.0f), Z(X.size(), 0.0f);

    for (size_t bandIndex = 0; bandIndex < bands.size(); ++bandIndex) {
        const float lambda = bands[bandIndex].first;
        const float weight = bands[bandIndex].second;
        Plane band(w, h);
        for (size_t i = 0; i < band.v.size(); ++i)
            band.v[i] = scale[i] * color::evalSpectrum(coeff[i], lambda);

        if (lp.doLateralCa)
            band = optics::warpPlane(band, bandDistortion, lp.lateralK * (lp.lambdaHat - lambda));

        if (lp.doPsf) {
            // psfAtField's samplesPerPixel means "grid samples spanned by one output
            // pixel" (psfrings.hpp: dx = (x - c) * samplesPerPixel indexes INTO the
            // ring grid as x steps by one output pixel) -- i.e. pixel size divided by
            // grid-sample size (verified against test_psfrings.cpp's "resampling to
            // coarser pixels keeps the energy": a LARGER samplesPerPixel is explicitly
            // the coarser-pixel case).
            //
            // The grid-sample size is psfSampleSpacingUm(lambda, fNumberWide) MULTIPLIED
            // by pupilFill: psfSampleSpacingUm is calibrated for a pupil that fills the
            // whole psfGrid FFT (apertureRadius == 1), and pupilFill (see the
            // apertureRadius line above) is exactly how much smaller than that our
            // actual working aperture is. This makes samplesPerPixel large (an
            // unaberrated PSF is genuinely sub-pixel at this project's default optics),
            // which used to break energy conservation for the *unaberrated* case only --
            // that was a resampling bug in psfAtField (point-sampling a sub-pixel PSF
            // aliases), fixed at the source in psfrings.hpp, not by changing this ratio.
            // Measured across bare-Airy through 20-wave defocus with that fix in place,
            // this (multiplying) is the formula that holds energy conservation in every
            // aberrated regime -- the regime this project actually renders -- not just
            // the unaberrated one; see psfrings.hpp's psfAtField for the measurements.
            const float gridSampleUm = optics::psfSampleSpacingUm(lambda, lp.fNumberWide) * lp.pupilFill;
            const float spp = lp.pixelPitchUm / gridSampleUm;

            // How wide a kernel this band's grid can actually cover. The
            // kernel's physical footprint (psfKernel pixels) cannot exceed the
            // ring's own grid footprint (psfGrid samples), or psfAtField's box
            // average draws from clamped edge samples for its outer pixels.
            //
            // This used to throw. That was right when the kernel carried its own
            // absolute scale, because truncating it silently lost light. The
            // convolution now normalises every kernel to unit sum before use, so
            // a narrower window costs a little of the PSF's outer tail and
            // nothing else -- and throwing here breaks legitimate requests: a
            // preview renders a proxy by telling the model its pixels are larger,
            // which shrinks this limit in exact proportion and made the whole
            // preview fail with "psfKernel (7) exceeds what psfGrid=256 covers at
            // 408nm". Shrinking the window is the correct response, since a PSF
            // that spans fewer pixels needs fewer pixels to hold it.
            int maxKernel = int(std::floor(lp.psfGrid * gridSampleUm / lp.pixelPitchUm));
            if (maxKernel % 2 == 0) --maxKernel;
            int bandKernel = std::min(lp.psfKernel, maxKernel);
            if (bandKernel % 2 == 0) --bandKernel;

            // Build the PSF rings ONCE per band, outside the convolution: the
            // effConvolve callback below only interpolates and rotates. spp
            // and bandKernel (this band's samplesPerPixel and output
            // window) are passed through so buildPsfRings measures axisEnergy
            // through the SAME window psfAtField actually resamples into --
            // see psfrings.hpp's Critical-2 fix: measuring it any other way
            // (e.g. over the whole ring) silently loses energy whenever the
            // window is narrower than the ring, which it always is at these
            // defaults (an unaberrated PSF is genuinely sub-pixel).
            // A grid that cannot span even a 3x3 kernel is telling us the PSF is
            // smaller than one pixel at this scale, and a sub-pixel PSF has
            // nothing to convolve: leave the band alone rather than fail.
            //
            // That is not a corner case, it is how a preview works. A proxy is
            // rendered by telling the model its pixels are larger, so a blur
            // covers the same FRACTION of the frame as it will at full size --
            // which shrinks this limit in exact proportion. Zoomed far enough
            // out, the blur genuinely is sub-pixel and drawing none is the
            // honest answer. Failing instead killed the whole preview with
            // "psfKernel (7) exceeds what psfGrid=256 covers at 408nm".
            if (bandKernel >= 3) {
            const optics::PsfRings rings = optics::buildPsfRings(
                    lp.pupil, [&](float t) { return wavefrontAt(lp, t, lambda); },
                    lambda, lp.lambdaHat, lp.psfRings, lp.psfGrid, spp, bandKernel);
                // Measure the resampler's own sub-pixel bias ON AXIS, where the PSF
                // is rotationally symmetric and any centroid offset is therefore
                // numerical rather than optical, and cancel it for every kernel.
                // Left in, it puts a bright rim along every high-contrast edge.
                float shiftX = 0.0f, shiftY = 0.0f;
                {
                    const Plane k0 = optics::psfAtField(rings, 0.0f, 0.0f, spp, bandKernel);
                    double sum = 0.0, mx = 0.0, my = 0.0;
                    for (int ky = 0; ky < k0.h; ++ky)
                        for (int kx = 0; kx < k0.w; ++kx) {
                            const double v = k0.at(kx, ky);
                            sum += v; mx += v * kx; my += v * ky;
                        }
                    if (sum > 1e-12) {
                        const double c = 0.5 * double(bandKernel - 1);
                        // Pixels to grid samples: the sampler steps spp grid samples
                        // for each output pixel.
                        shiftX = float((mx / sum - c) * spp);
                        shiftY = float((my / sum - c) * spp);
                    }
                }
    
                band = conv::effConvolve(band, lp.effPatch, [&](float cx, float cy) {
                    const float dx = (cx - frame.cx) / frame.halfDiag;
                    const float dy = (cy - frame.cy) / frame.halfDiag;
                    const float t  = std::clamp(std::sqrt(dx * dx + dy * dy), 0.0f, 1.0f);
                    Plane k = optics::psfAtField(rings, t, std::atan2(dy, dx), spp, bandKernel,
                                                 shiftX, shiftY, lp.psfSqueeze);
    
                    // Normalise each kernel to unit sum, then scale by the pupil's
                    // throughput at this field angle.
                    //
                    // psfAtField's energy depends on the ORIENTATION it is asked
                    // for, not just the field angle: measured at the frame corner on
                    // a perfect lens, the kernel sums 0.607 along the frame axes and
                    // 0.538 on the diagonals -- a 12% swing with the four-fold
                    // symmetry of the square resampling window. Since orientation
                    // varies with position, that painted a stationary cross of
                    // shading over the picture, clearly visible as a wedge through
                    // the centre of an otherwise flat frame.
                    //
                    // The RADIAL falloff, by contrast, was already right: kernel sum
                    // divided by pupil throughput held flat at 0.95 across the whole
                    // field, so the ring machinery was reproducing mechanical
                    // vignetting correctly and only its overall scale was off.
                    //
                    // So separate the two. Normalising here takes out both the 12%
                    // ripple and the 5% scale error and leaves the kernel as a pure
                    // redistribution of light. Mechanical vignetting is then applied
                    // per PIXEL by the vignette stage below, alongside the natural
                    // falloff it belongs with -- not per patch here, where the patch
                    // lattice is square and the falloff is radial, so blending
                    // neighbouring patches reintroduces an angular ripple of its own
                    // (measured 3.9% at mid-radius, against 12% for the original).
                    double sum = 0.0;
                    for (float v : k.v) sum += v;
                    if (sum > 1e-12) {
                        const float inv = float(1.0 / sum);
                        for (float& v : k.v) v *= inv;
                    }
    
                    return k;
                });
            }
        }

        // A plain linear combination -- which is exactly why it commutes with the
        // convolution above. Blur the bands and you get the blurred colour, with
        // no correction term to re-sharpen anything.
        const color::XYZ m = color::cmf(lambda);
        for (size_t i = 0; i < X.size(); ++i) {
            X[i] += weight * band.v[i] * m.x;
            Y[i] += weight * band.v[i] * m.y;
            Z[i] += weight * band.v[i] * m.z;
        }
        (void)bandIndex;
    }

    Image out(w, h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            const size_t i = size_t(y) * w + x;
            color::RGB c = color::xyzToRec2020(color::XYZ{
                X[i] / normY * wb[0], Y[i] / normY, Z[i] / normY * wb[1]});

            if (lp.doVignette) {
                const float dx = (float(x) - frame.cx) / frame.halfDiag;
                const float dy = (float(y) - frame.cy) / frame.halfDiag;
                const float t = std::clamp(std::sqrt(dx * dx + dy * dy), 0.0f, 1.0f);
                // Both physical terms, per pixel: natural cos^4 falloff and the
                // pupil's mechanical clipping. This stage is the ONLY thing in
                // the pipeline that darkens the frame.
                //
                // The mechanical term is close to LINEAR in field height -- two
                // circles sliding apart lose overlap at a near-constant rate --
                // and measured, the combined curve ran 0.94, 0.87, 0.80, 0.72,
                // 0.64, 0.57, 0.50 from a quarter of the way out to the corner.
                // A straight ramp does not look like a lens. Real vignetting is
                // flat across the middle and falls away sharply at the corner --
                // which is why it is conventionally modelled as a polynomial in
                // r SQUARED. So the loss is weighted by t^2, keeping the corner
                // where the physics puts it and flattening the middle.
                const float phys = optics::naturalFalloff(lp.vignette, t) * mechanicalAt(t);
                const float v = 1.0f - lp.vignetteAmount * (1.0f - phys) * t * t;
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
