#pragma once
#include "lenscore/geometry.hpp"
#include "lenscore/optics/psf.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace lens::optics {

// One PSF per field radius; psfAtField rotates the nearest ring(s) into place
// rather than recomputing a PSF per angle. This exploits Task 10's radial
// Zernike frame: the PSF's shape depends only on field radius, so at any
// other angle it is the same kernel, rotated.
//
// ring.size() must be >= 2: psfAtField always interpolates between a pair of
// adjacent rings (see i0/i0+1 below), so build via buildPsfRings, which
// enforces this, rather than assembling one by hand.
struct PsfRings {
    std::vector<Plane> ring;   // index i is field radius t = i/(rings-1)
    int gridN = 0;

    // Summed energy of ring 0 (the on-axis kernel), computed once when the
    // rings are built. psfFromPupil deliberately returns a raw, unnormalised
    // PSF -- its total energy is the optical vignetting -- but a convolution
    // kernel cannot carry that literally, or it multiplies every image by
    // N^2 * sum(A^2) (Parseval). Dividing every ring by axisEnergy instead
    // makes the on-axis kernel sum to 1 and every off-axis kernel sum to the
    // vignetting fraction *relative to the axis*, which is what a spatially
    // varying convolution actually wants. See axisEnergyOf for what happens
    // when ring 0 itself carries no usable energy.
    double axisEnergy = 1.0;
};

// Below this, a ring's total energy is not a safe normalisation reference.
// Exact zero happens for a fully occluded on-axis pupil (psfFromPupil writes
// literal 0.0f everywhere); a value that is technically positive but many
// orders of magnitude smaller -- denormal-range floating point noise, or a
// pathologically near-empty aperture -- is just as dangerous, because 1/s
// would still amplify every kernel by an arbitrary, physically meaningless
// factor (e.g. an axisEnergy of 1e-30 silently blows every kernel up by
// 1e30). The smallest energy a single realistically-illuminated pixel can
// contribute under this pupil model is many orders of magnitude above this
// floor, so it only ever fires for genuinely degenerate input.
inline constexpr double kMinAxisEnergy = 1e-6;

// Sum of a ring's raw intensity, floored against kMinAxisEnergy. Below the
// floor, the fallback is axisEnergy = 1.0: psfAtField's scale becomes just
// the resampling area factor, i.e. every ring keeps its raw psfFromPupil
// magnitude, unrescaled, instead of being divided by noise. Exposed as its
// own function (rather than inlined into buildPsfRings) so the floor is
// unit-testable directly, without needing a pupil configuration that happens
// to trigger it.
inline double axisEnergyOf(const Plane& axisRing) {
    double s = 0.0;
    for (float v : axisRing.v) s += v;
    return (s > kMinAxisEnergy) ? s : 1.0;
}

inline PsfRings buildPsfRings(const PupilParams& pp,
                              const std::function<Wavefront(float)>& wavefrontAtT,
                              float lambdaNm, float lambdaRefNm, int rings, int N) {
    // psfAtField always reads ring[i0] and ring[i0+1] together, so fewer than
    // two rings is a contract violation by the caller, not bad input data --
    // validated before the reserve below, which would otherwise turn a
    // negative rings into an enormous size_t and attempt a huge allocation.
    if (rings < 2) {
        throw std::invalid_argument(
            "buildPsfRings: rings must be >= 2 (psfAtField interpolates "
            "between a pair of rings); got " + std::to_string(rings));
    }

    PsfRings out;
    out.gridN = N;
    out.ring.reserve(size_t(rings));
    for (int i = 0; i < rings; ++i) {
        const float t = float(i) / float(rings - 1);
        out.ring.push_back(psfFromPupil(pp, wavefrontAtT(t), lambdaNm, lambdaRefNm, t, N));
    }
    out.axisEnergy = axisEnergyOf(out.ring.front());
    return out;
}

// Averages the ring's raw grid samples in a k x k box, roughly centred on
// (cx, cy) (floor-and-offset, not exact sub-sample centring -- "roughly" is
// the honest word for it), clamped at the ring's edge the same way
// sampleBilinear is. Used only when one output pixel covers MORE than one
// ring sample -- see psfAtField's oversampled branch below -- as a simple,
// deliberately unrigorous anti-alias: a point sample of a PSF narrower than
// the box would alias (miss it, or land squarely on its peak), and a box
// average is the cheapest fix that does not.
inline float sampleBoxAverage(const Plane& p, float cx, float cy, int k) {
    const int cxi = int(std::floor(cx));
    const int cyi = int(std::floor(cy));
    const int half = k / 2;
    double sum = 0.0;
    for (int j = 0; j < k; ++j) {
        const int yy = std::clamp(cyi - half + j, 0, p.h - 1);
        for (int i = 0; i < k; ++i) {
            const int xx = std::clamp(cxi - half + i, 0, p.w - 1);
            sum += double(p.at(xx, yy));
        }
    }
    return float(sum / double(k) / double(k));
}

// Interpolate between rings by radius, rotate into the radial frame, and
// resample from PSF grid samples to image pixels. Energy is conserved by
// scaling with the area ratio, then divided by axisEnergy so the on-axis
// kernel sums to 1 and off-axis kernels carry vignetting as a ratio to it.
//
// Which resampling is correct depends on which is finer, the output pixel or
// the ring's own grid sample:
//   samplesPerPixel <= 1: the output pixel is FINER than a ring sample (a
//     well-resolved PSF, or heavy oversampling), so bilinear point-sampling
//     the ring is the right thing to do -- there is no sub-pixel footprint to
//     average over.
//   samplesPerPixel >  1: the output pixel is COARSER than a ring sample --
//     at the defaults, an unaberrated diffraction spot is genuinely narrower
//     than one output pixel. Point-sampling a function narrower than the
//     sample spacing aliases: the single sample is whatever the ring happens
//     to be exactly at the pixel centre, and the areaFactor below then
//     multiplies THAT by the whole pixel's area, over- or under-counting by
//     however far that one sample was from the pixel's true average.
//     Averaging the ring over the pixel's own footprint (sampleBoxAverage)
//     instead of point-sampling it removes that aliasing; this was verified
//     by measuring the on-axis kernel's total energy (which must sum to 1)
//     across bare-Airy through 20-wave-defocus PSFs -- point-sampling
//     conserved energy only for the aberrated (several-pixel-wide) cases and
//     overshot the unaberrated, sub-pixel case by 16x; box-averaging holds
//     within a couple of percent at every one of them (see test_psfrings.cpp).
inline Plane psfAtField(const PsfRings& r, float t, float thetaRad,
                        float samplesPerPixel, int outSize) {
    // r was not necessarily built by buildPsfRings, so its invariant (at
    // least 2 rings) cannot be assumed to hold; without this check a bad
    // PsfRings drives i0 to -1 below and r.ring[-1] is out-of-bounds.
    if (r.ring.size() < 2) {
        throw std::invalid_argument(
            "psfAtField: PsfRings must contain at least 2 rings (got " +
            std::to_string(r.ring.size()) +
            "); build it with buildPsfRings(..., rings >= 2, ...)");
    }

    const int n = int(r.ring.size());
    const float ft = std::clamp(t, 0.0f, 1.0f) * float(n - 1);
    const int i0 = std::clamp(int(ft), 0, n - 2);
    const float f = ft - float(i0);

    const float ct = std::cos(thetaRad), st = std::sin(thetaRad);
    const float c = 0.5f * float(outSize - 1);
    // The ring's true centre is the FFT-shift centre psfFromPupil actually
    // produces (fftShift2d places DC at the integer index gridN/2, which is
    // where an unaberrated PSF peaks -- see test_psf.cpp), not the geometric
    // midpoint (gridN-1)/2. For even gridN those differ by half a sample;
    // sampling from the wrong one mis-registers the ring before it is even
    // rotated, and that mis-registration shows up as a spurious, rotation-
    // dependent asymmetry once the frame is rotated into place.
    const float gc = float(r.gridN / 2);

    const bool oversampled = samplesPerPixel > 1.0f;
    const int  boxK = oversampled ? std::max(1, int(std::ceil(samplesPerPixel))) : 1;

    Plane out(outSize, outSize);
    for (int y = 0; y < outSize; ++y) {
        for (int x = 0; x < outSize; ++x) {
            // Output pixel offset, rotated back into the ring's radial frame
            // -- the INVERSE of the forward rotation that would carry a point
            // from the radial frame out to angle thetaRad.
            const float dx = (float(x) - c) * samplesPerPixel;
            const float dy = (float(y) - c) * samplesPerPixel;
            const float ux =  ct * dx + st * dy;
            const float uy = -st * dx + ct * dy;
            const float sx = gc + ux, sy = gc + uy;
            float a, b;
            if (oversampled) {
                a = sampleBoxAverage(r.ring[i0],     sx, sy, boxK);
                b = sampleBoxAverage(r.ring[i0 + 1], sx, sy, boxK);
            } else {
                a = sampleBilinear(r.ring[i0],     sx, sy);
                b = sampleBilinear(r.ring[i0 + 1], sx, sy);
            }
            out.at(x, y) = (1.0f - f) * a + f * b;
        }
    }
    // samplesPerPixel <= 1: out.at() holds a point sample; one output pixel
    // covers samplesPerPixel^2 grid samples, so that is the area factor.
    // samplesPerPixel >  1: out.at() holds a MEAN over boxK^2 grid samples
    // (sampleBoxAverage above). Multiplying by samplesPerPixel^2 here as well
    // would double-count the averaging -- the mean already divided by boxK^2,
    // so the area factor that turns it back into the discrete SUM of ring
    // energy the pixel's box covers is boxK^2, not samplesPerPixel^2 (they
    // are close but not identical: boxK = ceil(samplesPerPixel)). That SUM,
    // as a fraction of axisEnergy (itself a discrete sum over the whole
    // ring), is what actually makes the on-axis kernel integrate to 1.
    const float areaFactor = oversampled ? float(boxK * boxK) : (samplesPerPixel * samplesPerPixel);
    const float scale = areaFactor / float(r.axisEnergy);
    for (float& v : out.v) v *= scale;
    return out;
}

}  // namespace lens::optics
