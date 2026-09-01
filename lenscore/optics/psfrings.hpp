#pragma once
#include "lenscore/geometry.hpp"
#include "lenscore/optics/psf.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace lens::optics {

// One PSF per field radius; psfAtField rotates the nearest ring(s) into place
// rather than recomputing a PSF per angle. This exploits Task 10's radial
// Zernike frame: the PSF's shape depends only on field radius, so at any
// other angle it is the same kernel, rotated.
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
    // varying convolution actually wants.
    double axisEnergy = 1.0;
};

inline PsfRings buildPsfRings(const PupilParams& pp,
                              const std::function<Wavefront(float)>& wavefrontAtT,
                              float lambdaNm, float lambdaRefNm, int rings, int N) {
    PsfRings out;
    out.gridN = N;
    out.ring.reserve(size_t(rings));
    for (int i = 0; i < rings; ++i) {
        const float t = (rings > 1) ? float(i) / float(rings - 1) : 0.0f;
        out.ring.push_back(psfFromPupil(pp, wavefrontAtT(t), lambdaNm, lambdaRefNm, t, N));
    }
    if (!out.ring.empty()) {
        double s = 0.0;
        for (float v : out.ring.front().v) s += v;
        out.axisEnergy = (s > 0.0) ? s : 1.0;
    }
    return out;
}

// Interpolate between rings by radius, rotate into the radial frame, and
// resample from PSF grid samples to image pixels. Energy is conserved by
// scaling with the area ratio, then divided by axisEnergy so the on-axis
// kernel sums to 1 and off-axis kernels carry vignetting as a ratio to it.
inline Plane psfAtField(const PsfRings& r, float t, float thetaRad,
                        float samplesPerPixel, int outSize) {
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
            const float a = sampleBilinear(r.ring[i0],     sx, sy);
            const float b = sampleBilinear(r.ring[i0 + 1], sx, sy);
            out.at(x, y) = (1.0f - f) * a + f * b;
        }
    }
    // One output pixel now covers samplesPerPixel^2 grid samples; dividing by
    // axisEnergy makes the on-axis kernel integrate to 1 (see PsfRings above).
    const float scale = (samplesPerPixel * samplesPerPixel) / float(r.axisEnergy);
    for (float& v : out.v) v *= scale;
    return out;
}

}  // namespace lens::optics
