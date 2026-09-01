#pragma once
#include "lenscore/geometry.hpp"
#include "lenscore/optics/distortion.hpp"
#include <cmath>

namespace lens::optics {

// m(lambda, t) = 1 + k_l * (lambda_hat - lambda) * t
inline float lateralMagnification(float k_l, float lambda_nm, float lambda_hat_nm, float t) {
    return 1.0f + k_l * (lambda_hat_nm - lambda_nm) * t;
}

// Solves t_out = t_src * (1 + K * t_src) for t_src.
//
// The textbook form (-1 + sqrt(1 + 4*K*tOut)) / (2*K) suffers catastrophic
// cancellation for small K: sqrt(1 + 4*K*tOut) sits close to 1, so the
// subtraction destroys most of its significant digits before the result is
// divided by the (also small) K, amplifying what little precision remains.
// Rationalising the numerator -- multiplying by (1 + sqrt(disc)) over itself --
// removes the subtraction entirely and is exact in the K -> 0 limit, so the
// K == 0 case falls out of the same expression with no branch needed for it.
// The computation runs in double so the near-cancellation inside the root
// itself keeps enough precision before the final cast back to float.
inline float inverseLateralRadius(float K, float tOut) {
    const double Kd = double(K), tOutd = double(tOut);
    const double disc = 1.0 + 4.0 * Kd * tOutd;
    if (disc <= 0.0) return tOut;                      // outside the invertible range
    return float(2.0 * tOutd / (1.0 + std::sqrt(disc)));
}

// One resample applies magnification and distortion together, never twice.
inline Plane warpPlane(const Plane& src, const Distortion& d, float K) {
    const Frame f = frameOf(src.w, src.h);
    Plane out(src.w, src.h);
    for (int y = 0; y < src.h; ++y) {
        for (int x = 0; x < src.w; ++x) {
            const float nx = (float(x) - f.cx) / f.halfDiag;
            const float ny = (float(y) - f.cy) / f.halfDiag;
            const float tOut = std::sqrt(nx * nx + ny * ny);
            float ux = nx, uy = ny;
            if (tOut > 1e-9f) {
                const float tSrc = inverseLateralRadius(K, tOut);
                const float s = tSrc / tOut;
                ux = nx * s; uy = ny * s;
            }
            float dx = 0, dy = 0;
            applyDistortion(d, ux, uy, dx, dy);
            out.at(x, y) = sampleBilinear(src, dx * f.halfDiag + f.cx, dy * f.halfDiag + f.cy);
        }
    }
    return out;
}

}  // namespace lens::optics
