#pragma once
#include "lenscore/plane.hpp"
#include <algorithm>
#include <cmath>

namespace lens::optics {

// Distance from centre to the aperture edge, as a fraction of the circumscribed
// radius. curvature 1 gives a circle, 0 a straight-sided polygon.
inline float apertureEdgeRadius(int blades, float curvature, float rotationRad, float theta) {
    if (blades < 3) return 1.0f;
    static constexpr float kPupilPi = 3.14159265358979323846f;
    const float seg = 2.0f * kPupilPi / float(blades);
    float a = std::fmod(theta + rotationRad, seg);
    if (a < 0.0f) a += seg;
    const float poly = std::cos(kPupilPi / float(blades)) / std::cos(a - kPupilPi / float(blades));
    const float c = std::clamp(curvature, 0.0f, 1.0f);
    return c * 1.0f + (1.0f - c) * poly;
}

struct PupilParams {
    int   blades           = 9;
    float curvature        = 0.15f;
    float rotationRad      = 0.0f;
    float apertureRadius   = 1.0f;   // in units of the wide-open radius
    float rEntrance        = 1.0f;
    float rExit            = 0.92f;
    float sepNorm          = 0.35f;
    float apodizationSlope = 0.0f;
};

// Amplitude on an N x N grid spanning [-1,1]^2. The radial direction is +x.
inline Plane rasterPupil(const PupilParams& p, float t, int N) {
    Plane out(N, N);
    if (out.w == 0 || out.h == 0) return out;
    const float d = p.sepNorm * t;
    for (int j = 0; j < N; ++j) {
        const float v = 2.0f * float(j) / float(N - 1) - 1.0f;
        for (int i = 0; i < N; ++i) {
            const float u = 2.0f * float(i) / float(N - 1) - 1.0f;
            const float rho = std::sqrt(u * u + v * v);
            if (rho > p.apertureRadius) continue;

            const float th = std::atan2(v, u);
            if (rho > p.apertureRadius * apertureEdgeRadius(p.blades, p.curvature, p.rotationRad, th))
                continue;

            // Front and rear barrel apertures displace in opposite directions.
            if (std::hypot(u - d, v) > p.rEntrance) continue;
            if (std::hypot(u + d, v) > p.rExit)     continue;

            out.at(i, j) = std::max(0.0f, 1.0f + p.apodizationSlope * t * u);
        }
    }
    return out;
}

// The optical vignetting term. Never normalised away.
inline float pupilEnergyFraction(const PupilParams& p, float t, int N) {
    const Plane amp = rasterPupil(p, t, N);
    if (amp.w == 0 || amp.h == 0) return 0.0f;
    double got = 0.0, full = 0.0;
    for (int j = 0; j < N; ++j) {
        const float v = 2.0f * float(j) / float(N - 1) - 1.0f;
        for (int i = 0; i < N; ++i) {
            const float u = 2.0f * float(i) / float(N - 1) - 1.0f;
            if (std::sqrt(u * u + v * v) <= p.apertureRadius) full += 1.0;
            got += amp.at(i, j);
        }
    }
    return (full > 0.0) ? float(got / full) : 0.0f;
}

}  // namespace lens::optics
