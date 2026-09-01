#pragma once
#include "lenscore/constants.hpp"
#include "lenscore/plane.hpp"
#include <algorithm>
#include <cmath>

namespace lens::optics {

// Distance from centre to the aperture edge, as a fraction of the circumscribed
// radius. curvature 1 gives a circle, 0 a straight-sided polygon.
inline float apertureEdgeRadius(int blades, float curvature, float rotationRad, float theta) {
    if (blades < 3) return 1.0f;
    const float seg = 2.0f * kPi / float(blades);
    float a = std::fmod(theta + rotationRad, seg);
    if (a < 0.0f) a += seg;
    const float poly = std::cos(kPi / float(blades)) / std::cos(a - kPi / float(blades));
    const float c = std::clamp(curvature, 0.0f, 1.0f);
    return c * 1.0f + (1.0f - c) * poly;
}

// Ceiling on |apodizationSlope|. The ramp 1 + slope*t*u is Aggarwal's *planar*
// approximation to off-axis pupil illumination; their measurement was 31%
// variation (slope 0.31) at a 10 degree field angle on a 16mm lens, which is
// the realistic magnitude this model targets. The ramp turns negative -- and
// its mean over the disc stops being exactly 1, breaking the energy-neutrality
// this model promises -- once slope*t*|u| >= 1; the worst case is t=1, u=-1,
// i.e. slope >= 1. Bounding at 0.9 keeps the ramp strictly positive (floor
// 1 - 0.9 = 0.1) everywhere inside the unit pupil for every t in [0,1], while
// still leaving nearly 3x the measured value as headroom for more extreme
// lenses. Past this bound the small-perturbation approximation itself no
// longer describes the pupil (it would require the pupil to go fully dark on
// one side), so the parameter is clamped to the bound rather than the result
// renormalised after the fact -- renormalising would silently pretend the
// model still applies in a regime it does not.
inline constexpr float kMaxApodizationSlope = 0.9f;

struct PupilParams {
    int   blades           = 9;
    float curvature        = 0.15f;
    float rotationRad      = 0.0f;
    float apertureRadius   = 1.0f;   // in units of the wide-open radius
    float rEntrance        = 1.0f;
    float rExit            = 0.92f;
    float sepNorm          = 0.35f;
    float apodizationSlope = 0.0f;   // clamped to +/- kMaxApodizationSlope, see above
};

// Amplitude on an N x N grid spanning [-1,1]^2. The radial direction is +x.
inline Plane rasterPupil(const PupilParams& p, float t, int N) {
    Plane out(N, N);
    if (out.w == 0 || out.h == 0) return out;
    const float d = p.sepNorm * t;
    const float slope = std::clamp(p.apodizationSlope, -kMaxApodizationSlope, kMaxApodizationSlope);
    for (int j = 0; j < N; ++j) {
        const float v = 2.0f * float(j) / float(N - 1) - 1.0f;
        for (int i = 0; i < N; ++i) {
            const float u = 2.0f * float(i) / float(N - 1) - 1.0f;
            const float rho = std::sqrt(u * u + v * v);
            if (rho > p.apertureRadius) continue;

            const float th = std::atan2(v, u);
            if (rho > p.apertureRadius * apertureEdgeRadius(p.blades, p.curvature, p.rotationRad, th))
                continue;

            // Front and rear barrel apertures displace in opposite directions
            // along the radial axis (the direction of the field point) -- that
            // opposition is what makes the clipped shape a lemon, symmetric
            // about the radial axis, rather than a crescent.
            if (std::hypot(u - d, v) > p.rEntrance) continue;
            if (std::hypot(u + d, v) > p.rExit)     continue;

            // std::max is a defensive floor, not the mechanism that keeps the
            // ramp energy-neutral -- the clamp on slope above already
            // guarantees 1 + slope*t*u > 0 everywhere in the unit pupil, so
            // this should never actually trigger.
            out.at(i, j) = std::max(0.0f, 1.0f + slope * t * u);
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
