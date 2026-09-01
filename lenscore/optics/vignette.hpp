#pragma once
#include "lenscore/constants.hpp"
#include <algorithm>
#include <cmath>

namespace lens::optics {

// Area of the lens-shaped intersection of two circles, centres d apart.
inline float circleOverlapArea(float d, float r1, float r2) {
    if (d >= r1 + r2) return 0.0f;
    const float rmin = std::min(r1, r2);
    if (d <= std::abs(r1 - r2)) return kPi * rmin * rmin;
    const float a1 = std::acos(std::clamp((d * d + r1 * r1 - r2 * r2) / (2.0f * d * r1), -1.0f, 1.0f));
    const float a2 = std::acos(std::clamp((d * d + r2 * r2 - r1 * r1) / (2.0f * d * r2), -1.0f, 1.0f));
    const float tri = 0.5f * std::sqrt(std::max(0.0f,
        (-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2)));
    return r1 * r1 * a1 + r2 * r2 * a2 - tri;
}

struct VignetteParams {
    float focal_mm          = 32.0f;
    float sensorHalfDiag_mm = 14.0f;
    float naturalExp        = 4.0f;   // free, not pinned at 4
    float tStopWide         = 2.0f;
    float tStop             = 2.0f;
    float rEntrance         = 1.0f;   // in units of the wide-open aperture radius
    float sepNorm           = 0.5f;   // pupil displacement per unit field radius
};

inline float naturalFalloff(const VignetteParams& p, float t) {
    // A non-positive focal length or sensor half-diagonal is a degenerate lens: there
    // is no field angle to fall off with. Without this guard, focal_mm == 0 with
    // t == 0 computes atan(0.0f / 0.0f), which is NaN and would silently propagate
    // through pow() and vignette() into every pixel. Treat it as "no natural
    // falloff" (1.0) instead, which is the least surprising value.
    if (p.focal_mm <= 0.0f || p.sensorHalfDiag_mm <= 0.0f) return 1.0f;
    const float theta = std::atan(t * p.sensorHalfDiag_mm / p.focal_mm);
    return std::pow(std::cos(theta), p.naturalExp);
}

// Aperture radius in units of the wide-open radius: 1.0 wide open, smaller stopped
// down. tStop <= 0 is invalid input; the division below then yields +inf (or, for
// negative tStop, a non-finite/negative value), which mechanicalFraction handles
// explicitly rather than relying on IEEE-754 inf/inf arithmetic.
inline float apertureRadius(const VignetteParams& p) { return p.tStopWide / p.tStop; }

// Sentinel meaning "mechanical vignetting never vanishes for this geometry" —
// reached when sepNorm >= rEntrance, so the aperture can never fit entirely inside
// the offset entrance circle at any tStop.
inline constexpr float kVanishStopNever = 1e9f;

inline float mechanicalVanishStop(const VignetteParams& p) {
    const float room = p.rEntrance - p.sepNorm;
    return (room > 0.0f) ? p.tStopWide / room : kVanishStopNever;
}

inline float mechanicalFraction(const VignetteParams& p, float t) {
    const float a = apertureRadius(p);
    if (a <= 0.0f) return 1.0f;
    // tStop <= 0 makes `a` non-finite (see apertureRadius above). Treat that
    // explicitly as full clipping (no light) instead of depending on inf/inf
    // arithmetic downstream, which a fast-math build could change.
    if (!std::isfinite(a)) return 0.0f;
    const float d = p.sepNorm * t;
    const float frac = circleOverlapArea(d, a, p.rEntrance) / (kPi * a * a);
    return std::clamp(frac, 0.0f, 1.0f);
}

inline float vignette(const VignetteParams& p, float t) {
    return naturalFalloff(p, t) * mechanicalFraction(p, t);
}

}  // namespace lens::optics
