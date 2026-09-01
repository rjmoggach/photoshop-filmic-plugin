#pragma once
#include <cassert>
#include <cmath>

namespace lens::color {

inline float srgbToLinear(float s) {
    return (s <= 0.04045f) ? s / 12.92f : std::pow((s + 0.055f) / 1.055f, 2.4f);
}

inline float linearToSrgb(float l) {
    return (l <= 0.0031308f) ? l * 12.92f : 1.055f * std::pow(l, 1.0f / 2.4f) - 0.055f;
}

struct Knee {
    float threshold = 0.85f;
    float peak      = 8.0f;
};

// Expands display-referred highlights back above 1.0 so that bloom, halation
// and bokeh have real energy to work with. C1 continuous at the threshold.
inline float expandHighlights(float x, const Knee& k) {
    assert(k.peak > 1.0f && k.threshold > 0.0f && k.threshold < 1.0f);
    if (x <= k.threshold) return x;
    const float A = 1.0f - k.threshold;
    const float B = k.peak - k.threshold;
    const float c = 1.0f - A / B;
    const float u = (x - k.threshold) / A;
    if (u <= 1.0f) return k.threshold + A * u / (1.0f - c * u);
    return k.peak + (x - 1.0f) * (B * B) / (A * A);   // linear continuation
}

inline float compressHighlights(float y, const Knee& k) {
    if (y <= k.threshold) return y;
    const float A = 1.0f - k.threshold;
    const float B = k.peak - k.threshold;
    const float c = 1.0f - A / B;
    if (y <= k.peak) {
        const float d = y - k.threshold;
        return k.threshold + A * (d / (A + c * d));
    }
    return 1.0f + (y - k.peak) * (A * A) / (B * B);
}

}  // namespace lens::color
