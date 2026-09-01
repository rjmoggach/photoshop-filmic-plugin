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

// Internal helper to validate Knee and compute coefficients.
// Returns true if the knee is valid and coefficients are safe to use.
// A degenerate knee (threshold out of bounds, peak <= 1, or peak <= threshold)
// returns false, and both functions treat it as identity.
struct KneeCoefficients {
    float A = 0.0f;
    float B = 0.0f;
    float c = 0.0f;
    bool valid = false;
};

inline KneeCoefficients validateKnee(const Knee& k) {
    KneeCoefficients result;

    // Check all degenerate conditions
    if (!(k.threshold > 0.0f && k.threshold < 1.0f)) return result;
    if (!(k.peak > 1.0f)) return result;
    if (!(k.peak > k.threshold)) return result;

    result.A = 1.0f - k.threshold;
    result.B = k.peak - k.threshold;
    result.c = 1.0f - result.A / result.B;
    result.valid = true;

    return result;
}

// Expands display-referred highlights back above 1.0 so that bloom, halation
// and bokeh have real energy to work with. C1 continuous at the threshold.
inline float expandHighlights(float x, const Knee& k) {
    KneeCoefficients coef = validateKnee(k);
    if (!coef.valid) return x;  // Degenerate knee: return input unchanged

    assert(k.peak > 1.0f && k.threshold > 0.0f && k.threshold < 1.0f);

    if (x <= k.threshold) return x;
    const float u = (x - k.threshold) / coef.A;
    if (u <= 1.0f) return k.threshold + coef.A * u / (1.0f - coef.c * u);
    return k.peak + (x - 1.0f) * (coef.B * coef.B) / (coef.A * coef.A);   // linear continuation
}

inline float compressHighlights(float y, const Knee& k) {
    KneeCoefficients coef = validateKnee(k);
    if (!coef.valid) return y;  // Degenerate knee: return input unchanged

    if (y <= k.threshold) return y;
    if (y <= k.peak) {
        const float d = y - k.threshold;
        return k.threshold + coef.A * (d / (coef.A + coef.c * d));
    }
    return 1.0f + (y - k.peak) * (coef.A * coef.A) / (coef.B * coef.B);
}

}  // namespace lens::color
