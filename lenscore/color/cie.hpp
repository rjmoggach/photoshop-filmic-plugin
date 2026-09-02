#pragma once
#include <cmath>

namespace lens::color {

struct XYZ { float x = 0, y = 0, z = 0; };
struct RGB { float r = 0, g = 0, b = 0; };

// Piecewise Gaussian: sigma differs either side of the mean.
inline float lobe(float x, float mu, float s1, float s2) {
    const float t = (x - mu) / (x < mu ? s1 : s2);
    return std::exp(-0.5f * t * t);
}

// Multi-lobe analytic fit to the CIE 1931 2-degree observer.
inline XYZ cmf(float lambda_nm) {
    XYZ o;
    o.x = 1.056f * lobe(lambda_nm, 599.8f, 37.9f, 31.0f)
        + 0.362f * lobe(lambda_nm, 442.0f, 16.0f, 26.7f)
        - 0.065f * lobe(lambda_nm, 501.1f, 20.4f, 26.2f);
    o.y = 0.821f * lobe(lambda_nm, 568.8f, 46.9f, 40.5f)
        + 0.286f * lobe(lambda_nm, 530.9f, 16.3f, 31.1f);
    o.z = 1.217f * lobe(lambda_nm, 437.0f, 11.8f, 36.0f)
        + 0.681f * lobe(lambda_nm, 459.0f, 26.0f, 13.8f);
    return o;
}

inline RGB xyzToRec2020(const XYZ& c) {
    return { 1.7166511880f * c.x - 0.3556707838f * c.y - 0.2533662814f * c.z,
            -0.6666843518f * c.x + 1.6164812366f * c.y + 0.0157685458f * c.z,
             0.0176398574f * c.x - 0.0427706133f * c.y + 0.9421031212f * c.z };
}

inline XYZ rec2020ToXyz(const RGB& c) {
    return { 0.6369580483f * c.r + 0.1446169036f * c.g + 0.1688809752f * c.b,
             0.2627002120f * c.r + 0.6779980715f * c.g + 0.0593017165f * c.b,
             0.0000000000f * c.r + 0.0280726930f * c.g + 1.0609850577f * c.b };
}

// Rec.2020 relative luminance: exactly the Y row of rec2020ToXyz above,
// exposed as its own function so a caller that only wants scalar luma (a
// test harness building a luminance image from a linear Rec.2020 render,
// say) shares the SAME three coefficients as the colour-management matrix
// instead of re-deriving a rounded copy that could silently drift from it.
inline float rec2020Luma(const RGB& c) {
    return 0.2627002120f * c.r + 0.6779980715f * c.g + 0.0593017165f * c.b;
}

// Divisor that keeps an equal-energy spectrum neutral under any band layout.
// Returns the sum of (weight * ybar) across all bands.
// Used as a divisor: an empty band list (n <= 0) or null pointers return 1.0f (neutral)
// to avoid inf/nan poisoning at the call site; an unnormalised spectrum is wrong but finite.
inline float cieYNormalisation(const float* lambdas, const float* weights, int n) {
    if (n <= 0 || !lambdas || !weights) return 1.0f;
    float s = 0.0f;
    for (int i = 0; i < n; ++i) s += weights[i] * cmf(lambdas[i]).y;
    return s;
}

}  // namespace lens::color
