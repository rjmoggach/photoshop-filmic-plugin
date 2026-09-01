#pragma once
#include "lenscore/color/cie.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

namespace lens::color {

inline constexpr float kLambdaMin = 380.0f;
inline constexpr float kLambdaMax = 780.0f;

struct Coeffs { float c0 = 0, c1 = 0, c2 = 0; };

inline float normaliseLambda(float l) { return (l - kLambdaMin) / (kLambdaMax - kLambdaMin); }

inline float sigmoid(float x) { return 0.5f + x / (2.0f * std::sqrt(1.0f + x * x)); }

inline float evalSpectrum(const Coeffs& c, float lambda_nm) {
    const float n = normaliseLambda(lambda_nm);
    return sigmoid(std::fma(std::fma(c.c0, n, c.c1), n, c.c2));
}

namespace detail {

// A flat (all-ones) reflectance spectrum integrated under illuminant E lands at a
// chromaticity near (1/3, 1/3) -- not the D65 white the Rec2020 matrices assume -- so
// its raw XYZ needs a diagonal (von-Kries style) adaptation onto the matrices' own
// reference white before conversion, or a perfectly white spectrum converts to visibly
// tinted RGB. Scale computed once, lazily, from cie.hpp's own cmf() and rec2020ToXyz().
inline std::array<float, 2> whiteBalanceScale() {
    static const std::array<float, 2> k = [] {
        XYZ acc{0, 0, 0};
        float norm = 0.0f;
        for (float l = kLambdaMin; l <= kLambdaMax; l += 5.0f) {
            const XYZ m = cmf(l);
            acc.x += m.x; acc.y += m.y; acc.z += m.z;
            norm  += m.y;
        }
        acc.x /= norm; acc.y /= norm; acc.z /= norm;
        const XYZ ref = rec2020ToXyz(RGB{1, 1, 1});
        return std::array<float, 2>{ref.x / acc.x, ref.z / acc.z};
    }();
    return k;
}

}  // namespace detail

// Integrates the spectrum against the CMFs under illuminant E, at 5nm, then adapts the
// result onto the Rec2020 matrices' reference white (see detail::whiteBalanceScale).
inline RGB spectrumToRec2020(const Coeffs& c) {
    XYZ acc{0, 0, 0};
    float norm = 0.0f;
    for (float l = kLambdaMin; l <= kLambdaMax; l += 5.0f) {
        const XYZ m = cmf(l);
        const float f = evalSpectrum(c, l);
        acc.x += f * m.x; acc.y += f * m.y; acc.z += f * m.z;
        norm  += m.y;
    }
    acc.x /= norm; acc.y /= norm; acc.z /= norm;
    const auto wb = detail::whiteBalanceScale();
    acc.x *= wb[0]; acc.z *= wb[1];
    return xyzToRec2020(acc);
}

// Levenberg-Marquardt on three parameters against three residuals.
inline Coeffs fitCoeffs(const RGB& target, Coeffs guess = {}) {
    assert(std::max({target.r, target.g, target.b}) <= 1.0001f &&
           "fitCoeffs expects a normalised colour; carry the scale separately");

    auto residual = [&](const Coeffs& c) {
        const RGB got = spectrumToRec2020(c);
        return std::array<float, 3>{got.r - target.r, got.g - target.g, got.b - target.b};
    };

    Coeffs c = guess;
    float damping = 1e-3f;

    for (int iter = 0; iter < 200; ++iter) {
        const auto r = residual(c);
        const float err = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        if (err < 1e-10f) break;

        // Central-difference Jacobian, 3x3.
        float J[3][3];
        float* p[3] = {&c.c0, &c.c1, &c.c2};
        for (int j = 0; j < 3; ++j) {
            const float save = *p[j];
            const float h = std::max(1e-3f, std::abs(save) * 1e-3f);
            *p[j] = save + h; const auto rp = residual(c);
            *p[j] = save - h; const auto rm = residual(c);
            *p[j] = save;
            for (int i = 0; i < 3; ++i) J[i][j] = (rp[i] - rm[i]) / (2.0f * h);
        }

        // Normal equations (J^T J + damping I) d = -J^T r, solved by elimination.
        double A[3][4] = {};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double s = 0; for (int k = 0; k < 3; ++k) s += double(J[k][i]) * J[k][j];
                A[i][j] = s + (i == j ? damping : 0.0);
            }
            double s = 0; for (int k = 0; k < 3; ++k) s += double(J[k][i]) * r[k];
            A[i][3] = -s;
        }
        for (int i = 0; i < 3; ++i) {
            int piv = i;
            for (int k = i + 1; k < 3; ++k) if (std::abs(A[k][i]) > std::abs(A[piv][i])) piv = k;
            if (std::abs(A[piv][i]) < 1e-18) { damping *= 10.0f; goto next; }
            for (int j = 0; j < 4; ++j) std::swap(A[i][j], A[piv][j]);
            for (int k = 0; k < 3; ++k) {
                if (k == i) continue;
                const double f = A[k][i] / A[i][i];
                for (int j = i; j < 4; ++j) A[k][j] -= f * A[i][j];
            }
        }
        {
            Coeffs trial{c.c0 + float(A[0][3] / A[0][0]),
                         c.c1 + float(A[1][3] / A[1][1]),
                         c.c2 + float(A[2][3] / A[2][2])};
            const auto rt = residual(trial);
            const float et = rt[0] * rt[0] + rt[1] * rt[1] + rt[2] * rt[2];
            if (et < err) { c = trial; damping = std::max(1e-7f, damping * 0.3f); }
            else          { damping *= 10.0f; }
        }
        next:;
        if (damping > 1e7f) break;
    }
    return c;
}

}  // namespace lens::color
