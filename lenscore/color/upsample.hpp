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

// Diagonal XYZ white-point scale that reconciles this file's equal-energy (illuminant
// E) spectral integration with the D65-referenced Rec.2020 matrices in cie.hpp. A flat
// (all-ones) reflectance spectrum integrated under illuminant E lands at a chromaticity
// near (1/3, 1/3); the Rec.2020 matrices assume a D65 white. Left uncorrected, a
// perfectly white spectrum converts to visibly tinted RGB -- and because evalSpectrum's
// sigmoid saturates at 1.0, no choice of Coeffs can compensate a red channel biased by
// more than 10%, so RGB{1,1,1} is not merely tinted, it becomes unreachable by
// fitCoeffs.
//
// This is NOT a von Kries transform: true von Kries adapts in cone-response (LMS)
// space. This is a plain diagonal scale applied directly in XYZ -- exact at the white
// point by construction, and only a linear extrapolation (not a physically modelled
// chromatic adaptation) away from it, so it is approximate for saturated colours. A
// physically stricter treatment would integrate the model spectrum against a real D65
// spectral power distribution instead of a flat illuminant E. That was a deliberate
// trade for this task, not an oversight, and is recorded here as a known limitation.
//
// Exported (not file-local) because any code that accumulates XYZ from cmf() and
// converts through xyzToRec2020 directly -- bypassing spectrumToRec2020 -- needs this
// same correction, or it will reproduce the same white-point bias. Computed once,
// lazily, from cie.hpp's own cmf() and rec2020ToXyz(); returns {scaleX, scaleZ} (Y
// needs no correction: both integrations already normalise it to 1).
inline std::array<float, 2> equalEnergyWhitePointScale() {
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

// Integrates the spectrum against the CMFs under illuminant E, at 5nm, then adapts the
// result onto the Rec2020 matrices' reference white (see equalEnergyWhitePointScale).
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
    const auto wb = equalEnergyWhitePointScale();
    acc.x *= wb[0]; acc.z *= wb[1];
    RGB out = xyzToRec2020(acc);
    // Highly saturated targets (near the Rec.2020 primaries) sit outside what a bounded
    // reflectance spectrum can reach -- the model gamut-maps them rather than
    // reproducing them exactly (accepted limitation, see the design spec), and the
    // XYZ->Rec2020 matrix can drive a channel slightly negative there. A negative
    // radiance is meaningless and would propagate as an artefact through the rest of
    // the pipeline, so clamp it at zero here, at the upsampling boundary.
    out.r = std::max(0.0f, out.r);
    out.g = std::max(0.0f, out.g);
    out.b = std::max(0.0f, out.b);
    return out;
}

// Levenberg-Marquardt on three parameters against three residuals.
//
// fitCoeffs expects a normalised target (max component <= 1); the caller carries any
// overall radiance scale (e.g. max(r,g,b) of the original pixel) separately, outside
// the model -- see the file-level note above. In debug builds that contract is
// checked by assert. An assert alone is not a contract: under NDEBUG it compiles away
// and a release build would silently fit an out-of-range target and hand back a poor,
// unsignalled fit. So the contract is also enforced on the release path: when the
// target's max component exceeds 1, it is normalised internally before fitting, and
// the coefficients for that normalised colour are returned -- the same semantics an
// in-range caller already gets, just recovered here rather than corrupted.
// The same Levenberg-Marquardt solve as fitCoeffs, against ANY forward model.
//
// fitCoeffs fits against the dense spectral integral, which is the right answer
// for building a table. A renderer integrating a handful of bands does not
// reproduce that integral, so a spectrum fitted to the dense model comes back as
// the wrong colour through a sparse one. Passing the renderer's own
// reconstruction in here fits the spectrum that its quadrature will actually
// turn back into the colour asked for.
template <typename Forward>
inline Coeffs fitCoeffsAgainst(const Forward& forward, const RGB& target, Coeffs guess) {
    auto residual = [&](const Coeffs& c) {
        const RGB got = forward(c);
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


inline Coeffs fitCoeffs(const RGB& target, Coeffs guess = {}) {
    assert(std::max({target.r, target.g, target.b}) <= 1.0001f &&
           "fitCoeffs expects a normalised colour; carry the scale separately");

    const float m = std::max({target.r, target.g, target.b, 1.0f});
    const RGB normalised = (m > 1.0f) ? RGB{target.r / m, target.g / m, target.b / m} : target;

    auto residual = [&](const Coeffs& c) {
        const RGB got = spectrumToRec2020(c);
        return std::array<float, 3>{got.r - normalised.r, got.g - normalised.g, got.b - normalised.b};
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
