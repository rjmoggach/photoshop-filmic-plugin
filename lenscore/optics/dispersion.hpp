#pragma once
#include <cmath>
#include <vector>

namespace lens::optics {

struct Dispersion {
    float B[3] = {1.03961212f, 0.231792344f, 1.01046945f};   // BK7
    float C[3] = {0.00600069867f, 0.0200179144f, 103.560653f};
    // 0 singlet, 2 achromat, 3 apochromat. Precondition: entries must be
    // distinct wavelengths (nanometres) -- a repeated (or too-close) entry
    // would divide by zero in the Lagrange interpolation. focusError guards
    // this: on a degenerate list it does not silently invent a correction,
    // it falls back to the raw (uncorrected) error instead.
    std::vector<float> correction_nm{};
    float residual = 1.0f;
};

inline Dispersion bk7() { return Dispersion{}; }

// Sellmeier equation. Note the constants above are defined for wavelength in
// MICROMETRES, while this function (like every function in this library)
// takes NANOMETRES -- the division by 1000 below is that conversion, not a
// simplification to remove.
inline float refractiveIndex(const Dispersion& d, float lambda_nm) {
    const double l = lambda_nm / 1000.0;      // micrometres
    const double l2 = l * l;
    double n2 = 1.0;
    for (int i = 0; i < 3; ++i) n2 += double(d.B[i]) * l2 / (l2 - double(d.C[i]));
    return float(std::sqrt(n2));
}

// Raw relative focal-length error before correction: F(lambda)/F_ref - 1.
// Focal length goes as 1/(n-1), so the ratio is (nRef-1)/(n-1). NOT (n-1)/(nRef-1),
// which is the relative optical power -- the negative of this -- and would flip the
// sign of chromatic defocus everywhere downstream.
inline float rawFocusError(const Dispersion& d, float lambda_nm, float lambda_ref_nm) {
    const float nRef = refractiveIndex(d, lambda_ref_nm);
    const float n    = refractiveIndex(d, lambda_nm);
    return (nRef - 1.0f) / (n - 1.0f) - 1.0f;
}

// Relative focal-length error, after achromatic/apochromatic correction.
// The number of entries in correction_nm selects
// the lens class with no branching: 0 = uncorrected singlet (raw Sellmeier
// error, monotonic), 2 = achromat (zero at two lines), 3 = apochromat (zero
// at three lines). What remains after subtracting the Lagrange interpolant
// through the correction wavelengths -- evaluated in x = 1/lambda^2, the
// variable in which chromatic focus is very nearly linear for normal glass
// -- is the secondary spectrum.
inline float focusError(const Dispersion& d, float lambda_nm, float lambda_ref_nm) {
    const float raw = rawFocusError(d, lambda_nm, lambda_ref_nm);
    const size_t k = d.correction_nm.size();
    if (k == 0) return raw * d.residual;

    // Lagrange interpolation of the raw error at the correction points,
    // in x = 1/lambda^2 with lambda in micrometres.
    auto xOf = [](float nm) { const double u = nm / 1000.0; return 1.0 / (u * u); };

    // Guard: a duplicate (or too-close) correction wavelength makes the
    // Lagrange denominator below exactly zero, propagating Inf/NaN into
    // every downstream wavefront. A malformed correction list should not
    // silently invent a correction, so fall back to the uncorrected error.
    for (size_t i = 0; i < k; ++i)
        for (size_t j = i + 1; j < k; ++j)
            if (std::fabs(xOf(d.correction_nm[i]) - xOf(d.correction_nm[j])) < 1e-9)
                return raw * d.residual;

    const double x = xOf(lambda_nm);
    double corr = 0.0;
    for (size_t i = 0; i < k; ++i) {
        double term = rawFocusError(d, d.correction_nm[i], lambda_ref_nm);
        for (size_t j = 0; j < k; ++j) {
            if (j == i) continue;
            term *= (x - xOf(d.correction_nm[j])) / (xOf(d.correction_nm[i]) - xOf(d.correction_nm[j]));
        }
        corr += term;
    }
    return float((raw - corr) * d.residual);
}

}  // namespace lens::optics
