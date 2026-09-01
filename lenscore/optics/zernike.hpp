#pragma once
#include <cmath>
#include <cstdlib>

namespace lens::optics {

inline double factorial(int n) {
    double f = 1.0;
    for (int i = 2; i <= n; ++i) f *= i;
    return f;
}

inline float zernikeRadial(int n, int m, float rho) {
    m = std::abs(m);
    if ((n - m) % 2 != 0) return 0.0f;
    double acc = 0.0;
    for (int k = 0; k <= (n - m) / 2; ++k) {
        const double num = ((k % 2) ? -1.0 : 1.0) * factorial(n - k);
        const double den = factorial(k) * factorial((n + m) / 2 - k) * factorial((n - m) / 2 - k);
        acc += num / den * std::pow(double(rho), n - 2 * k);
    }
    return float(acc);
}

// Orthonormal over the unit disc: the mean of Z^2 is exactly 1.
inline float zernike(int n, int m, float rho, float theta) {
    const float norm = std::sqrt(float(2 * (n + 1)) / (m == 0 ? 2.0f : 1.0f));
    const float r = zernikeRadial(n, m, rho);
    if (m == 0) return norm * r;
    return norm * r * (m > 0 ? std::cos(m * theta) : std::sin(-m * theta));
}

// Coefficients are in waves. theta is measured in the radial frame, from the
// direction pointing away from the image centre.
struct Wavefront {
    float defocus   = 0.0f;   // Z(2,0)
    float astig     = 0.0f;   // Z(2,2)
    float coma      = 0.0f;   // Z(3,1)
    float spherical = 0.0f;   // Z(4,0)
};

inline float wavefrontError(const Wavefront& w, float rho, float theta) {
    return w.defocus   * zernike(2, 0, rho, theta)
         + w.astig     * zernike(2, 2, rho, theta)
         + w.coma      * zernike(3, 1, rho, theta)
         + w.spherical * zernike(4, 0, rho, theta);
}

}  // namespace lens::optics
