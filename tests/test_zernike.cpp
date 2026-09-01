#include <doctest/doctest.h>
#include "lenscore/optics/zernike.hpp"
#include <cmath>

using namespace lens::optics;

// M_PI is a POSIX extension, not standard C++; this project builds with
// CXX_EXTENSIONS OFF and targets MSVC too, so use a local constant instead.
static constexpr double kPi = 3.14159265358979323846;

// Numerically integrate Z_a * Z_b over the unit disc, divided by the disc area.
static double innerProduct(int n1, int m1, int n2, int m2) {
    const int NR = 400, NT = 720;
    double acc = 0.0, norm = 0.0;
    for (int i = 0; i < NR; ++i) {
        const double rho = (i + 0.5) / NR;
        for (int j = 0; j < NT; ++j) {
            const double th = 2.0 * kPi * (j + 0.5) / NT;
            const double w = rho;                       // Jacobian
            acc  += w * zernike(n1, m1, float(rho), float(th)) * zernike(n2, m2, float(rho), float(th));
            norm += w;
        }
    }
    return acc / norm;
}

TEST_CASE("piston is unity everywhere") {
    CHECK(zernike(0, 0, 0.0f, 0.0f) == doctest::Approx(1.0f));
    CHECK(zernike(0, 0, 1.0f, 2.0f) == doctest::Approx(1.0f));
}

TEST_CASE("defocus has the known radial form") {
    CHECK(zernikeRadial(2, 0, 0.0f) == doctest::Approx(-1.0f));   // 2*rho^2 - 1
    CHECK(zernikeRadial(2, 0, 1.0f) == doctest::Approx(1.0f));
    CHECK(zernike(2, 0, 1.0f, 0.0f) == doctest::Approx(std::sqrt(3.0f)).epsilon(1e-4));
}

TEST_CASE("spherical has the known radial form") {
    // 6*rho^4 - 6*rho^2 + 1
    CHECK(zernikeRadial(4, 0, 0.0f) == doctest::Approx(1.0f));
    CHECK(zernikeRadial(4, 0, 1.0f) == doctest::Approx(1.0f));
    CHECK(zernikeRadial(4, 0, 0.5f) == doctest::Approx(6 * 0.0625f - 6 * 0.25f + 1.0f));
}

TEST_CASE("the basis is orthonormal over the unit disc") {
    const int modes[][2] = {{0,0}, {2,0}, {2,2}, {2,-2}, {3,1}, {3,-1}, {4,0}};
    for (const auto& a : modes)
        for (const auto& b : modes) {
            const double ip = innerProduct(a[0], a[1], b[0], b[1]);
            const double want = (a[0] == b[0] && a[1] == b[1]) ? 1.0 : 0.0;
            CAPTURE(a[0]); CAPTURE(a[1]); CAPTURE(b[0]); CAPTURE(b[1]);
            CHECK(ip == doctest::Approx(want).epsilon(0.02));
        }
}

TEST_CASE("astigmatism separates the sagittal and tangential axes") {
    CHECK(zernike(2, 2, 1.0f, 0.0f) > 0.0f);
    CHECK(zernike(2, 2, 1.0f, float(kPi) / 2.0f) < 0.0f);
}

// Coma sign, derived from the standard Zernike definition (not from the code
// under test):
//   R_3^1(rho) = 3*rho^3 - 2*rho          (textbook coma radial polynomial)
//   R_n^m(1)   = 1 for every valid (n, m) (standard Zernike edge property:
//                3*1 - 2*1 = 1, confirming the formula at rho = 1)
//   N_3^1      = sqrt(2*(3+1) / (1 + 0))  = sqrt(8)   (m != 0, no delta_m0 term)
// OSA/ANSI convention: for m >= 0, Z_n^m = N * R_n^m * cos(m*theta);
//                       for m <  0, Z_n^m = N * R_n^|m| * sin(|m|*theta).
// At rho = 1, theta = 0:      cos(0) = 1        -> Z(3, 1, 1, 0)        = +sqrt(8)
// At rho = 1, theta = pi/2:   sin(pi/2) = 1      -> Z(3,-1, 1, pi/2)    = +sqrt(8)
// At rho = 1, theta = -pi/2:  sin(-pi/2) = -1    -> Z(3,-1, 1, -pi/2)   = -sqrt(8)
// A flipped cosine/sine sign, or a flipped angular frame, changes these exact
// values -- unlike a bare antisymmetry check, which any such flip still passes.
TEST_CASE("coma pins the sign of the standard convention, cosine branch") {
    CHECK(zernike(3, 1, 1.0f, 0.0f) == doctest::Approx(std::sqrt(8.0f)));
}

TEST_CASE("coma pins the sign of the standard convention, sine branch (m < 0)") {
    CHECK(zernike(3, -1, 1.0f, float(kPi) / 2.0f) == doctest::Approx(std::sqrt(8.0f)));
    CHECK(zernike(3, -1, 1.0f, -float(kPi) / 2.0f) == doctest::Approx(-std::sqrt(8.0f)));
}

TEST_CASE("coma is antisymmetric about the tangential axis") {
    const float a = zernike(3, 1, 0.8f, 0.0f);
    const float b = zernike(3, 1, 0.8f, float(kPi));
    CHECK(a == doctest::Approx(-b).epsilon(1e-5));
}

TEST_CASE("a zero wavefront is flat and a defocused one is not") {
    CHECK(wavefrontError(Wavefront{}, 0.7f, 1.0f) == doctest::Approx(0.0f));
    Wavefront w; w.defocus = 0.5f;
    CHECK(wavefrontError(w, 1.0f, 0.0f) != doctest::Approx(0.0f));
}
