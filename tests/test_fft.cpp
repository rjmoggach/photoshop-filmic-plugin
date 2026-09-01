#include <doctest/doctest.h>
#include "lenscore/conv/fft.hpp"
#include <cmath>

using namespace lens::conv;

TEST_CASE("power of two detection") {
    CHECK(isPowerOfTwo(1));  CHECK(isPowerOfTwo(64));
    CHECK_FALSE(isPowerOfTwo(0)); CHECK_FALSE(isPowerOfTwo(48));
}

TEST_CASE("transform of a delta is flat") {
    std::vector<Cplx> v(16, Cplx(0, 0));
    v[0] = Cplx(1, 0);
    fft1d(v, false);
    for (const Cplx& c : v) {
        CHECK(c.real() == doctest::Approx(1.0f).epsilon(1e-5));
        CHECK(c.imag() == doctest::Approx(0.0f).epsilon(1e-5));
    }
}

TEST_CASE("transform of a constant is a delta at DC") {
    std::vector<Cplx> v(16, Cplx(1, 0));
    fft1d(v, false);
    CHECK(std::abs(v[0]) == doctest::Approx(16.0f).epsilon(1e-4));
    for (size_t i = 1; i < v.size(); ++i) CHECK(std::abs(v[i]) == doctest::Approx(0.0f).epsilon(1e-4));
}

TEST_CASE("forward then inverse is the identity in 1d") {
    std::vector<Cplx> v(32), original;
    for (size_t i = 0; i < v.size(); ++i) v[i] = Cplx(std::sin(0.3f * i), std::cos(0.7f * i));
    original = v;
    fft1d(v, false); fft1d(v, true);
    for (size_t i = 0; i < v.size(); ++i) CHECK(std::abs(v[i] - original[i]) < 1e-4f);
}

TEST_CASE("forward then inverse is the identity in 2d") {
    const int w = 16, h = 8;
    std::vector<Cplx> v(size_t(w) * h), original;
    for (size_t i = 0; i < v.size(); ++i) v[i] = Cplx(float(i % 7) - 3.0f, float(i % 5) - 2.0f);
    original = v;
    fft2d(v, w, h, false); fft2d(v, w, h, true);
    for (size_t i = 0; i < v.size(); ++i) CHECK(std::abs(v[i] - original[i]) < 1e-3f);
}

TEST_CASE("Parseval holds: energy is preserved up to the size factor") {
    const int w = 16, h = 16;
    std::vector<Cplx> v(size_t(w) * h);
    for (size_t i = 0; i < v.size(); ++i) v[i] = Cplx(float((i * 37) % 11) / 11.0f, 0.0f);
    double spatial = 0.0;
    for (const Cplx& c : v) spatial += std::norm(c);
    fft2d(v, w, h, false);
    double freq = 0.0;
    for (const Cplx& c : v) freq += std::norm(c);
    CHECK(freq == doctest::Approx(spatial * w * h).epsilon(1e-4));
}

TEST_CASE("fftShift moves DC to the centre and is its own inverse for even sizes") {
    const int w = 8, h = 8;
    std::vector<Cplx> v(size_t(w) * h, Cplx(0, 0));
    v[0] = Cplx(1, 0);
    fftShift2d(v, w, h);
    CHECK(std::abs(v[size_t(h / 2) * w + w / 2]) == doctest::Approx(1.0f));
    fftShift2d(v, w, h);
    CHECK(std::abs(v[0]) == doctest::Approx(1.0f));
}

TEST_CASE("fft1d rejects a non-power-of-two length") {
    std::vector<Cplx> v(48, Cplx(0, 0));
    CHECK_THROWS_AS(fft1d(v, false), std::invalid_argument);
}

TEST_CASE("fft2d rejects a non-power-of-two width or height") {
    std::vector<Cplx> v(size_t(48) * 16, Cplx(0, 0));
    CHECK_THROWS_AS(fft2d(v, 48, 16, false), std::invalid_argument);

    std::vector<Cplx> v2(size_t(16) * 48, Cplx(0, 0));
    CHECK_THROWS_AS(fft2d(v2, 16, 48, false), std::invalid_argument);
}

TEST_CASE("fft2d rejects a buffer whose size does not match w*h") {
    std::vector<Cplx> v(size_t(16) * 16 - 1, Cplx(0, 0));
    CHECK_THROWS_AS(fft2d(v, 16, 16, false), std::invalid_argument);
}

TEST_CASE("fftShift2d rejects odd dimensions") {
    std::vector<Cplx> v(size_t(9) * 8, Cplx(0, 0));
    CHECK_THROWS_AS(fftShift2d(v, 9, 8), std::invalid_argument);

    std::vector<Cplx> v2(size_t(8) * 9, Cplx(0, 0));
    CHECK_THROWS_AS(fftShift2d(v2, 8, 9), std::invalid_argument);
}
