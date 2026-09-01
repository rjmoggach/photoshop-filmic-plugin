#include <doctest/doctest.h>
#include "lenscore/optics/dispersion.hpp"

using namespace lens::optics;

TEST_CASE("BK7 hits its published index at the d line") {
    CHECK(refractiveIndex(bk7(), 587.6f) == doctest::Approx(1.5168f).epsilon(1e-3).scale(0));
}

TEST_CASE("dispersion is normal, index falls as wavelength rises") {
    const Dispersion d = bk7();
    float prev = refractiveIndex(d, 400.0f);
    for (float l = 420.0f; l <= 760.0f; l += 20.0f) {
        const float n = refractiveIndex(d, l);
        CHECK(n < prev);
        prev = n;
    }
}

TEST_CASE("an uncorrected singlet has monotonic focus error, rising with wavelength") {
    Dispersion d = bk7();
    d.correction_nm.clear();
    CHECK(focusError(d, 400.0f, 650.0f) < 0.0f);   // blue focuses closer
    CHECK(focusError(d, 760.0f, 650.0f) > 0.0f);   // red focuses further
    float prev = focusError(d, 400.0f, 650.0f);
    for (float l = 420.0f; l <= 760.0f; l += 20.0f) {
        const float e = focusError(d, l, 650.0f);
        CHECK(e > prev);
        prev = e;
    }
}

TEST_CASE("an achromat has exactly zero focus error at both correction lines") {
    Dispersion d = bk7();
    d.correction_nm = {486.1f, 656.3f};      // F and C lines
    CHECK(focusError(d, 486.1f, 650.0f) == doctest::Approx(0.0f).epsilon(1e-6));
    CHECK(focusError(d, 656.3f, 650.0f) == doctest::Approx(0.0f).epsilon(1e-6));
}

TEST_CASE("the achromat's secondary spectrum bulges with one sign between the lines") {
    Dispersion d = bk7();
    d.correction_nm = {486.1f, 656.3f};
    const float mid = focusError(d, 570.0f, 650.0f);
    CHECK(mid != doctest::Approx(0.0f).epsilon(1e-5));
    for (float l = 500.0f; l <= 640.0f; l += 20.0f)
        CHECK(focusError(d, l, 650.0f) * mid > 0.0f);   // no sign flip in between
}

TEST_CASE("an apochromat zeroes three lines") {
    Dispersion d = bk7();
    d.correction_nm = {436.0f, 546.1f, 656.3f};
    for (float l : d.correction_nm)
        CHECK(focusError(d, l, 650.0f) == doctest::Approx(0.0f).epsilon(1e-6));
}

TEST_CASE("residual scales the remaining error linearly") {
    Dispersion a = bk7(); a.correction_nm = {486.1f, 656.3f}; a.residual = 1.0f;
    Dispersion b = a;                                        b.residual = 0.25f;
    CHECK(focusError(b, 570.0f, 650.0f) == doctest::Approx(0.25f * focusError(a, 570.0f, 650.0f)));
}

TEST_CASE("a duplicated correction wavelength does not divide by zero") {
    Dispersion d = bk7();
    d.correction_nm = {486.1f, 486.1f};   // degenerate: same line twice
    for (float l : {400.0f, 486.1f, 550.0f, 656.3f, 760.0f})
        CHECK(std::isfinite(focusError(d, l, 650.0f)));
}

// A single index value (the d-line test above) only pins the curve's height
// at one point -- a wrong B/C pair that happens to land near 1.5168 there
// while distorting the curve's shape elsewhere would still pass. The Abbe
// number is built from three points (d, F, C), so it pins the shape of the
// dispersion curve, not just one height on it.
TEST_CASE("BK7 reproduces its published Abbe number, pinning the whole curve shape") {
    const Dispersion d = bk7();
    const float nd = refractiveIndex(d, 587.56f);
    const float nF = refractiveIndex(d, 486.13f);
    const float nC = refractiveIndex(d, 656.27f);
    CHECK(nd == doctest::Approx(1.51680f).epsilon(1e-4).scale(0));
    CHECK(nF == doctest::Approx(1.52238f).epsilon(1e-4).scale(0));
    CHECK(nC == doctest::Approx(1.51432f).epsilon(1e-4).scale(0));
    CHECK((nd - 1.0f) / (nF - nC) == doctest::Approx(64.17f).epsilon(2e-3).scale(0));
}
