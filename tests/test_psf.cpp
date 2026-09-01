#include <doctest/doctest.h>
#include "lenscore/optics/psf.hpp"
#include <cmath>
#include <vector>

using namespace lens;
using namespace lens::optics;

static PupilParams smallDisc() {
    PupilParams p;
    p.blades = 0; p.curvature = 1.0f; p.apertureRadius = 0.125f;
    p.rEntrance = 1e6f; p.rExit = 1e6f; p.sepNorm = 0.0f; p.apodizationSlope = 0.0f;
    return p;
}

// Mean intensity on a ring of the given radius about the PSF centre.
static double ringMean(const Plane& p, double r) {
    const double cx = p.w / 2.0, cy = p.h / 2.0;
    double acc = 0.0; int n = 0;
    for (int k = 0; k < 720; ++k) {
        const double a = 2.0 * kPi * k / 720.0;
        const int x = int(std::lround(cx + r * std::cos(a)));
        const int y = int(std::lround(cy + r * std::sin(a)));
        if (x >= 0 && y >= 0 && x < p.w && y < p.h) { acc += p.at(x, y); ++n; }
    }
    return n ? acc / n : 0.0;
}

TEST_CASE("sample spacing depends on the wide-open stop, not the working stop") {
    CHECK(psfSampleSpacingUm(550.0f, 2.0f) == doctest::Approx(1.1f).epsilon(1e-4));
    CHECK(psfSampleSpacingUm(550.0f, 4.0f) == doctest::Approx(2.2f).epsilon(1e-4));
}

TEST_CASE("an unaberrated PSF peaks dead centre") {
    const int N = 256;
    const Plane psf = psfFromPupil(smallDisc(), Wavefront{}, 550.0f, 550.0f, 0.0f, N);
    float best = -1.0f; int bx = 0, by = 0;
    for (int y = 0; y < N; ++y) for (int x = 0; x < N; ++x)
        if (psf.at(x, y) > best) { best = psf.at(x, y); bx = x; by = y; }
    CHECK(bx == N / 2);
    CHECK(by == N / 2);
}

TEST_CASE("the first Airy zero lands where diffraction says it should") {
    const int N = 256;
    const PupilParams p = smallDisc();
    const Plane psf = psfFromPupil(p, Wavefront{}, 550.0f, 550.0f, 0.0f, N);
    const double predicted = airyFirstZeroSamples(N, p.apertureRadius);   // 9.76
    CHECK(predicted == doctest::Approx(9.76).epsilon(0.01));

    double bestR = 0.0, bestV = 1e30;
    for (double r = 3.0; r < 18.0; r += 0.25) {
        const double v = ringMean(psf, r);
        if (v < bestV) { bestV = v; bestR = r; }
    }
    CHECK(bestR == doctest::Approx(predicted).epsilon(0.12));
}

TEST_CASE("stopping down widens the diffraction pattern") {
    const int N = 256;
    PupilParams wide = smallDisc();
    PupilParams stop = smallDisc(); stop.apertureRadius = 0.0625f;
    CHECK(airyFirstZeroSamples(N, stop.apertureRadius) >
          airyFirstZeroSamples(N, wide.apertureRadius));
}

TEST_CASE("total PSF energy tracks the pupil energy, so vignetting survives") {
    const int N = 256;
    PupilParams p = smallDisc();
    p.rEntrance = 0.13f; p.rExit = 1e6f; p.sepNorm = 0.10f;    // clips off axis

    auto energy = [&](float t) {
        const Plane psf = psfFromPupil(p, Wavefront{}, 550.0f, 550.0f, t, N);
        double s = 0.0; for (float v : psf.v) s += v; return s;
    };
    auto pupilSq = [&](float t) {
        const Plane a = rasterPupil(p, t, N);
        double s = 0.0; for (float v : a.v) s += double(v) * v; return s;
    };
    CHECK(energy(1.0f) / energy(0.0f) == doctest::Approx(pupilSq(1.0f) / pupilSq(0.0f)).epsilon(1e-3));
    CHECK(energy(1.0f) < energy(0.0f));      // it really did lose light
}

TEST_CASE("defocus broadens the PSF") {
    const int N = 256;
    auto secondMoment = [&](const Wavefront& w) {
        const Plane psf = psfFromPupil(smallDisc(), w, 550.0f, 550.0f, 0.0f, N);
        double s = 0.0, m = 0.0;
        for (int y = 0; y < N; ++y) for (int x = 0; x < N; ++x) {
            const double dx = x - N / 2.0, dy = y - N / 2.0;
            s += psf.at(x, y); m += psf.at(x, y) * (dx * dx + dy * dy);
        }
        return m / s;
    };
    Wavefront d; d.defocus = 1.5f;
    CHECK(secondMoment(d) > secondMoment(Wavefront{}));
}

TEST_CASE("astigmatism makes the PSF elliptical") {
    const int N = 256;
    Wavefront w; w.astig = 1.2f;
    const Plane psf = psfFromPupil(smallDisc(), w, 550.0f, 550.0f, 0.0f, N);
    double s = 0, mx = 0, my = 0;
    for (int y = 0; y < N; ++y) for (int x = 0; x < N; ++x) {
        const double dx = x - N / 2.0, dy = y - N / 2.0;
        s += psf.at(x, y); mx += psf.at(x, y) * dx * dx; my += psf.at(x, y) * dy * dy;
    }
    CHECK(std::abs(mx / s - my / s) / (mx / s) > 0.10);
}

TEST_CASE("coma throws the PSF off centre along the radial axis") {
    const int N = 256;
    Wavefront w; w.coma = 1.2f;
    const Plane psf = psfFromPupil(smallDisc(), w, 550.0f, 550.0f, 1.0f, N);
    double s = 0, cx = 0, cy = 0;
    for (int y = 0; y < N; ++y) for (int x = 0; x < N; ++x) {
        s += psf.at(x, y);
        cx += psf.at(x, y) * (x - N / 2.0);
        cy += psf.at(x, y) * (y - N / 2.0);
    }
    CHECK(std::abs(cx / s) > 0.5);                     // shifted along +/-x
    CHECK(std::abs(cy / s) < 0.1 * std::abs(cx / s));  // not sideways
}

TEST_CASE("an even blade count gives that many diffraction spikes") {
    const int N = 512;
    PupilParams p = smallDisc();
    p.blades = 6; p.curvature = 0.0f;
    const Plane psf = psfFromPupil(p, Wavefront{}, 550.0f, 550.0f, 0.0f, N);

    std::vector<double> ang(360);
    for (int k = 0; k < 360; ++k) {
        const double a = 2.0 * kPi * k / 360.0;
        double acc = 0.0;
        for (double r = 14.0; r < 34.0; r += 1.0) {
            const int x = int(std::lround(N / 2.0 + r * std::cos(a)));
            const int y = int(std::lround(N / 2.0 + r * std::sin(a)));
            if (x >= 0 && y >= 0 && x < N && y < N) acc += psf.at(x, y);
        }
        ang[k] = acc;
    }
    double mean = 0.0; for (double v : ang) mean += v; mean /= 360.0;
    int peaks = 0;
    for (int k = 0; k < 360; ++k) {
        const double p0 = ang[(k + 359) % 360], p1 = ang[k], p2 = ang[(k + 1) % 360];
        if (p1 > p0 && p1 >= p2 && p1 > 1.5 * mean) ++peaks;
    }
    CHECK(peaks == 6);
}
