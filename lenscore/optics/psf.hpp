#pragma once
#include "lenscore/conv/fft.hpp"
#include "lenscore/optics/pupil.hpp"
#include "lenscore/optics/zernike.hpp"
#include <cmath>

namespace lens::optics {

// Radius of the first Airy minimum, in output samples.
inline float airyFirstZeroSamples(int N, float apertureRadius) {
    const float Rs = apertureRadius * float(N) * 0.5f;
    return 1.22f * float(N) / (2.0f * Rs);
}

// Image-plane size of one output sample, in micrometres. The working stop
// cancels out: stopping down widens the pattern, it does not resample it.
inline float psfSampleSpacingUm(float lambdaNm, float fNumberWide) {
    return lambdaNm * 1e-3f * fNumberWide;
}

// PSF = |FFT{ A * exp(i 2 pi W lambdaRef/lambda) }|^2, fftshifted, unnormalised.
inline Plane psfFromPupil(const PupilParams& pp, const Wavefront& wf,
                          float lambdaNm, float lambdaRefNm, float t, int N) {
    const Plane amp = rasterPupil(pp, t, N);
    const float chroma = lambdaRefNm / lambdaNm;

    std::vector<conv::Cplx> field(size_t(N) * N, conv::Cplx(0.0f, 0.0f));
    for (int j = 0; j < N; ++j) {
        const float v = 2.0f * float(j) / float(N - 1) - 1.0f;
        for (int i = 0; i < N; ++i) {
            const float a = amp.at(i, j);
            if (a <= 0.0f) continue;
            const float u = 2.0f * float(i) / float(N - 1) - 1.0f;
            const float rho = std::sqrt(u * u + v * v) / pp.apertureRadius;
            const float th  = std::atan2(v, u);
            const float phase = kTwoPi * wavefrontError(wf, rho, th) * chroma;
            field[size_t(j) * N + i] = conv::Cplx(a * std::cos(phase), a * std::sin(phase));
        }
    }

    conv::fft2d(field, N, N, false);
    conv::fftShift2d(field, N, N);

    Plane out(N, N);
    for (size_t k = 0; k < field.size(); ++k) out.v[k] = std::norm(field[k]);
    return out;   // deliberately not normalised: the energy is the vignetting
}

}  // namespace lens::optics
