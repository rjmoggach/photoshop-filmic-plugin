#pragma once
#include "lenscore/conv/fft.hpp"
#include "lenscore/plane.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>

namespace lens::conv {

// Periodic Hann: at 50% overlap the shifted copies sum to exactly 1.
inline std::vector<float> hannWindow(int n) {
    std::vector<float> w(n);
    for (int i = 0; i < n; ++i) w[i] = 0.5f * (1.0f - std::cos(2.0f * kPi * float(i) / float(n)));
    return w;
}

inline Plane convolveDirect(const Plane& src, const Plane& k) {
    Plane out(src.w, src.h);
    const int kx = k.w / 2, ky = k.h / 2;
    for (int y = 0; y < src.h; ++y)
        for (int x = 0; x < src.w; ++x) {
            float acc = 0.0f;
            for (int j = 0; j < k.h; ++j)
                for (int i = 0; i < k.w; ++i) {
                    const int sx = std::clamp(x + i - kx, 0, src.w - 1);
                    const int sy = std::clamp(y + j - ky, 0, src.h - 1);
                    acc += src.at(sx, sy) * k.at(i, j);
                }
            out.at(x, y) = acc;
        }
    return out;
}

inline int nextPow2(int n) { int p = 1; while (p < n) p <<= 1; return p; }

// Spatially varying convolution ("Efficient Filter Flow"): the image is cut
// into overlapping square patches, each windowed by a separable Hann so the
// 50%-overlap patches blend back together without a seam, convolved with the
// kernel that applies at that patch's centre via FFT, and accumulated.
//
// Contract on psfAt: the FFT size S is fixed once, from a single probe call
// at the image centre, before any patch is processed. Every kernel psfAt
// returns for every other patch must therefore be no larger, in either
// width or height, than the probe kernel -- equal or smaller is fine, since
// that only grows the unused wraparound margin S was sized with. A larger
// kernel would make that patch's own wraparound region overlap real patch
// data; the modulo arithmetic below keeps every index in bounds regardless,
// so a violation would not crash or trip an assert, it would just silently
// add wrong values near that patch. That is checked and thrown on below
// rather than trusted.
inline Plane effConvolve(const Plane& src, int patch,
                         const std::function<Plane(float, float)>& psfAt) {
    const int hop = patch / 2;
    const Plane probe = psfAt(src.w * 0.5f, src.h * 0.5f);
    const int S = nextPow2(patch + std::max(probe.w, probe.h));

    const std::vector<float> win = hannWindow(patch);
    Plane out(src.w, src.h);
    std::vector<Cplx> buf(size_t(S) * S), ker(size_t(S) * S);

    // Origins start one hop before the image so every pixel is covered by full windows.
    for (int oy = -hop; oy < src.h; oy += hop) {
        for (int ox = -hop; ox < src.w; ox += hop) {
            const Plane k = psfAt(float(ox) + patch * 0.5f, float(oy) + patch * 0.5f);
            if (k.w > probe.w || k.h > probe.h) {
                throw std::invalid_argument(
                    "effConvolve: psfAt returned a " + std::to_string(k.w) + "x" +
                    std::to_string(k.h) + " kernel larger than the " + std::to_string(probe.w) +
                    "x" + std::to_string(probe.h) + " probe kernel the FFT size was sized from");
            }
            const int kh2 = k.h / 2, kw2 = k.w / 2;

            std::fill(buf.begin(), buf.end(), Cplx(0.0f, 0.0f));
            for (int j = 0; j < patch; ++j)
                for (int i = 0; i < patch; ++i) {
                    const int sx = std::clamp(ox + i, 0, src.w - 1);
                    const int sy = std::clamp(oy + j, 0, src.h - 1);
                    buf[size_t(j) * S + i] = Cplx(src.at(sx, sy) * win[i] * win[j], 0.0f);
                }

            // Kernel centred on the wraparound origin, so convolution does not shift.
            std::fill(ker.begin(), ker.end(), Cplx(0.0f, 0.0f));
            for (int j = 0; j < k.h; ++j)
                for (int i = 0; i < k.w; ++i) {
                    const int wx = ((i - kw2) % S + S) % S;
                    const int wy = ((j - kh2) % S + S) % S;
                    ker[size_t(wy) * S + wx] += Cplx(k.at(i, j), 0.0f);
                }

            fft2d(buf, S, S, false);
            fft2d(ker, S, S, false);
            for (size_t i = 0; i < buf.size(); ++i) buf[i] *= ker[i];
            fft2d(buf, S, S, true);

            // The kernel was centred on the wraparound origin, so circular
            // convolution parks the negative lags at the TOP of the FFT
            // output. Indices at or above S - k?2 are those negative lags
            // and must map back to negative offsets from (ox, oy), not to
            // the tail of the array. Since S >= patch + K, the wrapped
            // region never overlaps the real patch data, so this mapping
            // is unambiguous.
            for (int j = 0; j < S; ++j) {
                const int jj = (j >= S - kh2) ? j - S : j;
                const int dy = oy + jj;
                if (dy < 0 || dy >= src.h) continue;
                for (int i = 0; i < S; ++i) {
                    const int ii = (i >= S - kw2) ? i - S : i;
                    const int dx = ox + ii;
                    if (dx < 0 || dx >= src.w) continue;
                    out.at(dx, dy) += buf[size_t(j) * S + i].real();
                }
            }
        }
    }
    return out;
}

}  // namespace lens::conv
