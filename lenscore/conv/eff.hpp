#pragma once
#include "lenscore/conv/fft.hpp"
#include "lenscore/geometry.hpp"
#include "lenscore/plane.hpp"
#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

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
// Patches overlap by 50%, so their accumulations into `out` collide and cannot
// simply be threaded. Each patch's FFT work is independent though, and it is
// where essentially all the time goes -- measured, this stage costs 4.9 s per
// megapixel against 0.11 s for everything else in the pipeline combined. So the
// per-patch transforms run in parallel into private buffers and only the cheap
// accumulation is serialised. Peak extra memory is one S*S complex buffer per
// thread, not one per patch and not a copy of the image.
inline Plane effConvolve(const Plane& src, int patch,
                         const std::function<Plane(float, float)>& psfAt) {
    const int hop = patch / 2;
    const Plane probe = psfAt(src.w * 0.5f, src.h * 0.5f);
    const int S = nextPow2(patch + std::max(probe.w, probe.h));

    const std::vector<float> win = hannWindow(patch);
    Plane out(src.w, src.h);

    struct Origin { int ox, oy; };
    std::vector<Origin> origins;
    for (int oy = -hop; oy < src.h; oy += hop)
        for (int ox = -hop; ox < src.w; ox += hop) origins.push_back({ox, oy});

    // Everything one worker needs, reused across the chunks it processes.
    struct Scratch {
        std::vector<Cplx> buf, ker;
        int kw2 = 0, kh2 = 0;
        std::exception_ptr err;
    };

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    if (origins.size() < threads) threads = unsigned(origins.size());
    if (threads == 0) return out;

    // One scratch per in-flight patch: threads * kBatch of them (see the batch
    // loop below). Each holds two S*S complex buffers, so this is kilobytes per
    // slot, not a copy of the image.
    constexpr size_t kBatchSlots = 4;
    std::vector<Scratch> scratch(size_t(threads) * kBatchSlots);
    for (auto& sc : scratch) { sc.buf.resize(size_t(S) * S); sc.ker.resize(size_t(S) * S); }

    // One patch: window it, build its kernel, and convolve. Writes only into its
    // own scratch, so any number of these may run at once.
    auto convolveOne = [&](const Origin& o, Scratch& sc) {
        const Plane k = psfAt(float(o.ox) + patch * 0.5f, float(o.oy) + patch * 0.5f);
        if (k.w > probe.w || k.h > probe.h) {
            throw std::invalid_argument(
                "effConvolve: psfAt returned a " + std::to_string(k.w) + "x" +
                std::to_string(k.h) + " kernel larger than the " + std::to_string(probe.w) +
                "x" + std::to_string(probe.h) + " probe kernel the FFT size was sized from");
        }
        sc.kh2 = k.h / 2;
        sc.kw2 = k.w / 2;

        std::fill(sc.buf.begin(), sc.buf.end(), Cplx(0.0f, 0.0f));
        for (int j = 0; j < patch; ++j)
            for (int i = 0; i < patch; ++i) {
                // Mirrored rather than clamped. Patches start one hop outside the
                // image so every pixel is covered by full windows, so the gather
                // genuinely reads past the border; clamping replicates the edge
                // pixel and leaves a stripe down each side.
                const int sx = lens::mirrorIndex(o.ox + i, src.w);
                const int sy = lens::mirrorIndex(o.oy + j, src.h);
                sc.buf[size_t(j) * S + i] = Cplx(src.at(sx, sy) * win[i] * win[j], 0.0f);
            }

        // Kernel centred on the wraparound origin, so convolution does not shift.
        std::fill(sc.ker.begin(), sc.ker.end(), Cplx(0.0f, 0.0f));
        for (int j = 0; j < k.h; ++j)
            for (int i = 0; i < k.w; ++i) {
                const int wx = ((i - sc.kw2) % S + S) % S;
                const int wy = ((j - sc.kh2) % S + S) % S;
                sc.ker[size_t(wy) * S + wx] += Cplx(k.at(i, j), 0.0f);
            }

        fft2d(sc.buf, S, S, false);
        fft2d(sc.ker, S, S, false);
        for (size_t i = 0; i < sc.buf.size(); ++i) sc.buf[i] *= sc.ker[i];
        fft2d(sc.buf, S, S, true);
    };

    // The kernel was centred on the wraparound origin, so circular convolution
    // parks the negative lags at the TOP of the FFT output. Indices at or above
    // S - k?2 are those negative lags and must map back to negative offsets from
    // (ox, oy), not to the tail of the array. Since S >= patch + K, the wrapped
    // region never overlaps the real patch data, so this mapping is unambiguous.
    auto accumulate = [&](const Origin& o, const Scratch& sc) {
        for (int j = 0; j < S; ++j) {
            const int jj = (j >= S - sc.kh2) ? j - S : j;
            const int dy = o.oy + jj;
            if (dy < 0 || dy >= src.h) continue;
            for (int i = 0; i < S; ++i) {
                const int ii = (i >= S - sc.kw2) ? i - S : i;
                const int dx = o.ox + ii;
                if (dx < 0 || dx >= src.w) continue;
                out.at(dx, dy) += sc.buf[size_t(j) * S + i].real();
            }
        }
    };

    if (threads == 1) {
        for (const auto& o : origins) { convolveOne(o, scratch[0]); accumulate(o, scratch[0]); }
        return out;
    }

    // Each round has to join before the serial accumulation, and spawning a
    // thread costs about as much as one patch's transforms. Give every thread a
    // batch of patches per round so that barrier is paid once per batch rather
    // than once per patch.
    constexpr size_t kBatch = kBatchSlots;
    const size_t slots = size_t(threads) * kBatch;

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (size_t base = 0; base < origins.size(); base += slots) {
        const size_t n = std::min(slots, origins.size() - base);
        for (size_t t = 0; t < threads; ++t) {
            const size_t first = t * kBatch;
            if (first >= n) break;
            const size_t last = std::min(first + kBatch, n);
            for (size_t i = first; i < last; ++i) scratch[i].err = nullptr;
            pool.emplace_back([&, first, last, base] {
                // An exception must not escape a std::thread -- that calls
                // std::terminate. Carry it back and rethrow after the join.
                for (size_t i = first; i < last; ++i) {
                    try { convolveOne(origins[base + i], scratch[i]); }
                    catch (...) { scratch[i].err = std::current_exception(); }
                }
            });
        }
        for (auto& th : pool) th.join();
        pool.clear();
        for (size_t i = 0; i < n; ++i) if (scratch[i].err) std::rethrow_exception(scratch[i].err);
        for (size_t i = 0; i < n; ++i) accumulate(origins[base + i], scratch[i]);
    }
    return out;
}

}  // namespace lens::conv
