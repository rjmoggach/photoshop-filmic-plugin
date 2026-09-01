#pragma once
#include "lenscore/constants.hpp"
#include <cmath>
#include <complex>
#include <stdexcept>
#include <utility>
#include <vector>

// Radix-2 Cooley-Tukey FFT, in place, 1D and 2D, plus a quadrant shift.
//
// Every size this library ever hands in is a power of two (PSF grids and
// EFF patch sizes are both chosen that way), so radix-2 is sufficient.
// `assert` alone cannot enforce that: NDEBUG strips it from the shipping
// build, and a non-power-of-two length would then run the bit-reversal and
// butterfly stages over the wrong size, silently producing a plausible but
// wrong transform. Every optics result downstream (Task 13's PSF, Task 14's
// convolution) would then be quietly incorrect. So the size and shape
// invariants below are enforced with a real, always-present throw.

namespace lens::conv {

using Cplx = std::complex<float>;

// True for 1, 2, 4, 8, ... False for 0, negative numbers, and anything else.
inline bool isPowerOfTwo(int n) { return n > 0 && (n & (n - 1)) == 0; }

// In-place radix-2 FFT. `inverse` selects the +i twiddle direction and
// applies the 1/N scaling once, so that fft1d(forward) then fft1d(inverse)
// is the identity. a.size() must be a power of two (1 counts, and is a
// no-op transform); anything else throws rather than running UB-adjacent
// bit-reversal math over a size the algorithm was not designed for.
inline void fft1d(std::vector<Cplx>& a, bool inverse) {
    const int n = int(a.size());
    if (!isPowerOfTwo(n)) {
        throw std::invalid_argument("fft1d: length must be a power of two");
    }

    for (int i = 1, j = 0; i < n; ++i) {          // bit-reversal permutation
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        const double ang = 2.0 * double(kPi) / len * (inverse ? 1.0 : -1.0);
        const Cplx wl(float(std::cos(ang)), float(std::sin(ang)));
        for (int i = 0; i < n; i += len) {
            Cplx w(1.0f, 0.0f);
            for (int k = 0; k < len / 2; ++k) {
                const Cplx u = a[i + k];
                const Cplx v = a[i + k + len / 2] * w;
                a[i + k]           = u + v;
                a[i + k + len / 2] = u - v;
                w *= wl;
            }
        }
    }
    if (inverse) for (Cplx& c : a) c /= float(n);
}

// In-place 2D FFT: rows then columns, each a radix-2 fft1d. w and h must
// both be powers of two, and a.size() must equal w*h exactly — that check
// happens before row/col are sized from w and h, so those allocations are
// always bounded by a size the caller already committed to (a.size()),
// never blindly trusted from w/h alone.
inline void fft2d(std::vector<Cplx>& a, int w, int h, bool inverse) {
    if (!isPowerOfTwo(w) || !isPowerOfTwo(h)) {
        throw std::invalid_argument("fft2d: width and height must be powers of two");
    }
    if (a.size() != size_t(w) * size_t(h)) {
        throw std::invalid_argument("fft2d: buffer size must equal w*h");
    }

    std::vector<Cplx> row(w), col(h);
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) row[x] = a[size_t(y) * w + x];
        fft1d(row, inverse);
        for (int x = 0; x < w; ++x) a[size_t(y) * w + x] = row[x];
    }
    for (int x = 0; x < w; ++x) {
        for (int y = 0; y < h; ++y) col[y] = a[size_t(y) * w + x];
        fft1d(col, inverse);
        for (int y = 0; y < h; ++y) a[size_t(y) * w + x] = col[y];
    }
}

// Swaps quadrants so DC (index 0) lands at the centre. Built from disjoint
// pairwise swaps, so it is always its own inverse: applying it twice
// restores the original layout. w and h need not be powers of two, but both
// must be even — for an odd height the last row would never be touched by
// a pairwise swap, so DC would not land at the true centre, silently. That
// is a programmer error, not a data problem, so it throws rather than
// returning an asymmetric shift. a.size() must also equal w*h. Both checks
// run before any indexing, so a bad call throws instead of reading or
// writing out of bounds.
inline void fftShift2d(std::vector<Cplx>& a, int w, int h) {
    if (w <= 0 || h <= 0 || (w % 2) != 0 || (h % 2) != 0 ||
        a.size() != size_t(w) * size_t(h)) {
        throw std::invalid_argument("fftShift2d: w and h must be even and match buffer size");
    }

    const int hx = w / 2, hy = h / 2;
    for (int y = 0; y < hy; ++y)
        for (int x = 0; x < w; ++x) {
            const int sx = (x + hx) % w;
            std::swap(a[size_t(y) * w + x], a[size_t(y + hy) * w + sx]);
        }
}

}  // namespace lens::conv
