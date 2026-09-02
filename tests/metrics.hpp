#pragma once
#include "lenscore/conv/fft.hpp"
#include "lenscore/geometry.hpp"
#include "lenscore/image.hpp"
#include "lenscore/plane.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace lens::metrics {

inline Plane crop(const Plane& p, int x0, int y0, int w, int h) {
    Plane o(w, h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x)
            o.at(x, y) = p.at(std::clamp(x0 + x, 0, p.w - 1), std::clamp(y0 + y, 0, p.h - 1));
    return o;
}

// fitEdge/mtf50 scan ROWS and locate the edge along x -- they assume a
// near-VERTICAL edge (an angle near 0 degrees off vertical, in slantedEdge's
// convention). A near-horizontal edge (angle near 90) is not a rotated
// version of the same measurement to that code; transpose swaps x and y so a
// near-horizontal edge becomes near-vertical before mtf50 sees it, making the
// metric apply. Callers measuring resolution across a horizontal edge must
// go through this, not feed the untransposed ROI straight to mtf50.
inline Plane transpose(const Plane& p) {
    Plane o(p.h, p.w);
    for (int y = 0; y < p.h; ++y)
        for (int x = 0; x < p.w; ++x)
            o.at(y, x) = p.at(x, y);
    return o;
}

inline Plane gaussianBlur(const Plane& src, float sigma) {
    const int r = std::max(1, int(std::ceil(4.0f * sigma)));
    std::vector<float> k(size_t(2 * r + 1));
    float sum = 0.0f;
    for (int i = -r; i <= r; ++i) { k[i + r] = std::exp(-0.5f * i * i / (sigma * sigma)); sum += k[i + r]; }
    for (float& v : k) v /= sum;

    Plane tmp(src.w, src.h), out(src.w, src.h);
    for (int y = 0; y < src.h; ++y)
        for (int x = 0; x < src.w; ++x) {
            float a = 0.0f;
            for (int i = -r; i <= r; ++i) a += k[i + r] * src.at(std::clamp(x + i, 0, src.w - 1), y);
            tmp.at(x, y) = a;
        }
    for (int y = 0; y < src.h; ++y)
        for (int x = 0; x < src.w; ++x) {
            float a = 0.0f;
            for (int i = -r; i <= r; ++i) a += k[i + r] * tmp.at(x, std::clamp(y + i, 0, src.h - 1));
            out.at(x, y) = a;
        }
    return out;
}

// Fits the near-vertical edge; returns slope (px per row) and intercept at mid height.
inline void fitEdge(const Plane& roi, float& slope, float& intercept) {
    std::vector<float> cy(size_t(roi.h)), yy(size_t(roi.h));
    int n = 0;
    for (int y = 0; y < roi.h; ++y) {
        double num = 0.0, den = 0.0;
        for (int x = 0; x + 1 < roi.w; ++x) {
            const double d = std::abs(double(roi.at(x + 1, y)) - roi.at(x, y));
            num += d * (x + 0.5); den += d;
        }
        if (den > 1e-6) { cy[n] = float(num / den); yy[n] = float(y) - 0.5f * float(roi.h); ++n; }
    }
    double sx = 0, sy = 0, sxy = 0, sxx = 0;
    for (int i = 0; i < n; ++i) { sx += yy[i]; sy += cy[i]; sxy += yy[i] * cy[i]; sxx += yy[i] * yy[i]; }
    const double det = n * sxx - sx * sx;
    slope     = (std::abs(det) > 1e-9) ? float((n * sxy - sx * sy) / det) : 0.0f;
    intercept = float((sy - slope * sx) / std::max(1, n));
}

inline float edgePosition(const Plane& roi) {
    float s = 0, b = 0; fitEdge(roi, s, b); return b;
}

inline float mtf50(const Plane& roi) {
    float slope = 0, intercept = 0;
    fitEdge(roi, slope, intercept);

    // 4x oversampled edge spread function, projected onto the edge normal.
    const int OS = 4, NB = 256;
    std::vector<double> acc(NB, 0.0), cnt(NB, 0.0);
    for (int y = 0; y < roi.h; ++y) {
        const float ex = intercept + slope * (float(y) - 0.5f * float(roi.h));
        for (int x = 0; x < roi.w; ++x) {
            const float u = float(x) + 0.5f - ex;
            const int b = int(std::lround(u * OS)) + NB / 2;
            if (b >= 0 && b < NB) { acc[b] += roi.at(x, y); cnt[b] += 1.0; }
        }
    }
    std::vector<double> esf(NB, 0.0);
    // Seed the carry-forward value from the first real sample, not 0.0, so
    // any leading bins with no projected pixel are backfilled with the
    // edge's actual dark level instead of being pinned to zero. With a zero
    // dark level (every test target in this file) the two coincide, which
    // is why this only shows up once slantedEdge's lo != 0.
    double last = 0.0;
    for (int i = 0; i < NB; ++i) if (cnt[i] > 0) { last = acc[i] / cnt[i]; break; }
    for (int i = 0; i < NB; ++i) { if (cnt[i] > 0) last = acc[i] / cnt[i]; esf[i] = last; }

    // Line spread function, Hamming windowed to suppress ringing.
    std::vector<conv::Cplx> lsf(NB, conv::Cplx(0, 0));
    for (int i = 1; i < NB; ++i) {
        const double w = 0.54 - 0.46 * std::cos(2.0 * kPi * i / (NB - 1));
        lsf[i] = conv::Cplx(float((esf[i] - esf[i - 1]) * w), 0.0f);
    }
    conv::fft1d(lsf, false);

    const float dc = std::abs(lsf[0]);
    if (dc <= 0.0f) return 0.0f;
    // Bin width is 1/OS pixels, so bin k sits at k*OS/NB cycles per pixel.
    // The analytic-gaussian test's 8% tolerance catches a gross axis error
    // (e.g. k/NB instead of k*OS/NB, a factor of OS=4) loudly, but an
    // off-by-one of order (NB-1)/NB ~ 0.4% sits comfortably inside it. The
    // gate proves the axis is right to within a few percent, not exactly.
    float prevF = 0.0f, prevM = 1.0f;
    for (int k = 1; k < NB / 2; ++k) {
        const float f = float(k) * OS / float(NB);
        const float m = std::abs(lsf[k]) / dc;
        if (m <= 0.5f) return prevF + (prevM - 0.5f) / std::max(1e-6f, prevM - m) * (f - prevF);
        prevF = f; prevM = m;
    }
    return float(OS) * 0.5f;
}

// Bins by normalised field radius (frameOf: corner = 1.0), not pixels. The
// outermost bin is fed only by the four corner wedges of a rectangular
// image -- the inscribed circle of radius 1.0 has no pixels beyond it along
// the edge midpoints -- so it is a directionally biased sample, not a true
// azimuthal ring average. Acceptance-suite callers should treat the last
// bin as approximate.
inline std::vector<float> radialMean(const Plane& p, int bins) {
    const Frame f = frameOf(p.w, p.h);
    std::vector<double> acc(size_t(bins), 0.0), cnt(size_t(bins), 0.0);
    for (int y = 0; y < p.h; ++y)
        for (int x = 0; x < p.w; ++x) {
            const float dx = (float(x) - f.cx) / f.halfDiag, dy = (float(y) - f.cy) / f.halfDiag;
            const int b = std::min(bins - 1, int(std::sqrt(dx * dx + dy * dy) * bins));
            acc[b] += p.at(x, y); cnt[b] += 1.0;
        }
    std::vector<float> out(size_t(bins), 0.0f);
    for (int i = 0; i < bins; ++i) out[i] = cnt[i] > 0 ? float(acc[i] / cnt[i]) : 0.0f;
    return out;
}

inline Plane channel(const Image& im, int c) {
    Plane p(im.w, im.h);
    for (int y = 0; y < im.h; ++y) for (int x = 0; x < im.w; ++x) p.at(x, y) = im.at(x, y, c);
    return p;
}

inline float fringeWidthPx(const Image& roi) {
    return edgePosition(channel(roi, 0)) - edgePosition(channel(roi, 2));
}

inline double totalEnergy(const Plane& p) { double s = 0; for (float v : p.v) s += v; return s; }
inline double totalEnergy(const Image& im) { double s = 0; for (float v : im.px) s += v; return s; }

// Rotational (90-degree) asymmetry, normalised by peak magnitude. A
// non-square image has no well-defined 90-degree rotation, so that is a
// programmer error rather than a measurement -- it throws instead of
// returning a sentinel a caller could mistake for a real (extreme)
// asymmetry value.
inline float rot90Asymmetry(const Plane& p) {
    if (p.w != p.h) throw std::invalid_argument("rot90Asymmetry: image must be square");
    const int n = p.w;
    float peak = 0.0f, worst = 0.0f;
    for (float v : p.v) peak = std::max(peak, std::abs(v));
    if (peak <= 0.0f) return 0.0f;
    for (int y = 0; y < n; ++y)
        for (int x = 0; x < n; ++x)
            worst = std::max(worst, std::abs(p.at(x, y) - p.at(n - 1 - y, x)));
    return worst / peak;
}

}  // namespace lens::metrics
