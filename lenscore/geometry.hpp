#pragma once
#include "lenscore/plane.hpp"
#include <algorithm>
#include <cmath>

namespace lens {

// Normalised field coordinates: the image corner sits at t = 1.
struct Frame { float cx = 0, cy = 0, halfDiag = 1; };

inline Frame frameOf(int w, int h) {
    Frame f;
    f.cx = 0.5f * float(w - 1);
    f.cy = 0.5f * float(h - 1);
    f.halfDiag = std::sqrt(f.cx * f.cx + f.cy * f.cy);
    if (f.halfDiag <= 0.0f) f.halfDiag = 1.0f;
    return f;
}

// Reflect a coordinate back inside [0, n-1].
//
// Clamping is the obvious thing and it is visibly wrong here: every sample taken
// outside the frame returns the SAME edge pixel, so a warp that reaches past the
// border smears that one pixel into a long streak. Distortion reaches past the
// border by design, which is why the streaks appeared down the left and right of
// a wide picture. Reflecting continues the texture instead, so the border keeps
// the statistics of the picture rather than inventing a stripe.
inline float mirrorCoord(float x, int n) {
    if (n <= 1) return 0.0f;
    const float period = 2.0f * float(n - 1);
    x = std::fmod(x, period);
    if (x < 0.0f) x += period;
    return (x <= float(n - 1)) ? x : period - x;
}

inline int mirrorIndex(int i, int n) {
    if (n <= 1) return 0;
    const int period = 2 * n - 2;
    i %= period;
    if (i < 0) i += period;
    return (i < n) ? i : period - i;
}

// Mirrored at the border. Used by the convolution's patch gather, where the
// window genuinely overhangs the image and either clamping (a stripe) or black
// (a dark rim the lens never produced) would be an artefact of the method rather
// than of the optics.
inline float sampleBilinearMirror(const Plane& p, float x, float y) {
    x = mirrorCoord(x, p.w);
    y = mirrorCoord(y, p.h);
    const int x0 = std::min(int(x), p.w - 1), y0 = std::min(int(y), p.h - 1);
    const int x1 = mirrorIndex(x0 + 1, p.w), y1 = mirrorIndex(y0 + 1, p.h);
    const float fx = x - float(x0), fy = y - float(y0);
    return (1 - fy) * ((1 - fx) * p.at(x0, y0) + fx * p.at(x1, y0))
         +      fy  * ((1 - fx) * p.at(x0, y1) + fx * p.at(x1, y1));
}

// Black outside the frame. A geometric warp that reaches past the border has
// genuinely run out of picture, and saying so plainly is more useful than
// inventing content: it shows exactly where the frame edge went. Clamping smears
// one edge pixel into a stripe; mirroring invents plausible texture that was
// never photographed.
inline float sampleBilinearBorder(const Plane& p, float x, float y) {
    if (x < -1.0f || y < -1.0f || x > float(p.w) || y > float(p.h)) return 0.0f;
    const int x0 = int(std::floor(x)), y0 = int(std::floor(y));
    const float fx = x - float(x0), fy = y - float(y0);
    auto tap = [&](int xi, int yi) {
        return (xi < 0 || yi < 0 || xi >= p.w || yi >= p.h) ? 0.0f : p.at(xi, yi);
    };
    return (1 - fy) * ((1 - fx) * tap(x0, y0) + fx * tap(x0 + 1, y0))
         +      fy  * ((1 - fx) * tap(x0, y0 + 1) + fx * tap(x0 + 1, y0 + 1));
}

inline float sampleBilinear(const Plane& p, float x, float y) {
    x = std::clamp(x, 0.0f, float(p.w - 1));
    y = std::clamp(y, 0.0f, float(p.h - 1));
    const int x0 = std::min(int(x), p.w - 1), y0 = std::min(int(y), p.h - 1);
    const int x1 = std::min(x0 + 1, p.w - 1), y1 = std::min(y0 + 1, p.h - 1);
    const float fx = x - float(x0), fy = y - float(y0);
    return (1 - fy) * ((1 - fx) * p.at(x0, y0) + fx * p.at(x1, y0))
         +      fy  * ((1 - fx) * p.at(x0, y1) + fx * p.at(x1, y1));
}

}  // namespace lens
