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
