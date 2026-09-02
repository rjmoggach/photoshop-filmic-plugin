#pragma once
#include "lenscore/color/cie.hpp"
#include "lenscore/constants.hpp"
#include "lenscore/geometry.hpp"
#include "lenscore/image.hpp"
#include "lenscore/plane.hpp"
#include <algorithm>
#include <cmath>

namespace lens::targets {

inline Plane flatField(int w, int h, float value) {
    Plane p(w, h);
    std::fill(p.v.begin(), p.v.end(), value);
    return p;
}

inline Plane pointGrid(int w, int h, int spacing) {
    Plane p(w, h);
    for (int y = spacing; y < h; y += spacing)
        for (int x = spacing; x < w; x += spacing)
            if (x < w && y < h) p.at(x, y) = 1.0f;
    return p;
}

// Area-exact antialiasing: coverage of the pixel by the half-plane, computed
// analytically, so the generator adds no blur of its own.
inline Plane slantedEdge(int w, int h, float angleDeg, float lo, float hi) {
    Plane p(w, h);
    const float a = angleDeg * kPi / 180.0f;
    const float tan_a = std::tan(a);
    const float xMid = 0.5f * float(w);
    for (int y = 0; y < h; ++y) {
        const float edgeX = xMid + tan_a * (float(y) - 0.5f * float(h));
        for (int x = 0; x < w; ++x) {
            const float d = float(x) + 0.5f - edgeX;      // signed distance in pixels
            const float cov = std::clamp(d + 0.5f, 0.0f, 1.0f);
            p.at(x, y) = lo + (hi - lo) * cov;
        }
    }
    return p;
}

inline Plane siemensStar(int w, int h, int spokes) {
    Plane p(w, h);
    const Frame f = frameOf(w, h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            const float dx = float(x) - f.cx, dy = float(y) - f.cy;
            const float th = std::atan2(dy, dx);
            p.at(x, y) = 0.5f + 0.5f * ((std::cos(th * float(spokes)) > 0.0f) ? 1.0f : -1.0f);
        }
    return p;
}

inline Plane singlePoint(int w, int h) {
    Plane p(w, h);
    p.at(w / 2, h / 2) = 1.0f;
    return p;
}

inline Image toImage(const Plane& p) {
    Image im(p.w, p.h);
    for (int y = 0; y < p.h; ++y)
        for (int x = 0; x < p.w; ++x)
            for (int c = 0; c < 3; ++c) im.at(x, y, c) = p.at(x, y);
    return im;
}

inline Plane luminance(const Image& im) {
    Plane p(im.w, im.h);
    for (int y = 0; y < im.h; ++y)
        for (int x = 0; x < im.w; ++x)
            p.at(x, y) = color::rec2020Luma(color::RGB{im.at(x, y, 0), im.at(x, y, 1), im.at(x, y, 2)});
    return p;
}

}  // namespace lens::targets
