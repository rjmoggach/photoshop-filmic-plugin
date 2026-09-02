#pragma once
#include "lenscore/raster_bounds.hpp"
#include <vector>
#include <cstddef>

namespace lens {

struct Image {
    int w = 0, h = 0;
    std::vector<float> px;   // interleaved RGB

    Image() = default;
    // Same bound and same empty-on-out-of-range behaviour as Plane (see
    // raster_bounds.hpp) -- one shared guard rather than two that could drift.
    Image(int width, int height) {
        const size_t n = validatedPixelCount(width, height);
        if (n > 0) { w = width; h = height; }
        px.assign(n * 3, 0.0f);
    }

    float&       at(int x, int y, int c)       { return px[(size_t(y) * w + x) * 3 + c]; }
    const float& at(int x, int y, int c) const { return px[(size_t(y) * w + x) * 3 + c]; }

    bool sameSize(const Image& o) const { return w == o.w && h == o.h; }
};

}  // namespace lens
