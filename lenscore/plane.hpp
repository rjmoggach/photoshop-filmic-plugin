#pragma once
#include "lenscore/raster_bounds.hpp"
#include <vector>

namespace lens {

struct Plane {
    int w = 0, h = 0;
    std::vector<float> v;
    Plane() = default;
    // Dimensions outside validatedPixelCount's bounds (non-positive, or a
    // product exceeding kMaxRasterPixels) construct an empty Plane (w = h = 0,
    // v empty) -- the same zero-size state the default constructor already
    // produces, rather than a half-set width/height with a mismatched buffer.
    Plane(int width, int height) {
        const size_t n = validatedPixelCount(width, height);
        if (n > 0) { w = width; h = height; }
        v.assign(n, 0.0f);
    }
    float&       at(int x, int y)       { return v[size_t(y) * w + x]; }
    const float& at(int x, int y) const { return v[size_t(y) * w + x]; }
};

}  // namespace lens
