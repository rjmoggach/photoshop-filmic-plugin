#pragma once
#include <vector>

namespace lens {

struct Plane {
    int w = 0, h = 0;
    std::vector<float> v;
    Plane() = default;
    Plane(int width, int height) : w(width), h(height), v(size_t(width) * height, 0.0f) {}
    float&       at(int x, int y)       { return v[size_t(y) * w + x]; }
    const float& at(int x, int y) const { return v[size_t(y) * w + x]; }
};

}  // namespace lens
