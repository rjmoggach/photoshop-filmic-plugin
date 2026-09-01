#pragma once
#include <vector>
#include <cstddef>

namespace lens {

struct Image {
    int w = 0, h = 0;
    std::vector<float> px;   // interleaved RGB

    Image() = default;
    Image(int width, int height) : w(width), h(height), px(size_t(width) * height * 3, 0.0f) {}

    float&       at(int x, int y, int c)       { return px[(size_t(y) * w + x) * 3 + c]; }
    const float& at(int x, int y, int c) const { return px[(size_t(y) * w + x) * 3 + c]; }

    bool sameSize(const Image& o) const { return w == o.w && h == o.h; }
};

}  // namespace lens
