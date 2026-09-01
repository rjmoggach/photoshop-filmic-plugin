#pragma once
#include "lenscore/image.hpp"
#include <cstdio>
#include <optional>
#include <string>

namespace lens::pfm {

// PFM stores rows bottom-to-top. Negative scale means little-endian.
inline bool write(const std::string& path, const Image& im) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::fprintf(f, "PF\n%d %d\n-1.0\n", im.w, im.h);
    for (int y = im.h - 1; y >= 0; --y)
        std::fwrite(&im.px[size_t(y) * im.w * 3], sizeof(float), size_t(im.w) * 3, f);
    std::fclose(f);
    return true;
}

inline std::optional<Image> read(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return std::nullopt;
    char magic[3] = {};
    int w = 0, h = 0; double scale = 0;
    if (std::fscanf(f, "%2s %d %d %lf", magic, &w, &h, &scale) != 4 ||
        std::string(magic) != "PF" || w <= 0 || h <= 0) { std::fclose(f); return std::nullopt; }
    std::fgetc(f);  // single whitespace byte after the scale
    Image im(w, h);
    for (int y = h - 1; y >= 0; --y)
        if (std::fread(&im.px[size_t(y) * w * 3], sizeof(float), size_t(w) * 3, f) != size_t(w) * 3)
            { std::fclose(f); return std::nullopt; }
    std::fclose(f);
    return im;
}

}  // namespace lens::pfm
