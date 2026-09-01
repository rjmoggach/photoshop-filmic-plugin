#pragma once
#include "lenscore/color/upsample.hpp"
#include <algorithm>
#include <cstdio>
#include <optional>
#include <string>
#include <vector>

namespace lens::color {

struct SpecTable {
    int res = 0;
    std::vector<float> data;   // [i][z][y][x][0..2], i is the largest-component index
};

inline float smoothstep(float x) { return x * x * (3.0f - 2.0f * x); }

// Warped brightness axis: resolution concentrated near black and near white.
inline float axisScale(int k, int res) { return smoothstep(smoothstep(float(k) / float(res - 1))); }

inline size_t tableIndex(int res, int i, int z, int y, int x) {
    return ((((size_t(i) * res + z) * res + y) * res + x)) * 3;
}

inline SpecTable buildTable(int res) {
    SpecTable t; t.res = res; t.data.assign(size_t(3) * res * res * res * 3, 0.0f);
    for (int i = 0; i < 3; ++i) {
        for (int z = 0; z < res; ++z) {
            const float scale = std::max(1e-4f, axisScale(z, res));
            Coeffs warm{};                                   // warm start down each row
            for (int y = 0; y < res; ++y) {
                for (int x = 0; x < res; ++x) {
                    const float a = float(x) / float(res - 1);
                    const float b = float(y) / float(res - 1);
                    float rgb[3];
                    rgb[i]           = scale;
                    rgb[(i + 1) % 3] = a * scale;
                    rgb[(i + 2) % 3] = b * scale;
                    // rgb[] already has max component == scale <= 1, satisfying fitCoeffs's
                    // normalised-input contract directly -- do NOT divide by scale here, or
                    // every z-slice collapses onto the same amplitude-1 target and the
                    // brightness axis stops encoding brightness (see lookup's z-bracket,
                    // which compares against the true unnormalised target value).
                    const RGB target{rgb[0], rgb[1], rgb[2]};
                    Coeffs c = fitCoeffs(target, warm);
                    // Near-degenerate targets (very low chroma on one axis) occasionally
                    // strand the warm-started solve in a bad basin: the fit still "completes"
                    // but lands far from target with wildly larger coefficients than its
                    // neighbours. Left alone this poisons every later warm start down the row
                    // AND breaks lookup's bilinear blend, since interpolating two very
                    // different coefficient triples does not track a blend of the spectra
                    // they produce. Detect it and retry cold; keep whichever solve is closer.
                    auto sqErr = [&](const Coeffs& cand) {
                        const RGB got = spectrumToRec2020(cand);
                        const float dr = got.r - target.r, dg = got.g - target.g, db = got.b - target.b;
                        return dr * dr + dg * dg + db * db;
                    };
                    float e = sqErr(c);
                    if (e > 1e-4f) {
                        const Coeffs cold = fitCoeffs(target, Coeffs{});
                        const float ec = sqErr(cold);
                        if (ec < e) { c = cold; e = ec; }
                    }
                    warm = c;
                    const size_t o = tableIndex(res, i, z, y, x);
                    t.data[o + 0] = c.c0; t.data[o + 1] = c.c1; t.data[o + 2] = c.c2;
                }
            }
        }
    }
    return t;
}

inline Coeffs lookup(const SpecTable& t, const RGB& colour) {
    const float v[3] = {colour.r, colour.g, colour.b};
    int i = 0;
    if (v[1] > v[i]) i = 1;
    if (v[2] > v[i]) i = 2;
    const float mx = v[i];
    if (mx <= 0.0f) return Coeffs{0.0f, 0.0f, -1e4f};   // black: flat zero spectrum

    const float a = v[(i + 1) % 3] / mx;
    const float b = v[(i + 2) % 3] / mx;

    // Invert the warped brightness axis by search; res is small so this is cheap.
    const float target = std::min(mx, 1.0f);
    int z0 = 0;
    while (z0 + 2 < t.res && axisScale(z0 + 1, t.res) < target) ++z0;
    const float s0 = axisScale(z0, t.res), s1 = axisScale(z0 + 1, t.res);
    const float fz = (s1 > s0) ? std::clamp((target - s0) / (s1 - s0), 0.0f, 1.0f) : 0.0f;

    const float fx = a * (t.res - 1), fy = b * (t.res - 1);
    const int x0 = std::clamp(int(fx), 0, t.res - 2), y0 = std::clamp(int(fy), 0, t.res - 2);
    const float dx = fx - x0, dy = fy - y0;

    Coeffs out{};
    float* o[3] = {&out.c0, &out.c1, &out.c2};
    for (int k = 0; k < 3; ++k) {
        float acc = 0.0f;
        for (int zz = 0; zz < 2; ++zz)
            for (int yy = 0; yy < 2; ++yy)
                for (int xx = 0; xx < 2; ++xx) {
                    const float w = (zz ? fz : 1 - fz) * (yy ? dy : 1 - dy) * (xx ? dx : 1 - dx);
                    acc += w * t.data[tableIndex(t.res, i, z0 + zz, y0 + yy, x0 + xx) + k];
                }
        *o[k] = acc;
    }
    return out;
}

inline bool writeTable(const std::string& path, const SpecTable& t) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::fwrite("LSPT", 1, 4, f);
    std::fwrite(&t.res, sizeof(int), 1, f);
    std::fwrite(t.data.data(), sizeof(float), t.data.size(), f);
    std::fclose(f);
    return true;
}

inline std::optional<SpecTable> readTable(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return std::nullopt;
    char magic[4];
    SpecTable t;
    if (std::fread(magic, 1, 4, f) != 4 || std::string(magic, 4) != "LSPT" ||
        std::fread(&t.res, sizeof(int), 1, f) != 1 || t.res < 2) { std::fclose(f); return std::nullopt; }
    t.data.resize(size_t(3) * t.res * t.res * t.res * 3);
    const bool ok = std::fread(t.data.data(), sizeof(float), t.data.size(), f) == t.data.size();
    std::fclose(f);
    return ok ? std::optional<SpecTable>(t) : std::nullopt;
}

}  // namespace lens::color
