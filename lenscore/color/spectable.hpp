#pragma once
#include "lenscore/color/upsample.hpp"
#include <algorithm>
#include <cstdio>
#include <optional>
#include <stdexcept>
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

// z=0 is an exact-black target, which is degenerate for fitCoeffs (there is no finite
// (c0,c1,c2) that reaches zero at every wavelength; the sigmoid only approaches it as
// c2 -> -inf). A tiny floor gives it a well-posed, fittable target instead. This must
// be the ONLY place that floors the brightness axis: buildTable calls it to choose the
// fit target for each slice, and lookup calls it to find that slice's bracket position.
// If the two ever used different floors -- or one used the raw axisScale -- the slices
// they'd store data at and the slices lookup would search for would silently disagree,
// which is exactly the class of bug this function exists to rule out by construction.
inline constexpr float kMinBrightness = 1e-6f;
inline float axisBrightness(int k, int res) { return std::max(kMinBrightness, axisScale(k, res)); }

inline size_t tableIndex(int res, int i, int z, int y, int x) {
    return ((((size_t(i) * res + z) * res + y) * res + x)) * 3;
}

// Real tables ship at res=64. This bounds a corrupted or hostile header well above any
// table this project would plausibly generate, while keeping the res^3 arithmetic below
// -- and the allocation it drives -- far from a size_t overflow or a surprise
// multi-gigabyte allocation. Shared by readTable (below) and buildTable (immediately
// below this), so a generated table and a loaded table can never disagree about what
// counts as a valid resolution.
inline constexpr int kMaxTableRes = 256;

// buildTable fits every entry against the dense spectral integral. A renderer
// that integrates a handful of bands does not reproduce that integral, so a
// table fitted to the dense model hands it the wrong colour. buildTableAgainst
// fits against whatever reconstruction the caller actually uses.
template <typename Forward>
inline SpecTable buildTableAgainst(int res, const Forward& forward) {
    if (res < 2 || res > kMaxTableRes) {
        throw std::invalid_argument("buildTableAgainst: res out of range");
    }
    SpecTable t; t.res = res; t.data.assign(size_t(3) * res * res * res * 3, 0.0f);
    for (int i = 0; i < 3; ++i)
        for (int z = 0; z < res; ++z) {
            const float scale = axisBrightness(z, res);
            Coeffs warm{};
            for (int y = 0; y < res; ++y)
                for (int x = 0; x < res; ++x) {
                    const float a = float(x) / float(res - 1);
                    const float b = float(y) / float(res - 1);
                    float rgb[3];
                    rgb[i]           = scale;
                    rgb[(i + 1) % 3] = a * scale;
                    rgb[(i + 2) % 3] = b * scale;
                    const RGB target{rgb[0], rgb[1], rgb[2]};
                    Coeffs c = fitCoeffsAgainst(forward, target, warm);
                    auto sqErr = [&](const Coeffs& cand) {
                        const RGB got = forward(cand);
                        const float dr = got.r - target.r, dg = got.g - target.g, db = got.b - target.b;
                        return dr * dr + dg * dg + db * db;
                    };
                    if (sqErr(c) > 1e-4f) {
                        const Coeffs cold = fitCoeffsAgainst(forward, target, Coeffs{});
                        if (sqErr(cold) < sqErr(c)) c = cold;
                    }
                    warm = c;
                    const size_t o = tableIndex(res, i, z, y, x);
                    t.data[o + 0] = c.c0; t.data[o + 1] = c.c1; t.data[o + 2] = c.c2;
                }
        }
    return t;
}

inline SpecTable buildTable(int res) {
    // res is a caller-controlled dimension that directly sizes an allocation
    // (3*res^3*3 floats below) -- readTable applies the identical [2, kMaxTableRes]
    // bound to the same field read from a file header (see below); buildTable's res
    // reaches this same arithmetic from rgb2spec/main.cpp's atoi(argv[1]), completely
    // unvalidated before this fix, and res <= 1 also divides by (res - 1) in
    // axisScale/lookup. This is a programmer/caller error, not untrusted file bytes,
    // so it throws rather than returning an empty table a caller might not check.
    if (res < 2 || res > kMaxTableRes) {
        throw std::invalid_argument(
            "buildTable: res must be in [2, " + std::to_string(kMaxTableRes) +
            "]; got " + std::to_string(res));
    }
    SpecTable t; t.res = res; t.data.assign(size_t(3) * res * res * res * 3, 0.0f);
    for (int i = 0; i < 3; ++i) {
        for (int z = 0; z < res; ++z) {
            const float scale = axisBrightness(z, res);
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
    // Must use axisBrightness (the same floored convention buildTable fit against), not
    // the raw axisScale, or the bracket this finds disagrees with what's actually stored.
    const float target = std::min(mx, 1.0f);
    int z0 = 0;
    while (z0 + 2 < t.res && axisBrightness(z0 + 1, t.res) < target) ++z0;
    const float s0 = axisBrightness(z0, t.res), s1 = axisBrightness(z0 + 1, t.res);
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
        std::fread(&t.res, sizeof(int), 1, f) != 1 ||
        t.res < 2 || t.res > kMaxTableRes) { std::fclose(f); return std::nullopt; }

    const size_t floatCount = size_t(3) * t.res * t.res * t.res * 3;
    const size_t expectedBytes = floatCount * sizeof(float);

    // Cross-check the header's claimed size against the bytes actually left in the file
    // before trusting it enough to size an allocation. An assert would compile away
    // under NDEBUG and leave this exact allocation sized from an unverified header in
    // a release build -- the one build where a corrupt file is most likely to appear.
    if (std::fseek(f, 0, SEEK_END) != 0) { std::fclose(f); return std::nullopt; }
    const long end = std::ftell(f);
    if (end < 0 || size_t(end) != 8 + expectedBytes) { std::fclose(f); return std::nullopt; }
    if (std::fseek(f, 8, SEEK_SET) != 0) { std::fclose(f); return std::nullopt; }

    t.data.resize(floatCount);
    const bool ok = std::fread(t.data.data(), sizeof(float), t.data.size(), f) == t.data.size();
    std::fclose(f);
    return ok ? std::optional<SpecTable>(t) : std::nullopt;
}

}  // namespace lens::color
