#pragma once
#include "lenscore/image.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace lens::film {

// Film grain.
//
// Grain is not white noise laid over a picture. On film it is silver halide
// crystals developing or not developing, so the visible fluctuation is largest
// in the mid densities and vanishes at both ends: clear film has nothing to
// develop, and fully exposed film has nothing left undeveloped. The variance of
// a binomial process, v(1-v), captures that shape, and its square root is the
// amplitude. That is why grain reads as "in the mids" rather than as an even
// dusting, and why simply adding uniform noise never looks like film.
//
// Grain also has a SIZE, independent of image resolution: the same emulsion
// scanned larger gives larger clumps in pixels. So the noise is generated on a
// coarser lattice and interpolated up, rather than being per-pixel.

struct GrainParams {
    float amount = 0.0f;    // 0 = off, 1 = strong
    float size   = 1.5f;    // clump size in pixels
    float colour = 0.35f;   // 0 = all three channels identical, 1 = independent
    uint32_t seed = 1;
};

namespace detail {

// A cheap integer hash. Deterministic, so a preview and the final render agree
// and dragging a slider does not make the grain crawl.
inline float hashNoise(uint32_t x, uint32_t y, uint32_t c, uint32_t seed) {
    uint32_t h = x * 374761393u + y * 668265263u + c * 2246822519u + seed * 3266489917u;
    h = (h ^ (h >> 13)) * 1274126177u;
    h ^= h >> 16;
    // Two hashes averaged: a single uniform sample has heavy tails compared to
    // the roughly normal fluctuation of a developed emulsion.
    uint32_t g = h * 2654435761u;
    g ^= g >> 15;
    const float a = float(h) * (1.0f / 4294967296.0f);
    const float b = float(g) * (1.0f / 4294967296.0f);
    return a + b - 1.0f;    // -1 .. 1, centre-weighted
}

inline float lattice(float fx, float fy, uint32_t c, uint32_t seed) {
    const int x0 = int(std::floor(fx)), y0 = int(std::floor(fy));
    const float tx = fx - float(x0), ty = fy - float(y0);
    // Smoothstep, so clumps blend without visible lattice lines.
    const float sx = tx * tx * (3.0f - 2.0f * tx);
    const float sy = ty * ty * (3.0f - 2.0f * ty);
    const float n00 = hashNoise(uint32_t(x0),     uint32_t(y0),     c, seed);
    const float n10 = hashNoise(uint32_t(x0 + 1), uint32_t(y0),     c, seed);
    const float n01 = hashNoise(uint32_t(x0),     uint32_t(y0 + 1), c, seed);
    const float n11 = hashNoise(uint32_t(x0 + 1), uint32_t(y0 + 1), c, seed);
    return (n00 * (1 - sx) + n10 * sx) * (1 - sy) + (n01 * (1 - sx) + n11 * sx) * sy;
}

}  // namespace detail

// Applies grain in place. pixelScale is how many full-size pixels one pixel of
// this buffer covers, so a preview proxy gets clumps of the same apparent size
// as the final render rather than a much finer sprinkle.
inline void applyGrain(Image& img, const GrainParams& g, double pixelScale = 1.0) {
    if (g.amount <= 0.0f || img.w <= 0 || img.h <= 0) return;

    const float size = std::max(0.25f, float(g.size / std::max(0.0001, pixelScale)));
    const float inv = 1.0f / size;
    const float colour = std::clamp(g.colour, 0.0f, 1.0f);

    for (int y = 0; y < img.h; ++y) {
        for (int x = 0; x < img.w; ++x) {
            // One shared sample plus a per-channel one, mixed by `colour`: real
            // emulsions are correlated across layers, not independent.
            const float shared = detail::lattice(float(x) * inv, float(y) * inv, 0, g.seed);
            for (int c = 0; c < 3; ++c) {
                float& v = img.at(x, y, c);
                const float own = colour > 0.0f
                    ? detail::lattice(float(x) * inv, float(y) * inv, uint32_t(c + 1), g.seed)
                    : 0.0f;
                const float n = shared * (1.0f - colour) + own * colour;

                // Amplitude follows sqrt(v(1-v)): nothing in clear film, nothing
                // in solid black, most in the mids.
                const float d = std::clamp(v, 0.0f, 1.0f);
                // 0.16, not 0.5. At full travel the old constant moved a mid grey
                // by a tenth of full scale, which reads as sand rather than as
                // grain; this lands nearer 3% and leaves headroom above 100.
                v += n * g.amount * 0.16f * std::sqrt(std::max(0.0f, d * (1.0f - d)));
            }
        }
    }
}

}  // namespace lens::film
