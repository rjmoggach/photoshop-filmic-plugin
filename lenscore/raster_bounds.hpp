#pragma once
#include <cstddef>

namespace lens {

// Shared ceiling for any raster's pixel count (width * height), used by both
// Plane (single channel) and Image (interleaved RGB, 3x this many floats).
//
// 1,000,000,000 (1e9) pixels is comfortably larger than any real Photoshop
// document this plugin will encounter: a 20,000 x 20,000px canvas -- already
// far beyond a typical print or web asset -- is 4e8 pixels, under half this
// ceiling. At the same time a single Plane at the ceiling is 4GB (1e9 floats),
// so the bound also stops a bad width/height pair (e.g. from corrupted input
// or an upstream overflow) from silently trying to allocate an amount of
// memory that would thrash or crash the host process before the allocation
// even starts.
inline constexpr size_t kMaxRasterPixels = 1'000'000'000;

// Validates width/height and returns the pixel count (width * height), or 0
// if the dimensions are invalid: non-positive, or a product that would
// exceed kMaxRasterPixels. Never computes width * height before checking it
// against the ceiling, so a pair of dimensions large enough to overflow
// size_t is rejected rather than silently wrapping into a small, wrong count.
inline constexpr size_t validatedPixelCount(int width, int height) {
    if (width <= 0 || height <= 0) return 0;
    const size_t w = size_t(width), h = size_t(height);
    if (w > kMaxRasterPixels / h) return 0;   // w * h would exceed the ceiling (or overflow)
    return w * h;
}

}  // namespace lens
