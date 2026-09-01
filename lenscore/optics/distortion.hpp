#pragma once

namespace lens::optics {

struct Distortion { float k1 = 0, k2 = 0, k3 = 0, p1 = 0, p2 = 0; };

// Maps ideal normalised coordinates to the distorted position they came from.
inline void applyDistortion(const Distortion& d, float ux, float uy, float& dx, float& dy) {
    const float r2 = ux * ux + uy * uy;
    const float radial = 1.0f + d.k1 * r2 + d.k2 * r2 * r2 + d.k3 * r2 * r2 * r2;
    dx = ux * radial + 2.0f * d.p1 * ux * uy + d.p2 * (r2 + 2.0f * ux * ux);
    dy = uy * radial + d.p1 * (r2 + 2.0f * uy * uy) + 2.0f * d.p2 * ux * uy;
}

}  // namespace lens::optics
