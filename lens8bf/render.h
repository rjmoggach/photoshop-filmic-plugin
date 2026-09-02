#pragma once
#include "LensParams.h"
#include <cstddef>
#include <string>

// One entry point for both the dialog's proxy and the final full-size pass, so
// the preview cannot drift from what Apply produces.
//
// rgb is interleaved float RGB, w*h*3, nominally 0..1. It is rendered in place.
//
// pixelScale is how many full-image pixels one pixel of THIS buffer covers: 1.0
// for the final render, and the proxy's scale factor for the preview. Blur radii
// are in pixels, so without it a proxy would show a proportionally smaller blur
// than the real thing and the preview would lie.
//
// Returns an empty string on success, or a message to show the user.
std::string lensRender(float* rgb, int w, int h, const LensControls& c,
                       double pixelScale, const std::string& tablePath);
