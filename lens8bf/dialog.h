#pragma once
#include "LensParams.h"
#include <string>
#include <vector>

// Everything the dialog needs from Photoshop, with no Photoshop headers in it.
struct ProxySource {
    virtual ~ProxySource() = default;

    // Fetch a view of the image being filtered.
    //   maxW,maxH      the size of the well it has to land in. Both are needed:
    //                  clamping width and height to the same number crops a
    //                  SQUARE out of a wide frame the moment you zoom in.
    //   scale          how many image pixels one preview pixel covers; 1 is 100%,
    //                  0 asks for whatever fits the whole area
    //   cxFrac,cyFrac  where to centre the view, 0..1 across the filter area
    // Reports the scale actually used through pixelScale.
    virtual bool fetchProxy(int maxW, int maxH, int scale, double cxFrac, double cyFrac,
                            std::vector<float>& rgb, int& w, int& h, double& pixelScale) = 0;

    // Size of the whole filter area, in image pixels.
    virtual void areaSize(int& w, int& h) const = 0;

    virtual std::string tablePath() const = 0;

    // Filled in when fetchProxy fails, so the failure can be shown rather than
    // looking like a plug-in that does nothing.
    std::string lastError;
};

// Runs the modal dialog. Returns true if the user pressed OK.
bool showLensDialog(LensControls& controls, ProxySource& source);
