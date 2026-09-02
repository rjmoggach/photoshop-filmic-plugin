#pragma once
#include "LensParams.h"
#include <string>
#include <vector>

// Starting points, not destinations. Each is a plausible combination that gets
// you past a blank slate; every control stays live afterwards.
struct LensPreset {
    const char*  name;
    const char*  group;
    const char*  desc;    // one line, shown under the menu
    LensControls c;
};

inline std::vector<LensPreset> lensPresets() {
    std::vector<LensPreset> p;

    p.push_back({"Default", "", "Everything off. A clean starting point.", LensControls{}});

    // ---- cinema ----
    LensControls vintage{};
    vintage.distortion = -18; vintage.chromaticAberration = 28; vintage.vignette = 42;
    vintage.edgeBlur = 46;    vintage.astigmatism = 24;         vintage.coma = 18;
    vintage.tStop = 2.0;      vintage.grain = 26;               vintage.grainSize = 1.6;
    p.push_back({"Vintage 35mm", "Cinema", "Soft, warm-edged and grainy. Barrel bow and heavy corners.", vintage});

    LensControls scope{};
    scope.squeeze = 2.0;      scope.chromaticAberration = 20;   scope.vignette = 38;
    scope.edgeBlur = 34;      scope.coma = 26;                  scope.tStop = 2.0;
    scope.grain = 18;         scope.grainSize = 1.8;
    p.push_back({"Anamorphic scope", "Cinema", "Oval bokeh and sideways streaks, the classic scope look.", scope});

    LensControls toy{};
    toy.distortion = -62;     toy.chromaticAberration = 55;     toy.vignette = 78;
    toy.edgeBlur = 70;        toy.astigmatism = 45;             toy.coma = 40;
    toy.tStop = 1.4;          toy.grain = 45;                   toy.grainSize = 2.2;
    p.push_back({"Toy camera", "Cinema", "Cheap plastic optics: heavy bulge, colour fringes, dark corners.", toy});

    LensControls portrait{};
    portrait.vignette = 26;   portrait.edgeBlur = 30;           portrait.astigmatism = 14;
    portrait.tStop = 1.4;     portrait.grain = 12;              portrait.grainSize = 1.2;
    p.push_back({"Soft portrait", "Cinema", "Gentle edge softness and a light vignette to hold the face.", portrait});

    LensControls stock{};
    stock.grain = 40;         stock.grainSize = 1.4;            stock.grainColour = 20;
    p.push_back({"Film stock only", "Cinema", "No lens character at all -- just emulsion.", stock});

    // ---- photographic ----
    // A zoom at its wide end is the classic barrel case, and the corners are
    // where its colour error shows.
    LensControls wide{};
    wide.distortion = -46;    wide.chromaticAberration = 30;    wide.vignette = 44;
    wide.edgeBlur = 26;       wide.astigmatism = 16;            wide.tStop = 4.0;
    p.push_back({"Wide zoom, wide end", "Photographic", "The classic barrel and corner colour of a zoom at 24mm.", wide});

    // Wide open, a fast normal prime is soft and dark in the corners and nearly
    // straight; that combination is what makes it read as "fast fifty".
    LensControls fifty{};
    fifty.distortion = -8;    fifty.chromaticAberration = 16;   fifty.vignette = 62;
    fifty.edgeBlur = 34;      fifty.astigmatism = 20;           fifty.coma = 22;
    fifty.tStop = 1.4;
    p.push_back({"Fast 50mm wide open", "Photographic", "Nearly straight, dark in the corners, soft off-centre.", fifty});

    LensControls eightyfive{};
    eightyfive.distortion = 6; eightyfive.chromaticAberration = 8; eightyfive.vignette = 30;
    eightyfive.edgeBlur = 16;  eightyfive.tStop = 2.0;
    p.push_back({"85mm portrait", "Photographic", "Well corrected, a hint of pincushion, mild falloff.", eightyfive});

    // Long lenses pinch rather than bulge, and correct well.
    LensControls tele{};
    tele.distortion = 18;     tele.chromaticAberration = 6;     tele.vignette = 22;
    tele.edgeBlur = 10;       tele.tStop = 4.0;
    p.push_back({"Telephoto 200mm", "Photographic", "Long glass: slight pinch, very clean colour, little falloff.", tele});

    LensControls compact{};
    compact.distortion = -56; compact.chromaticAberration = 46; compact.vignette = 54;
    compact.edgeBlur = 44;    compact.astigmatism = 30;         compact.coma = 24;
    compact.tStop = 2.8;      compact.grain = 10;
    p.push_back({"Compact camera zoom", "Photographic", "Small-sensor compromises: strong bulge and visible fringing.", compact});

    LensControls fish{};
    fish.distortion = -94;    fish.chromaticAberration = 34;    fish.vignette = 48;
    fish.edgeBlur = 38;       fish.astigmatism = 22;            fish.tStop = 5.6;
    p.push_back({"Fisheye", "Photographic", "An effect, not a lens. Extreme bulge well past any real optic.", fish});

    // ---- characters drawn from families of real glass ----
    //
    // The distortion figures are deliberately modest. Almost every real lens sits
    // inside a Brown-Conrady k1 of about +-0.05, which is +-20 on this control;
    // the "fisheye" and "toy" entries above sit far outside that on purpose,
    // being effects rather than lenses.

    // Rangefinder glass is nearly rectilinear and very well colour-corrected,
    // and pays for its speed with corners that fall away hard wide open.
    LensControls rf{};
    rf.distortion = -6;       rf.chromaticAberration = 10;      rf.vignette = 72;
    rf.edgeBlur = 22;         rf.astigmatism = 14;              rf.coma = 10;
    rf.tStop = 1.4;           rf.grain = 14;                    rf.grainSize = 1.2;
    p.push_back({"Rangefinder 35mm f/1.4", "Inspired by", "Nearly rectilinear and very clean, corners fall hard wide open.", rf});

    // A fast reporter's normal lens: a touch of barrel, visible colour at the
    // edges, and coma on point highlights until it is stopped down.
    LensControls rep{};
    rep.distortion = -11;     rep.chromaticAberration = 20;     rep.vignette = 58;
    rep.edgeBlur = 26;        rep.astigmatism = 16;             rep.coma = 22;
    rep.tStop = 1.4;          rep.grain = 16;                   rep.grainSize = 1.4;
    p.push_back({"Reportage 50mm f/1.4", "Inspired by", "A touch of barrel, edge colour, coma on point highlights.", rep});

    // Pre-war double-Gauss descendants swirl: heavy coma and astigmatism at the
    // edge of the field drag the out-of-focus background around the subject.
    LensControls swirl{};
    swirl.distortion = -14;   swirl.chromaticAberration = 28;   swirl.vignette = 66;
    swirl.edgeBlur = 40;      swirl.astigmatism = 42;           swirl.coma = 50;
    swirl.tStop = 2.0;        swirl.grain = 22;                 swirl.grainSize = 1.7;
    p.push_back({"Swirly 58mm f/2 (old glass)", "Inspired by", "Pre-war descent: coma and astigmatism swirl the background.", swirl});

    // Projector glass was built to be flat and fast, never to be sharp at the
    // edge of a still frame: almost no distortion, and a field that curves away.
    LensControls proj{};
    proj.distortion = 4;      proj.chromaticAberration = 12;    proj.vignette = 52;
    proj.edgeBlur = 54;       proj.astigmatism = 28;            proj.coma = 14;
    proj.tStop = 2.0;         proj.grain = 18;                  proj.grainSize = 1.5;
    p.push_back({"Projection 100mm f/2", "Inspired by", "Projector glass: almost no distortion, strongly curved field.", proj});

    return p;
}
