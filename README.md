# Filmic Lens Plugin

Physically accurate lens optics for Photoshop: chromatic aberration, vignetting,
field curvature and edge blur, computed from optical models whose parameters can
be measured from real glass.

**Status: loadable Photoshop plugin, plus the optics core it runs on.** The
panel gives you sliders for chromatic aberration, vignetting, edge blur,
astigmatism, coma and aperture, and applies them to the active layer. The
physics runs natively — the panel is UI only.

## What's here

| Component | What it is |
|---|---|
| `lenscore/` | Header-only C++20 optics library. Links nothing, includes only the standard library. |
| `plugin/` | The UXP panel: manifest, UI, and the JavaScript that moves pixels. |
| `lensaddon/` | The native `.uxpaddon` — a shared library that hands pixels to `lenscore`. |
| `lenscli/` | Terminal renderer, for measurement rather than viewing. |
| `lensdata/` | `.lens` parameter files and the 64³ spectral upsampling table. |
| `tests/` | 166 cases, 29,467 assertions. |
| `docs/superpowers/specs/` | The design spec, including every accepted limitation with its measurements. |

The split is deliberate: `lenscore` has no host dependency, so the physics is
unaffected by how it is eventually driven. That paid for itself once already —
the host route changed from an 8BF filter to a UXP Hybrid plugin partway through
and not one line of optics moved.

## Build

Requires CMake 3.24+ and a C++20 compiler. No other dependencies; doctest is
fetched automatically for the tests.

```
cmake -B build -S .
cmake --build build
ctest --test-dir build --output-on-failure
```

Expect 166 cases and 29,467 assertions, all passing, in about 20 seconds.

## Install the plugin

Requires Photoshop 24.2 or newer (you are on 2026) and Adobe's UXP Developer
Tools.

```
# 1. Build the native addon. Point at your unpacked UXP Hybrid SDK.
cmake -B build -S . -DUXP_ADDON_SDK_ROOT=/path/to/uxp-hybrid-plugin-sdk-main
cmake --build build --target lensaddon
```

That stages `lens.uxpaddon` and the spectral table into `plugin/`.

```
# 2. In UXP Developer Tools: Add Plugin -> select plugin/manifest.json
# 3. Click ... next to the entry -> Load
# 4. In Photoshop: Plugins -> Filmic Lens
```

Open an image, select a layer, move the sliders, click **Apply to layer**. The
panel reports the render time; the result is a single undoable history step.

If the panel says the native addon failed to load, the addon was not built or
not staged — re-run step 1 and reload in UDT.

## What you can do with this today

**Run the test suite.** This is the meaningful thing to do right now, and it is
not a formality — the suite is the project's evidence. Five acceptance tests
measure physical predictions rather than merely executing:

- lateral chromatic aberration fringe width grows 0.0006 → 0.88 px across the field
- corner MTF50 falls to 24% of centre under field curvature
- astigmatism separates sagittal from tangential resolution by 83%
- an achromat's focus error is exactly zero at both correction wavelengths, so it
  fringes green and magenta rather than red and blue
- stopping down removes mechanical vignetting, leaving the natural falloff

**Known rough edges.** Renders are synchronous, so a large layer will hold the
UI while it works — the panel reports elapsed time when it finishes. Start with
a modest layer, and with **Spectral bands** at 3; raise it to 7 or 11 for
smoother fringes once you have a look you like. `lenscli` still reads and writes
PFM only, which is useful for measurement, not for viewing.

## Accuracy, and its limits

Every effect derives from an optical model with parameters that can be measured,
and every model is checked against a numeric metric in CI. The slanted-edge MTF
verifies itself against the analytic Gaussian result `0.1874/σ` before it is used
to judge any lens.

Three limitations are measured and recorded in the spec rather than left implicit:

- **Saturated colour.** The working space is Rec.2020, which is wider than a
  bounded reflectance spectrum can represent. Colours near the primaries render
  as their nearest representable metamer. More table resolution does not help.
- **White point.** The equal-energy to D65 adaptation is a diagonal scale, exact
  at white and approximate elsewhere.
- **Corner vignetting.** The spatially varying convolution blends overlapping
  patches, so a patch straddling the corner mixes in a less-vignetted interior
  kernel. The rendered corner sits up to ~17% above the closed-form prediction.

## Where the physics comes from

Jakob & Hanika (spectral upsampling), Jeong et al. (chromatic aberration),
Cholewiak et al. (PSF from the pupil function), Aggarwal, Hua & Ahuja (vignetting
beyond `cos⁴`, and pupil apodization), Schuler et al. (efficient filter flow),
Newson, Delon & Galerne (grain), Geigel & Musgrave (film), Hullin et al.
(polynomial optics). The papers are not redistributed here.

## Next

The film stage, the UXP Hybrid host, and `lensfit` — the calibration tool that
turns chart photographs into a `.lens` file — each get their own plan. `lensfit`
is what makes "accurate" a measurement rather than a claim.
