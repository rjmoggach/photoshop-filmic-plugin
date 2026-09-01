# Filmic Lens Plugin

Physically accurate lens optics for Photoshop: chromatic aberration, vignetting,
field curvature and edge blur, computed from optical models whose parameters can
be measured from real glass.

**Status: the optics core is complete and tested. There is no Photoshop plugin
yet.** This release is a host-free C++20 library and a terminal harness. See
[What you can do with this today](#what-you-can-do-with-this-today) before you
plan around it.

## What's here

| Component | What it is |
|---|---|
| `lenscore/` | Header-only C++20 optics library. Links nothing, includes only the standard library. |
| `lenscli/` | Terminal renderer. |
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

**What you cannot do yet:** point this at a photograph. `lenscli render` reads and
writes PFM, and nothing in the repo generates a PFM to feed it or converts one to
something you can look at. A synthetic-target generator and a viewable output
format are the obvious next step and are small; they are simply not in this
release.

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
