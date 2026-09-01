# Filmic Lens Plugin for Photoshop — Design

Date: 2026-08-31
Status: approved design, pre-implementation

## 1. Goal

A Photoshop filter plugin that applies **physically accurate** lens and film
artefacts to a still image: chromatic aberration, vignetting, edge blur and
field curvature, halation and bloom, plus a film-stock stage.

"Accurate" has a specific meaning here. Every effect derives from an optical
or photochemical model with named parameters that can be **measured from real
glass or real stock**, and every model is verified against a golden-image
metric in CI. Art-directable controls exist, but they sit on top of the
physics as explicit deviations, never in place of it.

## 2. Non-goals

- Video, batch processing, or Lightroom/Camera Raw integration.
- Lens flare ghosts (inter-element reflections). Deferred to a later release;
  the data format reserves a `flare.ghosts` field for it.
- Lens *correction* (removing aberration). We only add.
- Automatic lens detection from EXIF.

## 3. Architecture

Four components. Only one is coupled to Adobe.

| Component | Language | Depends on | Purpose |
|---|---|---|---|
| `lenscore` | C++20 | nothing | All optics and film physics. Linear float image in, image out. |
| `lenscli` | C++20 | `lenscore` | Terminal renderer and test harness. |
| `lens8bf` | C++20 | `lenscore`, Photoshop SDK | The `.8bf` filter. Tiling, bit depth, colour, scripting. |
| `lensui` | C++20 | `lenscore`, Dear ImGui | Cross-platform dialog. Metal on macOS, D3D11 on Windows. |

### 3.1 Why this split

`lenscore` has no Adobe headers, no UI, and no I/O. It is a pure function of
`(image, params) -> image`. Three consequences:

1. Optics can be developed and tuned from a terminal, with numeric assertions,
   without launching Photoshop. Tuning optics through a GUI is not viable.
2. The Photoshop-specific risk (SDK gating, tiling, signing, ImGui on two
   platforms) is confined to `lens8bf` and `lensui`. If the C++ host path
   proves too costly, `lenscore` compiles to WebAssembly and drops into a UXP
   panel with no change to the physics.
3. `lenscore` is the only thing that needs to be fast, so optimisation effort
   has one target.

### 3.2 File layout

```
lenscore/
  color/     transfer.hpp     EOTF / inverse EOTF, highlight recovery
             cie.hpp          CIE 1931 CMFs, spectrum -> XYZ -> linear RGB
             upsample.hpp     RGB -> reflectance spectrum (Jakob-Hanika)
  optics/    dispersion.hpp   Sellmeier / Abbe -> n(lambda) -> F(lambda)
             distortion.hpp   Brown-Conrady radial + tangential
             lateralca.hpp    per-wavelength magnification
             pupil.hpp        entrance/exit pupil clip, cat's-eye, blade polygon
             zernike.hpp      Zernike basis for defocus/astig/coma/spherical
             psf.hpp          assembles the PSF at a given field radius
             vignette.hpp     natural + mechanical + optical falloff
             glare.hpp        veiling glare PSF, diffraction star
  film/      halation.hpp     per-channel back-scatter re-exposure
             grain.hpp        stochastic Boolean-disc grain
             printvig.hpp     print/gate vignette, gate weave
  conv/      fft.hpp          real FFT
             eff.hpp          Efficient Filter Flow: patched varying convolution
  pipeline.hpp                stage ordering, the one public entry point
  params.hpp                  versioned POD parameter struct
lensdata/    *.lens  *.stock  JSON parameter sets
lenscli/     main.cpp
lens8bf/     pipl.r  main.cpp  tiles.cpp  descriptors.cpp  colour.cpp
lensui/      app.cpp  preview.cpp  backends/{metal,d3d11}
lensfit/     calibration tool: chart photos -> .lens file
tests/       golden/  metrics/{mtf.cpp,vignette.cpp,fringe.cpp,energy.cpp}
```

A file that grows past roughly 400 lines is a signal it holds two
responsibilities and should split.

## 4. Pipeline

The stage order is physical, not arbitrary. Lens effects form the image;
film effects happen in and behind the emulsion; grain is the last thing
that touches the signal.

```
document pixels
  |
  1. decode          -> linear scene-referred RGB, highlight recovery
  2. spectral uplift -> per-pixel reflectance spectrum (CA stage only)
  3. lateral CA      -> per-wavelength radial magnification
  4. varying PSF     -> axial CA + field curvature + astigmatism + coma
                        + cat's-eye pupil clip, convolved per wavelength
  5. integrate       -> spectrum back to linear RGB through CIE CMFs
  6. vignette        -> natural cos^4 x mechanical pupil clip
  7. glare / bloom   -> veiling glare halo + aperture diffraction star
  --- lens ends, film begins ---
  8. halation        -> per-channel back-scatter, red widest
  9. print vignette  -> print stage falloff
 10. grain           -> signal-dependent stochastic grain
 11. encode          -> back to document space, dither if 8-bit
```

Stages 2-5 are one fused loop over wavelength bands, not four passes. Each
band is magnified, convolved, and accumulated into the XYZ sum in a single
traversal, so the intermediate spectral image never exists in memory.

### 4.1 Stage 1 — decode and highlight recovery

Convert the document to linear scene-referred RGB using the document's ICC
profile transfer function.

**Highlight recovery is not optional.** Halation, bloom, and defocus bokeh are
all driven by values far above diffuse white. An 8-bit or display-referred
image clips everything to 1.0, and blooming a clipped image produces the flat,
grey, obviously-fake glow that gives away every cheap lens plugin. We apply an
inverse-knee expansion above a threshold (default 0.85) that lifts near-white
pixels to a configurable peak (default 8.0 linear). Users working in 32-bit
scene-referred float can switch it off, because their highlights are real.

### 4.2 Stage 2 — spectral uplift

The optics are spectral but the input is RGB. We reconstruct a plausible
smooth spectrum per pixel using the Jakob-Hanika sigmoid-polynomial model:
three coefficients per pixel, evaluated analytically at each wavelength, with
a precomputed coefficient lookup table. It is smooth, energy-sane, and cheap.

Two quality tiers:

- **Spectral** (default, N = 11 bands over 380-780nm): full dispersion, real
  purple and green fringes.
- **RGB** (preview, N = 3 at 650 / 510 / 475nm, per Jeong): three times
  faster, visibly banded on strong dispersion. Preview only.

Band placement uses **spectral importance sampling** through the inverse CDF
of the spectral equalizer weights (Jeong Sec. 5.2). Uniform sampling at N = 11
bands produces visible spectral banding; importance sampling does not.

### 4.3 Stage 3 — lateral chromatic aberration

Per-wavelength magnification about the optical centre (Jeong Eq. 8):

```
m(lambda, t) = 1 + k_l * (lambda_hat - lambda) * L(t, theta)
```

`t` is normalised field radius (0 at centre, 1 at the corner), `lambda_hat`
= 650nm is the reference so magnification is never negative, `k_l` comes from
the `.lens` file. `L(t, theta) = t` is the physical baseline; the expressive
layer replaces `L` with a composite mapping (Sec. 6).

Brown-Conrady radial and tangential distortion applies in the same resample,
so the image is warped once, not twice.

### 4.4 Stage 4 — the spatially varying PSF

This is the hard stage and the one that makes or breaks the look.

**Kernel shape.** The base defocus kernel is a **disc**, not a Gaussian
(ChromaBlur Alg. 1: defocus of a planar subject is convolution with a cylinder
kernel). A Gaussian defocus is the single most common tell of a fake lens
effect. On top of the disc:

- **Blade polygon.** At apertures below wide open, the disc becomes an
  N-sided polygon with per-blade curvature, rotated by the blade angle.
- **Cat's-eye clip.** Off axis, the pupil is the intersection of the aperture
  disc with the entrance and exit pupil circles, displaced by field angle.
  The kernel becomes lemon-shaped and tilts tangentially. This produces swirl
  as an emergent property, not as a rotation filter.
- **Astigmatism.** Sagittal and tangential focus differ off axis, so the
  kernel stretches into an ellipse oriented radially or tangentially. Modelled
  as Zernike `Z(2,+/-2)` scaled by `t^2`.
- **Coma.** Zernike `Z(3,1)`, giving the kernel a comet tail pointing away
  from centre, scaled by `t^3`.
- **Diffraction.** At small apertures, convolve the pupil with the Airy
  pattern. This also produces the aperture star for point highlights.

**Defocus radius.** Two contributions sum:

1. **Axial CA.** Focal length varies with wavelength via Sellmeier (Jeong
   Eq. 4) and the lens-maker's equation (Eq. 3), giving a per-wavelength focal
   depth `d_f(lambda)` (Eq. 5) and a circle of confusion `C = ((d_f - d)/d_f) * E`
   (Eq. 2). BK7 defaults: `B = [1.03961, 0.23179, 1.01047]`,
   `C = [0.00600, 0.02002, 103.56065]` in micrometres squared.
2. **Petzval field curvature.** The focal surface is a sphere, so defocus
   grows with `t^2` even for a flat subject.

**The depth problem, and the decision.** A 2D photograph has no depth map, but
axial CA needs an object distance `d`. Resolved as follows:

- **Default**: the image is taken to lie on the focal surface, so `d` comes
  entirely from field curvature. Axial CA then varies with `t^2`, which is the
  physically correct behaviour for a flat subject and needs no extra input.
- **Optional**: the user designates a layer or channel as a depth map. `d` is
  then read per pixel and true focus-plane axial CA appears. Exposed as a
  "Depth source" dropdown, default None.

This is stated explicitly because it is the one place the model is
under-determined by the input.

**Fast evaluation.** A naive per-pixel elliptical gather at a 20px radius over
50 megapixels is roughly 2e10 operations per wavelength band. Not viable.
We use **Efficient Filter Flow** (Hirsch/Schuler): partition into overlapping
patches, assign one PSF per patch, FFT-convolve each, and recombine with a
partition-of-unity window. Cost is `O(N log P)`.

One refinement that matters: our PSF is **rotationally covariant**. Its shape
depends only on field radius `t`; at a different angle `theta` it is the same
kernel rotated. So the PSF table collapses from 2D to a **1D function of `t`**,
patches are laid out in polar rings, and each ring's kernel is computed once
and rotated per patch. This cuts PSF construction cost by roughly the number
of angular divisions and guarantees exact rotational symmetry, which is
directly asserted in the test suite.

### 4.5 Stage 6 — vignetting

Three multiplicative terms, following Aggarwal/Hua/Ahuja on why `cos^4` alone
is insufficient for real lenses:

1. **Natural falloff.** `cos^n(theta)` with `theta = atan(r / f)`. `n = 4` is
   the ideal case; the exponent is exposed because real lenses with
   non-coincident pupils deviate from it.
2. **Mechanical vignetting.** The entrance pupil is clipped by the barrel and
   by the rear element. Computed as the overlap area of the aperture circle
   with two circles offset by field angle, divided by the unclipped area. This
   yields the correct hard shoulder wide open, and correctly rounds toward
   pure `cos^n` as the lens stops down.
3. **Optical vignetting.** Already implicit: stage 4 clips the kernel, so
   energy is lost off axis automatically. We only normalise here, and assert
   the two paths agree.

### 4.6 Stage 7 — glare and bloom

Veiling glare is scattering inside the barrel, and it is **not Gaussian**. Its
PSF has a heavy power-law tail spanning most of the frame at very low
amplitude. We model it as a sum of three Gaussians fitted to a power-law
profile, applied to the highlight-thresholded image in linear light.

The diffraction star comes free from stage 4's Airy convolution with the blade
polygon; blade count sets the point count (even blades give N points, odd
blades give 2N).

### 4.7 Stages 8-10 — the film layer

- **Halation.** Light passes through the emulsion, reflects off the film base,
  and re-exposes the emulsion from behind. Red penetrates deepest, so the red
  radius is the largest and the halo is orange-red. Modelled per channel:
  threshold in linear, convolve with a sum-of-Gaussians scatter profile, scale,
  add back. Anti-halation backing strength is a stock parameter.
- **Print vignette and gate weave.** A second, softer, non-optical falloff
  from the printing stage, plus a sub-pixel positional jitter.
- **Grain.** Uniform noise is wrong: real grain is signal-dependent, has a
  size distribution, and is spatially correlated. We use the Newson/Delon/
  Galerne stochastic Boolean-disc model — a Poisson field of discs with a
  log-normal radius distribution, per channel, with density driven by local
  exposure. It is resolution-independent, which matters because Photoshop
  documents vary enormously in size and grain must not change character when
  the same stock is applied at a different resolution.

### 4.8 Stage 11 — encode

Inverse of stage 1. When the document is 8-bit, apply triangular-PDF dither
before quantising; the wide, low-amplitude gradients produced by vignette and
glare band severely without it.

## 5. Lens and stock data

Effects are driven by JSON data files, so new glass needs no rebuild.

```jsonc
// lensdata/cooke-s4-32mm.lens
{
  "schema": 1,
  "name": "Cooke S4/i 32mm",
  "focal_mm": 32.0,
  "aperture": { "t_stop": 2.0, "blades": 9, "curvature": 0.15, "rotation_deg": 0 },
  "dispersion": { "model": "sellmeier", "B": [1.03961, 0.23179, 1.01047],
                                        "C": [0.00600, 0.02002, 103.56065] },
  "distortion": { "k1": -0.008, "k2": 0.001, "k3": 0.0, "p1": 0.0, "p2": 0.0 },
  "lateral_ca": { "k_l": 0.0021 },
  "vignette": { "natural_exp": 3.6,
                "pupil": { "r_entrance": 1.0, "r_exit": 0.92, "sep_norm": 0.35 } },
  "field": { "petzval": 0.42, "astig_sag": 0.18, "astig_tan": -0.11, "coma": 0.06,
             "spherical": 0.04 },
  "glare": { "amplitude": 0.004, "alpha": 2.1 },
  "flare": { "ghosts": [] }
}
```

```jsonc
// lensdata/kodak-5219.stock
{
  "schema": 1,
  "name": "Kodak Vision3 500T 5219",
  "halation": { "threshold": 1.0, "strength": 0.22,
                "radius_px_at_2k": { "r": 34.0, "g": 12.0, "b": 5.0 } },
  "grain": { "density": { "r": 0.9, "g": 1.0, "b": 1.35 },
             "radius_um": { "mu": 0.55, "sigma": 0.22 },
             "signal_exponent": 0.5 },
  "print_vignette": { "strength": 0.06, "exponent": 2.2 },
  "gate_weave": { "amplitude_px": 0.0 }
}
```

Radii are specified at a 2K reference width and scale with document width, so
a stock looks the same on a 2K and an 8K document.

### 5.1 `lensfit` — how a preset becomes real

A preset is only "accurate" if its numbers came from measurement. `lensfit`
takes chart photographs (or digitised manufacturer MTF and vignette charts)
and fits the parameters:

| Measurement | Input | Fits |
|---|---|---|
| Distortion | grid chart | `k1, k2, k3, p1, p2` |
| Vignette | flat evenly-lit field, per stop | `natural_exp`, pupil geometry |
| Lateral CA | high-contrast point grid | `k_l` |
| Field aberrations | slanted edges at 5+ field positions, sagittal and tangential | `petzval`, `astig_*`, `coma` |
| Glare | point light in a dark frame, HDR bracket | `amplitude`, `alpha` |

This is what separates the plugin from a slider pack.

## 6. Art-directable controls

Physics sets the default; the artist deviates on purpose. Jeong's composite
mapping supplies the mechanism: the field-dependence function `L(t, theta)`
and the axial function `A(t, theta)` are replaced by
`M(t, theta) = g(t, theta) . f(t, theta)`, built from atomic maps
(`t^2`, `t^3`, `exp(t)`, `1 - t`, `1/t`, `cos(2*pi*t)`).

In the UI this is not exposed as function composition. It appears as:

- **Amount** — global scale per effect, 0-200%, where 100% is the measured
  physical value. Above 100% is explicitly a deviation.
- **Falloff shape** — a small curve widget selecting the mapping.
- **Spectral equalizer** — a per-wavelength weight curve (Jeong Sec. 5.2),
  for pushing fringes toward magenta or cyan.

Every control's neutral position reproduces the measured lens exactly.

## 7. Photoshop adapter (`lens8bf`)

### 7.1 Selectors

`filterSelectorAbout`, `Parameters`, `Prepare`, `Start`, `Continue`, `Finish`.
`Parameters` shows the ImGui dialog; `Start`/`Continue` drive rendering.

### 7.2 Whole-image read, not tile-local

Halation and veiling glare have support spanning hundreds of pixels, so a
tile-local filter with a modest `inRect` overlap cannot produce a correct
result. The adapter therefore reads the **whole layer** through
`sPSChannelProcs->ReadPixelsFromLevel()` on `filterRecord->documentInfo`,
processes it once in `lenscore`, and streams the result out through the normal
tiled `outData` path. Memory is the cost: a 50MP RGBA float image is 800MB.
Mitigations: process in float16 where precision allows, and fall back to a
banded whole-width strip decomposition with PSF-radius overlap for documents
above a configurable size.

Photoshop delivers channels in **planar** order (all red, then all green);
`lenscore` wants interleaved. Conversion happens at the adapter boundary only.

### 7.3 Bit depth

`filterRecord->depth` is 8, 16, or 32.

- 8-bit: `0..255`.
- 16-bit: **`0..32768`, not `0..65535`.** This is the classic Photoshop trap
  and it is asserted in a unit test.
- 32-bit: float, already scene-referred, nominally `0..1` but legally above.

### 7.4 Colour management

Read the document profile through the ICC suite. Convert to a linear working
space (Rec.2020 linear, chosen for gamut headroom under saturated spectral
fringes), operate, and convert back. Refuse to guess: if no profile is
present, assume sRGB and say so in the dialog.

### 7.5 Selection masks

Respect `haveMask`. Leave `autoMask` at its default so Photoshop composites the
selection, except when feathering interacts badly with wide-support stages, in
which case we set `autoMask = false` and composite ourselves.

### 7.6 Scripting

Register parameters through `PIDescriptorParameters` with a unique event ID so
the filter is recordable in Actions and callable from UXP `batchPlay`. Parameter
keys mirror `params.hpp` field names one-to-one.

### 7.7 Threading

The Photoshop API is **not thread-safe**. All SDK calls happen on the calling
thread. Parallelism exists only inside `lenscore`, over a buffer the adapter
already owns.

## 8. UI (`lensui`)

Dear ImGui, one codebase, two rendering backends (Metal, D3D11). Layout: four
collapsible sections matching the effect groups, a lens and stock preset
picker, and a preview.

**Preview strategy.** Full-resolution spatially varying convolution is seconds,
not milliseconds. The preview renders a proxy: Photoshop supplies a
downsampled proxy through `advanceState()`, and all PSF radii scale by the
proxy ratio so the look is faithful at reduced resolution. Full quality runs
once, on commit. A "Preview quality" control trades proxy resolution against
latency, and the RGB three-band spectral tier is the preview default.

## 9. Testing

Correctness is asserted numerically against synthetic targets, in `lenscli`,
with no Photoshop involved.

**Targets:** flat field, slanted-edge grid, point grid, Siemens star, single
point in black.

**Assertions:**

| Test | Method | Passes when |
|---|---|---|
| Vignette curve | render flat field, measure radial mean | matches `.lens` model within 1% |
| MTF by field | slanted-edge MTF50 at 5 radii, sagittal and tangential | matches predicted defocus within 5% |
| Lateral CA | point grid, measure R-B fringe width vs radius | matches Eq. 8 within 0.5px |
| Energy | total linear luminance before and after stages 3-5 | conserved within 0.1% |
| Rotational symmetry | rotationally symmetric input | output symmetric within 1e-4 |
| Bit depth | 16-bit round trip | `0..32768` handled, no clipping |
| Resolution invariance | same stock at 2K and 8K | grain and halation scale-match |

Energy conservation and rotational symmetry are the two cheapest tests and
catch the most bugs: they fail immediately on kernel normalisation errors and
on any accidental axis-aligned assumption.

Golden images are stored as EXR with a perceptual-difference tolerance, not
byte equality.

## 10. Build and distribution

- CMake, C++20, one tree, `LENS_BUILD_8BF` off by default so `lenscore` and
  `lenscli` build without the gated SDK.
- macOS: universal binary (arm64 + x86_64), hardened runtime, codesign and
  notarize. PiPL via `Rez`.
- Windows: MSVC, `/arch:AVX2`. PiPL via `Cvntpipl.exe` then `rc.exe`.
- Install to the Photoshop `Plug-ins` folder; development uses a symlink.

**Prerequisite the project cannot satisfy itself:** the Photoshop Plug-In and
Connection SDK is gated behind an Adobe developer account and licence
acceptance at `developer.adobe.com/console` (Downloads, Photoshop). A human
must download it and set `PHOTOSHOP_SDK_ROOT`.

## 11. Phases

| # | Deliverable | Done when |
|---|---|---|
| 0 | Toolchain, CMake, PiPL, hello-world `.8bf` that inverts pixels | It loads and runs in Photoshop on both platforms |
| 1 | `lenscore` skeleton, colour stages, `lenscli`, golden harness | Energy and symmetry tests pass on a null pipeline |
| 2 | Lateral + axial CA, vignetting | Fringe-width and vignette-curve tests pass |
| 3 | PSF assembly, Zernike, pupil clip, EFF convolver | MTF-by-field test passes |
| 4 | `lens8bf` real adapter: whole-image read, bit depths, colour, mask | Full-res render matches `lenscli` output |
| 5 | ImGui dialog, proxy preview, descriptor scripting | Recordable in an Action, replayable |
| 6 | Film stage: halation, grain, print vignette | Resolution-invariance test passes |
| 7 | `lensfit` plus three measured presets | A fitted preset round-trips within tolerance |
| 8 | Universal binary, signing, notarization, installer | Installs clean on a machine with no dev tools |

Phases 1-3 are the accuracy work and carry the most risk of being wrong.
Phases 4-5 are the Photoshop work and carry the most risk of being slow.

## 12. Risks

| Risk | Impact | Mitigation |
|---|---|---|
| Full-res render takes seconds | Poor feel | Proxy preview; commit-time full quality; optional Metal/D3D compute path |
| 800MB whole-image buffer | Fails on large docs | float16 intermediates; banded strip fallback with PSF overlap |
| SDK gating blocks phase 0 | Total block | Phases 1-3 need no SDK; sequence them in parallel |
| ImGui on two platforms is more work than expected | Schedule | Ship macOS first; `lenscore` is unaffected |
| C++ host path proves wrong overall | Rework | `lenscore` compiles to WASM for a UXP panel; physics survives intact |
| Spectral uplift is inaccurate for saturated inputs | Wrong fringe colour | Rec.2020 linear working space; assert gamut in tests |

## 13. References

Held in `ref/`:

- Jeong, Lee, Kwon, Lee. *Expressive Chromatic Accumulation Buffering for
  Defocus Blur.* The Visual Computer, 2016. Sellmeier dispersion, axial and
  lateral CA, spectral importance sampling, expressive mapping functions.
- Cholewiak, Love, Srinivasan, Ng, Banks. *ChromaBlur.* ACM TOG 36(6), 2017 —
  **supplement only**. Disc/cylinder defocus kernel, per-channel focus shift.

Wanted, in priority order:

1. Aggarwal, Hua, Ahuja. *On cosine-fourth and vignetting effects in real
   lenses.* ICCV 2001. — stage 6.
2. Newson, Delon, Galerne. *A Stochastic Film Grain Model for
   Resolution-Independent Rendering.* CGF 2017. — stage 10.
3. Jakob, Hanika. *A Low-Dimensional Function Space for Efficient Spectral
   Upsampling.* Eurographics 2019. — stage 2.
4. Talvala, Adams, Horowitz, Levoy. *Veiling Glare in High Dynamic Range
   Imaging.* SIGGRAPH 2007. — stage 7.
5. Geigel, Musgrave. *A Model for Simulating the Photographic Development
   Process on Digital Images.* SIGGRAPH 1997. — stages 8-10.
6. Hullin, Hanika, Heidrich. *Polynomial Optics.* EGSR 2012. — stage 4.
7. Schuler, Hirsch, Harmeling, Schölkopf. *Non-stationary correction of
   optical aberrations.* ICCV 2011. — EFF convolver, Sec. 4.4.
8. Ritschel et al. *Temporal Glare.* Eurographics 2009. — diffraction star.
9. Cholewiak et al. *ChromaBlur* main paper. — completes the supplement.
10. Wu, Zheng, Hu, Xu. *Rendering realistic spectral bokeh due to lens stops
    and aberrations.* The Visual Computer, 2013. — cat's-eye validation.
