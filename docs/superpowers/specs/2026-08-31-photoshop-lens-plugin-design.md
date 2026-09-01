# Filmic Lens Plugin for Photoshop — Design

Date: 2026-08-31
Status: approved design, revised after reading the reference set

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
- Lens flare ghosts (inter-element reflections). Deferred; the data format
  reserves `flare.ghosts`, and Hullin's two-reflection path enumeration is the
  route in when we want it.
- Lens *correction* (removing aberration). We only add.
- Automatic lens detection from EXIF.
- Azimuthal asymmetry of vignetting from iris-blade decentring. Real but small
  (Aggarwal measured ~1.2% over the full 2*pi), and it costs rotational
  symmetry, which is one of our cheapest correctness tests.

## 3. Architecture

Four components. Only one is coupled to Adobe.

| Component | Language | Depends on | Purpose |
|---|---|---|---|
| `lenscore` | C++20 | nothing | All optics and film physics. Linear float image in, image out. |
| `lenscli` | C++20 | `lenscore` | Terminal renderer and test harness. |
| `lensaddon` | C++20 | `lenscore`, UXP Hybrid SDK | The `.uxpaddon` native module. Pixel buffers in and out, parameter marshalling. |
| `lenspanel` | JS/HTML | UXP, Spectrum | The panel UI. Controls, presets, preview, document I/O via the Imaging API. |

### 3.1 Why this split

`lenscore` has no Adobe headers, no UI, and no I/O. It is a pure function of
`(image, params) -> image`. Three consequences:

1. Optics can be developed and tuned from a terminal, with numeric assertions,
   without launching Photoshop. Tuning optics through a GUI is not viable.
2. The host-specific risk — SDK versions, pixel marshalling, signing — is
   confined to `lensaddon` and `lenspanel`. This paid for itself on
   2026-08-31, when the host route changed from an 8BF filter to a UXP Hybrid
   plugin and not one line of physics moved.
3. `lenscore` is the only thing that needs to be fast, so optimisation effort
   has one target.

### 3.2 File layout

```
lenscore/
  color/     transfer.hpp     EOTF / inverse EOTF, highlight recovery
             cie.hpp          CIE 1931 CMFs, spectrum -> XYZ -> linear RGB
             upsample.hpp     RGB -> spectrum (Jakob-Hanika sigmoid model)
  optics/    dispersion.hpp   Sellmeier n(lambda) + secondary-spectrum residual
             distortion.hpp   Brown-Conrady radial + tangential
             lateralca.hpp    per-wavelength magnification
             pupil.hpp        aperture polygon, cat's-eye clip, apodization ramp
             zernike.hpp      Zernike basis for the wavefront error
             psf.hpp          generalized pupil function -> PSF via one FFT
             vignette.hpp     natural falloff + mechanical clip
             glare.hpp        veiling glare spread function
  film/      sensitivity.hpp  per-layer spectral sensitivity
             mtf.hpp          emulsion and print internal scattering
             halation.hpp     back-reflection re-exposure
             hd.hpp           H&D characteristic curve, negative and print
             grain.hpp        Newson stochastic Boolean-disc grain
             printvig.hpp     print falloff, gate weave
  conv/      fft.hpp          real FFT
             eff.hpp          Efficient Filter Flow: patched varying convolution
  pipeline.hpp                stage ordering, the one public entry point
  params.hpp                  versioned POD parameter struct
lensdata/    *.lens  *.stock  JSON parameter sets
             rgb2spec/        precomputed spectral upsampling tables
lenscli/     main.cpp
lensaddon/   module.cpp  marshal.cpp  colour.cpp    the .uxpaddon native module
lenspanel/   manifest.json  index.html  main.js  ui/  the UXP panel
lensfit/     calibration tool: chart photos or a prescription -> .lens file
rgb2spec/    one-off generator for the spectral upsampling tables
tests/       golden/  metrics/{mtf.cpp,vignette.cpp,fringe.cpp,energy.cpp}
```

A file that grows past roughly 400 lines is a signal it holds two
responsibilities and should split.

## 4. Pipeline

The stage order is physical, not arbitrary. The lens forms an image; the
emulsion records it; the print reproduces it.

```
document pixels
  |
  1. decode          -> linear scene-referred RGB, highlight recovery
  2. spectral uplift -> per-pixel smooth spectrum (Jakob-Hanika)
  --- lens ---
  3. lateral CA      -> per-wavelength radial magnification + distortion
  4. varying PSF     -> one generalized-pupil FFT per (field radius, lambda):
                        defocus, secondary spectrum, field curvature,
                        astigmatism, coma, spherical, cat's-eye clip,
                        pupil apodization, blade polygon, diffraction
  5. integrate       -> spectrum back to linear RGB through CIE CMFs
  6. vignette        -> natural falloff x mechanical pupil clip
  7. veiling glare   -> broad low-amplitude scatter halo
  --- film ---
  8. sensitivity     -> per-layer (R/G/B dye) spectral response
  9. emulsion MTF    -> internal scattering within the emulsion
 10. halation        -> back-reflection re-exposure, red widest
 11. H&D (negative)  -> characteristic curve to density
 12. grain           -> stochastic Boolean-disc, density-dependent
 13. print           -> paper H&D cascade, print MTF, print vignette
 14. encode          -> density to transmission, back to document space
```

Stages 2-5 are one fused loop over wavelength bands, not four passes. Each
band is magnified, convolved, and accumulated into the running XYZ sum in a
single traversal, so the intermediate spectral image never exists in memory.

Every stage is individually switchable. A user who wants only chromatic
aberration pays for stages 1-5 and 14.

### 4.1 Stage 1 — decode and highlight recovery

Convert the document to linear scene-referred RGB using the document's ICC
profile transfer function.

**Highlight recovery is not optional.** Halation, glare, and defocus bokeh are
all driven by values far above diffuse white. A display-referred image clips
everything to 1.0, and blooming a clipped image produces the flat grey glow
that gives away every cheap lens plugin. We apply an inverse-knee expansion
above a threshold (default 0.85) that lifts near-white pixels to a
configurable peak (default 8.0 linear). Users in 32-bit scene-referred float
switch it off, because their highlights are real.

### 4.2 Stage 2 — spectral uplift

The optics are spectral but the input is RGB. We use the Jakob-Hanika model:

```
f(lambda) = S(c0*lambda^2 + c1*lambda + c2)
S(x)      = 1/2 + x / (2*sqrt(1 + x^2))
```

Six floating-point operations per wavelength (two FMAs and an rsqrt),
trivially vectorised. Coefficients `(c0,c1,c2)` come from a trilinear lookup
into a `3 x 64^3` table: find the largest RGB component `i`, normalise the
other two by it, and index `table[i]` by `(c[(i+1)%3]/max, c[(i+2)%3]/max,
max)`. The `alpha` axis is discretised as `smoothstep(smoothstep(i/63))`,
which concentrates resolution where the coefficients move fastest.

**Radiance, not reflectance.** The model produces bounded reflectance spectra.
Our pixels are radiance and legally exceed 1.0. We therefore upsample the
*normalised* colour and carry `max(r,g,b)` as a separate scalar multiplier.
HDR values pass through untouched.

**Table generation.** The tables are ours to build: `rgb2spec` runs a
Gauss-Newton fit per grid cell against the CIE objective, ~10-60s per colour
space, and the result is committed as a binary asset. We ship tables for
sRGB/Rec.709 and Rec.2020 linear.

**Two caveats to record.** The model minimises error projected onto the CIE
colour matching functions. It is valid because we integrate back through those
same CMFs. Retargeting to a camera's response curves would require refitting,
not merely a different `3x3` matrix.

Second, the reconstruction treats the spectrum as a reflectance under an
**equal-energy illuminant**, while the Rec.2020 matrices are referenced to
**D65**. We reconcile the two with a diagonal white-point scale in XYZ, applied
once inside `spectrumToRec2020`. It is exact at the white point and an
approximation elsewhere — a true chromatic adaptation would work in
cone-response space, and a more physical treatment would integrate against a
real D65 spectral power distribution so that reflectance times illuminant is
radiance.

This cannot disturb round-trip identity, because the fit optimises against the
same function that applies the scale. What it does change is *which metamer* we
reconstruct for a given RGB — and since different metamers disperse
differently, it is a second-order limitation on fringe colour for saturated
subjects. Accepted deliberately: a D65 table costs roughly eighty more
constants, and no test or downstream consumer currently distinguishes the two.
Revisit if a later stage needs physically accurate spectral radiance rather
than round-trip fidelity.

Two quality tiers:

- **Spectral** (default, N = 11 bands over 380-780nm) with **spectral
  importance sampling** through the inverse CDF of the equalizer weights.
  Uniform sampling at this count bands visibly; importance sampling does not.
- **RGB** (preview, N = 3 at 650 / 510 / 475nm). Three times faster, banded on
  strong dispersion. Preview only.

### 4.3 Stage 3 — lateral chromatic aberration

Per-wavelength magnification about the optical centre:

```
m(lambda, t) = 1 + k_l * (lambda_hat - lambda) * L(t, theta)
```

`t` is normalised field radius, `lambda_hat` = 650nm is the reference so
magnification never goes negative, `k_l` comes from the `.lens` file.
`L(t, theta) = t` is the physical baseline; the expressive layer replaces `L`
with a composite mapping (Sec. 6).

Brown-Conrady radial and tangential distortion applies in the same resample,
so the image is warped once, not twice. Aggarwal's point-grid photographs
confirm the two are entangled in reality: PSF centres are measurably further
apart at the periphery than at the centre.

### 4.4 Stage 4 — the PSF, from one Fourier transform

This is the hard stage and the one that makes or breaks the look.

**The unifying formulation.** Rather than assembling a kernel from separate
disc, ellipse, star and clip terms, we compute the PSF the way physical optics
defines it — as the squared modulus of the Fourier transform of the
generalized pupil function:

```
PSF(u, v; lambda, t) = | FFT{ A(x, y; t) * exp( i * 2*pi/lambda * W(x, y; t) ) } |^2
```

- `A` is the **pupil amplitude**: the aperture polygon (blade count, per-blade
  curvature, rotation), intersected with the cat's-eye clip, multiplied by the
  apodization ramp.
- `W` is the **wavefront error** in a Zernike basis: defocus `Z(2,0)`,
  astigmatism `Z(2,+/-2)`, coma `Z(3,+/-1)`, spherical `Z(4,0)`.

Everything falls out of this single calculation. Diffraction, the aperture
star, the soft-edged disc, the lemon-shaped off-axis kernel, the comet tail,
the intensity gradient across the pupil — all of it, with no special cases and
no risk of the parts contradicting each other. It is also the formulation
ChromaBlur uses.

**Pupil apodization — the term that is easy to miss.** An off-axis point does
*not* illuminate the pupil uniformly. Nonlinear refraction through the
elements skews the distribution, an effect Aggarwal calls pupil aberration and
measured at up to **31% variation across the pupil at a 10 degree field angle
on a 16mm lens**. It shifts the pupil centroid, which makes the fall-off
surface depend on aperture in ways `cos^4` and mechanical vignetting together
cannot explain. We model it as Aggarwal's own planar approximation: a linear
intensity ramp across the pupil, radially oriented, slope proportional to
field angle. One extra parameter, real measured effect.

**Defocus radius.** Two contributions sum:

1. **Longitudinal chromatic aberration**, below.
2. **Petzval field curvature**: the focal surface is a sphere, so defocus grows
   with `t^2` even for a flat subject.

**Secondary spectrum — the correction that matters most.** Modelling
dispersion with a bare Sellmeier `n(lambda)` describes a *single uncorrected
element*. Its focal length varies monotonically with wavelength, producing a
red-to-blue rainbow. **Real photographic lenses are achromats or apochromats**
— corrected so that focus coincides at two (or three) wavelengths. Their
residual longitudinal error is the *secondary spectrum*: a curve that crosses
zero at the correction wavelengths and bulges between and beyond them. It is
why real lenses fringe green and magenta rather than red and blue.

So: take Sellmeier for the shape, then subtract the linear (achromat) or
quadratic (apochromat) fit in `1/lambda^2` that zeroes the error at the
correction wavelengths, and scale by a measured residual magnitude.

```
d_f(lambda) = d_f_ref * ( 1 + secondary(lambda; correction_nm[], residual) )
```

The `.lens` file names the correction wavelengths, so a cheap achromatic
doublet, a modern apochromatic cine prime, and an uncorrected vintage singlet
are three settings of one model rather than three code paths.

**Depth, and the decision.** A 2D photograph has no depth map, but defocus
needs an object distance. Resolved as follows:

- **Default**: the subject lies on the focal surface, so defocus comes entirely
  from field curvature and grows with `t^2`. Physically correct for a flat
  subject, and needs no extra input.
- **Optional**: the user designates a layer or channel as a depth map, and true
  focus-plane defocus appears. A "Depth source" dropdown, default None.

Stated explicitly because it is the one place the model is under-determined by
the input.

**Fast evaluation.** A naive per-pixel elliptical gather at a 20px radius over
50 megapixels is roughly 2e10 operations per wavelength band. Not viable.
We use **Efficient Filter Flow**: cover the image with overlapping patches,
assign one PSF per patch, and evaluate

```
y = sum_r  L_r^T F^H Diag(F P f^(r)) F K_r Diag(w^(r)) x
```

where `K_r` crops patch `r`, `Diag(w)` applies its window, `F` is the DFT, `P`
zero-pads the kernel, and `L_r^T` accumulates the result. The windows form a
partition of unity, so patch PSFs interpolate smoothly. Cost is `O(N log P)`.
Schuler used an 18x27 support grid on a 21MP image and found bilinear
interpolation between measured PSFs sufficient; our PSFs are analytic and
smoother still.

**Do not normalise the PSF.** Its total energy *is* the optical vignetting —
the cat's-eye clip removes light exactly as the barrel does. Schuler makes the
same point from the correcting side: relaxing the filters from `sum = 1` to
`sum <= 1` folds vignetting into the same operator for free. Stage 6 therefore
adds only the terms the pupil does not already contain.

**The rotational-covariance collapse.** The PSF's shape depends only on field
radius `t`; at a different angle it is the same kernel, rotated. So the PSF
table is a **1D function of `t`**, patches are laid out in polar rings, one
kernel is computed per ring and rotated per patch. This cuts PSF construction
by the number of angular divisions and makes exact rotational symmetry a
structural guarantee, which the test suite asserts directly.

**Alternate derivation path.** When the input is a lens *prescription* (radii,
thicknesses, glasses) rather than photographs, Hullin's polynomial optics gives
the coefficients. A degree-3 Taylor expansion of the ray mapping decomposes
into exactly the classical groups, in reduced ray coordinates
`(r_x, r_y, d_x, d_y)`:

```
distortion:      r'_x = -r_x  +/-  a_radial * r_x (r_x^2 + r_y^2)
coma:            r'_x = -r_x  +/-  a_coma   * (3 r_x d_x^2 + 2 r_y d_x d_y + r_x d_y^2)
field curvature: r'_x = -r_x  +/-  a_curv   * (3 r_x^2 d_x + 2 r_x r_y d_y + r_y^2 d_x)
spherical:       r'_x = -r_x  +/-  a_sph    * d_x (d_x^2 + d_y^2)
```

Degree 3 suffices; degree 5 changes the result marginally at real cost. Note
two cautions Hullin records: the expansion is centred on the axis and degrades
at the extreme periphery, which is another argument for our per-ring PSFs; and
Seidel's aberration names do not map exactly onto this basis, though as they
put it, these phenomena do not occur in isolated form in real optics either.
`lensfit` uses this path to produce Zernike coefficients; the renderer stays
on the Fourier formulation.

### 4.5 Stage 6 — vignetting

The pupil clip in stage 4 already removed the optical component. Two terms
remain:

1. **Natural falloff**, `cos^n(theta)` with `theta = atan(r/f)`. `n = 4` is the
   textbook ideal. **The exponent must be free**: Aggarwal measured three real
   lenses across six apertures and the fall-off curves *"do not coincide with
   the cosine-fourth curves for any of the aperture settings"* — not even
   stopped fully down, where vignetting is absent and `cos^4` should be exact.
2. **Mechanical vignetting**, the overlap area of the aperture circle with the
   entrance and exit pupil circles displaced by field angle. It must **vanish
   above roughly f/4** and dominate at f/2.8 and wider; that threshold is a
   measured result, and the model reproduces it rather than assuming it.

A consistency test asserts that the energy lost by the stage-4 pupil clip
matches the mechanical term computed independently here.

### 4.6 Stage 7 — veiling glare

Scattering inside the barrel and off the sensor. Talvala's measurements give
us both the shape and the magnitude:

- The glare spread function is `gsf = delta * (1-k) + k * b`, energy-conserving,
  where `b` is the broad component. Applying it forward is one convolution.
- `b` is **not** a single Gaussian. Glare has a steep high-frequency falloff
  near a bright point *and* a smooth low-frequency floor across the whole
  frame. Talvala fit theirs as a combination of Gaussians of different widths;
  we do the same, with three components.
- **Magnitude is specified in f-stops below the source**, which is how it was
  measured and how a lens person reasons: a good DSLR prime sits about **20
  stops down**, a compact camera about **16**. Ghost count and glare level both
  scale with the number of air-glass surfaces.

Talvala notes the GSF is not truly shift-invariant — its shape, size and even
colour vary across the field. That defeats their *inverse* problem. For
forward synthesis a shift-invariant GSF is a defensible approximation, and we
expose an optional field-radius scale for users who want the corner falloff.

The aperture diffraction star is not modelled here. It already came out of
stage 4's pupil transform, at the correct brightness relative to the halo.

### 4.7 Stages 8-13 — the film

Not a look-up curve. Following Geigel and Musgrave, film is a cascade of
measurable sensitometric stages, run **twice** — once for the negative and
once for the print. That double cascade is where film's characteristic
contrast actually comes from.

- **Spectral sensitivity (8).** Each dye layer has its own response curve.
  Exposure per layer is the spectrum integrated against that curve. This is
  what makes a tungsten stock render daylight the way it does.
- **Emulsion MTF (9).** Light scatters *between grains inside* the emulsion,
  degrading resolution independently of grain itself. A typical emulsion MTF
  falls to 0.1 around 100-150 cycles/mm. Modelled as a per-layer low-pass.
- **Halation (10).** Light crosses the emulsion, reflects off the film base,
  and re-exposes from behind. Red penetrates deepest, so the red radius is
  largest and the halo reads orange-red. Threshold in linear, convolve with a
  per-layer sum-of-Gaussians scatter profile, add back. Anti-halation backing
  strength is a stock parameter; stocks with effective backing get a small
  radius rather than a special case.
- **H&D characteristic curve (11).** The density response: toe, straight-line
  section whose slope is gamma, shoulder, and the solarisation rollover beyond
  it. Per layer. Specified as sampled `(log E, D)` points read straight off a
  datasheet, linearly interpolated. **This stage was missing from the first
  draft and it is roughly half of what people mean by a filmic look.**
- **Grain (12).** Newson's stochastic Boolean model. Grains are discs of
  log-normally distributed radius placed by an inhomogeneous Poisson process
  whose local intensity is derived from the image itself:

  ```
  lambda(y) = 1 / (pi * (mu_r^2 + sigma_r^2)) * log( 1 / (1 - u_tilde(floor(y))) )
  ```

  This choice makes the mean area covered equal the local grey level, so
  **contrast is preserved exactly** and grain is signal-dependent by
  construction rather than by a fudge. Rendering is Monte Carlo: `N` sample
  offsets drawn from `N(0, sigma^2 I)`, averaged; `sigma` is the perceptual
  low-pass.

  **Use the pixel-wise algorithm, not the grain-wise one.** Storing grain
  positions costs 35GB at 2048x2048 with a small radius; a 50MP Photoshop
  document is out of the question. The pixel-wise variant generates grains
  on the fly from a per-cell PRNG keyed by cell coordinates, so memory depends
  only on output size and `N`, and it parallelises cleanly.

  Datasheets quote RMS granularity, not disc radii. `lensfit` converts via
  Selwyn granularity `G = sqrt(2A) * sigma`, with granularity rising roughly as
  the cube root of density.
- **Print (13).** The negative is printed onto paper: a second H&D cascade with
  the paper's own gamma and speed, a print MTF, and the print/gate falloff.
  Speed follows `SP = K/E_m` with `K = 0.8, m = 0.1` for films and
  `K = 1000, m = 0.6` for papers.

### 4.8 Stage 14 — encode

Density to transmission (`D = log10(1/T)`), then the inverse of stage 1. When
the document is 8-bit, apply triangular-PDF dither before quantising; the
wide, low-amplitude gradients from vignette and glare band badly without it.

## 5. Lens and stock data

Effects are driven by JSON, so new glass needs no rebuild.

```jsonc
// lensdata/cooke-s4-32mm.lens
{
  "schema": 2,
  "name": "Cooke S4/i 32mm",
  "focal_mm": 32.0,
  "aperture": { "t_stop": 2.0, "blades": 9, "curvature": 0.15, "rotation_deg": 0 },

  // Sellmeier gives the shape; correction_nm and residual make it a real lens.
  "dispersion": { "model": "sellmeier",
                  "B": [1.03961, 0.23179, 1.01047],
                  "C": [0.00600, 0.02002, 103.56065],
                  "correction_nm": [486.1, 656.3],   // achromat: F and C lines
                  "residual": 0.00042 },

  "distortion":  { "k1": -0.008, "k2": 0.001, "k3": 0.0, "p1": 0.0, "p2": 0.0 },
  "lateral_ca":  { "k_l": 0.0021 },

  "pupil": { "r_entrance": 1.0, "r_exit": 0.92, "sep_norm": 0.35,
             "apodization_slope": 0.31 },          // Aggarwal ramp, per radian

  "vignette":  { "natural_exp": 3.6, "mech_vanish_fstop": 4.0 },

  // Zernike wavefront coefficients, in waves, as functions of field radius t
  "wavefront": { "petzval": 0.42, "astig": 0.18, "coma": 0.06, "spherical": 0.04 },

  "glare": { "floor_stops": 20.0,
             "gaussians": [ {"sigma_norm": 0.004, "weight": 0.55},
                            {"sigma_norm": 0.03,  "weight": 0.30},
                            {"sigma_norm": 0.25,  "weight": 0.15} ] },

  "flare": { "ghosts": [] }
}
```

```jsonc
// lensdata/kodak-5219.stock
{
  "schema": 2,
  "name": "Kodak Vision3 500T 5219",
  "layers": {
    "r": { "sensitivity_nm": [[580,0.1],[620,0.9],[660,1.0],[700,0.4]],
           "hd_curve": [[-3.0,0.08],[-1.5,0.35],[0.0,1.10],[1.2,1.95],[2.0,2.20]],
           "mtf_cycles_per_mm_at_50": 78,
           "halation_radius_px_at_2k": 34.0 },
    "g": { "...": "..." },
    "b": { "...": "..." }
  },
  "halation":  { "threshold": 1.0, "strength": 0.22, "backing": 0.6 },
  "grain":     { "rms_granularity": 11.0, "radius_um": {"mu": 0.55, "sigma": 0.22},
                 "density": { "r": 0.9, "g": 1.0, "b": 1.35 } },
  "print":     { "paper_gamma": 1.67, "paper_speed": 250,
                 "vignette_strength": 0.06, "vignette_exponent": 2.2 },
  "gate_weave": { "amplitude_px": 0.0 }
}
```

Radii are given at a 2K reference width and scale with document width, so a
stock looks the same on a 2K and an 8K document.

### 5.1 `lensfit` — how a preset becomes real

A preset is only "accurate" if its numbers came from measurement. `lensfit`
takes either chart photographs or a lens prescription and fits the parameters:

| Measurement | Input | Fits |
|---|---|---|
| Distortion | grid chart | `k1, k2, k3, p1, p2` |
| Vignette | flat field, **per stop** | `natural_exp`, pupil geometry, `mech_vanish_fstop` |
| Lateral CA | high-contrast point grid | `k_l` |
| Secondary spectrum | through-focus series on a narrowband target | `correction_nm`, `residual` |
| Wavefront | slanted edges at 5+ radii, sagittal and tangential | `petzval`, `astig`, `coma`, `spherical` |
| Pupil apodization | defocused point grid, measure the intensity ramp across each blur disc | `apodization_slope` |
| Glare | point light in a dark frame, HDR bracket | `floor_stops`, Gaussian widths and weights |
| From a prescription | radii, thicknesses, glasses | all of the above via Hullin degree-3 |
| Stock | datasheet H&D, MTF, RMS granularity | layer curves, grain radii |

Aggarwal's point-grid protocol is directly reusable: photograph a grid of
point sources, and each blur patch *is* the PSF at that field position.

## 6. Art-directable controls

Physics sets the default; the artist deviates on purpose. The field-dependence
function `L(t, theta)` and the axial function `A(t, theta)` may be replaced by
a composite `M = g . f` built from atomic maps (`t^2`, `t^3`, `exp t`, `1-t`,
`1/t`, `cos 2*pi*t`).

In the UI this is not exposed as function composition. It appears as:

- **Amount** — per-effect scale, 0-200%, where 100% is the measured physical
  value. Above 100% is explicitly a deviation.
- **Falloff shape** — a small curve widget selecting the mapping.
- **Spectral equalizer** — a per-wavelength weight curve, for pushing fringes
  toward magenta or cyan. It also drives the importance sampling in stage 2,
  so emphasising a band costs no extra noise.

Every control's neutral position reproduces the measured lens exactly.

## 7. Photoshop host (`lensaddon` + `lenspanel`)

**Decision, 2026-08-31: UXP Hybrid, not an 8BF filter.** Adobe shipped UXP
Hybrid plugins for Photoshop in v24.2 — a UXP panel front-end plus a native
C++ `.uxpaddon` back-end, with a Node-API-style surface. This is the supported
answer to the gap recorded earlier in this project: that there is no supported
UXP UI for a C++ filter plugin. Taking it deletes four things the 8BF route
demanded — the Dear ImGui dialog, its Metal and D3D11 backends, the PiPL
resource chain (`Rez`, `Cvntpipl.exe`, `rc.exe`), and the per-platform IDE
projects. `lenscore` is unchanged either way, which is why the decision was
cheap to take this late.

### 7.1 What a `.uxpaddon` actually is

A plain dynamic library — `.dylib` on macOS, `.dll` on Windows — with the file
extension renamed to `.uxpaddon`. JavaScript loads it with
`require("lens.uxpaddon")`. Because it is an ordinary shared library, **one
CMake target builds it on both platforms**; the SDK's Xcode and Visual Studio
projects are samples, not a requirement.

Minimum host version is Photoshop **24.2.0**, declared in `manifest.json`.

### 7.2 Division of labour

- **`lenspanel`** (UXP, JS/HTML/Spectrum): all UI, preset loading, document
  selection, progress, and reading and writing pixels through the Imaging API
  (`imaging.getPixels` / `imaging.putPixels`).
- **`lensaddon`** (C++): receives a linear float buffer plus a parameter
  object, calls `lens::render`, returns the result. It holds no UI and no
  Photoshop state.

The JS-to-native boundary is cheap to cross, so the panel may call the addon
per preview at proxy resolution and once more on commit.

### 7.3 What this route costs

Two things the 8BF filter would have given us for free, both accepted
deliberately:

1. **No Filter menu entry and no Smart Filter re-editability.** The plugin is
   a panel. Re-editability is recovered by writing results to a new layer and
   storing the parameter set in that layer's metadata, so the panel can reload
   and re-render — an approximation of a Smart Filter, not the real thing.
2. **No host-provided tiling or proxy.** The panel asks the Imaging API for
   the pixels it wants, so proxy generation and any banding for large
   documents are ours to implement rather than the host's.

If Smart Filter support later proves necessary, an 8BF filter can be added
beside the panel: it would link the same `lenscore` and need only the adapter
and dialog work this decision defers. The SDK's `PIUXPSuite.h` suggests the
two can be bridged, but that is unverified and not relied upon.

### 7.4 Colour and precision

The panel requests pixels in a known working space and bit depth rather than
inheriting whatever the document happens to be. Convert to **Rec.2020 linear**
at the boundary, operate, convert back. If the document has no profile, assume
sRGB and say so in the panel rather than guessing silently.

Photoshop's 16-bit values run **`0..32768`, not `0..65535`** — the classic
trap, asserted in a unit test wherever the addon touches 16-bit data.

### 7.5 Threading

`lenscore` may use all cores over a buffer the addon owns. UXP addon entry
points follow the SDK's threading contract (`UxpTask`), and no Photoshop or
UXP object is touched from a worker thread.

## 8. Panel UI (`lenspanel`)

Spectrum web components, matching Photoshop's own look with no styling work.
Collapsible sections per effect group, a lens and stock preset picker, and a
preview.

**Preview strategy.** Full-resolution spatially varying convolution plus Monte
Carlo grain is seconds, not milliseconds. The panel fetches a downsampled
proxy through the Imaging API, scales every PSF and halation radius by the
proxy ratio, drops grain `N` to 8, and uses the three-band spectral tier. Full
quality runs once, on commit.

## 9. Testing

Correctness is asserted numerically against synthetic targets, in `lenscli`,
with no Photoshop involved.

**Targets:** flat field, slanted-edge grid, point grid, Siemens star, single
point in black, uniform grey patch ladder.

**Assertions:**

| Test | Method | Passes when |
|---|---|---|
| Vignette curve | flat field, radial mean, per stop | matches `.lens` model within 1% |
| Vignette consistency | stage-4 pupil energy vs stage-6 mechanical term | agree within 1% |
| Mechanical vanishing | flat field at f/8 vs f/2 | mechanical term ~0 above `mech_vanish_fstop` |
| MTF by field | slanted-edge MTF50 at 5 radii, sagittal and tangential | matches predicted defocus within 5% |
| Lateral CA | point grid, R-B fringe width vs radius | matches the model within 0.5px |
| Secondary spectrum | through-focus on narrowband targets | focus error zero at `correction_nm` |
| Diffraction star | point source, small aperture | point count = blades (even) or 2x (odd) |
| Energy | total linear luminance, stages 3-5, PSF normalised | conserved within 0.1% |
| Rotational symmetry | rotationally symmetric input | output symmetric within 1e-4 |
| Grain contrast | grey ladder before and after grain | patch means unchanged within 0.5% |
| H&D round trip | grey ladder through negative and print | matches the datasheet curve within 2% |
| Bit depth | 16-bit round trip | `0..32768` handled, no clipping |
| Resolution invariance | same stock at 2K and 8K | grain and halation scale-match |

Energy conservation, rotational symmetry, and grain contrast are the three
cheapest tests and catch the most bugs: they fail immediately on kernel
normalisation errors, on any accidental axis-aligned assumption, and on a
grain model that has drifted off Newson's density relation.

Golden images are stored as PFM (plain float raster, ~30 lines of reader and
writer) with a perceptual-difference tolerance, not byte equality. The format
is chosen so that `lenscore`, `lenscli` and the whole test suite carry **zero
third-party code** — the accuracy work has no dependency to blame. `lensdata`
JSON parsing is the sole exception and lives outside `lenscore`.

## 10. Build and distribution

- CMake, C++20, one tree. `LENS_BUILD_ADDON` off by default so `lenscore` and
  `lenscli` build with no SDK present at all.
- The addon is one `SHARED` CMake target whose output extension is set to
  `.uxpaddon`. The same target definition serves both platforms.
- macOS: universal binary (arm64 + x86_64), hardened runtime, codesign and
  notarize. Windows: MSVC, `/arch:AVX2`.
- Development loads the plugin through UXP Developer Tools (UDT) by pointing
  it at `manifest.json`; no install step and no symlink into the app bundle.

**Prerequisite already satisfied:** the UXP Hybrid Plugin SDK (v6.5.0) and the
Photoshop Plug-In and Connection SDK (2026, mac, v2) are both unpacked under
`ref/`, git-ignored. Only the hybrid SDK is on the critical path; the C++ SDK
is retained for its documentation and for the deferred 8BF option.

## 11. Phases

| # | Deliverable | Done when |
|---|---|---|
| 0 | Toolchain, CMake, hello-world `.uxpaddon` that inverts pixels, loaded via UDT | The panel calls it and Photoshop shows the result |
| 1 | `lenscore` skeleton, colour stages, `rgb2spec`, `lenscli`, golden harness | Energy and symmetry tests pass on a null pipeline |
| 2 | Dispersion with secondary spectrum, lateral CA, vignetting | Fringe-width, vignette-curve and mechanical-vanishing tests pass |
| 3 | Zernike wavefront, pupil model, PSF-by-FFT, EFF convolver | MTF-by-field and diffraction-star tests pass |
| 4 | `lensaddon` real bridge: Imaging API buffers, bit depths, colour, selection | Full-res render matches `lenscli` output |
| 5 | `lenspanel`: Spectrum UI, proxy preview, preset picker, layer-metadata round trip | Re-opening a rendered layer restores its exact settings |
| 6 | Film stage: sensitivity, MTF, halation, H&D, grain, print | Grain-contrast, H&D and resolution-invariance tests pass |
| 7 | `lensfit` plus three measured lens presets and two stocks | A fitted preset round-trips within tolerance |
| 8 | Universal binary, signing, notarization, installer | Installs clean on a machine with no dev tools |

Phases 1-3 and 6 are the accuracy work and carry the most risk of being wrong.
Phases 4-5 are the Photoshop work and carry the most risk of being slow.

The host decision is recorded in section 7; phases 0, 4, 5 and 8 changed with it,
and phases 1-3, 6 and 7 — the accuracy work — did not move at all.

## 12. Risks

| Risk | Impact | Mitigation |
|---|---|---|
| Monte Carlo grain is slow. Newson report 37.8s for 2048x2048 at N=800; a 50MP document extrapolates to minutes | Unusable commit times | `N` is the quality knob, not a constant: 8 for preview, 64 for commit. Embarrassingly parallel. Pixel-wise algorithm keeps memory flat |
| Full-res PSF render takes seconds | Poor feel | Proxy preview; commit-time full quality; optional Metal/D3D compute path |
| 800MB whole-image buffer | Fails on large docs | float16 intermediates; banded strip processing with PSF-radius overlap |
| Panel and addon version drift | Silent breakage | Addon exports a schema version; the panel refuses to load a mismatch |
| Imaging API round trip is slow on large documents | Poor feel | Proxy for preview; band the commit pass; the JS-to-native boundary itself is cheap |
| Losing Smart Filter re-editability proves unacceptable | Rework | Layer-metadata round trip first; an 8BF filter can be added later against the same `lenscore` |
| Spectral uplift inaccurate for saturated inputs | Wrong fringe colour | Rec.2020 linear working space; round-trip `deltaE` asserted in tests |
| Stock datasheet curves are hard to source | Weak presets | H&D and MTF are digitisable from published datasheet plots; `lensfit` accepts sampled points |

## 13. References

All held in `ref/`.

| Paper | Used for |
|---|---|
| Jeong, Lee, Kwon, Lee. *Expressive Chromatic Accumulation Buffering for Defocus Blur.* The Visual Computer, 2016. | Sellmeier dispersion, axial and lateral CA, spectral importance sampling, expressive mapping functions (Secs. 4.2, 4.3, 4.4, 6) |
| Cholewiak, Love, Srinivasan, Ng, Banks. *ChromaBlur.* ACM TOG 36(6), 2017, plus supplement. | PSF as the squared modulus of the transformed complex aperture function; disc not Gaussian; per-channel focus shift (Sec. 4.4) |
| Aggarwal, Hua, Ahuja. *On cosine-fourth and vignetting effects in real lenses.* ICCV 2001. | Pupil aberration and the apodization ramp; free `cos^n` exponent; mechanical vignetting vanishing above f/4; the point-grid PSF protocol (Secs. 4.4, 4.5, 5.1) |
| Jakob, Hanika. *A Low-Dimensional Function Space for Efficient Spectral Upsampling.* Eurographics 2019. | RGB to spectrum: the sigmoid-polynomial model, table layout and lookup (Sec. 4.2) |
| Schuler, Hirsch, Harmeling, Schölkopf. *Non-stationary correction of optical aberrations.* ICCV 2011. | Efficient Filter Flow forward operator; vignetting as relaxed filter normalisation; support-grid density (Sec. 4.4) |
| Talvala, Adams, Horowitz, Levoy. *Veiling Glare in High Dynamic Range Imaging.* SIGGRAPH 2007. | Glare spread function shape and magnitude in f-stops; multi-width Gaussian fit; shift-variance caveat (Sec. 4.6) |
| Newson, Delon, Galerne. *A Stochastic Film Grain Model for Resolution-Independent Rendering.* CGF 2017. | Inhomogeneous Boolean grain model, the density relation, the pixel-wise algorithm (Sec. 4.7) |
| Geigel, Musgrave. *A Model for Simulating the Photographic Development Process on Digital Images.* SIGGRAPH 1997. | The whole film stage: sensitivity, emulsion MTF, H&D curves, negative-print cascade, Selwyn granularity (Sec. 4.7) |
| Hullin, Hanika, Heidrich. *Polynomial Optics.* EGSR 2012. | Degree-3 aberration terms in ray space; prescription-to-coefficient path for `lensfit`; ghost enumeration for future flare (Secs. 4.4, 5.1) |
