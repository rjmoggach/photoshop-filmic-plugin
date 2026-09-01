# Lenscore Optics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `lenscore` (host-free optics library) and `lenscli` (terminal renderer and test harness) through Phase 3 of the spec, so that physically accurate chromatic aberration, vignetting, and spatially varying edge blur can be rendered and numerically verified without Photoshop.

**Architecture:** A pure C++20 static library with no third-party code, exposing one entry point `lens::render(Image, Params) -> Image`. Optics are spectral: RGB is upsampled to a smooth spectrum, each wavelength band is magnified and convolved with a point spread function computed as the squared modulus of the Fourier transform of a generalized pupil function, and the bands are integrated back through the CIE colour matching functions. The spatially varying convolution uses Efficient Filter Flow over polar patch rings, exploiting the fact that the PSF depends only on field radius.

**Tech Stack:** C++20, CMake 3.24+, doctest (tests only, via FetchContent), nlohmann/json (in `lensdata` only, via FetchContent). `lenscore`, `lenscli` and every test carry zero third-party code.

**Spec:** `docs/superpowers/specs/2026-08-31-photoshop-lens-plugin-design.md`

## Global Constraints

- C++20. No compiler extensions (`CXX_EXTENSIONS OFF`).
- `lenscore` links **nothing**. No Adobe headers, no JSON, no image libraries, no test framework. JSON loading lives in `lensdata`, outside the core.
- Working colour space is **Rec.2020 linear**. All internal maths is scene-referred linear float, legally above 1.0.
- Spectral range **380-780nm**. Default band count **N = 11**, importance-sampled. Preview tier is **N = 3** at 650 / 510 / 475nm.
- Reference wavelength for lateral CA is **`lambda_hat = 650nm`**, so magnification never goes negative.
- **The PSF is never normalised.** Its total energy is the optical vignetting term. Any code that divides a PSF by its sum is a bug.
- Any source file that exceeds **400 lines** must be split by responsibility.
- Golden images are **PFM**. Comparison is by metric tolerance, never byte equality.
- Every public function that takes an angle takes **radians**. Every function that takes a wavelength takes **nanometres**. Every field radius `t` is **normalised so that the image corner is 1.0**.
- Tests run under `ctest`. A task is not done until `ctest --output-on-failure` is green.
- **Never use `M_PI`.** It is a POSIX extension, not standard C++, and MSVC omits it
  unless `_USE_MATH_DEFINES` is set — while this project builds with `CXX_EXTENSIONS OFF`
  and targets Windows as well as macOS. Any test needing pi declares it locally:
  `static constexpr double kPi = 3.14159265358979323846;` in tests, and in library code
  include `lenscore/constants.hpp`, which defines `kPi` and `kTwoPi` ONCE in
  `lens::` — never redefine either in a header, or two headers in the same namespace
  collide at the point some third file includes both.

## File Structure

| File | Responsibility |
|---|---|
| `CMakeLists.txt` | Top-level build, options, FetchContent for doctest |
| `lenscore/image.hpp` | `Image` value type: interleaved float RGB, width, height |
| `lenscore/constants.hpp` | `kPi`, `kTwoPi`. ONE definition, shared. Six local copies caused an ODR collision. |
| `lenscore/pfm.hpp` | PFM read and write. Test-harness I/O, header-only |
| `lenscore/color/transfer.hpp` | sRGB EOTF, inverse EOTF, highlight-recovery knee |
| `lenscore/color/cie.hpp` | CIE 1931 CMF table, spectrum to XYZ, XYZ to Rec.2020 linear |
| `lenscore/color/upsample.hpp` | Jakob-Hanika sigmoid model evaluation and table lookup |
| `lenscore/optics/dispersion.hpp` | Sellmeier `n(lambda)`, focal length, secondary-spectrum residual |
| `lenscore/optics/distortion.hpp` | Brown-Conrady radial and tangential distortion |
| `lenscore/optics/lateralca.hpp` | Per-wavelength magnification |
| `lenscore/optics/zernike.hpp` | Zernike polynomial basis over the unit disc |
| `lenscore/optics/pupil.hpp` | Aperture polygon, cat's-eye clip, apodization ramp |
| `lenscore/optics/psf.hpp` | Generalized pupil function to PSF via one FFT |
| `lenscore/optics/vignette.hpp` | Natural `cos^n` falloff and mechanical pupil overlap |
| `lenscore/conv/fft.hpp` | Radix-2 complex FFT, 1D and 2D |
| `lenscore/conv/eff.hpp` | Efficient Filter Flow spatially varying convolution |
| `lenscore/params.hpp` | Versioned POD parameter struct |
| `lenscore/pipeline.hpp` | Stage ordering, the one public entry point |
| `lensdata/lensfile.hpp` | `.lens` JSON to `Params`. The only file that touches JSON |
| `lenscli/main.cpp` | CLI: render, generate targets, measure metrics |
| `tests/targets.hpp` | Synthetic target generators |
| `tests/metrics.hpp` | MTF50, radial mean, fringe width, energy, symmetry |
| `tests/test_*.cpp` | One test file per module |

---
### Task 1: Build scaffold and test harness

**Files:**
- Create: `CMakeLists.txt`, `lenscore/CMakeLists.txt`, `tests/CMakeLists.txt`, `tests/test_smoke.cpp`

**Interfaces:**
- Consumes: nothing
- Produces: CMake targets `lenscore` (INTERFACE library, header-only for now), `lens_tests` (executable registered with `ctest`)

- [ ] **Step 1: Write the failing test**

`tests/test_smoke.cpp`:
```cpp
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "lenscore/version.hpp"

TEST_CASE("version string is present") {
    CHECK(std::string(lens::kVersion) == "0.1.0");
}
```

- [ ] **Step 2: Write the build files**

`CMakeLists.txt`:
```cmake
cmake_minimum_required(VERSION 3.24)
project(lens LANGUAGES CXX VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LENS_BUILD_TESTS "Build the test suite" ON)

add_subdirectory(lenscore)

if(LENS_BUILD_TESTS)
  include(FetchContent)
  FetchContent_Declare(doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.11)
  FetchContent_MakeAvailable(doctest)
  enable_testing()
  add_subdirectory(tests)
endif()
```

`lenscore/CMakeLists.txt`:
```cmake
add_library(lenscore INTERFACE)
target_include_directories(lenscore INTERFACE ${CMAKE_SOURCE_DIR})
target_compile_features(lenscore INTERFACE cxx_std_20)
```

`tests/CMakeLists.txt`:
```cmake
add_executable(lens_tests test_smoke.cpp)
target_link_libraries(lens_tests PRIVATE lenscore doctest::doctest)
add_test(NAME lens_tests COMMAND lens_tests)
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `cmake -B build -S . && cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/version.hpp` not found.

- [ ] **Step 4: Write the minimal implementation**

`lenscore/version.hpp`:
```cpp
#pragma once
namespace lens { inline constexpr const char* kVersion = "0.1.0"; }
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: `1/1 Test #1: lens_tests ... Passed`

- [ ] **Step 6: Commit**

```bash
git add CMakeLists.txt lenscore tests
git commit -m "build: CMake scaffold with doctest"
```

---

### Task 2: Image type and PFM I/O

**Files:**
- Create: `lenscore/image.hpp`, `lenscore/pfm.hpp`, `tests/test_image.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `struct lens::Image { int w, h; std::vector<float> px; }` with interleaved RGB, `px.size() == w*h*3`
  - `float& Image::at(int x, int y, int c)` and its const overload
  - `bool lens::pfm::write(const std::string& path, const Image&)`
  - `std::optional<Image> lens::pfm::read(const std::string& path)`

- [ ] **Step 1: Write the failing test**

`tests/test_image.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/image.hpp"
#include "lenscore/pfm.hpp"
#include <cstdio>

TEST_CASE("image indexing is interleaved RGB") {
    lens::Image im(4, 3);
    CHECK(im.px.size() == 4u * 3u * 3u);
    im.at(2, 1, 1) = 0.5f;
    CHECK(im.px[(1 * 4 + 2) * 3 + 1] == doctest::Approx(0.5f));
}

TEST_CASE("pfm round trip preserves float values exactly") {
    lens::Image a(5, 2);
    for (size_t i = 0; i < a.px.size(); ++i) a.px[i] = float(i) * 0.125f - 1.0f;
    const std::string p = "roundtrip.pfm";
    REQUIRE(lens::pfm::write(p, a));
    auto b = lens::pfm::read(p);
    REQUIRE(b.has_value());
    CHECK(b->w == a.w);
    CHECK(b->h == a.h);
    for (size_t i = 0; i < a.px.size(); ++i) CHECK(b->px[i] == a.px[i]);
    std::remove(p.c_str());
}
```

Add `test_image.cpp` to the `add_executable` list in `tests/CMakeLists.txt`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/image.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/image.hpp`:
```cpp
#pragma once
#include <vector>
#include <cassert>

namespace lens {

struct Image {
    int w = 0, h = 0;
    std::vector<float> px;   // interleaved RGB

    Image() = default;
    Image(int width, int height) : w(width), h(height), px(size_t(width) * height * 3, 0.0f) {}

    float&       at(int x, int y, int c)       { return px[(size_t(y) * w + x) * 3 + c]; }
    const float& at(int x, int y, int c) const { return px[(size_t(y) * w + x) * 3 + c]; }

    bool sameSize(const Image& o) const { return w == o.w && h == o.h; }
};

}  // namespace lens
```

`lenscore/pfm.hpp`:
```cpp
#pragma once
#include "lenscore/image.hpp"
#include <cstdio>
#include <optional>
#include <string>

namespace lens::pfm {

// PFM stores rows bottom-to-top. Negative scale means little-endian.
inline bool write(const std::string& path, const Image& im) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::fprintf(f, "PF\n%d %d\n-1.0\n", im.w, im.h);
    for (int y = im.h - 1; y >= 0; --y)
        std::fwrite(&im.px[size_t(y) * im.w * 3], sizeof(float), size_t(im.w) * 3, f);
    std::fclose(f);
    return true;
}

inline std::optional<Image> read(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return std::nullopt;
    char magic[3] = {};
    int w = 0, h = 0; double scale = 0;
    if (std::fscanf(f, "%2s %d %d %lf", magic, &w, &h, &scale) != 4 ||
        std::string(magic) != "PF" || w <= 0 || h <= 0) { std::fclose(f); return std::nullopt; }
    std::fgetc(f);  // single whitespace byte after the scale
    Image im(w, h);
    for (int y = h - 1; y >= 0; --y)
        if (std::fread(&im.px[size_t(y) * w * 3], sizeof(float), size_t(w) * 3, f) != size_t(w) * 3)
            { std::fclose(f); return std::nullopt; }
    std::fclose(f);
    return im;
}

}  // namespace lens::pfm
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all assertions pass.

- [ ] **Step 5: Commit**

```bash
git add lenscore/image.hpp lenscore/pfm.hpp tests/test_image.cpp tests/CMakeLists.txt
git commit -m "feat: Image type and PFM round trip"
```

---
### Task 3: Colour transfer and highlight recovery

**Files:**
- Create: `lenscore/color/transfer.hpp`, `tests/test_transfer.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `float lens::color::srgbToLinear(float)` and `float lens::color::linearToSrgb(float)`
  - `struct lens::color::Knee { float threshold = 0.85f; float peak = 8.0f; }`
  - `float lens::color::expandHighlights(float x, const Knee&)`
  - `float lens::color::compressHighlights(float y, const Knee&)`

The knee is a hyperbola chosen so the curve is continuous **and has slope exactly
1.0 at the threshold**, so no visible crease appears where recovery starts:

```
A = 1 - threshold          B = peak - threshold          c = 1 - A/B
u = (x - threshold) / A
y = threshold + A*u / (1 - c*u)
```

At `u = 0` the slope is 1; at `u = 1` the output is exactly `peak`. Above
`x = 1` the curve continues linearly with the slope it had at `u = 1`, which is
`B^2 / A^2`. It requires `peak > 1`, which the constructor asserts.

- [ ] **Step 1: Write the failing test**

`tests/test_transfer.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/color/transfer.hpp"

using namespace lens::color;

TEST_CASE("sRGB transfer round trips") {
    for (float v : {0.0f, 0.001f, 0.04f, 0.5f, 1.0f})
        CHECK(linearToSrgb(srgbToLinear(v)) == doctest::Approx(v).epsilon(1e-5));
}

TEST_CASE("sRGB transfer hits the known anchor") {
    CHECK(srgbToLinear(0.5f) == doctest::Approx(0.2140f).epsilon(1e-3));
}

TEST_CASE("knee is identity below the threshold") {
    Knee k;
    CHECK(expandHighlights(0.5f, k) == doctest::Approx(0.5f));
    CHECK(expandHighlights(k.threshold, k) == doctest::Approx(k.threshold));
}

TEST_CASE("knee lifts white to the peak") {
    Knee k;
    CHECK(expandHighlights(1.0f, k) == doctest::Approx(k.peak).epsilon(1e-4));
}

TEST_CASE("knee has unit slope at the threshold so there is no crease") {
    Knee k;
    const float h = 1e-4f;
    const float slope = (expandHighlights(k.threshold + h, k) - k.threshold) / h;
    CHECK(slope == doctest::Approx(1.0f).epsilon(1e-2));
}

TEST_CASE("knee is strictly increasing") {
    Knee k;
    float prev = -1.0f;
    for (int i = 0; i <= 200; ++i) {
        const float y = expandHighlights(float(i) / 100.0f, k);  // spans 0..2
        CHECK(y > prev);
        prev = y;
    }
}

TEST_CASE("knee inverts exactly") {
    Knee k;
    for (float v : {0.1f, 0.85f, 0.9f, 0.99f, 1.0f, 1.5f})
        CHECK(compressHighlights(expandHighlights(v, k), k) == doctest::Approx(v).epsilon(1e-4));
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/color/transfer.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/color/transfer.hpp`:
```cpp
#pragma once
#include <cassert>
#include <cmath>

namespace lens::color {

inline float srgbToLinear(float s) {
    return (s <= 0.04045f) ? s / 12.92f : std::pow((s + 0.055f) / 1.055f, 2.4f);
}

inline float linearToSrgb(float l) {
    return (l <= 0.0031308f) ? l * 12.92f : 1.055f * std::pow(l, 1.0f / 2.4f) - 0.055f;
}

struct Knee {
    float threshold = 0.85f;
    float peak      = 8.0f;
};

// Expands display-referred highlights back above 1.0 so that bloom, halation
// and bokeh have real energy to work with. C1 continuous at the threshold.
inline float expandHighlights(float x, const Knee& k) {
    assert(k.peak > 1.0f && k.threshold > 0.0f && k.threshold < 1.0f);
    if (x <= k.threshold) return x;
    const float A = 1.0f - k.threshold;
    const float B = k.peak - k.threshold;
    const float c = 1.0f - A / B;
    const float u = (x - k.threshold) / A;
    if (u <= 1.0f) return k.threshold + A * u / (1.0f - c * u);
    return k.peak + (x - 1.0f) * (B * B) / (A * A);   // linear continuation
}

inline float compressHighlights(float y, const Knee& k) {
    if (y <= k.threshold) return y;
    const float A = 1.0f - k.threshold;
    const float B = k.peak - k.threshold;
    const float c = 1.0f - A / B;
    if (y <= k.peak) {
        const float d = y - k.threshold;
        return k.threshold + A * (d / (A + c * d));
    }
    return 1.0f + (y - k.peak) * (A * A) / (B * B);
}

}  // namespace lens::color
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all seven cases pass.

- [ ] **Step 5: Commit**

```bash
git add lenscore/color/transfer.hpp tests/test_transfer.cpp tests/CMakeLists.txt
git commit -m "feat: sRGB transfer and C1 highlight recovery knee"
```

---

### Task 4: CIE colour matching and spectrum integration

**Files:**
- Create: `lenscore/color/cie.hpp`, `tests/test_cie.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `struct lens::color::XYZ { float x, y, z; }`
  - `struct lens::color::RGB { float r, g, b; }`
  - `XYZ lens::color::cmf(float lambda_nm)` — the CIE 1931 2-degree observer
  - `RGB lens::color::xyzToRec2020(const XYZ&)` and `XYZ lens::color::rec2020ToXyz(const RGB&)`
  - `float lens::color::cieYNormalisation(const float* lambdas, const float* weights, int n)` — the `sum(w * ybar)` divisor used to keep an equal-energy spectrum neutral

The colour matching functions use the multi-lobe piecewise-Gaussian analytic fit
(Wyman, Sloan and Shirley). Ten lines instead of a 243-number table, accurate to
about 1% of peak, which is far below our other error sources. If a future test
ever needs better, swap in the tabulated 5nm data behind the same `cmf()`
signature; nothing else changes.

- [ ] **Step 1: Write the failing test**

`tests/test_cie.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/color/cie.hpp"
#include <vector>

using namespace lens::color;

TEST_CASE("ybar peaks near 555nm") {
    float best = 0.0f, arg = 0.0f;
    for (float l = 380.0f; l <= 780.0f; l += 0.5f)
        if (cmf(l).y > best) { best = cmf(l).y; arg = l; }
    CHECK(arg == doctest::Approx(555.0f).epsilon(0.02));
    CHECK(best == doctest::Approx(1.0f).epsilon(0.05));
}

TEST_CASE("equal energy spectrum is illuminant E, chromaticity one third") {
    XYZ sum{0, 0, 0};
    for (float l = 380.0f; l <= 780.0f; l += 1.0f) {
        const XYZ c = cmf(l);
        sum.x += c.x; sum.y += c.y; sum.z += c.z;
    }
    const float t = sum.x + sum.y + sum.z;
    CHECK(sum.x / t == doctest::Approx(1.0f / 3.0f).epsilon(0.02));
    CHECK(sum.y / t == doctest::Approx(1.0f / 3.0f).epsilon(0.02));
}

TEST_CASE("Rec2020 matrices are mutual inverses") {
    const RGB in{0.3f, 0.7f, 0.2f};
    const RGB out = xyzToRec2020(rec2020ToXyz(in));
    CHECK(out.r == doctest::Approx(in.r).epsilon(1e-4));
    CHECK(out.g == doctest::Approx(in.g).epsilon(1e-4));
    CHECK(out.b == doctest::Approx(in.b).epsilon(1e-4));
}

TEST_CASE("Rec2020 white maps to a neutral chromaticity") {
    const XYZ w = rec2020ToXyz(RGB{1.0f, 1.0f, 1.0f});
    const float t = w.x + w.y + w.z;
    CHECK(w.x / t == doctest::Approx(0.3127f).epsilon(0.01));  // D65
    CHECK(w.y / t == doctest::Approx(0.3290f).epsilon(0.01));
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/color/cie.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/color/cie.hpp`:
```cpp
#pragma once
#include <cmath>

namespace lens::color {

struct XYZ { float x = 0, y = 0, z = 0; };
struct RGB { float r = 0, g = 0, b = 0; };

// Piecewise Gaussian: sigma differs either side of the mean.
inline float lobe(float x, float mu, float s1, float s2) {
    const float t = (x - mu) / (x < mu ? s1 : s2);
    return std::exp(-0.5f * t * t);
}

// Multi-lobe analytic fit to the CIE 1931 2-degree observer.
inline XYZ cmf(float l) {
    XYZ o;
    o.x = 1.056f * lobe(l, 599.8f, 37.9f, 31.0f)
        + 0.362f * lobe(l, 442.0f, 16.0f, 26.7f)
        - 0.065f * lobe(l, 501.1f, 20.4f, 26.2f);
    o.y = 0.821f * lobe(l, 568.8f, 46.9f, 40.5f)
        + 0.286f * lobe(l, 530.9f, 16.3f, 31.1f);
    o.z = 1.217f * lobe(l, 437.0f, 11.8f, 36.0f)
        + 0.681f * lobe(l, 459.0f, 26.0f, 13.8f);
    return o;
}

inline RGB xyzToRec2020(const XYZ& c) {
    return { 1.7166511880f * c.x - 0.3556707838f * c.y - 0.2533662814f * c.z,
            -0.6666843518f * c.x + 1.6164812366f * c.y + 0.0157685458f * c.z,
             0.0176398574f * c.x - 0.0427706133f * c.y + 0.9421031212f * c.z };
}

inline XYZ rec2020ToXyz(const RGB& c) {
    return { 0.6369580483f * c.r + 0.1446169036f * c.g + 0.1688809752f * c.b,
             0.2627002120f * c.r + 0.6779980715f * c.g + 0.0593017165f * c.b,
             0.0000000000f * c.r + 0.0280726930f * c.g + 1.0609850577f * c.b };
}

// Divisor that keeps an equal-energy spectrum neutral under any band layout.
inline float cieYNormalisation(const float* lambdas, const float* weights, int n) {
    float s = 0.0f;
    for (int i = 0; i < n; ++i) s += weights[i] * cmf(lambdas[i]).y;
    return s;
}

}  // namespace lens::color
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all four cases pass.

- [ ] **Step 5: Commit**

```bash
git add lenscore/color/cie.hpp tests/test_cie.cpp tests/CMakeLists.txt
git commit -m "feat: CIE colour matching and Rec2020 conversion"
```

---
### Task 5: Spectral upsampling model and per-colour fit

**Files:**
- Create: `lenscore/color/upsample.hpp`, `tests/test_upsample.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `lens::color::cmf`, `xyzToRec2020` (Task 4)
- Produces:
  - `struct lens::color::Coeffs { float c0, c1, c2; }`
  - `float lens::color::evalSpectrum(const Coeffs&, float lambda_nm)`
  - `RGB lens::color::spectrumToRec2020(const Coeffs&)`
  - `Coeffs lens::color::fitCoeffs(const RGB& target, Coeffs guess = {0,0,0})`

The model is Jakob and Hanika's sigmoid polynomial. Six floating-point
operations per wavelength:

```
f(lambda) = S(c0*n^2 + c1*n + c2),   S(x) = 1/2 + x / (2*sqrt(1 + x^2))
```

`n` is the wavelength rescaled to `[0,1]` over 380-780nm, which keeps the
coefficients order-one and the fit well conditioned. `S` is bounded to `[0,1]`
by construction, so the spectrum is always physically realisable — that is the
whole point of the sigmoid and the reason no clamping is ever needed.

**The model produces reflectance, our pixels are radiance.** Callers fit the
*normalised* colour and carry `max(r,g,b)` separately as a scalar multiplier.
`fitCoeffs` therefore asserts its target is normalised.

- [ ] **Step 1: Write the failing test**

`tests/test_upsample.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/color/upsample.hpp"
#include <cmath>

using namespace lens::color;

TEST_CASE("sigmoid is bounded and centred") {
    CHECK(evalSpectrum(Coeffs{0, 0, 0}, 550.0f) == doctest::Approx(0.5f));
    for (float c2 : {-50.0f, -1.0f, 0.0f, 1.0f, 50.0f}) {
        const float v = evalSpectrum(Coeffs{0, 0, c2}, 550.0f);
        CHECK(v >= 0.0f);
        CHECK(v <= 1.0f);
    }
}

TEST_CASE("large positive constant gives a flat white spectrum") {
    const Coeffs c{0.0f, 0.0f, 1e4f};
    for (float l = 400.0f; l <= 700.0f; l += 50.0f)
        CHECK(evalSpectrum(c, l) == doctest::Approx(1.0f).epsilon(1e-3));
    const RGB rgb = spectrumToRec2020(c);
    CHECK(rgb.r == doctest::Approx(1.0f).epsilon(0.02));
    CHECK(rgb.g == doctest::Approx(1.0f).epsilon(0.02));
    CHECK(rgb.b == doctest::Approx(1.0f).epsilon(0.02));
}

TEST_CASE("fit round trips a range of normalised colours") {
    const RGB targets[] = {
        {1.0f, 1.0f, 1.0f}, {1.0f, 0.5f, 0.2f}, {0.2f, 1.0f, 0.3f},
        {0.1f, 0.3f, 1.0f}, {1.0f, 1.0f, 0.1f}, {0.6f, 0.6f, 0.6f},
    };
    for (const RGB& t : targets) {
        const Coeffs c = fitCoeffs(t);
        const RGB got = spectrumToRec2020(c);
        CAPTURE(t.r); CAPTURE(t.g); CAPTURE(t.b);
        CHECK(got.r == doctest::Approx(t.r).epsilon(0.03));
        CHECK(got.g == doctest::Approx(t.g).epsilon(0.03));
        CHECK(got.b == doctest::Approx(t.b).epsilon(0.03));
    }
}

TEST_CASE("fitted spectra stay smooth, no oscillation") {
    const Coeffs c = fitCoeffs(RGB{1.0f, 0.35f, 0.1f});
    int reversals = 0;
    float prev = evalSpectrum(c, 380.0f), prevSlope = 0.0f;
    for (float l = 385.0f; l <= 780.0f; l += 5.0f) {
        const float v = evalSpectrum(c, l);
        const float slope = v - prev;
        if (slope * prevSlope < 0.0f) ++reversals;
        prevSlope = slope; prev = v;
    }
    CHECK(reversals <= 2);   // a quadratic through a monotone sigmoid cannot wiggle more
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/color/upsample.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/color/upsample.hpp`:
```cpp
#pragma once
#include "lenscore/color/cie.hpp"
#include <cassert>
#include <cmath>

namespace lens::color {

inline constexpr float kLambdaMin = 380.0f;
inline constexpr float kLambdaMax = 780.0f;

struct Coeffs { float c0 = 0, c1 = 0, c2 = 0; };

inline float normaliseLambda(float l) { return (l - kLambdaMin) / (kLambdaMax - kLambdaMin); }

inline float sigmoid(float x) { return 0.5f + x / (2.0f * std::sqrt(1.0f + x * x)); }

inline float evalSpectrum(const Coeffs& c, float lambda_nm) {
    const float n = normaliseLambda(lambda_nm);
    return sigmoid(std::fma(std::fma(c.c0, n, c.c1), n, c.c2));
}

// Integrates the spectrum against the CMFs under illuminant E, at 5nm.
inline RGB spectrumToRec2020(const Coeffs& c) {
    XYZ acc{0, 0, 0};
    float norm = 0.0f;
    for (float l = kLambdaMin; l <= kLambdaMax; l += 5.0f) {
        const XYZ m = cmf(l);
        const float f = evalSpectrum(c, l);
        acc.x += f * m.x; acc.y += f * m.y; acc.z += f * m.z;
        norm  += m.y;
    }
    acc.x /= norm; acc.y /= norm; acc.z /= norm;
    return xyzToRec2020(acc);
}

// Levenberg-Marquardt on three parameters against three residuals.
inline Coeffs fitCoeffs(const RGB& target, Coeffs guess = {}) {
    assert(std::max({target.r, target.g, target.b}) <= 1.0001f &&
           "fitCoeffs expects a normalised colour; carry the scale separately");

    auto residual = [&](const Coeffs& c) {
        const RGB got = spectrumToRec2020(c);
        return std::array<float, 3>{got.r - target.r, got.g - target.g, got.b - target.b};
    };

    Coeffs c = guess;
    float damping = 1e-3f;

    for (int iter = 0; iter < 200; ++iter) {
        const auto r = residual(c);
        const float err = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        if (err < 1e-10f) break;

        // Central-difference Jacobian, 3x3.
        float J[3][3];
        float* p[3] = {&c.c0, &c.c1, &c.c2};
        for (int j = 0; j < 3; ++j) {
            const float save = *p[j];
            const float h = std::max(1e-3f, std::abs(save) * 1e-3f);
            *p[j] = save + h; const auto rp = residual(c);
            *p[j] = save - h; const auto rm = residual(c);
            *p[j] = save;
            for (int i = 0; i < 3; ++i) J[i][j] = (rp[i] - rm[i]) / (2.0f * h);
        }

        // Normal equations (J^T J + damping I) d = -J^T r, solved by elimination.
        double A[3][4] = {};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double s = 0; for (int k = 0; k < 3; ++k) s += double(J[k][i]) * J[k][j];
                A[i][j] = s + (i == j ? damping : 0.0);
            }
            double s = 0; for (int k = 0; k < 3; ++k) s += double(J[k][i]) * r[k];
            A[i][3] = -s;
        }
        for (int i = 0; i < 3; ++i) {
            int piv = i;
            for (int k = i + 1; k < 3; ++k) if (std::abs(A[k][i]) > std::abs(A[piv][i])) piv = k;
            if (std::abs(A[piv][i]) < 1e-18) { damping *= 10.0f; goto next; }
            for (int j = 0; j < 4; ++j) std::swap(A[i][j], A[piv][j]);
            for (int k = 0; k < 3; ++k) {
                if (k == i) continue;
                const double f = A[k][i] / A[i][i];
                for (int j = i; j < 4; ++j) A[k][j] -= f * A[i][j];
            }
        }
        {
            Coeffs trial{c.c0 + float(A[0][3] / A[0][0]),
                         c.c1 + float(A[1][3] / A[1][1]),
                         c.c2 + float(A[2][3] / A[2][2])};
            const auto rt = residual(trial);
            const float et = rt[0] * rt[0] + rt[1] * rt[1] + rt[2] * rt[2];
            if (et < err) { c = trial; damping = std::max(1e-7f, damping * 0.3f); }
            else          { damping *= 10.0f; }
        }
        next:;
        if (damping > 1e7f) break;
    }
    return c;
}

}  // namespace lens::color
```

Add `#include <array>` and `#include <algorithm>` at the top alongside the
existing includes.

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all four cases pass. If the round-trip case fails only on the most
saturated targets, raise the iteration cap before loosening the tolerance — the
tolerance is the specification, the iteration count is not.

- [ ] **Step 5: Commit**

```bash
git add lenscore/color/upsample.hpp tests/test_upsample.cpp tests/CMakeLists.txt
git commit -m "feat: Jakob-Hanika spectral upsampling model and fit"
```

---
### Task 6: Spectral upsampling table and fast lookup

**Files:**
- Create: `lenscore/color/spectable.hpp`, `rgb2spec/main.cpp`, `rgb2spec/CMakeLists.txt`, `tests/test_spectable.cpp`
- Modify: `CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `lens::color::fitCoeffs`, `Coeffs`, `RGB` (Task 5)
- Produces:
  - `struct lens::color::SpecTable { int res; std::vector<float> data; }` holding `3 * res^3 * 3` floats
  - `SpecTable lens::color::buildTable(int res)`
  - `bool lens::color::writeTable(const std::string&, const SpecTable&)` / `std::optional<SpecTable> lens::color::readTable(const std::string&)`
  - `Coeffs lens::color::lookup(const SpecTable&, const RGB& colour)` — accepts unnormalised colour, handles the scale internally

Layout follows Jakob Algorithm 2. Pick the largest component `i`. The two
remaining components divided by the maximum give two axes; the maximum itself
gives the third. That third axis is warped by `smoothstep(smoothstep(k/(res-1)))`
because the coefficients change fastest near black and near white, and a linear
axis wastes resolution in the middle where nothing happens.

- [ ] **Step 1: Write the failing test**

`tests/test_spectable.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/color/spectable.hpp"
#include <cstdio>

using namespace lens::color;

static const SpecTable& smallTable() {
    static SpecTable t = buildTable(12);   // small, so the suite stays fast
    return t;
}

TEST_CASE("table has the documented size") {
    const SpecTable& t = smallTable();
    CHECK(t.res == 12);
    CHECK(t.data.size() == size_t(3) * 12 * 12 * 12 * 3);
}

TEST_CASE("lookup round trips grey at several brightnesses") {
    for (float v : {0.15f, 0.4f, 0.75f, 1.0f}) {
        const RGB in{v, v, v};
        const RGB got = spectrumToRec2020(lookup(smallTable(), in));
        CAPTURE(v);
        CHECK(got.r / v == doctest::Approx(1.0f).epsilon(0.06));
        CHECK(got.g / v == doctest::Approx(1.0f).epsilon(0.06));
        CHECK(got.b / v == doctest::Approx(1.0f).epsilon(0.06));
    }
}

TEST_CASE("lookup is scale equivariant at and beyond the gamut boundary") {
    // The table spans [0,1]^3 only, so everything at or past the boundary along a ray
    // clamps to the same cell and must agree exactly. Below the boundary the model is
    // deliberately NOT scale-invariant: that is the entire reason the table has a third,
    // brightness axis. A muted colour legitimately needs a different spectral shape than
    // the same hue at full brightness. Comparing an in-gamut point against an out-of-gamut
    // one would assert a property this design does not have.
    const Coeffs a = lookup(smallTable(), RGB{1.0f, 0.5f, 0.25f});
    const Coeffs b = lookup(smallTable(), RGB{10.0f, 5.0f, 2.5f});
    CHECK(a.c0 == doctest::Approx(b.c0).epsilon(1e-4));
    CHECK(a.c1 == doctest::Approx(b.c1).epsilon(1e-4));
    CHECK(a.c2 == doctest::Approx(b.c2).epsilon(1e-4));
}

TEST_CASE("lookup round trips saturated hues within tolerance") {
    const RGB targets[] = {{1.0f, 0.2f, 0.1f}, {0.1f, 1.0f, 0.2f}, {0.15f, 0.2f, 1.0f}};
    for (const RGB& t : targets) {
        const RGB got = spectrumToRec2020(lookup(smallTable(), t));
        CAPTURE(t.r); CAPTURE(t.g); CAPTURE(t.b);
        CHECK(got.r == doctest::Approx(t.r).epsilon(0.10));
        CHECK(got.g == doctest::Approx(t.g).epsilon(0.10));
        CHECK(got.b == doctest::Approx(t.b).epsilon(0.10));
    }
}

TEST_CASE("table serialises and reloads bit exactly") {
    const std::string p = "tbl.bin";
    REQUIRE(writeTable(p, smallTable()));
    auto back = readTable(p);
    REQUIRE(back.has_value());
    CHECK(back->res == smallTable().res);
    CHECK(back->data == smallTable().data);
    std::remove(p.c_str());
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/color/spectable.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/color/spectable.hpp`:
```cpp
#pragma once
#include "lenscore/color/upsample.hpp"
#include <algorithm>
#include <cstdio>
#include <optional>
#include <string>
#include <vector>

namespace lens::color {

struct SpecTable {
    int res = 0;
    std::vector<float> data;   // [i][z][y][x][0..2], i is the largest-component index
};

inline float smoothstep(float x) { return x * x * (3.0f - 2.0f * x); }

// Warped brightness axis: resolution concentrated near black and near white.
inline float axisScale(int k, int res) { return smoothstep(smoothstep(float(k) / float(res - 1))); }

inline size_t tableIndex(int res, int i, int z, int y, int x) {
    return ((((size_t(i) * res + z) * res + y) * res + x)) * 3;
}

inline SpecTable buildTable(int res) {
    SpecTable t; t.res = res; t.data.assign(size_t(3) * res * res * res * 3, 0.0f);
    for (int i = 0; i < 3; ++i) {
        for (int z = 0; z < res; ++z) {
            const float scale = std::max(1e-4f, axisScale(z, res));
            Coeffs warm{};                                   // warm start down each row
            for (int y = 0; y < res; ++y) {
                for (int x = 0; x < res; ++x) {
                    const float a = float(x) / float(res - 1);
                    const float b = float(y) / float(res - 1);
                    float rgb[3];
                    rgb[i]           = scale;
                    rgb[(i + 1) % 3] = a * scale;
                    rgb[(i + 2) % 3] = b * scale;
                    // Fit the target AS IS. Its max component is already `scale` <= 1, which
                    // satisfies fitCoeffs's normalised-input contract. Do NOT divide by
                    // `scale`: that collapses every z-slice onto the same amplitude-1 target,
                    // so all slices fit identically and the brightness axis stops encoding
                    // brightness -- defeating the point of a 3D table.
                    const Coeffs c = fitCoeffs(RGB{rgb[0], rgb[1], rgb[2]}, warm);
                    warm = c;
                    const size_t o = tableIndex(res, i, z, y, x);
                    t.data[o + 0] = c.c0; t.data[o + 1] = c.c1; t.data[o + 2] = c.c2;
                }
            }
        }
    }
    return t;
}

inline Coeffs lookup(const SpecTable& t, const RGB& colour) {
    const float v[3] = {colour.r, colour.g, colour.b};
    int i = 0;
    if (v[1] > v[i]) i = 1;
    if (v[2] > v[i]) i = 2;
    const float mx = v[i];
    if (mx <= 0.0f) return Coeffs{0.0f, 0.0f, -1e4f};   // black: flat zero spectrum

    const float a = v[(i + 1) % 3] / mx;
    const float b = v[(i + 2) % 3] / mx;

    // Invert the warped brightness axis by search; res is small so this is cheap.
    const float target = std::min(mx, 1.0f);
    int z0 = 0;
    while (z0 + 2 < t.res && axisScale(z0 + 1, t.res) < target) ++z0;
    const float s0 = axisScale(z0, t.res), s1 = axisScale(z0 + 1, t.res);
    const float fz = (s1 > s0) ? std::clamp((target - s0) / (s1 - s0), 0.0f, 1.0f) : 0.0f;

    const float fx = a * (t.res - 1), fy = b * (t.res - 1);
    const int x0 = std::clamp(int(fx), 0, t.res - 2), y0 = std::clamp(int(fy), 0, t.res - 2);
    const float dx = fx - x0, dy = fy - y0;

    Coeffs out{};
    float* o[3] = {&out.c0, &out.c1, &out.c2};
    for (int k = 0; k < 3; ++k) {
        float acc = 0.0f;
        for (int zz = 0; zz < 2; ++zz)
            for (int yy = 0; yy < 2; ++yy)
                for (int xx = 0; xx < 2; ++xx) {
                    const float w = (zz ? fz : 1 - fz) * (yy ? dy : 1 - dy) * (xx ? dx : 1 - dx);
                    acc += w * t.data[tableIndex(t.res, i, z0 + zz, y0 + yy, x0 + xx) + k];
                }
        *o[k] = acc;
    }
    return out;
}

inline bool writeTable(const std::string& path, const SpecTable& t) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::fwrite("LSPT", 1, 4, f);
    std::fwrite(&t.res, sizeof(int), 1, f);
    std::fwrite(t.data.data(), sizeof(float), t.data.size(), f);
    std::fclose(f);
    return true;
}

inline std::optional<SpecTable> readTable(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return std::nullopt;
    char magic[4];
    SpecTable t;
    if (std::fread(magic, 1, 4, f) != 4 || std::string(magic, 4) != "LSPT" ||
        std::fread(&t.res, sizeof(int), 1, f) != 1 || t.res < 2) { std::fclose(f); return std::nullopt; }
    t.data.resize(size_t(3) * t.res * t.res * t.res * 3);
    const bool ok = std::fread(t.data.data(), sizeof(float), t.data.size(), f) == t.data.size();
    std::fclose(f);
    return ok ? std::optional<SpecTable>(t) : std::nullopt;
}

}  // namespace lens::color
```

`rgb2spec/main.cpp`:
```cpp
#include "lenscore/color/spectable.hpp"
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
    const int res = (argc > 1) ? std::atoi(argv[1]) : 64;
    const char* out = (argc > 2) ? argv[2] : "lensdata/rgb2spec/rec2020.bin";
    std::printf("building %d^3 table...\n", res);
    const lens::color::SpecTable t = lens::color::buildTable(res);
    if (!lens::color::writeTable(out, t)) { std::fprintf(stderr, "write failed: %s\n", out); return 1; }
    std::printf("wrote %s (%zu floats)\n", out, t.data.size());
    return 0;
}
```

`rgb2spec/CMakeLists.txt`:
```cmake
add_executable(rgb2spec main.cpp)
target_link_libraries(rgb2spec PRIVATE lenscore)
```

Add `add_subdirectory(rgb2spec)` to the top-level `CMakeLists.txt`.

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all five cases pass.

- [ ] **Step 5: Generate the shipping table**

```bash
mkdir -p lensdata/rgb2spec
./build/rgb2spec/rgb2spec 64 lensdata/rgb2spec/rec2020.bin
```
Expected: `wrote lensdata/rgb2spec/rec2020.bin (2359296 floats)`, about 9.4MB, one to two minutes.

- [ ] **Step 6: Commit**

```bash
git add lenscore/color/spectable.hpp rgb2spec tests/test_spectable.cpp \
        tests/CMakeLists.txt CMakeLists.txt lensdata/rgb2spec/rec2020.bin
git commit -m "feat: spectral upsampling table, generator and lookup"
```

---
### Task 7: Dispersion and secondary spectrum

**Files:**
- Create: `lenscore/optics/dispersion.hpp`, `tests/test_dispersion.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `struct lens::optics::Dispersion { float B[3]; float C[3]; std::vector<float> correction_nm; float residual; }`
  - `Dispersion lens::optics::bk7()` — the reference glass
  - `float lens::optics::refractiveIndex(const Dispersion&, float lambda_nm)`
  - `float lens::optics::focusError(const Dispersion&, float lambda_nm, float lambda_ref_nm)` — relative focal-length error `F(lambda)/F_ref - 1`, after correction

Sellmeier, with wavelength in micrometres:

```
n^2(lambda) = 1 + sum_i  B_i * lambda^2 / (lambda^2 - C_i)
```

Focal length goes as `1/(n-1)`, so the raw relative error is
`(n_ref - 1)/(n(lambda) - 1) - 1`.

**This is where a singlet becomes a real lens.** Raw Sellmeier error is
monotonic in wavelength: that is an uncorrected element, and it fringes red to
blue. A real photographic lens is an achromat or apochromat, corrected so focus
coincides at two or three wavelengths, and what remains is the *secondary
spectrum* — the reason real glass fringes green and magenta. We model it by
Lagrange-interpolating the raw error at the correction wavelengths as a
polynomial in `x = 1/lambda^2` (the variable in which chromatic focus is very
nearly linear for normal glass), subtracting it, and scaling by `residual`.

The number of correction wavelengths selects the lens class with no branching:
zero for an uncorrected singlet, one for a simple lens focused at one colour,
two for an achromat, three for an apochromat.

- [ ] **Step 1: Write the failing test**

`tests/test_dispersion.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/optics/dispersion.hpp"

using namespace lens::optics;

TEST_CASE("BK7 hits its published index at the d line") {
    CHECK(refractiveIndex(bk7(), 587.6f) == doctest::Approx(1.5168f).epsilon(1e-3));
}

TEST_CASE("BK7 reproduces its published Abbe number, pinning the whole curve shape") {
    // One index value would tolerate a wrong B/C pair that happens to land near 1.5168
    // at the d-line while distorting dispersion elsewhere. The Abbe number is built from
    // three points, so it pins the SHAPE of the curve, not just one height on it.
    const Dispersion d = bk7();
    const float nd = refractiveIndex(d, 587.56f);
    const float nF = refractiveIndex(d, 486.13f);
    const float nC = refractiveIndex(d, 656.27f);
    CHECK(nd == doctest::Approx(1.51680f).epsilon(1e-4));
    CHECK(nF == doctest::Approx(1.52238f).epsilon(1e-4));
    CHECK(nC == doctest::Approx(1.51432f).epsilon(1e-4));
    CHECK((nd - 1.0f) / (nF - nC) == doctest::Approx(64.17f).epsilon(2e-3));
}

TEST_CASE("dispersion is normal, index falls as wavelength rises") {
    const Dispersion d = bk7();
    float prev = refractiveIndex(d, 400.0f);
    for (float l = 420.0f; l <= 760.0f; l += 20.0f) {
        const float n = refractiveIndex(d, l);
        CHECK(n < prev);
        prev = n;
    }
}

TEST_CASE("an uncorrected singlet has monotonic focus error, rising with wavelength") {
    // Sign convention: focusError returns the relative focal LENGTH error,
    // F(lambda)/F_ref - 1. Blue has a higher refractive index, so a shorter focal
    // length, so it focuses CLOSER than the reference -- a negative error. Red
    // focuses further away -- positive. The error therefore RISES with wavelength.
    // Getting this backwards returns relative optical power instead, which is the
    // exact negative, and would make chromatic defocus oppose field curvature
    // downstream rather than add to it.
    Dispersion d = bk7();
    d.correction_nm.clear();
    CHECK(focusError(d, 400.0f, 650.0f) < 0.0f);   // blue focuses closer
    CHECK(focusError(d, 760.0f, 650.0f) > 0.0f);   // red focuses further
    float prev = focusError(d, 400.0f, 650.0f);
    for (float l = 420.0f; l <= 760.0f; l += 20.0f) {
        const float e = focusError(d, l, 650.0f);
        CHECK(e > prev);
        prev = e;
    }
}

TEST_CASE("an achromat has exactly zero focus error at both correction lines") {
    Dispersion d = bk7();
    d.correction_nm = {486.1f, 656.3f};      // F and C lines
    CHECK(focusError(d, 486.1f, 650.0f) == doctest::Approx(0.0f).epsilon(1e-6));
    CHECK(focusError(d, 656.3f, 650.0f) == doctest::Approx(0.0f).epsilon(1e-6));
}

TEST_CASE("the achromat's secondary spectrum bulges with one sign between the lines") {
    Dispersion d = bk7();
    d.correction_nm = {486.1f, 656.3f};
    const float mid = focusError(d, 570.0f, 650.0f);
    CHECK(mid != doctest::Approx(0.0f).epsilon(1e-5));
    for (float l = 500.0f; l <= 640.0f; l += 20.0f)
        CHECK(focusError(d, l, 650.0f) * mid > 0.0f);   // no sign flip in between
}

TEST_CASE("an apochromat zeroes three lines") {
    Dispersion d = bk7();
    d.correction_nm = {436.0f, 546.1f, 656.3f};
    for (float l : d.correction_nm)
        CHECK(focusError(d, l, 650.0f) == doctest::Approx(0.0f).epsilon(1e-6));
}

TEST_CASE("residual scales the remaining error linearly") {
    Dispersion a = bk7(); a.correction_nm = {486.1f, 656.3f}; a.residual = 1.0f;
    Dispersion b = a;                                        b.residual = 0.25f;
    CHECK(focusError(b, 570.0f, 650.0f) == doctest::Approx(0.25f * focusError(a, 570.0f, 650.0f)));
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/optics/dispersion.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/optics/dispersion.hpp`:
```cpp
#pragma once
#include <cmath>
#include <vector>

namespace lens::optics {

struct Dispersion {
    float B[3] = {1.03961212f, 0.231792344f, 1.01046945f};   // BK7
    float C[3] = {0.00600069867f, 0.0200179144f, 103.560653f};
    std::vector<float> correction_nm{};   // 0 singlet, 2 achromat, 3 apochromat
    float residual = 1.0f;
};

inline Dispersion bk7() { return Dispersion{}; }

inline float refractiveIndex(const Dispersion& d, float lambda_nm) {
    const double l = lambda_nm / 1000.0;      // micrometres
    const double l2 = l * l;
    double n2 = 1.0;
    for (int i = 0; i < 3; ++i) n2 += double(d.B[i]) * l2 / (l2 - double(d.C[i]));
    return float(std::sqrt(n2));
}

// Raw relative focal-length error before correction: F(lambda)/F_ref - 1.
// Focal length goes as 1/(n-1), so the ratio is (nRef-1)/(n-1). NOT (n-1)/(nRef-1),
// which is the relative optical power -- the negative of this -- and would flip the
// sign of chromatic defocus everywhere downstream.
inline float rawFocusError(const Dispersion& d, float lambda_nm, float lambda_ref_nm) {
    const float nRef = refractiveIndex(d, lambda_ref_nm);
    const float n    = refractiveIndex(d, lambda_nm);
    return (nRef - 1.0f) / (n - 1.0f) - 1.0f;
}

inline float focusError(const Dispersion& d, float lambda_nm, float lambda_ref_nm) {
    const float raw = rawFocusError(d, lambda_nm, lambda_ref_nm);
    const size_t k = d.correction_nm.size();
    if (k == 0) return raw * d.residual;

    // Lagrange interpolation of the raw error at the correction points,
    // in x = 1/lambda^2 with lambda in micrometres.
    auto xOf = [](float nm) { const double u = nm / 1000.0; return 1.0 / (u * u); };
    const double x = xOf(lambda_nm);
    double corr = 0.0;
    for (size_t i = 0; i < k; ++i) {
        double term = rawFocusError(d, d.correction_nm[i], lambda_ref_nm);
        for (size_t j = 0; j < k; ++j) {
            if (j == i) continue;
            term *= (x - xOf(d.correction_nm[j])) / (xOf(d.correction_nm[i]) - xOf(d.correction_nm[j]));
        }
        corr += term;
    }
    return float((raw - corr) * d.residual);
}

}  // namespace lens::optics
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all seven cases pass.

- [ ] **Step 5: Commit**

```bash
git add lenscore/optics/dispersion.hpp tests/test_dispersion.cpp tests/CMakeLists.txt
git commit -m "feat: Sellmeier dispersion with achromat secondary spectrum"
```

---
### Task 8: Plane, geometry, distortion and lateral chromatic aberration

**Files:**
- Create: `lenscore/plane.hpp`, `lenscore/geometry.hpp`, `lenscore/optics/distortion.hpp`, `lenscore/optics/lateralca.hpp`, `tests/test_warp.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `struct lens::Plane { int w, h; std::vector<float> v; float& at(int,int); }`
  - `struct lens::Frame { float cx, cy, halfDiag; }` and `Frame lens::frameOf(int w, int h)` — the normalised field coordinate system, `t = 1` at the image corner
  - `float lens::sampleBilinear(const Plane&, float x, float y)` — clamped at the edges
  - `struct lens::optics::Distortion { float k1, k2, k3, p1, p2; }`
  - `void lens::optics::applyDistortion(const Distortion&, float ux, float uy, float& dx, float& dy)`
  - `float lens::optics::lateralMagnification(float k_l, float lambda_nm, float lambda_hat_nm, float t)`
  - `float lens::optics::inverseLateralRadius(float K, float tOut)`
  - `Plane lens::optics::warpPlane(const Plane& src, const Distortion&, float K)`

**The inverse-radius detail.** Magnification depends on radius, so undoing it is
not a division. The forward map is `t_out = t_src * (1 + K * t_src)`. Resampling
needs the source radius for a known output radius, so we solve the quadratic
exactly:

```
t_src = ( -1 + sqrt(1 + 4*K*t_out) ) / (2*K)      for K != 0
t_src = t_out                                      for K == 0
```

Getting this wrong produces a magnification error that grows toward the corners
— precisely where lateral chromatic aberration is supposed to be measured.

- [ ] **Step 1: Write the failing test**

`tests/test_warp.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/optics/lateralca.hpp"
#include <cmath>

using namespace lens;
using namespace lens::optics;

static Plane dot(int w, int h, int px, int py) {
    Plane p(w, h);
    p.at(px, py) = 1.0f;
    return p;
}

TEST_CASE("frame puts t = 1 at the corner and 0 at the centre") {
    const Frame f = frameOf(101, 51);
    CHECK(f.cx == doctest::Approx(50.0f));
    CHECK(f.cy == doctest::Approx(25.0f));
    const float t = std::sqrt(f.cx * f.cx + f.cy * f.cy) / f.halfDiag;
    CHECK(t == doctest::Approx(1.0f));
}

TEST_CASE("inverse radius exactly inverts the forward magnification") {
    for (float K : {-0.05f, -0.001f, 0.0f, 0.002f, 0.08f}) {
        for (float tOut : {0.0f, 0.25f, 0.5f, 1.0f}) {
            const float ts = inverseLateralRadius(K, tOut);
            CHECK(ts * (1.0f + K * ts) == doctest::Approx(tOut).epsilon(1e-5));
        }
    }
}

TEST_CASE("magnification is unity at the reference wavelength") {
    CHECK(lateralMagnification(0.02f, 650.0f, 650.0f, 1.0f) == doctest::Approx(1.0f));
}

TEST_CASE("blue is magnified differently from red, and only off axis") {
    const float mBlue = lateralMagnification(1e-4f, 450.0f, 650.0f, 1.0f);
    CHECK(mBlue != doctest::Approx(1.0f));
    CHECK(lateralMagnification(1e-4f, 450.0f, 650.0f, 0.0f) == doctest::Approx(1.0f));
}

TEST_CASE("identity warp reproduces the source") {
    Plane src(33, 21);
    for (int y = 0; y < 21; ++y)
        for (int x = 0; x < 33; ++x) src.at(x, y) = float((x * 7 + y * 3) % 11) / 11.0f;
    const Plane out = warpPlane(src, Distortion{}, 0.0f);
    for (int y = 0; y < 21; ++y)
        for (int x = 0; x < 33; ++x) CHECK(out.at(x, y) == doctest::Approx(src.at(x, y)).epsilon(1e-5));
}

TEST_CASE("positive K pushes a feature outward from the centre") {
    // K must be large enough to move a DISCRETE pixel argmax. At K = 0.01 the feature at
    // this radius shifts well under half a pixel, so the assertion below would be
    // mathematically unsatisfiable no matter how correct the code is.
    const int w = 129, h = 129;
    const Plane src = dot(w, h, 100, 64);          // right of centre
    const Plane out = warpPlane(src, Distortion{}, 0.08f);
    float best = -1.0f; int bx = 0;
    for (int x = 0; x < w; ++x) if (out.at(x, 64) > best) { best = out.at(x, 64); bx = x; }
    CHECK(bx > 100);
}

TEST_CASE("barrel distortion pulls the corners inward") {
    Distortion d; d.k1 = -0.10f;
    float dx = 0, dy = 0;
    applyDistortion(d, 0.7071f, 0.7071f, dx, dy);
    CHECK(std::sqrt(dx * dx + dy * dy) < 1.0f);
}

TEST_CASE("warp conserves total energy to within edge effects") {
    Plane src(65, 65);
    for (int y = 20; y < 45; ++y) for (int x = 20; x < 45; ++x) src.at(x, y) = 1.0f;
    const Plane out = warpPlane(src, Distortion{}, 0.005f);
    double a = 0, b = 0;
    for (float v : src.v) a += v;
    for (float v : out.v) b += v;
    CHECK(b / a == doctest::Approx(1.0).epsilon(0.03));
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/optics/lateralca.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/plane.hpp`:
```cpp
#pragma once
#include <vector>

namespace lens {

struct Plane {
    int w = 0, h = 0;
    std::vector<float> v;
    Plane() = default;
    Plane(int width, int height) : w(width), h(height), v(size_t(width) * height, 0.0f) {}
    float&       at(int x, int y)       { return v[size_t(y) * w + x]; }
    const float& at(int x, int y) const { return v[size_t(y) * w + x]; }
};

}  // namespace lens
```

`lenscore/geometry.hpp`:
```cpp
#pragma once
#include "lenscore/plane.hpp"
#include <algorithm>
#include <cmath>

namespace lens {

// Normalised field coordinates: the image corner sits at t = 1.
struct Frame { float cx = 0, cy = 0, halfDiag = 1; };

inline Frame frameOf(int w, int h) {
    Frame f;
    f.cx = 0.5f * float(w - 1);
    f.cy = 0.5f * float(h - 1);
    f.halfDiag = std::sqrt(f.cx * f.cx + f.cy * f.cy);
    if (f.halfDiag <= 0.0f) f.halfDiag = 1.0f;
    return f;
}

inline float sampleBilinear(const Plane& p, float x, float y) {
    x = std::clamp(x, 0.0f, float(p.w - 1));
    y = std::clamp(y, 0.0f, float(p.h - 1));
    const int x0 = std::min(int(x), p.w - 1), y0 = std::min(int(y), p.h - 1);
    const int x1 = std::min(x0 + 1, p.w - 1), y1 = std::min(y0 + 1, p.h - 1);
    const float fx = x - float(x0), fy = y - float(y0);
    return (1 - fy) * ((1 - fx) * p.at(x0, y0) + fx * p.at(x1, y0))
         +      fy  * ((1 - fx) * p.at(x0, y1) + fx * p.at(x1, y1));
}

}  // namespace lens
```

`lenscore/optics/distortion.hpp`:
```cpp
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
```

`lenscore/optics/lateralca.hpp`:
```cpp
#pragma once
#include "lenscore/geometry.hpp"
#include "lenscore/optics/distortion.hpp"
#include <cmath>

namespace lens::optics {

// m(lambda, t) = 1 + k_l * (lambda_hat - lambda) * t
inline float lateralMagnification(float k_l, float lambda_nm, float lambda_hat_nm, float t) {
    return 1.0f + k_l * (lambda_hat_nm - lambda_nm) * t;
}

// Solves t_out = t_src * (1 + K * t_src) for t_src.
//
// The textbook root, (-1 + sqrt(disc)) / (2K), is NUMERICALLY WRONG here. For small K
// it subtracts two nearly-equal numbers and then divides by a small denominator --
// catastrophic cancellation. Measured in float32 its worst relative error is 1.4e-3,
// over a hundred times this function's own test threshold, and it grows as K shrinks.
// Small K is a WELL-CORRECTED lens, so the naive form is least accurate on the best glass.
//
// Rationalising by (1 + sqrt(disc)) gives an algebraically identical expression with no
// subtraction of like quantities and no division by K, accurate to ~1e-7 and continuous
// as K approaches zero:
//     (-1 + sqrt(d)) / 2K  ==  (d - 1) / (2K (1 + sqrt(d)))  ==  2*tOut / (1 + sqrt(d))
inline float inverseLateralRadius(float K, float tOut) {
    const float disc = 1.0f + 4.0f * K * tOut;
    if (disc <= 0.0f) return tOut;                    // outside the invertible range
    return 2.0f * tOut / (1.0f + std::sqrt(disc));    // K == 0 falls out as tOut, no branch
}

// One resample applies magnification and distortion together, never twice.
inline Plane warpPlane(const Plane& src, const Distortion& d, float K) {
    const Frame f = frameOf(src.w, src.h);
    Plane out(src.w, src.h);
    for (int y = 0; y < src.h; ++y) {
        for (int x = 0; x < src.w; ++x) {
            const float nx = (float(x) - f.cx) / f.halfDiag;
            const float ny = (float(y) - f.cy) / f.halfDiag;
            const float tOut = std::sqrt(nx * nx + ny * ny);
            float ux = nx, uy = ny;
            if (tOut > 1e-9f) {
                const float tSrc = inverseLateralRadius(K, tOut);
                const float s = tSrc / tOut;
                ux = nx * s; uy = ny * s;
            }
            float dx = 0, dy = 0;
            applyDistortion(d, ux, uy, dx, dy);
            out.at(x, y) = sampleBilinear(src, dx * f.halfDiag + f.cx, dy * f.halfDiag + f.cy);
        }
    }
    return out;
}

}  // namespace lens::optics
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all eight cases pass.

- [ ] **Step 5: Commit**

```bash
git add lenscore/plane.hpp lenscore/geometry.hpp lenscore/optics/distortion.hpp \
        lenscore/optics/lateralca.hpp tests/test_warp.cpp tests/CMakeLists.txt
git commit -m "feat: field geometry, distortion and lateral CA resampling"
```

---
### Task 9: Vignetting, natural and mechanical

**Files:**
- Create: `lenscore/optics/vignette.hpp`, `tests/test_vignette.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `float lens::optics::circleOverlapArea(float d, float r1, float r2)`
  - `struct lens::optics::VignetteParams { float focal_mm, sensorHalfDiag_mm, naturalExp, tStopWide, tStop, rEntrance, sepNorm; }`
  - `float lens::optics::naturalFalloff(const VignetteParams&, float t)`
  - `float lens::optics::mechanicalFraction(const VignetteParams&, float t)`
  - `float lens::optics::mechanicalVanishStop(const VignetteParams&)`
  - `float lens::optics::vignette(const VignetteParams&, float t)`

**Natural falloff** is `cos^n(theta)` with `theta = atan(r/f)`. The exponent is a
parameter, not the constant 4. Aggarwal measured three real lenses across six
apertures and found the fall-off curves *do not* coincide with `cos^4` at any
setting — not even fully stopped down, where mechanical vignetting is absent and
the textbook says `cos^4` should be exact.

**Mechanical vignetting** is the aperture disc clipped by the barrel. In units
of the wide-open aperture radius, the aperture has radius `a = tStopWide/tStop`
and the clipping circle has radius `rEntrance`, displaced by `sepNorm * t`. The
overlap fraction is a closed-form circular-lens area.

The vanishing threshold falls out rather than being asserted: clipping stops
when the aperture fits entirely inside the offset circle, `a <= rEntrance -
sepNorm`, which happens at `tStop >= tStopWide / (rEntrance - sepNorm)`. With
the default geometry that is **f/4**, matching Aggarwal's measurement. Anyone
tuning `sepNorm` is tuning a number they can verify against a flat-field
photograph.

- [ ] **Step 1: Write the failing test**

`tests/test_vignette.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/optics/vignette.hpp"
#include <cmath>

using namespace lens::optics;

static VignetteParams defaults() {
    VignetteParams p;
    p.focal_mm = 32.0f; p.sensorHalfDiag_mm = 14.0f; p.naturalExp = 4.0f;
    p.tStopWide = 2.0f; p.tStop = 2.0f; p.rEntrance = 1.0f; p.sepNorm = 0.5f;
    return p;
}

TEST_CASE("overlap area handles containment and disjointness") {
    CHECK(circleOverlapArea(5.0f, 1.0f, 1.0f) == doctest::Approx(0.0f));            // disjoint
    CHECK(circleOverlapArea(0.0f, 1.0f, 2.0f) == doctest::Approx(pi));              // contained
    CHECK(circleOverlapArea(0.0f, 1.0f, 1.0f) == doctest::Approx(pi));              // coincident
}

TEST_CASE("overlap area of two unit circles at unit separation is the known value") {
    // 2*acos(1/2) - sqrt(3)/2 = 1.22837
    CHECK(circleOverlapArea(1.0f, 1.0f, 1.0f) == doctest::Approx(1.22837f).epsilon(1e-4));
}

TEST_CASE("everything is exactly unity at the optical centre") {
    const VignetteParams p = defaults();
    CHECK(naturalFalloff(p, 0.0f) == doctest::Approx(1.0f));
    CHECK(mechanicalFraction(p, 0.0f) == doctest::Approx(1.0f));
    CHECK(vignette(p, 0.0f) == doctest::Approx(1.0f));
}

TEST_CASE("natural falloff matches the cosine law analytically") {
    VignetteParams p = defaults();
    const float r = 1.0f * p.sensorHalfDiag_mm;
    const float expected = std::pow(std::cos(std::atan(r / p.focal_mm)), 4.0f);
    CHECK(naturalFalloff(p, 1.0f) == doctest::Approx(expected).epsilon(1e-5));
}

TEST_CASE("the exponent is free, not pinned at four") {
    VignetteParams a = defaults(); a.naturalExp = 4.0f;
    VignetteParams b = defaults(); b.naturalExp = 3.2f;
    CHECK(naturalFalloff(b, 1.0f) > naturalFalloff(a, 1.0f));
}

TEST_CASE("natural falloff decreases monotonically toward the corner") {
    const VignetteParams p = defaults();
    float prev = 1.0f;
    for (float t = 0.1f; t <= 1.0f; t += 0.1f) {
        const float v = naturalFalloff(p, t);
        CHECK(v < prev);
        prev = v;
    }
}

TEST_CASE("mechanical vignetting bites wide open") {
    const VignetteParams p = defaults();
    CHECK(mechanicalFraction(p, 1.0f) < 0.99f);
    CHECK(mechanicalFraction(p, 1.0f) > 0.0f);
}

TEST_CASE("mechanical vignetting vanishes above f/4, as measured") {
    VignetteParams p = defaults();
    CHECK(mechanicalVanishStop(p) == doctest::Approx(4.0f));
    p.tStop = 4.0f;  CHECK(mechanicalFraction(p, 1.0f) == doctest::Approx(1.0f).epsilon(1e-4));
    p.tStop = 8.0f;  CHECK(mechanicalFraction(p, 1.0f) == doctest::Approx(1.0f).epsilon(1e-4));
    p.tStop = 2.8f;  CHECK(mechanicalFraction(p, 1.0f) < 1.0f);
}

TEST_CASE("stopping down monotonically reduces mechanical clipping") {
    VignetteParams p = defaults();
    float prev = 0.0f;
    for (float s : {2.0f, 2.4f, 2.8f, 3.4f, 4.0f}) {
        p.tStop = s;
        const float m = mechanicalFraction(p, 1.0f);
        CHECK(m >= prev - 1e-6f);
        prev = m;
    }
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/optics/vignette.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/optics/vignette.hpp`:
```cpp
#pragma once
#include <algorithm>
#include <cmath>

namespace lens::optics {

// Area of the lens-shaped intersection of two circles, centres d apart.
inline float circleOverlapArea(float d, float r1, float r2) {
    if (d >= r1 + r2) return 0.0f;
    const float rmin = std::min(r1, r2);
    if (d <= std::abs(r1 - r2)) return kPi * rmin * rmin;
    const float a1 = std::acos(std::clamp((d * d + r1 * r1 - r2 * r2) / (2.0f * d * r1), -1.0f, 1.0f));
    const float a2 = std::acos(std::clamp((d * d + r2 * r2 - r1 * r1) / (2.0f * d * r2), -1.0f, 1.0f));
    const float tri = 0.5f * std::sqrt(std::max(0.0f,
        (-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2)));
    return r1 * r1 * a1 + r2 * r2 * a2 - tri;
}

struct VignetteParams {
    float focal_mm          = 32.0f;
    float sensorHalfDiag_mm = 14.0f;
    float naturalExp        = 4.0f;   // free, not pinned at 4
    float tStopWide         = 2.0f;
    float tStop             = 2.0f;
    float rEntrance         = 1.0f;   // in units of the wide-open aperture radius
    float sepNorm           = 0.5f;   // pupil displacement per unit field radius
};

inline float naturalFalloff(const VignetteParams& p, float t) {
    const float theta = std::atan(t * p.sensorHalfDiag_mm / p.focal_mm);
    return std::pow(std::cos(theta), p.naturalExp);
}

// Aperture radius in units of the wide-open radius: 1.0 wide open, smaller stopped down.
inline float apertureRadius(const VignetteParams& p) { return p.tStopWide / p.tStop; }

inline float mechanicalVanishStop(const VignetteParams& p) {
    const float room = p.rEntrance - p.sepNorm;
    return (room > 0.0f) ? p.tStopWide / room : 1e9f;
}

inline float mechanicalFraction(const VignetteParams& p, float t) {
    const float a = apertureRadius(p);
    if (a <= 0.0f) return 1.0f;
    const float d = p.sepNorm * t;
    const float frac = circleOverlapArea(d, a, p.rEntrance) / (kPi * a * a);
    return std::clamp(frac, 0.0f, 1.0f);
}

inline float vignette(const VignetteParams& p, float t) {
    return naturalFalloff(p, t) * mechanicalFraction(p, t);
}

}  // namespace lens::optics
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all nine cases pass.

- [ ] **Step 5: Commit**

```bash
git add lenscore/optics/vignette.hpp tests/test_vignette.cpp tests/CMakeLists.txt
git commit -m "feat: natural and mechanical vignetting with measured f/4 vanishing"
```

---
### Task 10: Zernike wavefront basis

**Files:**
- Create: `lenscore/optics/zernike.hpp`, `tests/test_zernike.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `float lens::optics::zernikeRadial(int n, int m, float rho)`
  - `float lens::optics::zernike(int n, int m, float rho, float theta)` — orthonormal over the unit disc, `m < 0` selects the sine term
  - `struct lens::optics::Wavefront { float defocus, astig, coma, spherical; }`
  - `float lens::optics::wavefrontError(const Wavefront&, float rho, float theta)` — total optical path error in waves

The wavefront is what makes the PSF an aberrated one. In the Zernike basis the
four aberrations we model are `Z(2,0)` defocus, `Z(2,2)` astigmatism, `Z(3,1)`
coma and `Z(4,0)` spherical. Orientation is fixed in the **radial frame**:
`theta` is measured from the direction pointing away from the image centre, so
astigmatism separates sagittal from tangential and coma's tail points outward.
Fixing the frame this way is what makes the PSF depend on field radius alone,
which Task 15 exploits.

- [ ] **Step 1: Write the failing test**

`tests/test_zernike.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/optics/zernike.hpp"
#include <cmath>

using namespace lens::optics;

// Numerically integrate Z_a * Z_b over the unit disc, divided by the disc area.
static double innerProduct(int n1, int m1, int n2, int m2) {
    const int NR = 400, NT = 720;
    double acc = 0.0, norm = 0.0;
    for (int i = 0; i < NR; ++i) {
        const double rho = (i + 0.5) / NR;
        for (int j = 0; j < NT; ++j) {
            const double th = 2.0 * kPi * (j + 0.5) / NT;
            const double w = rho;                       // Jacobian
            acc  += w * zernike(n1, m1, float(rho), float(th)) * zernike(n2, m2, float(rho), float(th));
            norm += w;
        }
    }
    return acc / norm;
}

TEST_CASE("piston is unity everywhere") {
    CHECK(zernike(0, 0, 0.0f, 0.0f) == doctest::Approx(1.0f));
    CHECK(zernike(0, 0, 1.0f, 2.0f) == doctest::Approx(1.0f));
}

TEST_CASE("defocus has the known radial form") {
    CHECK(zernikeRadial(2, 0, 0.0f) == doctest::Approx(-1.0f));   // 2*rho^2 - 1
    CHECK(zernikeRadial(2, 0, 1.0f) == doctest::Approx(1.0f));
    CHECK(zernike(2, 0, 1.0f, 0.0f) == doctest::Approx(std::sqrt(3.0f)).epsilon(1e-4));
}

TEST_CASE("spherical has the known radial form") {
    // 6*rho^4 - 6*rho^2 + 1
    CHECK(zernikeRadial(4, 0, 0.0f) == doctest::Approx(1.0f));
    CHECK(zernikeRadial(4, 0, 1.0f) == doctest::Approx(1.0f));
    CHECK(zernikeRadial(4, 0, 0.5f) == doctest::Approx(6 * 0.0625f - 6 * 0.25f + 1.0f));
}

TEST_CASE("the basis is orthonormal over the unit disc") {
    const int modes[][2] = {{0,0}, {2,0}, {2,2}, {2,-2}, {3,1}, {3,-1}, {4,0}};
    for (const auto& a : modes)
        for (const auto& b : modes) {
            const double ip = innerProduct(a[0], a[1], b[0], b[1]);
            const double want = (a[0] == b[0] && a[1] == b[1]) ? 1.0 : 0.0;
            CAPTURE(a[0]); CAPTURE(a[1]); CAPTURE(b[0]); CAPTURE(b[1]);
            CHECK(ip == doctest::Approx(want).epsilon(0.02));
        }
}

TEST_CASE("astigmatism separates the sagittal and tangential axes") {
    CHECK(zernike(2, 2, 1.0f, 0.0f) > 0.0f);
    CHECK(zernike(2, 2, 1.0f, float(kPi) / 2.0f) < 0.0f);
}

TEST_CASE("coma is antisymmetric about the tangential axis") {
    const float a = zernike(3, 1, 0.8f, 0.0f);
    const float b = zernike(3, 1, 0.8f, float(kPi));
    CHECK(a == doctest::Approx(-b).epsilon(1e-5));
}

TEST_CASE("a zero wavefront is flat and a defocused one is not") {
    CHECK(wavefrontError(Wavefront{}, 0.7f, 1.0f) == doctest::Approx(0.0f));
    Wavefront w; w.defocus = 0.5f;
    CHECK(wavefrontError(w, 1.0f, 0.0f) != doctest::Approx(0.0f));
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/optics/zernike.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/optics/zernike.hpp`:
```cpp
#pragma once
#include <cmath>
#include <cstdlib>

namespace lens::optics {

inline double factorial(int n) {
    double f = 1.0;
    for (int i = 2; i <= n; ++i) f *= i;
    return f;
}

inline float zernikeRadial(int n, int m, float rho) {
    m = std::abs(m);
    if ((n - m) % 2 != 0) return 0.0f;
    double acc = 0.0;
    for (int k = 0; k <= (n - m) / 2; ++k) {
        const double num = ((k % 2) ? -1.0 : 1.0) * factorial(n - k);
        const double den = factorial(k) * factorial((n + m) / 2 - k) * factorial((n - m) / 2 - k);
        acc += num / den * std::pow(double(rho), n - 2 * k);
    }
    return float(acc);
}

// Orthonormal over the unit disc: the mean of Z^2 is exactly 1.
inline float zernike(int n, int m, float rho, float theta) {
    const float norm = std::sqrt(float(2 * (n + 1)) / (m == 0 ? 2.0f : 1.0f));
    const float r = zernikeRadial(n, m, rho);
    if (m == 0) return norm * r;
    return norm * r * (m > 0 ? std::cos(m * theta) : std::sin(-m * theta));
}

// Coefficients are in waves. theta is measured in the radial frame, from the
// direction pointing away from the image centre.
struct Wavefront {
    float defocus   = 0.0f;   // Z(2,0)
    float astig     = 0.0f;   // Z(2,2)
    float coma      = 0.0f;   // Z(3,1)
    float spherical = 0.0f;   // Z(4,0)
};

inline float wavefrontError(const Wavefront& w, float rho, float theta) {
    return w.defocus   * zernike(2, 0, rho, theta)
         + w.astig     * zernike(2, 2, rho, theta)
         + w.coma      * zernike(3, 1, rho, theta)
         + w.spherical * zernike(4, 0, rho, theta);
}

}  // namespace lens::optics
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all seven cases pass. The orthonormality case is the slow one, around
a second; that is acceptable for the guarantee it buys.

- [ ] **Step 5: Commit**

```bash
git add lenscore/optics/zernike.hpp tests/test_zernike.cpp tests/CMakeLists.txt
git commit -m "feat: orthonormal Zernike wavefront basis"
```

---
### Task 11: Pupil amplitude — aperture polygon, cat's-eye clip, apodization

**Files:**
- Create: `lenscore/optics/pupil.hpp`, `tests/test_pupil.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `lens::Plane` (Task 8), `lens::optics::VignetteParams`, `mechanicalFraction` (Task 9, tests only)
- Produces:
  - `float lens::optics::apertureEdgeRadius(int blades, float curvature, float rotationRad, float theta)`
  - `struct lens::optics::PupilParams { int blades; float curvature; float rotationRad; float apertureRadius; float rEntrance; float rExit; float sepNorm; float apodizationSlope; }`
  - `lens::Plane lens::optics::rasterPupil(const PupilParams&, float t, int N)` — amplitude on an `N x N` grid spanning `[-1,1]^2`
  - `float lens::optics::pupilEnergyFraction(const PupilParams&, float t, int N)`

Three things multiply into the pupil amplitude:

1. **The aperture polygon.** `blades` sides with per-blade `curvature`, where 1
   is a perfect circle and 0 is a straight-sided polygon.
2. **The cat's-eye clip.** Off axis, the front and rear barrel apertures appear
   displaced in *opposite* directions, by `sepNorm * t`. Their intersection with
   the aperture is the lemon shape that produces swirl. This is one geometric
   fact, not a swirl filter.
3. **The apodization ramp.** An off-axis point does not illuminate the pupil
   uniformly — Aggarwal measured 31% variation at a 10 degree field angle on a
   16mm lens. Modelled as their own planar approximation, `1 + slope * t * x`,
   with `x` the pupil coordinate along the radial direction. Its mean over the
   disc is exactly 1, so it redistributes energy without inventing any.

**`pupilEnergyFraction` is the optical vignetting term.** Nothing normalises it
away. Its agreement with Task 9's independent closed form is asserted below.

- [ ] **Step 1: Write the failing test**

`tests/test_pupil.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/optics/pupil.hpp"
#include "lenscore/optics/vignette.hpp"
#include <cmath>

using namespace lens;
using namespace lens::optics;

static PupilParams circular() {
    PupilParams p;
    p.blades = 0; p.curvature = 1.0f; p.rotationRad = 0.0f;
    p.apertureRadius = 1.0f; p.rEntrance = 1e6f; p.rExit = 1e6f;
    p.sepNorm = 0.0f; p.apodizationSlope = 0.0f;
    return p;
}

TEST_CASE("a circular aperture has unit edge radius in every direction") {
    for (float th = 0.0f; th < 6.28f; th += 0.3f)
        CHECK(apertureEdgeRadius(0, 1.0f, 0.0f, th) == doctest::Approx(1.0f));
}

TEST_CASE("a straight-sided hexagon touches 1 at vertices and cos(pi/6) at flats") {
    float lo = 2.0f, hi = 0.0f;
    for (float th = 0.0f; th < 6.283f; th += 0.001f) {
        const float r = apertureEdgeRadius(6, 0.0f, 0.0f, th);
        lo = std::min(lo, r); hi = std::max(hi, r);
    }
    CHECK(hi == doctest::Approx(1.0f).epsilon(1e-3));
    CHECK(lo == doctest::Approx(std::cos(3.14159265f / 6.0f)).epsilon(1e-3));
}

TEST_CASE("curvature interpolates the polygon back to a circle") {
    const float straight = apertureEdgeRadius(6, 0.0f, 0.0f, 3.14159265f / 6.0f);
    const float round    = apertureEdgeRadius(6, 1.0f, 0.0f, 3.14159265f / 6.0f);
    const float half     = apertureEdgeRadius(6, 0.5f, 0.0f, 3.14159265f / 6.0f);
    CHECK(round == doctest::Approx(1.0f));
    CHECK(half > straight);
    CHECK(half < round);
}

TEST_CASE("an unclipped circular pupil has full energy on axis") {
    CHECK(pupilEnergyFraction(circular(), 0.0f, 256) == doctest::Approx(1.0f).epsilon(0.01));
}

TEST_CASE("a hexagonal aperture passes the known polygon area ratio") {
    PupilParams p = circular();
    p.blades = 6; p.curvature = 0.0f;
    // regular hexagon in a unit circle: area = 3*sin(60) = 2.598, over pi = 0.827
    CHECK(pupilEnergyFraction(p, 0.0f, 512) == doctest::Approx(0.827f).epsilon(0.02));
}

TEST_CASE("the cat's-eye clip removes energy off axis and not on axis") {
    PupilParams p = circular();
    p.rEntrance = 1.0f; p.rExit = 1.0f; p.sepNorm = 0.5f;
    CHECK(pupilEnergyFraction(p, 0.0f, 256) == doctest::Approx(1.0f).epsilon(0.01));
    CHECK(pupilEnergyFraction(p, 1.0f, 256) < 0.95f);
}

TEST_CASE("stopping down removes the cat's-eye clip entirely") {
    PupilParams p = circular();
    p.rEntrance = 1.0f; p.rExit = 1.0f; p.sepNorm = 0.5f;
    p.apertureRadius = 0.5f;                       // f/4 with a wide-open f/2
    CHECK(pupilEnergyFraction(p, 1.0f, 256) == doctest::Approx(1.0f).epsilon(0.01));
}

TEST_CASE("rasterised pupil energy agrees with the independent closed form") {
    // One clipping circle only, so the two models describe the same geometry.
    VignetteParams v;
    v.tStopWide = 2.0f; v.rEntrance = 1.0f; v.sepNorm = 0.5f;
    for (float stop : {2.0f, 2.4f, 2.8f, 3.2f, 4.0f}) {
        v.tStop = stop;
        PupilParams p = circular();
        p.apertureRadius = v.tStopWide / v.tStop;
        p.rEntrance = v.rEntrance; p.rExit = 1e6f; p.sepNorm = v.sepNorm;
        CAPTURE(stop);
        CHECK(pupilEnergyFraction(p, 1.0f, 512) ==
              doctest::Approx(mechanicalFraction(v, 1.0f)).epsilon(0.02));
    }
}

TEST_CASE("apodization shifts the pupil centroid outward without adding energy") {
    PupilParams flat = circular();
    PupilParams ramp = circular(); ramp.apodizationSlope = 0.31f;
    const Plane a = rasterPupil(flat, 1.0f, 256);
    const Plane b = rasterPupil(ramp, 1.0f, 256);
    double sa = 0, sb = 0, cx = 0;
    for (int y = 0; y < 256; ++y)
        for (int x = 0; x < 256; ++x) {
            const double u = (2.0 * x) / 255.0 - 1.0;
            sa += a.at(x, y); sb += b.at(x, y); cx += b.at(x, y) * u;
        }
    CHECK(sb / sa == doctest::Approx(1.0).epsilon(0.01));   // energy preserved
    CHECK(cx / sb > 0.02);                                   // centroid moved outward
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/optics/pupil.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/optics/pupil.hpp`:
```cpp
#pragma once
#include "lenscore/plane.hpp"
#include <algorithm>
#include <cmath>

namespace lens::optics {

// Distance from centre to the aperture edge, as a fraction of the circumscribed
// radius. curvature 1 gives a circle, 0 a straight-sided polygon.
inline float apertureEdgeRadius(int blades, float curvature, float rotationRad, float theta) {
    if (blades < 3) return 1.0f;
    const float seg = 2.0f * kPi / float(blades);
    float a = std::fmod(theta + rotationRad, seg);
    if (a < 0.0f) a += seg;
    const float poly = std::cos(kPi / float(blades)) / std::cos(a - kPi / float(blades));
    return std::clamp(curvature, 0.0f, 1.0f) * 1.0f + (1.0f - std::clamp(curvature, 0.0f, 1.0f)) * poly;
}

struct PupilParams {
    int   blades           = 9;
    float curvature        = 0.15f;
    float rotationRad      = 0.0f;
    float apertureRadius   = 1.0f;   // in units of the wide-open radius
    float rEntrance        = 1.0f;
    float rExit            = 0.92f;
    float sepNorm          = 0.35f;
    float apodizationSlope = 0.0f;
};

// Amplitude on an N x N grid spanning [-1,1]^2. The radial direction is +x.
inline Plane rasterPupil(const PupilParams& p, float t, int N) {
    Plane out(N, N);
    const float d = p.sepNorm * t;
    for (int j = 0; j < N; ++j) {
        const float v = 2.0f * float(j) / float(N - 1) - 1.0f;
        for (int i = 0; i < N; ++i) {
            const float u = 2.0f * float(i) / float(N - 1) - 1.0f;
            const float rho = std::sqrt(u * u + v * v);
            if (rho > p.apertureRadius) continue;

            const float th = std::atan2(v, u);
            if (rho > p.apertureRadius * apertureEdgeRadius(p.blades, p.curvature, p.rotationRad, th))
                continue;

            // Front and rear barrel apertures displace in opposite directions.
            if (std::hypot(u - d, v) > p.rEntrance) continue;
            if (std::hypot(u + d, v) > p.rExit)     continue;

            out.at(i, j) = std::max(0.0f, 1.0f + p.apodizationSlope * t * u);
        }
    }
    return out;
}

// The optical vignetting term. Never normalised away.
inline float pupilEnergyFraction(const PupilParams& p, float t, int N) {
    const Plane amp = rasterPupil(p, t, N);
    double got = 0.0, full = 0.0;
    for (int j = 0; j < N; ++j) {
        const float v = 2.0f * float(j) / float(N - 1) - 1.0f;
        for (int i = 0; i < N; ++i) {
            const float u = 2.0f * float(i) / float(N - 1) - 1.0f;
            if (std::sqrt(u * u + v * v) <= p.apertureRadius) full += 1.0;
            got += amp.at(i, j);
        }
    }
    return (full > 0.0) ? float(got / full) : 0.0f;
}

}  // namespace lens::optics
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all nine cases pass. The consistency case is the important one: it
proves the rasterised pupil and the closed-form vignette describe the same lens.

- [ ] **Step 5: Commit**

```bash
git add lenscore/optics/pupil.hpp tests/test_pupil.cpp tests/CMakeLists.txt
git commit -m "feat: pupil amplitude with cat's-eye clip and Aggarwal apodization"
```

---
### Task 12: Complex FFT

**Files:**
- Create: `lenscore/conv/fft.hpp`, `tests/test_fft.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `using lens::conv::Cplx = std::complex<float>`
  - `void lens::conv::fft1d(std::vector<Cplx>&, bool inverse)`
  - `void lens::conv::fft2d(std::vector<Cplx>&, int w, int h, bool inverse)`
  - `void lens::conv::fftShift2d(std::vector<Cplx>&, int w, int h)`
  - `bool lens::conv::isPowerOfTwo(int)`

Radix-2 Cooley-Tukey, in place. Every size this library uses is a power of two —
PSF grids and EFF patch sizes are both chosen that way — so higher radices buy
nothing. The inverse divides by `N`, so `fft2d(inverse)` after `fft2d(forward)`
is the identity.

- [ ] **Step 1: Write the failing test**

`tests/test_fft.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/conv/fft.hpp"
#include <cmath>

using namespace lens::conv;

TEST_CASE("power of two detection") {
    CHECK(isPowerOfTwo(1));  CHECK(isPowerOfTwo(64));
    CHECK_FALSE(isPowerOfTwo(0)); CHECK_FALSE(isPowerOfTwo(48));
}

TEST_CASE("transform of a delta is flat") {
    std::vector<Cplx> v(16, Cplx(0, 0));
    v[0] = Cplx(1, 0);
    fft1d(v, false);
    for (const Cplx& c : v) {
        CHECK(c.real() == doctest::Approx(1.0f).epsilon(1e-5));
        CHECK(c.imag() == doctest::Approx(0.0f).epsilon(1e-5));
    }
}

TEST_CASE("transform of a constant is a delta at DC") {
    std::vector<Cplx> v(16, Cplx(1, 0));
    fft1d(v, false);
    CHECK(std::abs(v[0]) == doctest::Approx(16.0f).epsilon(1e-4));
    for (size_t i = 1; i < v.size(); ++i) CHECK(std::abs(v[i]) == doctest::Approx(0.0f).epsilon(1e-4));
}

TEST_CASE("forward then inverse is the identity in 1d") {
    std::vector<Cplx> v(32), original;
    for (size_t i = 0; i < v.size(); ++i) v[i] = Cplx(std::sin(0.3f * i), std::cos(0.7f * i));
    original = v;
    fft1d(v, false); fft1d(v, true);
    for (size_t i = 0; i < v.size(); ++i) CHECK(std::abs(v[i] - original[i]) < 1e-4f);
}

TEST_CASE("forward then inverse is the identity in 2d") {
    const int w = 16, h = 8;
    std::vector<Cplx> v(size_t(w) * h), original;
    for (size_t i = 0; i < v.size(); ++i) v[i] = Cplx(float(i % 7) - 3.0f, float(i % 5) - 2.0f);
    original = v;
    fft2d(v, w, h, false); fft2d(v, w, h, true);
    for (size_t i = 0; i < v.size(); ++i) CHECK(std::abs(v[i] - original[i]) < 1e-3f);
}

TEST_CASE("Parseval holds: energy is preserved up to the size factor") {
    const int w = 16, h = 16;
    std::vector<Cplx> v(size_t(w) * h);
    for (size_t i = 0; i < v.size(); ++i) v[i] = Cplx(float((i * 37) % 11) / 11.0f, 0.0f);
    double spatial = 0.0;
    for (const Cplx& c : v) spatial += std::norm(c);
    fft2d(v, w, h, false);
    double freq = 0.0;
    for (const Cplx& c : v) freq += std::norm(c);
    CHECK(freq == doctest::Approx(spatial * w * h).epsilon(1e-4));
}

TEST_CASE("fftShift moves DC to the centre and is its own inverse for even sizes") {
    const int w = 8, h = 8;
    std::vector<Cplx> v(size_t(w) * h, Cplx(0, 0));
    v[0] = Cplx(1, 0);
    fftShift2d(v, w, h);
    CHECK(std::abs(v[size_t(h / 2) * w + w / 2]) == doctest::Approx(1.0f));
    fftShift2d(v, w, h);
    CHECK(std::abs(v[0]) == doctest::Approx(1.0f));
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/conv/fft.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/conv/fft.hpp`:
```cpp
#pragma once
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

namespace lens::conv {

using Cplx = std::complex<float>;

inline bool isPowerOfTwo(int n) { return n > 0 && (n & (n - 1)) == 0; }

inline void fft1d(std::vector<Cplx>& a, bool inverse) {
    const int n = int(a.size());
    assert(isPowerOfTwo(n));

    for (int i = 1, j = 0; i < n; ++i) {          // bit-reversal permutation
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        const double ang = 2.0 * double(kPi) / len * (inverse ? 1.0 : -1.0);
        const Cplx wl(float(std::cos(ang)), float(std::sin(ang)));
        for (int i = 0; i < n; i += len) {
            Cplx w(1.0f, 0.0f);
            for (int k = 0; k < len / 2; ++k) {
                const Cplx u = a[i + k];
                const Cplx v = a[i + k + len / 2] * w;
                a[i + k]           = u + v;
                a[i + k + len / 2] = u - v;
                w *= wl;
            }
        }
    }
    if (inverse) for (Cplx& c : a) c /= float(n);
}

inline void fft2d(std::vector<Cplx>& a, int w, int h, bool inverse) {
    assert(isPowerOfTwo(w) && isPowerOfTwo(h) && a.size() == size_t(w) * h);
    std::vector<Cplx> row(w), col(h);
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) row[x] = a[size_t(y) * w + x];
        fft1d(row, inverse);
        for (int x = 0; x < w; ++x) a[size_t(y) * w + x] = row[x];
    }
    for (int x = 0; x < w; ++x) {
        for (int y = 0; y < h; ++y) col[y] = a[size_t(y) * w + x];
        fft1d(col, inverse);
        for (int y = 0; y < h; ++y) a[size_t(y) * w + x] = col[y];
    }
}

inline void fftShift2d(std::vector<Cplx>& a, int w, int h) {
    const int hx = w / 2, hy = h / 2;
    for (int y = 0; y < hy; ++y)
        for (int x = 0; x < w; ++x) {
            const int sx = (x + hx) % w;
            std::swap(a[size_t(y) * w + x], a[size_t(y + hy) * w + sx]);
        }
}

}  // namespace lens::conv
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all seven cases pass.

- [ ] **Step 5: Commit**

```bash
git add lenscore/conv/fft.hpp tests/test_fft.cpp tests/CMakeLists.txt
git commit -m "feat: radix-2 FFT, 1d and 2d, with shift"
```

---
### Task 13: The PSF, from one Fourier transform

**Files:**
- Create: `lenscore/optics/psf.hpp`, `tests/test_psf.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `rasterPupil`, `PupilParams` (Task 11), `Wavefront`, `wavefrontError` (Task 10), `fft2d`, `fftShift2d` (Task 12)
- Produces:
  - `lens::Plane lens::optics::psfFromPupil(const PupilParams&, const Wavefront&, float lambdaNm, float lambdaRefNm, float t, int N)`
  - `float lens::optics::airyFirstZeroSamples(int N, float apertureRadius)`
  - `float lens::optics::psfSampleSpacingUm(float lambdaNm, float fNumberWide)`

The whole aberration model collapses into one line of physics:

```
PSF = | FFT{ A(x,y) * exp(i * 2*pi * W(rho,theta) * lambdaRef/lambda) } |^2
```

`A` is the pupil amplitude from Task 11 — polygon, cat's-eye clip, apodization.
`W` is the Zernike wavefront from Task 10, in waves at the reference wavelength;
the `lambdaRef/lambda` factor is why an error fixed in *length* costs more waves
at short wavelengths, which is most of what makes blue behave differently from
red. Defocus, astigmatism, coma, spherical, the aperture star, the soft disc
edge, the lemon shape and diffraction are all consequences. There are no special
cases, and no way for the parts to contradict each other.

**Two scale facts worth stating, because they are easy to get wrong.**

The first Airy zero sits at `1.22 * N / (2 * R_s)` output samples, where `R_s =
apertureRadius * N/2` is the pupil radius in grid samples. Stopping down shrinks
`R_s`, which widens the PSF in samples — diffraction, for free, in the right
direction.

The image-plane spacing of one output sample is `lambda * fNumberWide` — note
what is *not* in it. Because `apertureRadius = fNumberWide/fNumber`, the stop
cancels: the sampling is fixed while the pattern grows. Getting this backwards
makes stopping down sharpen the diffraction blur, which is the opposite of
physics.

- [ ] **Step 1: Write the failing test**

`tests/test_psf.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/optics/psf.hpp"
#include <cmath>
#include <vector>

using namespace lens;
using namespace lens::optics;

static PupilParams smallDisc() {
    PupilParams p;
    p.blades = 0; p.curvature = 1.0f; p.apertureRadius = 0.125f;
    p.rEntrance = 1e6f; p.rExit = 1e6f; p.sepNorm = 0.0f; p.apodizationSlope = 0.0f;
    return p;
}

// Mean intensity on a ring of the given radius about the PSF centre.
static double ringMean(const Plane& p, double r) {
    const double cx = p.w / 2.0, cy = p.h / 2.0;
    double acc = 0.0; int n = 0;
    for (int k = 0; k < 720; ++k) {
        const double a = 2.0 * kPi * k / 720.0;
        const int x = int(std::lround(cx + r * std::cos(a)));
        const int y = int(std::lround(cy + r * std::sin(a)));
        if (x >= 0 && y >= 0 && x < p.w && y < p.h) { acc += p.at(x, y); ++n; }
    }
    return n ? acc / n : 0.0;
}

TEST_CASE("sample spacing depends on the wide-open stop, not the working stop") {
    CHECK(psfSampleSpacingUm(550.0f, 2.0f) == doctest::Approx(1.1f).epsilon(1e-4));
    CHECK(psfSampleSpacingUm(550.0f, 4.0f) == doctest::Approx(2.2f).epsilon(1e-4));
}

TEST_CASE("an unaberrated PSF peaks dead centre") {
    const int N = 256;
    const Plane psf = psfFromPupil(smallDisc(), Wavefront{}, 550.0f, 550.0f, 0.0f, N);
    float best = -1.0f; int bx = 0, by = 0;
    for (int y = 0; y < N; ++y) for (int x = 0; x < N; ++x)
        if (psf.at(x, y) > best) { best = psf.at(x, y); bx = x; by = y; }
    CHECK(bx == N / 2);
    CHECK(by == N / 2);
}

TEST_CASE("the first Airy zero lands where diffraction says it should") {
    const int N = 256;
    const PupilParams p = smallDisc();
    const Plane psf = psfFromPupil(p, Wavefront{}, 550.0f, 550.0f, 0.0f, N);
    const double predicted = airyFirstZeroSamples(N, p.apertureRadius);   // 9.76
    CHECK(predicted == doctest::Approx(9.76).epsilon(0.01));

    double bestR = 0.0, bestV = 1e30;
    for (double r = 3.0; r < 18.0; r += 0.25) {
        const double v = ringMean(psf, r);
        if (v < bestV) { bestV = v; bestR = r; }
    }
    CHECK(bestR == doctest::Approx(predicted).epsilon(0.12));
}

TEST_CASE("stopping down widens the diffraction pattern") {
    const int N = 256;
    PupilParams wide = smallDisc();
    PupilParams stop = smallDisc(); stop.apertureRadius = 0.0625f;
    CHECK(airyFirstZeroSamples(N, stop.apertureRadius) >
          airyFirstZeroSamples(N, wide.apertureRadius));
}

TEST_CASE("total PSF energy tracks the pupil energy, so vignetting survives") {
    const int N = 256;
    PupilParams p = smallDisc();
    p.rEntrance = 0.13f; p.rExit = 1e6f; p.sepNorm = 0.10f;    // clips off axis

    auto energy = [&](float t) {
        const Plane psf = psfFromPupil(p, Wavefront{}, 550.0f, 550.0f, t, N);
        double s = 0.0; for (float v : psf.v) s += v; return s;
    };
    auto pupilSq = [&](float t) {
        const Plane a = rasterPupil(p, t, N);
        double s = 0.0; for (float v : a.v) s += double(v) * v; return s;
    };
    CHECK(energy(1.0f) / energy(0.0f) == doctest::Approx(pupilSq(1.0f) / pupilSq(0.0f)).epsilon(1e-3));
    CHECK(energy(1.0f) < energy(0.0f));      // it really did lose light
}

TEST_CASE("defocus broadens the PSF") {
    const int N = 256;
    auto secondMoment = [&](const Wavefront& w) {
        const Plane psf = psfFromPupil(smallDisc(), w, 550.0f, 550.0f, 0.0f, N);
        double s = 0.0, m = 0.0;
        for (int y = 0; y < N; ++y) for (int x = 0; x < N; ++x) {
            const double dx = x - N / 2.0, dy = y - N / 2.0;
            s += psf.at(x, y); m += psf.at(x, y) * (dx * dx + dy * dy);
        }
        return m / s;
    };
    Wavefront d; d.defocus = 1.5f;
    CHECK(secondMoment(d) > secondMoment(Wavefront{}));
}

TEST_CASE("astigmatism makes the PSF elliptical") {
    const int N = 256;
    Wavefront w; w.astig = 1.2f;
    const Plane psf = psfFromPupil(smallDisc(), w, 550.0f, 550.0f, 0.0f, N);
    double s = 0, mx = 0, my = 0;
    for (int y = 0; y < N; ++y) for (int x = 0; x < N; ++x) {
        const double dx = x - N / 2.0, dy = y - N / 2.0;
        s += psf.at(x, y); mx += psf.at(x, y) * dx * dx; my += psf.at(x, y) * dy * dy;
    }
    CHECK(std::abs(mx / s - my / s) / (mx / s) > 0.10);
}

TEST_CASE("coma throws the PSF off centre along the radial axis") {
    const int N = 256;
    Wavefront w; w.coma = 1.2f;
    const Plane psf = psfFromPupil(smallDisc(), w, 550.0f, 550.0f, 1.0f, N);
    double s = 0, cx = 0, cy = 0;
    for (int y = 0; y < N; ++y) for (int x = 0; x < N; ++x) {
        s += psf.at(x, y);
        cx += psf.at(x, y) * (x - N / 2.0);
        cy += psf.at(x, y) * (y - N / 2.0);
    }
    CHECK(std::abs(cx / s) > 0.5);                     // shifted along +/-x
    CHECK(std::abs(cy / s) < 0.1 * std::abs(cx / s));  // not sideways
}

TEST_CASE("an even blade count gives that many diffraction spikes") {
    const int N = 512;
    PupilParams p = smallDisc();
    p.blades = 6; p.curvature = 0.0f;
    const Plane psf = psfFromPupil(p, Wavefront{}, 550.0f, 550.0f, 0.0f, N);

    std::vector<double> ang(360);
    for (int k = 0; k < 360; ++k) {
        const double a = 2.0 * kPi * k / 360.0;
        double acc = 0.0;
        for (double r = 14.0; r < 34.0; r += 1.0) {
            const int x = int(std::lround(N / 2.0 + r * std::cos(a)));
            const int y = int(std::lround(N / 2.0 + r * std::sin(a)));
            if (x >= 0 && y >= 0 && x < N && y < N) acc += psf.at(x, y);
        }
        ang[k] = acc;
    }
    double mean = 0.0; for (double v : ang) mean += v; mean /= 360.0;
    int peaks = 0;
    for (int k = 0; k < 360; ++k) {
        const double p0 = ang[(k + 359) % 360], p1 = ang[k], p2 = ang[(k + 1) % 360];
        if (p1 > p0 && p1 >= p2 && p1 > 1.5 * mean) ++peaks;
    }
    CHECK(peaks == 6);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/optics/psf.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/optics/psf.hpp`:
```cpp
#pragma once
#include "lenscore/conv/fft.hpp"
#include "lenscore/optics/pupil.hpp"
#include "lenscore/optics/zernike.hpp"
#include <cmath>

namespace lens::optics {

// Radius of the first Airy minimum, in output samples.
inline float airyFirstZeroSamples(int N, float apertureRadius) {
    const float Rs = apertureRadius * float(N) * 0.5f;
    return 1.22f * float(N) / (2.0f * Rs);
}

// Image-plane size of one output sample, in micrometres. The working stop
// cancels out: stopping down widens the pattern, it does not resample it.
inline float psfSampleSpacingUm(float lambdaNm, float fNumberWide) {
    return lambdaNm * 1e-3f * fNumberWide;
}

// PSF = |FFT{ A * exp(i 2 pi W lambdaRef/lambda) }|^2, fftshifted, unnormalised.
inline Plane psfFromPupil(const PupilParams& pp, const Wavefront& wf,
                          float lambdaNm, float lambdaRefNm, float t, int N) {
    const Plane amp = rasterPupil(pp, t, N);
    const float chroma = lambdaRefNm / lambdaNm;

    std::vector<conv::Cplx> field(size_t(N) * N, conv::Cplx(0.0f, 0.0f));
    for (int j = 0; j < N; ++j) {
        const float v = 2.0f * float(j) / float(N - 1) - 1.0f;
        for (int i = 0; i < N; ++i) {
            const float a = amp.at(i, j);
            if (a <= 0.0f) continue;
            const float u = 2.0f * float(i) / float(N - 1) - 1.0f;
            const float rho = std::sqrt(u * u + v * v) / pp.apertureRadius;
            const float th  = std::atan2(v, u);
            const float phase = kTwoPi * wavefrontError(wf, rho, th) * chroma;
            field[size_t(j) * N + i] = conv::Cplx(a * std::cos(phase), a * std::sin(phase));
        }
    }

    conv::fft2d(field, N, N, false);
    conv::fftShift2d(field, N, N);

    Plane out(N, N);
    for (size_t k = 0; k < field.size(); ++k) out.v[k] = std::norm(field[k]);
    return out;   // deliberately not normalised: the energy is the vignetting
}

}  // namespace lens::optics
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all nine cases pass. If the blade-spike count comes out as 12 rather
than 6, the peak threshold is picking up the secondary lobes between spikes —
raise `1.5 * mean`, do not change the physics.

- [ ] **Step 5: Commit**

```bash
git add lenscore/optics/psf.hpp tests/test_psf.cpp tests/CMakeLists.txt
git commit -m "feat: PSF from the generalized pupil function"
```

---
### Task 14: Efficient Filter Flow — spatially varying convolution

**Files:**
- Create: `lenscore/conv/eff.hpp`, `tests/test_eff.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `lens::Plane` (Task 8), `fft2d` (Task 12)
- Produces:
  - `std::vector<float> lens::conv::hannWindow(int n)`
  - `lens::Plane lens::conv::convolveDirect(const Plane& src, const Plane& kernel)` — reference implementation, tests only
  - `lens::Plane lens::conv::effConvolve(const Plane& src, int patch, const std::function<Plane(float, float)>& psfAt)`

`psfAt(cx, cy)` returns the kernel for a patch centred on that pixel. Patches
are square, `patch` wide, stepped by `patch/2`. Each is windowed by a separable
Hann, transformed, multiplied by the transformed kernel, inverse transformed and
accumulated. Hann windows at 50% overlap **sum to exactly 1**, which is what
makes the patches blend without a seam and makes the whole thing reduce to plain
convolution when every kernel is the same. That reduction is the first test,
because if it does not hold nothing downstream can be trusted.

Cost is `O(N log S)` rather than `O(N K^2)`.

- [ ] **Step 1: Write the failing test**

`tests/test_eff.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/conv/eff.hpp"
#include <cmath>

using namespace lens;
using namespace lens::conv;

static Plane noise(int w, int h) {
    Plane p(w, h);
    unsigned s = 12345u;
    for (auto& v : p.v) { s = s * 1664525u + 1013904223u; v = float(s >> 8 & 0xFFFF) / 65535.0f; }
    return p;
}

static Plane deltaKernel(int k) { Plane p(k, k); p.at(k / 2, k / 2) = 1.0f; return p; }

static Plane boxKernel(int k) {
    Plane p(k, k);
    for (auto& v : p.v) v = 1.0f / float(k * k);
    return p;
}

TEST_CASE("Hann windows at fifty percent overlap sum to one") {
    const int P = 32;
    const std::vector<float> w = hannWindow(P);
    for (int i = 0; i < P / 2; ++i)
        CHECK(w[i] + w[i + P / 2] == doctest::Approx(1.0f).epsilon(1e-5));
}

TEST_CASE("a delta kernel everywhere is the identity") {
    const Plane src = noise(64, 64);
    const Plane out = effConvolve(src, 32, [](float, float) { return deltaKernel(9); });
    for (int y = 8; y < 56; ++y)
        for (int x = 8; x < 56; ++x)
            CHECK(out.at(x, y) == doctest::Approx(src.at(x, y)).epsilon(1e-3));
}

TEST_CASE("a uniform kernel reduces exactly to plain convolution") {
    const Plane src = noise(64, 64);
    const Plane k = boxKernel(7);
    const Plane want = convolveDirect(src, k);
    const Plane got  = effConvolve(src, 32, [&](float, float) { return k; });
    for (int y = 12; y < 52; ++y)
        for (int x = 12; x < 52; ++x)
            CHECK(got.at(x, y) == doctest::Approx(want.at(x, y)).epsilon(2e-3));
}

TEST_CASE("a normalised kernel preserves energy in the interior") {
    Plane src(64, 64);
    for (int y = 16; y < 48; ++y) for (int x = 16; x < 48; ++x) src.at(x, y) = 1.0f;
    const Plane out = effConvolve(src, 32, [](float, float) { return boxKernel(5); });
    double a = 0, b = 0;
    for (float v : src.v) a += v;
    for (float v : out.v) b += v;
    CHECK(b / a == doctest::Approx(1.0).epsilon(0.01));
}

TEST_CASE("a varying kernel really does vary across the frame") {
    Plane src(128, 128);
    for (int y = 0; y < 128; ++y)
        for (int x = 0; x < 128; ++x) src.at(x, y) = ((x / 4) % 2) ? 1.0f : 0.0f;   // vertical bars

    // Sharp on the left, blurry on the right.
    const Plane out = effConvolve(src, 32, [](float cx, float) {
        return (cx < 64.0f) ? deltaKernel(11) : boxKernel(11);
    });

    auto contrast = [&](const Plane& p, int x0, int x1) {
        double lo = 1e9, hi = -1e9;
        for (int y = 40; y < 88; ++y) for (int x = x0; x < x1; ++x) {
            lo = std::min(lo, double(p.at(x, y))); hi = std::max(hi, double(p.at(x, y)));
        }
        return hi - lo;
    };
    CHECK(contrast(out, 10, 40) > 0.9);
    CHECK(contrast(out, 90, 120) < 0.5);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/conv/eff.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/conv/eff.hpp`:
```cpp
#pragma once
#include "lenscore/conv/fft.hpp"
#include "lenscore/geometry.hpp"
#include "lenscore/plane.hpp"
#include <cmath>
#include <functional>

namespace lens::conv {

// Periodic Hann: at 50% overlap the shifted copies sum to exactly 1.
inline std::vector<float> hannWindow(int n) {
    std::vector<float> w(n);
    for (int i = 0; i < n; ++i) w[i] = 0.5f * (1.0f - std::cos(2.0f * kPi * float(i) / float(n)));
    return w;
}

inline Plane convolveDirect(const Plane& src, const Plane& k) {
    Plane out(src.w, src.h);
    const int kx = k.w / 2, ky = k.h / 2;
    for (int y = 0; y < src.h; ++y)
        for (int x = 0; x < src.w; ++x) {
            float acc = 0.0f;
            for (int j = 0; j < k.h; ++j)
                for (int i = 0; i < k.w; ++i) {
                    const int sx = std::clamp(x + i - kx, 0, src.w - 1);
                    const int sy = std::clamp(y + j - ky, 0, src.h - 1);
                    acc += src.at(sx, sy) * k.at(i, j);
                }
            out.at(x, y) = acc;
        }
    return out;
}

inline int nextPow2(int n) { int p = 1; while (p < n) p <<= 1; return p; }

inline Plane effConvolve(const Plane& src, int patch,
                         const std::function<Plane(float, float)>& psfAt) {
    const int hop = patch / 2;
    const Plane probe = psfAt(src.w * 0.5f, src.h * 0.5f);
    const int S = nextPow2(patch + std::max(probe.w, probe.h));

    const std::vector<float> win = hannWindow(patch);
    Plane out(src.w, src.h);
    std::vector<Cplx> buf(size_t(S) * S), ker(size_t(S) * S);

    // Origins start one hop before the image so every pixel is covered by full windows.
    for (int oy = -hop; oy < src.h; oy += hop) {
        for (int ox = -hop; ox < src.w; ox += hop) {
            const Plane k = psfAt(float(ox) + patch * 0.5f, float(oy) + patch * 0.5f);

            std::fill(buf.begin(), buf.end(), Cplx(0.0f, 0.0f));
            for (int j = 0; j < patch; ++j)
                for (int i = 0; i < patch; ++i) {
                    const int sx = std::clamp(ox + i, 0, src.w - 1);
                    const int sy = std::clamp(oy + j, 0, src.h - 1);
                    buf[size_t(j) * S + i] = Cplx(src.at(sx, sy) * win[i] * win[j], 0.0f);
                }

            // Kernel centred on the wraparound origin, so convolution does not shift.
            std::fill(ker.begin(), ker.end(), Cplx(0.0f, 0.0f));
            for (int j = 0; j < k.h; ++j)
                for (int i = 0; i < k.w; ++i) {
                    const int wx = ((i - k.w / 2) % S + S) % S;
                    const int wy = ((j - k.h / 2) % S + S) % S;
                    ker[size_t(wy) * S + wx] += Cplx(k.at(i, j), 0.0f);
                }

            fft2d(buf, S, S, false);
            fft2d(ker, S, S, false);
            for (size_t i = 0; i < buf.size(); ++i) buf[i] *= ker[i];
            fft2d(buf, S, S, true);

            for (int j = 0; j < S; ++j) {
                const int dy = oy + j;
                if (dy < 0 || dy >= src.h) continue;
                for (int i = 0; i < S; ++i) {
                    const int dx = ox + i;
                    if (dx < 0 || dx >= src.w) continue;
                    out.at(dx, dy) += buf[size_t(j) * S + i].real();
                }
            }
        }
    }
    return out;
}

}  // namespace lens::conv
```

Add `#include <algorithm>` alongside the existing includes.

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all five cases pass. The uniform-kernel case is the load-bearing one.

- [ ] **Step 5: Commit**

```bash
git add lenscore/conv/eff.hpp tests/test_eff.cpp tests/CMakeLists.txt
git commit -m "feat: Efficient Filter Flow spatially varying convolution"
```

---
### Task 15: Polar PSF rings and rotational covariance

**Files:**
- Create: `lenscore/optics/psfrings.hpp`, `tests/test_psfrings.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `psfFromPupil`, `PupilParams`, `Wavefront` (Tasks 11, 13)
- Produces:
  - `struct lens::optics::PsfRings { std::vector<Plane> ring; int gridN; }`
  - `PsfRings lens::optics::buildPsfRings(const PupilParams&, const std::function<Wavefront(float t)>&, float lambdaNm, float lambdaRefNm, int rings, int N)`
  - `Plane lens::optics::psfAtField(const PsfRings&, float t, float thetaRad, float samplesPerPixel, int outSize)`

**The collapse that makes this affordable.** Because Task 10 fixes the Zernike
frame to the radial direction, the PSF's *shape* depends only on field radius.
At any other angle it is the same kernel, rotated. So we compute one PSF per
radius ring — not per patch — and rotate it into place. With 16 rings and 24
angular divisions that is 16 Fourier transforms instead of 384, and exact
rotational symmetry stops being something to test for and becomes something the
structure cannot violate.

`samplesPerPixel` bridges the PSF grid to image pixels: it is
`psfSampleSpacingUm(lambda, fNumberWide) / pixelPitchUm`.

- [ ] **Step 1: Write the failing test**

`tests/test_psfrings.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/optics/psfrings.hpp"
#include <cmath>

using namespace lens;
using namespace lens::optics;

static PupilParams disc() {
    PupilParams p;
    p.blades = 0; p.curvature = 1.0f; p.apertureRadius = 0.125f;
    p.rEntrance = 1e6f; p.rExit = 1e6f; p.sepNorm = 0.0f; p.apodizationSlope = 0.0f;
    return p;
}

static void centroid(const Plane& p, double& cx, double& cy) {
    double s = 0; cx = 0; cy = 0;
    for (int y = 0; y < p.h; ++y) for (int x = 0; x < p.w; ++x) {
        s  += p.at(x, y);
        cx += p.at(x, y) * (x - p.w / 2.0);
        cy += p.at(x, y) * (y - p.h / 2.0);
    }
    if (s > 0) { cx /= s; cy /= s; }
}

TEST_CASE("rings are built at the requested count") {
    const PsfRings r = buildPsfRings(disc(), [](float) { return Wavefront{}; },
                                     550.0f, 550.0f, 8, 128);
    CHECK(r.ring.size() == 8u);
    CHECK(r.gridN == 128);
}

TEST_CASE("a rotationally symmetric aberration gives the same PSF at every angle") {
    const PsfRings r = buildPsfRings(disc(), [](float t) {
        Wavefront w; w.defocus = 1.0f * t * t; return w; }, 550.0f, 550.0f, 8, 128);
    const Plane a = psfAtField(r, 0.8f, 0.0f,             1.0f, 33);
    const Plane b = psfAtField(r, 0.8f, 1.0471975f,       1.0f, 33);   // 60 degrees
    const Plane c = psfAtField(r, 0.8f, 2.7f,             1.0f, 33);
    for (int i = 0; i < 33 * 33; ++i) {
        CHECK(b.v[i] == doctest::Approx(a.v[i]).epsilon(0.02));
        CHECK(c.v[i] == doctest::Approx(a.v[i]).epsilon(0.02));
    }
}

TEST_CASE("coma's tail follows the radial direction round the frame") {
    const PsfRings r = buildPsfRings(disc(), [](float t) {
        Wavefront w; w.coma = 1.5f * t; return w; }, 550.0f, 550.0f, 8, 256);

    double cx = 0, cy = 0;
    centroid(psfAtField(r, 1.0f, 0.0f, 1.0f, 65), cx, cy);
    const double along = std::abs(cx);
    CHECK(along > 0.3);
    CHECK(std::abs(cy) < 0.25 * along);

    centroid(psfAtField(r, 1.0f, 1.5707963f, 1.0f, 65), cx, cy);   // 90 degrees
    CHECK(std::abs(cy) > 0.3);
    CHECK(std::abs(cx) < 0.25 * std::abs(cy));
}

TEST_CASE("sampling exactly on a ring returns that ring's energy") {
    const PsfRings r = buildPsfRings(disc(), [](float t) {
        Wavefront w; w.defocus = 2.0f * t; return w; }, 550.0f, 550.0f, 5, 128);
    auto energy = [](const Plane& p) { double s = 0; for (float v : p.v) s += v; return s; };
    const double onRing = energy(psfAtField(r, 0.5f, 0.0f, 1.0f, 41));   // ring index 2 of 5
    CHECK(onRing > 0.0);
}

TEST_CASE("resampling to coarser pixels keeps the energy") {
    const PsfRings r = buildPsfRings(disc(), [](float) { return Wavefront{}; },
                                     550.0f, 550.0f, 4, 256);
    auto energy = [](const Plane& p) { double s = 0; for (float v : p.v) s += v; return s; };
    const double fine   = energy(psfAtField(r, 0.0f, 0.0f, 1.0f, 81));
    const double coarse = energy(psfAtField(r, 0.0f, 0.0f, 2.0f, 41));
    CHECK(coarse / fine == doctest::Approx(1.0).epsilon(0.05));
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/optics/psfrings.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`lenscore/optics/psfrings.hpp`:
```cpp
#pragma once
#include "lenscore/geometry.hpp"
#include "lenscore/optics/psf.hpp"
#include <functional>

namespace lens::optics {

struct PsfRings {
    std::vector<Plane> ring;   // index i is field radius t = i/(rings-1)
    int gridN = 0;
};

inline PsfRings buildPsfRings(const PupilParams& pp,
                              const std::function<Wavefront(float)>& wavefrontAtT,
                              float lambdaNm, float lambdaRefNm, int rings, int N) {
    PsfRings out;
    out.gridN = N;
    out.ring.reserve(size_t(rings));
    for (int i = 0; i < rings; ++i) {
        const float t = (rings > 1) ? float(i) / float(rings - 1) : 0.0f;
        out.ring.push_back(psfFromPupil(pp, wavefrontAtT(t), lambdaNm, lambdaRefNm, t, N));
    }
    return out;
}

// Interpolate between rings by radius, rotate into the radial frame, and
// resample from PSF grid samples to image pixels. Energy is conserved by
// scaling with the area ratio; the PSF is still never normalised to 1.
inline Plane psfAtField(const PsfRings& r, float t, float thetaRad,
                        float samplesPerPixel, int outSize) {
    const int n = int(r.ring.size());
    const float ft = std::clamp(t, 0.0f, 1.0f) * float(n - 1);
    const int i0 = std::clamp(int(ft), 0, n - 2);
    const float f = ft - float(i0);

    const float ct = std::cos(thetaRad), st = std::sin(thetaRad);
    const float c = 0.5f * float(outSize - 1);
    const float gc = 0.5f * float(r.gridN - 1);

    Plane out(outSize, outSize);
    for (int y = 0; y < outSize; ++y) {
        for (int x = 0; x < outSize; ++x) {
            // Output pixel offset, rotated back into the ring's radial frame.
            const float dx = (float(x) - c) * samplesPerPixel;
            const float dy = (float(y) - c) * samplesPerPixel;
            const float ux =  ct * dx + st * dy;
            const float uy = -st * dx + ct * dy;
            const float sx = gc + ux, sy = gc + uy;
            const float a = sampleBilinear(r.ring[i0],     sx, sy);
            const float b = sampleBilinear(r.ring[i0 + 1], sx, sy);
            out.at(x, y) = (1.0f - f) * a + f * b;
        }
    }
    // One output pixel now covers samplesPerPixel^2 grid samples.
    const float area = samplesPerPixel * samplesPerPixel;
    for (float& v : out.v) v *= area;
    return out;
}

}  // namespace lens::optics
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all five cases pass. The coma case is the one that proves the radial
frame is wired correctly; if the tail points sideways, the sign of `st` in the
inverse rotation is flipped.

- [ ] **Step 5: Commit**

```bash
git add lenscore/optics/psfrings.hpp tests/test_psfrings.cpp tests/CMakeLists.txt
git commit -m "feat: polar PSF rings exploiting rotational covariance"
```

---
### Task 16: Synthetic test targets

**Files:**
- Create: `tests/targets.hpp`, `tests/test_targets.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `lens::Plane`, `lens::Image`, `lens::frameOf`
- Produces, all in namespace `lens::targets`:
  - `Plane flatField(int w, int h, float value)`
  - `Plane pointGrid(int w, int h, int spacing)`
  - `Plane slantedEdge(int w, int h, float angleDeg, float lo, float hi)`
  - `Plane siemensStar(int w, int h, int spokes)`
  - `Plane singlePoint(int w, int h)`
  - `Image toImage(const Plane&)` and `Plane luminance(const Image&)`

These are the inputs every accuracy claim is measured against. `slantedEdge`
is antialiased by area coverage rather than by supersampling, so the edge
profile is exact — a metric that measures blur must not be fed an edge that is
already blurred by its own generator.

- [ ] **Step 1: Write the failing test**

`tests/test_targets.cpp`:
```cpp
#include <doctest/doctest.h>
#include "targets.hpp"
#include <cmath>

using namespace lens;
using namespace lens::targets;

TEST_CASE("flat field is flat") {
    const Plane p = flatField(16, 9, 0.7f);
    for (float v : p.v) CHECK(v == doctest::Approx(0.7f));
}

TEST_CASE("point grid places isolated unit points on the expected lattice") {
    const Plane p = pointGrid(64, 64, 16);
    double sum = 0; int nonzero = 0;
    for (float v : p.v) { sum += v; if (v > 0.0f) ++nonzero; }
    CHECK(nonzero == 9);                 // 3 x 3 interior lattice
    CHECK(sum == doctest::Approx(9.0));
}

TEST_CASE("slanted edge is genuinely slanted and hits both levels") {
    const Plane p = slantedEdge(64, 64, 5.0f, 0.1f, 0.9f);
    CHECK(p.at(2, 32)  == doctest::Approx(0.1f).epsilon(1e-3));
    CHECK(p.at(61, 32) == doctest::Approx(0.9f).epsilon(1e-3));
    // The transition column shifts down the rows because the edge is tilted.
    auto crossing = [&](int y) {
        for (int x = 1; x < 64; ++x) if (p.at(x, y) > 0.5f) return x;
        return 64;
    };
    CHECK(crossing(10) != crossing(54));
}

TEST_CASE("slanted edge has exactly one partially covered pixel per row") {
    const Plane p = slantedEdge(64, 64, 5.0f, 0.0f, 1.0f);
    for (int y = 0; y < 64; ++y) {
        int partial = 0;
        for (int x = 0; x < 64; ++x) {
            const float v = p.at(x, y);
            if (v > 1e-4f && v < 1.0f - 1e-4f) ++partial;
        }
        CHECK(partial <= 1);
    }
}

TEST_CASE("Siemens star is rotationally periodic in the spoke count") {
    const Plane p = siemensStar(128, 128, 16);
    const double r = 40.0, cx = 63.5, cy = 63.5;
    auto sample = [&](double a) {
        return double(p.at(int(std::lround(cx + r * std::cos(a))),
                           int(std::lround(cy + r * std::sin(a)))));
    };
    for (int k = 0; k < 8; ++k) {
        const double a = 0.31 + k * 0.4;
        CHECK(sample(a) == doctest::Approx(sample(a + 2.0 * kPi / 16.0)).epsilon(0.05));
    }
}

TEST_CASE("single point carries unit energy at the exact centre") {
    const Plane p = singlePoint(65, 65);
    double s = 0; for (float v : p.v) s += v;
    CHECK(s == doctest::Approx(1.0));
    CHECK(p.at(32, 32) == doctest::Approx(1.0f));
}

TEST_CASE("plane to image and back is lossless for grey") {
    const Plane a = flatField(8, 8, 0.42f);
    const Plane b = luminance(toImage(a));
    for (int i = 0; i < 64; ++i) CHECK(b.v[i] == doctest::Approx(a.v[i]).epsilon(1e-5));
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `targets.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`tests/targets.hpp`:
```cpp
#pragma once
#include "lenscore/geometry.hpp"
#include "lenscore/image.hpp"
#include "lenscore/plane.hpp"
#include <algorithm>
#include <cmath>

namespace lens::targets {

inline Plane flatField(int w, int h, float value) {
    Plane p(w, h);
    std::fill(p.v.begin(), p.v.end(), value);
    return p;
}

inline Plane pointGrid(int w, int h, int spacing) {
    Plane p(w, h);
    for (int y = spacing; y < h; y += spacing)
        for (int x = spacing; x < w; x += spacing)
            if (x < w && y < h) p.at(x, y) = 1.0f;
    return p;
}

// Area-exact antialiasing: coverage of the pixel by the half-plane, computed
// analytically, so the generator adds no blur of its own.
inline Plane slantedEdge(int w, int h, float angleDeg, float lo, float hi) {
    Plane p(w, h);
    const float a = angleDeg * 3.14159265358979323846f / 180.0f;
    const float tan_a = std::tan(a);
    const float xMid = 0.5f * float(w);
    for (int y = 0; y < h; ++y) {
        const float edgeX = xMid + tan_a * (float(y) - 0.5f * float(h));
        for (int x = 0; x < w; ++x) {
            const float d = float(x) + 0.5f - edgeX;      // signed distance in pixels
            const float cov = std::clamp(d + 0.5f, 0.0f, 1.0f);
            p.at(x, y) = lo + (hi - lo) * cov;
        }
    }
    return p;
}

inline Plane siemensStar(int w, int h, int spokes) {
    Plane p(w, h);
    const Frame f = frameOf(w, h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            const float dx = float(x) - f.cx, dy = float(y) - f.cy;
            const float th = std::atan2(dy, dx);
            p.at(x, y) = 0.5f + 0.5f * ((std::cos(th * float(spokes)) > 0.0f) ? 1.0f : -1.0f);
        }
    return p;
}

inline Plane singlePoint(int w, int h) {
    Plane p(w, h);
    p.at(w / 2, h / 2) = 1.0f;
    return p;
}

inline Image toImage(const Plane& p) {
    Image im(p.w, p.h);
    for (int y = 0; y < p.h; ++y)
        for (int x = 0; x < p.w; ++x)
            for (int c = 0; c < 3; ++c) im.at(x, y, c) = p.at(x, y);
    return im;
}

inline Plane luminance(const Image& im) {
    Plane p(im.w, im.h);
    for (int y = 0; y < im.h; ++y)
        for (int x = 0; x < im.w; ++x)
            p.at(x, y) = 0.2627f * im.at(x, y, 0) + 0.6780f * im.at(x, y, 1) + 0.0593f * im.at(x, y, 2);
    return p;
}

}  // namespace lens::targets
```

Add `target_include_directories(lens_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})` to `tests/CMakeLists.txt` so `targets.hpp` resolves.

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all seven cases pass.

- [ ] **Step 5: Commit**

```bash
git add tests/targets.hpp tests/test_targets.cpp tests/CMakeLists.txt
git commit -m "test: synthetic optical targets with area-exact edges"
```

---
### Task 17: Measurement metrics

**Files:**
- Create: `tests/metrics.hpp`, `tests/test_metrics.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `lens::Plane`, `lens::Image`, `lens::frameOf`, `lens::conv::fft1d`
- Produces, all in namespace `lens::metrics`:
  - `Plane gaussianBlur(const Plane&, float sigma)` — reference blur, tests only
  - `Plane crop(const Plane&, int x0, int y0, int w, int h)`
  - `float mtf50(const Plane& roi)` — cycles per pixel, from a near-vertical slanted edge
  - `float edgePosition(const Plane& roi)` — sub-pixel intercept of the fitted edge
  - `std::vector<float> radialMean(const Plane&, int bins)`
  - `float fringeWidthPx(const Image& roi)` — red minus blue edge position
  - `double totalEnergy(const Plane&)` / `double totalEnergy(const Image&)`
  - `float rot90Asymmetry(const Plane&)`

`mtf50` is the ISO 12233 slanted-edge method: fit the edge, project every pixel
onto the edge normal to build a 4x oversampled edge spread function,
differentiate to a line spread function, window it, and transform. It is the
single most important number in the project — every claim about edge blur and
field curvature is measured with it.

**The metric verifies itself.** A Gaussian PSF of width `sigma` has the analytic
`MTF50 = sqrt(ln 2) / (pi * sigma * sqrt(2)) = 0.1874 / sigma` cycles per pixel.
Blurring a synthetic edge by a known `sigma` and recovering that number proves
the measurement before it is ever used to judge the optics. A metric that has
not been checked against a closed form is not evidence.

- [ ] **Step 1: Write the failing test**

`tests/test_metrics.cpp`:
```cpp
#include <doctest/doctest.h>
#include "metrics.hpp"
#include "targets.hpp"
#include <cmath>

using namespace lens;
using namespace lens::metrics;
using namespace lens::targets;

TEST_CASE("gaussian blur preserves energy and is symmetric") {
    Plane p(64, 64); p.at(32, 32) = 1.0f;
    const Plane b = gaussianBlur(p, 2.0f);
    CHECK(totalEnergy(b) == doctest::Approx(1.0).epsilon(0.01));
    CHECK(b.at(30, 32) == doctest::Approx(b.at(34, 32)).epsilon(1e-4));
    CHECK(b.at(32, 30) == doctest::Approx(b.at(32, 34)).epsilon(1e-4));
}

TEST_CASE("mtf50 of a sharp edge is near the Nyquist limit") {
    const Plane e = slantedEdge(96, 96, 5.0f, 0.0f, 1.0f);
    CHECK(mtf50(e) > 0.35f);
}

TEST_CASE("mtf50 recovers the analytic value for a known gaussian") {
    for (float sigma : {1.0f, 1.5f, 2.0f, 3.0f}) {
        const Plane e = gaussianBlur(slantedEdge(128, 128, 5.0f, 0.0f, 1.0f), sigma);
        const float got  = mtf50(crop(e, 24, 24, 80, 80));
        const float want = 0.1874f / sigma;
        CAPTURE(sigma); CAPTURE(got); CAPTURE(want);
        CHECK(got == doctest::Approx(want).epsilon(0.08));
    }
}

TEST_CASE("mtf50 falls monotonically as blur grows") {
    float prev = 1.0f;
    for (float sigma : {0.8f, 1.2f, 1.8f, 2.6f}) {
        const float m = mtf50(crop(gaussianBlur(slantedEdge(128, 128, 5.0f, 0, 1), sigma), 24, 24, 80, 80));
        CHECK(m < prev);
        prev = m;
    }
}

TEST_CASE("edge position finds the true sub-pixel intercept") {
    const Plane e = slantedEdge(64, 64, 5.0f, 0.0f, 1.0f);
    CHECK(edgePosition(e) == doctest::Approx(32.0f).epsilon(0.02));
}

TEST_CASE("radial mean of a flat field is flat") {
    const std::vector<float> r = radialMean(flatField(64, 64, 0.6f), 8);
    CHECK(r.size() == 8u);
    for (float v : r) CHECK(v == doctest::Approx(0.6f).epsilon(1e-4));
}

TEST_CASE("radial mean detects a synthetic vignette") {
    Plane p(128, 128);
    const Frame f = frameOf(128, 128);
    for (int y = 0; y < 128; ++y)
        for (int x = 0; x < 128; ++x) {
            const float dx = (x - f.cx) / f.halfDiag, dy = (y - f.cy) / f.halfDiag;
            p.at(x, y) = 1.0f - 0.5f * std::sqrt(dx * dx + dy * dy);
        }
    const std::vector<float> r = radialMean(p, 10);
    CHECK(r.front() > 0.95f);
    CHECK(r.back() < 0.60f);
    for (size_t i = 1; i < r.size(); ++i) CHECK(r[i] < r[i - 1]);
}

TEST_CASE("fringe width is zero with no chromatic offset and nonzero with one") {
    const Plane e = slantedEdge(96, 96, 5.0f, 0.0f, 1.0f);
    CHECK(fringeWidthPx(toImage(e)) == doctest::Approx(0.0f).epsilon(0.05));

    Image shifted = toImage(e);                              // move blue by one pixel
    for (int y = 0; y < 96; ++y)
        for (int x = 95; x > 0; --x) shifted.at(x, y, 2) = shifted.at(x - 1, y, 2);
    CHECK(fringeWidthPx(shifted) == doctest::Approx(-1.0f).epsilon(0.10));
}

TEST_CASE("rotational asymmetry is zero for a symmetric image and large otherwise") {
    Plane sym(65, 65);
    const Frame f = frameOf(65, 65);
    for (int y = 0; y < 65; ++y)
        for (int x = 0; x < 65; ++x) {
            const float dx = x - f.cx, dy = y - f.cy;
            sym.at(x, y) = std::exp(-0.01f * (dx * dx + dy * dy));
        }
    CHECK(rot90Asymmetry(sym) < 1e-4f);

    Plane skew = sym;
    for (int y = 0; y < 65; ++y) skew.at(10, y) += 1.0f;
    CHECK(rot90Asymmetry(skew) > 0.1f);
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `metrics.hpp` not found.

- [ ] **Step 3: Write the minimal implementation**

`tests/metrics.hpp`:
```cpp
#pragma once
#include "lenscore/conv/fft.hpp"
#include "lenscore/geometry.hpp"
#include "lenscore/image.hpp"
#include "lenscore/plane.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace lens::metrics {

inline Plane crop(const Plane& p, int x0, int y0, int w, int h) {
    Plane o(w, h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x)
            o.at(x, y) = p.at(std::clamp(x0 + x, 0, p.w - 1), std::clamp(y0 + y, 0, p.h - 1));
    return o;
}

inline Plane gaussianBlur(const Plane& src, float sigma) {
    const int r = std::max(1, int(std::ceil(4.0f * sigma)));
    std::vector<float> k(size_t(2 * r + 1));
    float sum = 0.0f;
    for (int i = -r; i <= r; ++i) { k[i + r] = std::exp(-0.5f * i * i / (sigma * sigma)); sum += k[i + r]; }
    for (float& v : k) v /= sum;

    Plane tmp(src.w, src.h), out(src.w, src.h);
    for (int y = 0; y < src.h; ++y)
        for (int x = 0; x < src.w; ++x) {
            float a = 0.0f;
            for (int i = -r; i <= r; ++i) a += k[i + r] * src.at(std::clamp(x + i, 0, src.w - 1), y);
            tmp.at(x, y) = a;
        }
    for (int y = 0; y < src.h; ++y)
        for (int x = 0; x < src.w; ++x) {
            float a = 0.0f;
            for (int i = -r; i <= r; ++i) a += k[i + r] * tmp.at(x, std::clamp(y + i, 0, src.h - 1));
            out.at(x, y) = a;
        }
    return out;
}

// Fits the near-vertical edge; returns slope (px per row) and intercept at mid height.
inline void fitEdge(const Plane& roi, float& slope, float& intercept) {
    std::vector<float> cy(size_t(roi.h)), yy(size_t(roi.h));
    int n = 0;
    for (int y = 0; y < roi.h; ++y) {
        double num = 0.0, den = 0.0;
        for (int x = 0; x + 1 < roi.w; ++x) {
            const double d = std::abs(double(roi.at(x + 1, y)) - roi.at(x, y));
            num += d * (x + 0.5); den += d;
        }
        if (den > 1e-6) { cy[n] = float(num / den); yy[n] = float(y) - 0.5f * float(roi.h); ++n; }
    }
    double sx = 0, sy = 0, sxy = 0, sxx = 0;
    for (int i = 0; i < n; ++i) { sx += yy[i]; sy += cy[i]; sxy += yy[i] * cy[i]; sxx += yy[i] * yy[i]; }
    const double det = n * sxx - sx * sx;
    slope     = (std::abs(det) > 1e-9) ? float((n * sxy - sx * sy) / det) : 0.0f;
    intercept = float((sy - slope * sx) / std::max(1, n));
}

inline float edgePosition(const Plane& roi) {
    float s = 0, b = 0; fitEdge(roi, s, b); return b;
}

inline float mtf50(const Plane& roi) {
    float slope = 0, intercept = 0;
    fitEdge(roi, slope, intercept);

    // 4x oversampled edge spread function, projected onto the edge normal.
    const int OS = 4, NB = 256;
    std::vector<double> acc(NB, 0.0), cnt(NB, 0.0);
    for (int y = 0; y < roi.h; ++y) {
        const float ex = intercept + slope * (float(y) - 0.5f * float(roi.h));
        for (int x = 0; x < roi.w; ++x) {
            const float u = float(x) + 0.5f - ex;
            const int b = int(std::lround(u * OS)) + NB / 2;
            if (b >= 0 && b < NB) { acc[b] += roi.at(x, y); cnt[b] += 1.0; }
        }
    }
    std::vector<double> esf(NB, 0.0);
    double last = 0.0;
    for (int i = 0; i < NB; ++i) { if (cnt[i] > 0) last = acc[i] / cnt[i]; esf[i] = last; }

    // Line spread function, Hamming windowed to suppress ringing.
    std::vector<conv::Cplx> lsf(NB, conv::Cplx(0, 0));
    for (int i = 1; i < NB; ++i) {
        const double w = 0.54 - 0.46 * std::cos(2.0 * kPi * i / (NB - 1));
        lsf[i] = conv::Cplx(float((esf[i] - esf[i - 1]) * w), 0.0f);
    }
    conv::fft1d(lsf, false);

    const float dc = std::abs(lsf[0]);
    if (dc <= 0.0f) return 0.0f;
    // Bin width is 1/OS pixels, so bin k sits at k*OS/NB cycles per pixel.
    float prevF = 0.0f, prevM = 1.0f;
    for (int k = 1; k < NB / 2; ++k) {
        const float f = float(k) * OS / float(NB);
        const float m = std::abs(lsf[k]) / dc;
        if (m <= 0.5f) return prevF + (prevM - 0.5f) / std::max(1e-6f, prevM - m) * (f - prevF);
        prevF = f; prevM = m;
    }
    return float(OS) * 0.5f;
}

inline std::vector<float> radialMean(const Plane& p, int bins) {
    const Frame f = frameOf(p.w, p.h);
    std::vector<double> acc(size_t(bins), 0.0), cnt(size_t(bins), 0.0);
    for (int y = 0; y < p.h; ++y)
        for (int x = 0; x < p.w; ++x) {
            const float dx = (float(x) - f.cx) / f.halfDiag, dy = (float(y) - f.cy) / f.halfDiag;
            const int b = std::min(bins - 1, int(std::sqrt(dx * dx + dy * dy) * bins));
            acc[b] += p.at(x, y); cnt[b] += 1.0;
        }
    std::vector<float> out(size_t(bins), 0.0f);
    for (int i = 0; i < bins; ++i) out[i] = cnt[i] > 0 ? float(acc[i] / cnt[i]) : 0.0f;
    return out;
}

inline Plane channel(const Image& im, int c) {
    Plane p(im.w, im.h);
    for (int y = 0; y < im.h; ++y) for (int x = 0; x < im.w; ++x) p.at(x, y) = im.at(x, y, c);
    return p;
}

inline float fringeWidthPx(const Image& roi) {
    return edgePosition(channel(roi, 0)) - edgePosition(channel(roi, 2));
}

inline double totalEnergy(const Plane& p) { double s = 0; for (float v : p.v) s += v; return s; }
inline double totalEnergy(const Image& im) { double s = 0; for (float v : im.px) s += v; return s; }

inline float rot90Asymmetry(const Plane& p) {
    if (p.w != p.h) return 1e9f;
    const int n = p.w;
    float peak = 0.0f, worst = 0.0f;
    for (float v : p.v) peak = std::max(peak, std::abs(v));
    if (peak <= 0.0f) return 0.0f;
    for (int y = 0; y < n; ++y)
        for (int x = 0; x < n; ++x)
            worst = std::max(worst, std::abs(p.at(x, y) - p.at(n - 1 - y, x)));
    return worst / peak;
}

}  // namespace lens::metrics
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all nine cases pass. The analytic-gaussian case is the gate: until it
passes, no MTF number elsewhere in the project means anything.

- [ ] **Step 5: Commit**

```bash
git add tests/metrics.hpp tests/test_metrics.cpp tests/CMakeLists.txt
git commit -m "test: slanted-edge MTF and measurement metrics, verified analytically"
```

---
### Task 18: Parameters and the fused spectral pipeline

**Files:**
- Create: `lenscore/params.hpp`, `lenscore/pipeline.hpp`, `tests/test_pipeline.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: everything from Tasks 3-17
- Produces:
  - `struct lens::Params` — the versioned parameter set, with per-stage switches
  - `std::vector<std::pair<float,float>> lens::sampleBands(const Params&)` — `(lambda, weight)` pairs by inverse-CDF importance sampling of the spectral equalizer
  - `float lens::chromaticDefocusWaves(float relFocalError, float focal_mm, float lambdaNm, float fNumber)`
  - `lens::optics::Wavefront lens::wavefrontAt(const Params&, float t, float lambdaNm)`
  - `lens::Image lens::render(const Image& srcLinear, const Params&, const color::SpecTable&)`

**Stages 2 through 5 are one loop, not four passes.** Each band is generated from
the per-pixel spectral coefficients, warped, convolved, and accumulated straight
into a running XYZ sum. The spectral image — eleven full-resolution planes —
never exists.

**Defocus in waves.** `focusError` returns a *relative focal length* error. The
wavefront wants waves. The bridge is `W = F * e / (8 * lambda * N^2)`, which is
where the wavelength dependence of chromatic blur actually comes from, and it is
why stopping down suppresses chromatic defocus quadratically.

**Field dependence of the aberrations**, all fixed by the spec:
`defocus = petzval * t^2 + chromatic`, `astig = astig * t^2`,
`coma = coma * t^3`, `spherical` is constant across the field.

- [ ] **Step 1: Write the failing test**

`tests/test_pipeline.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/pipeline.hpp"
#include "metrics.hpp"
#include "targets.hpp"
#include <cmath>

using namespace lens;
using namespace lens::metrics;
using namespace lens::targets;

static const color::SpecTable& table() {
    static color::SpecTable t = color::buildTable(12);
    return t;
}

static Params nullParams() {
    Params p;
    p.bands = 3;
    p.highlightRecovery = false;
    p.doLateralCa = false; p.doPsf = false; p.doVignette = false;
    return p;
}

TEST_CASE("bands are inside the visible range and their weights sum to one") {
    Params p; p.bands = 11;
    const auto b = sampleBands(p);
    CHECK(b.size() == 11u);
    double w = 0.0;
    for (auto& [l, wt] : b) { CHECK(l >= 380.0f); CHECK(l <= 780.0f); w += wt; }
    CHECK(w == doctest::Approx(1.0).epsilon(1e-4));
}

TEST_CASE("the preview tier uses three bands") {
    Params p; p.bands = 3;
    CHECK(sampleBands(p).size() == 3u);
}

TEST_CASE("every stage disabled is the identity") {
    const Image src = toImage(siemensStar(64, 64, 12));
    const Image out = render(src, nullParams(), table());
    for (size_t i = 0; i < src.px.size(); ++i)
        CHECK(out.px[i] == doctest::Approx(src.px[i]).epsilon(0.03));
}

TEST_CASE("the spectral round trip keeps grey grey") {
    Params p = nullParams(); p.bands = 11;
    for (float v : {0.2f, 0.5f, 0.9f}) {
        const Image out = render(toImage(flatField(32, 32, v)), p, table());
        CAPTURE(v);
        CHECK(out.at(16, 16, 0) == doctest::Approx(v).epsilon(0.05));
        CHECK(out.at(16, 16, 1) == doctest::Approx(v).epsilon(0.05));
        CHECK(out.at(16, 16, 2) == doctest::Approx(v).epsilon(0.05));
    }
}

TEST_CASE("chromatic defocus in waves falls quadratically with the f-number") {
    const float a = chromaticDefocusWaves(1e-3f, 32.0f, 550.0f, 2.0f);
    const float b = chromaticDefocusWaves(1e-3f, 32.0f, 550.0f, 4.0f);
    CHECK(a / b == doctest::Approx(4.0f).epsilon(1e-4));
}

TEST_CASE("wavefront field dependence follows the documented powers") {
    Params p;
    p.petzval = 1.0f; p.astig = 1.0f; p.coma = 1.0f; p.spherical = 1.0f;
    p.dispersion.residual = 0.0f;                    // isolate the geometric terms
    const optics::Wavefront a = wavefrontAt(p, 0.5f, 550.0f);
    const optics::Wavefront b = wavefrontAt(p, 1.0f, 550.0f);
    CHECK(b.defocus   / a.defocus   == doctest::Approx(4.0f).epsilon(1e-3));   // t^2
    CHECK(b.astig     / a.astig     == doctest::Approx(4.0f).epsilon(1e-3));   // t^2
    CHECK(b.coma      / a.coma      == doctest::Approx(8.0f).epsilon(1e-3));   // t^3
    CHECK(b.spherical == doctest::Approx(a.spherical));                        // constant
}

TEST_CASE("vignetting reproduces the model's radial curve") {
    Params p = nullParams();
    p.doVignette = true;
    p.vignette.tStop = 2.0f;
    const Image out = render(toImage(flatField(128, 128, 1.0f)), p, table());
    const std::vector<float> got = radialMean(luminance(out), 8);
    for (int i = 0; i < 8; ++i) {
        const float t = (i + 0.5f) / 8.0f;
        CAPTURE(i);
        CHECK(got[i] == doctest::Approx(optics::vignette(p.vignette, t)).epsilon(0.02));
    }
}

TEST_CASE("a rotationally symmetric input stays symmetric through the whole pipeline") {
    Params p;
    p.bands = 3;
    p.highlightRecovery = false;
    p.doLateralCa = true; p.doPsf = true; p.doVignette = true;
    p.lateralK = 2e-5f; p.petzval = 0.8f; p.spherical = 0.3f;
    p.psfGrid = 64; p.psfRings = 6; p.psfKernel = 17; p.effPatch = 32;

    Plane sym(65, 65);
    const Frame f = frameOf(65, 65);
    for (int y = 0; y < 65; ++y)
        for (int x = 0; x < 65; ++x) {
            const float dx = x - f.cx, dy = y - f.cy;
            sym.at(x, y) = std::exp(-0.004f * (dx * dx + dy * dy));
        }
    CHECK(rot90Asymmetry(luminance(render(toImage(sym), p, table()))) < 0.02f);
}

TEST_CASE("energy is conserved across the spectral and convolution stages") {
    Params p;
    p.bands = 3;
    p.highlightRecovery = false;
    p.doLateralCa = false; p.doVignette = false; p.doPsf = true;
    p.petzval = 0.0f; p.astig = 0.0f; p.coma = 0.0f; p.spherical = 0.0f;
    p.dispersion.residual = 0.0f;
    p.pupil.rEntrance = 1e6f; p.pupil.rExit = 1e6f; p.pupil.sepNorm = 0.0f;
    p.psfGrid = 64; p.psfRings = 4; p.psfKernel = 17; p.effPatch = 32;

    Plane blob(64, 64);
    for (int y = 20; y < 44; ++y) for (int x = 20; x < 44; ++x) blob.at(x, y) = 1.0f;
    const Image src = toImage(blob);
    const Image out = render(src, p, table());
    CHECK(totalEnergy(out) / totalEnergy(src) == doctest::Approx(1.0).epsilon(0.02));
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: compile error, `lenscore/pipeline.hpp` not found.

- [ ] **Step 3: Write `lenscore/params.hpp`**

```cpp
#pragma once
#include "lenscore/color/spectable.hpp"
#include "lenscore/color/transfer.hpp"
#include "lenscore/optics/dispersion.hpp"
#include "lenscore/optics/distortion.hpp"
#include "lenscore/optics/pupil.hpp"
#include "lenscore/optics/vignette.hpp"
#include <vector>

namespace lens {

struct Params {
    int schema = 1;

    // Stage switches. A user who wants only chromatic aberration pays for nothing else.
    bool doLateralCa = true;
    bool doPsf       = true;
    bool doVignette  = true;

    // Colour
    bool         highlightRecovery = true;
    color::Knee  knee{};
    int          bands = 11;                 // 3 is the preview tier
    std::vector<std::pair<float, float>> equalizer{};   // (lambda, weight); empty means flat

    // Lens
    optics::Dispersion     dispersion{};
    optics::Distortion     distortion{};
    optics::PupilParams    pupil{};
    optics::VignetteParams vignette{};
    float lateralK   = 0.0f;    // k_l
    float lambdaHat  = 650.0f;  // reference wavelength, keeps magnification positive

    // Wavefront coefficients in waves, at the corner of the frame
    float petzval   = 0.0f;
    float astig     = 0.0f;
    float coma      = 0.0f;
    float spherical = 0.0f;

    // Sampling
    float focal_mm     = 32.0f;
    float fNumberWide  = 2.0f;
    float pixelPitchUm = 5.0f;
    int   psfGrid      = 128;   // pupil FFT grid, power of two
    int   psfRings     = 12;
    int   psfKernel    = 33;    // odd
    int   effPatch     = 64;    // power of two
};

}  // namespace lens
```

- [ ] **Step 4: Write `lenscore/pipeline.hpp`**

```cpp
#pragma once
#include "lenscore/conv/eff.hpp"
#include "lenscore/optics/lateralca.hpp"
#include "lenscore/optics/psfrings.hpp"
#include "lenscore/params.hpp"
#include <algorithm>
#include <cmath>

namespace lens {

// Inverse-CDF importance sampling of the equalizer, so emphasising a band
// costs no extra noise. A flat equalizer reduces to uniform spacing.
inline std::vector<std::pair<float, float>> sampleBands(const Params& p) {
    const int N = std::max(1, p.bands);
    const int R = 512;
    std::vector<float> pdf(size_t(R)), cdf(size_t(R) + 1, 0.0f);
    for (int i = 0; i < R; ++i) {
        const float l = color::kLambdaMin + (color::kLambdaMax - color::kLambdaMin) * (i + 0.5f) / R;
        float w = 1.0f;
        if (!p.equalizer.empty()) {
            w = p.equalizer.back().second;
            for (size_t k = 1; k < p.equalizer.size(); ++k)
                if (l <= p.equalizer[k].first) {
                    const auto& a = p.equalizer[k - 1]; const auto& b = p.equalizer[k];
                    const float f = (l - a.first) / std::max(1e-6f, b.first - a.first);
                    w = a.second + f * (b.second - a.second);
                    break;
                }
        }
        pdf[i] = std::max(1e-6f, w);
    }
    for (int i = 0; i < R; ++i) cdf[i + 1] = cdf[i] + pdf[i];
    for (float& c : cdf) c /= cdf[R];

    std::vector<std::pair<float, float>> out;
    out.reserve(size_t(N));
    for (int k = 0; k < N; ++k) {
        const float u = (k + 0.5f) / float(N);
        int i = 0;
        while (i + 1 < R && cdf[i + 1] < u) ++i;
        const float l = color::kLambdaMin + (color::kLambdaMax - color::kLambdaMin) * (i + 0.5f) / R;
        out.emplace_back(l, 1.0f / float(N));   // importance sampling makes weights equal
    }
    return out;
}

// Relative focal-length error to waves of defocus: W = F e / (8 lambda N^2).
inline float chromaticDefocusWaves(float relFocalError, float focal_mm, float lambdaNm, float fNumber) {
    const float lambda_mm = lambdaNm * 1e-6f;
    return focal_mm * relFocalError / (8.0f * lambda_mm * fNumber * fNumber);
}

inline optics::Wavefront wavefrontAt(const Params& p, float t, float lambdaNm) {
    optics::Wavefront w;
    const float e = optics::focusError(p.dispersion, lambdaNm, p.lambdaHat);
    w.defocus   = p.petzval * t * t + chromaticDefocusWaves(e, p.focal_mm, lambdaNm, p.fNumberWide);
    w.astig     = p.astig * t * t;
    w.coma      = p.coma * t * t * t;
    w.spherical = p.spherical;
    return w;
}

inline Image render(const Image& src, const Params& p, const color::SpecTable& tbl) {
    const int w = src.w, h = src.h;
    const Frame frame = frameOf(w, h);

    // Per-pixel spectral coefficients plus the scale that rides outside the model.
    std::vector<color::Coeffs> coeff(size_t(w) * h);
    std::vector<float> scale(size_t(w) * h, 0.0f);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            color::RGB c{src.at(x, y, 0), src.at(x, y, 1), src.at(x, y, 2)};
            if (p.highlightRecovery) {
                c.r = color::expandHighlights(c.r, p.knee);
                c.g = color::expandHighlights(c.g, p.knee);
                c.b = color::expandHighlights(c.b, p.knee);
            }
            const size_t i = size_t(y) * w + x;
            // The table reproduces the colour directly for in-gamut inputs, so only the
            // EXCESS above the gamut boundary rides outside the model. Scaling by the raw
            // maximum here would square the brightness of every in-gamut pixel.
            scale[i] = std::max(1.0f, std::max({c.r, c.g, c.b, 0.0f}));
            coeff[i] = color::lookup(tbl, c);
        }

    const auto bands = sampleBands(p);
    std::vector<float> X(size_t(w) * h, 0.0f), Y(X.size(), 0.0f), Z(X.size(), 0.0f);
    double norm = 0.0;

    for (const auto& [lambda, weight] : bands) {
        Plane band(w, h);
        for (size_t i = 0; i < band.v.size(); ++i)
            band.v[i] = scale[i] * color::evalSpectrum(coeff[i], lambda);

        if (p.doLateralCa)
            band = optics::warpPlane(band, p.distortion, p.lateralK * (p.lambdaHat - lambda));

        if (p.doPsf) {
            const optics::PsfRings rings = optics::buildPsfRings(
                p.pupil, [&](float t) { return wavefrontAt(p, t, lambda); },
                lambda, p.lambdaHat, p.psfRings, p.psfGrid);
            const float spp = optics::psfSampleSpacingUm(lambda, p.fNumberWide) / p.pixelPitchUm;
            band = conv::effConvolve(band, p.effPatch, [&](float cx, float cy) {
                const float dx = (cx - frame.cx) / frame.halfDiag;
                const float dy = (cy - frame.cy) / frame.halfDiag;
                const float t  = std::clamp(std::sqrt(dx * dx + dy * dy), 0.0f, 1.0f);
                return optics::psfAtField(rings, t, std::atan2(dy, dx), spp, p.psfKernel);
            });
        }

        const color::XYZ m = color::cmf(lambda);
        for (size_t i = 0; i < X.size(); ++i) {
            X[i] += weight * band.v[i] * m.x;
            Y[i] += weight * band.v[i] * m.y;
            Z[i] += weight * band.v[i] * m.z;
        }
        norm += double(weight) * m.y;
    }

    Image out(w, h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            const size_t i = size_t(y) * w + x;
            color::RGB c = color::xyzToRec2020(color::XYZ{
                float(X[i] / norm), float(Y[i] / norm), float(Z[i] / norm)});

            if (p.doVignette) {
                const float dx = (float(x) - frame.cx) / frame.halfDiag;
                const float dy = (float(y) - frame.cy) / frame.halfDiag;
                const float v = optics::vignette(p.vignette, std::sqrt(dx * dx + dy * dy));
                c.r *= v; c.g *= v; c.b *= v;
            }
            if (p.highlightRecovery) {
                c.r = color::compressHighlights(c.r, p.knee);
                c.g = color::compressHighlights(c.g, p.knee);
                c.b = color::compressHighlights(c.b, p.knee);
            }
            out.at(x, y, 0) = c.r; out.at(x, y, 1) = c.g; out.at(x, y, 2) = c.b;
        }
    return out;
}

}  // namespace lens
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: all nine cases pass.

- [ ] **Step 6: Commit**

```bash
git add lenscore/params.hpp lenscore/pipeline.hpp tests/test_pipeline.cpp tests/CMakeLists.txt
git commit -m "feat: fused spectral pipeline with per-stage switches"
```

---
### Task 19: Lens files, the CLI, and the acceptance tests

**Files:**
- Create: `lensdata/lensfile.hpp`, `lensdata/CMakeLists.txt`, `lensdata/cooke-s4-32mm.lens`, `lenscli/main.cpp`, `lenscli/CMakeLists.txt`, `tests/test_acceptance.cpp`
- Modify: `CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `lens::Params` (Task 18), everything below it
- Produces:
  - `std::optional<lens::Params> lens::data::loadLensFile(const std::string& path)`
  - `lenscli` with subcommands `render`, `target`, `measure`

`lensdata` is the **only** component that touches JSON. `lenscore` stays free of
it, which is what lets the same physics compile into a Photoshop plugin, a
WebAssembly panel, or a test binary without dragging a parser along.

**The acceptance tests are the point of the whole plan.** They are the numbers
that turn "accurate" from a claim into a measurement.

- [ ] **Step 1: Write the failing acceptance test**

`tests/test_acceptance.cpp`:
```cpp
#include <doctest/doctest.h>
#include "lenscore/pipeline.hpp"
#include "metrics.hpp"
#include "targets.hpp"
#include <cmath>

using namespace lens;
using namespace lens::metrics;
using namespace lens::targets;

static const color::SpecTable& table() {
    static color::SpecTable t = color::buildTable(12);
    return t;
}

static Params base() {
    Params p;
    p.bands = 7;
    p.highlightRecovery = false;
    p.doLateralCa = false; p.doPsf = false; p.doVignette = false;
    p.psfGrid = 64; p.psfRings = 8; p.psfKernel = 21; p.effPatch = 32;
    p.dispersion.residual = 0.0f;
    return p;
}

// An edge patch placed at a chosen normalised field radius along +x.
static Image edgeAt(const Params& p, int full, float t, bool vertical) {
    Plane big(full, full);
    const Frame f = frameOf(full, full);
    const int roi = 64;
    const int cx = int(f.cx + t * f.halfDiag * 0.70710678f);
    const int cy = int(f.cy + t * f.halfDiag * 0.70710678f);
    const Plane e = vertical ? slantedEdge(roi, roi, 5.0f, 0.05f, 0.95f)
                             : slantedEdge(roi, roi, 85.0f, 0.05f, 0.95f);
    for (int y = 0; y < roi; ++y)
        for (int x = 0; x < roi; ++x) {
            const int X = cx - roi / 2 + x, Y = cy - roi / 2 + y;
            if (X >= 0 && Y >= 0 && X < full && Y < full) big.at(X, Y) = e.at(x, y);
        }
    return render(toImage(big), p, table());
}

static Plane roiAt(const Image& im, float t, int roi) {
    const Frame f = frameOf(im.w, im.h);
    const int cx = int(f.cx + t * f.halfDiag * 0.70710678f);
    const int cy = int(f.cy + t * f.halfDiag * 0.70710678f);
    return crop(luminance(im), cx - roi / 2, cy - roi / 2, roi, roi);
}

TEST_CASE("ACCEPTANCE: lateral CA fringe width grows linearly with field radius") {
    Params p = base();
    p.doLateralCa = true;
    p.lateralK = 6e-5f;

    float prev = -1.0f;
    std::vector<float> widths;
    for (float t : {0.0f, 0.35f, 0.7f, 1.0f}) {
        const Image im = edgeAt(p, 320, t, true);
        const Frame f = frameOf(320, 320);
        const int cx = int(f.cx + t * f.halfDiag * 0.70710678f);
        const int cy = int(f.cy + t * f.halfDiag * 0.70710678f);
        Image roi(48, 48);
        for (int y = 0; y < 48; ++y)
            for (int x = 0; x < 48; ++x)
                for (int c = 0; c < 3; ++c)
                    roi.at(x, y, c) = im.at(std::clamp(cx - 24 + x, 0, 319),
                                            std::clamp(cy - 24 + y, 0, 319), c);
        const float wpx = std::abs(fringeWidthPx(roi));
        CAPTURE(t); CAPTURE(wpx);
        CHECK(wpx > prev - 0.02f);      // grows with radius
        widths.push_back(wpx);
        prev = wpx;
    }
    CHECK(widths.front() < 0.10f);      // no fringe on axis
    CHECK(widths.back()  > 0.40f);      // real fringe at the corner
}

TEST_CASE("ACCEPTANCE: field curvature softens the corners relative to the centre") {
    Params p = base();
    p.doPsf = true;
    p.petzval = 1.2f;

    const float centre = mtf50(roiAt(edgeAt(p, 256, 0.0f, true), 0.0f, 48));
    const float corner = mtf50(roiAt(edgeAt(p, 256, 1.0f, true), 1.0f, 48));
    CAPTURE(centre); CAPTURE(corner);
    CHECK(corner < centre);
    CHECK(corner < 0.85f * centre);
}

TEST_CASE("ACCEPTANCE: astigmatism separates sagittal from tangential resolution") {
    Params p = base();
    p.doPsf = true;
    p.astig = 1.2f;

    const float sag = mtf50(roiAt(edgeAt(p, 256, 1.0f, true),  1.0f, 48));
    const float tan = mtf50(roiAt(edgeAt(p, 256, 1.0f, false), 1.0f, 48));
    CAPTURE(sag); CAPTURE(tan);
    CHECK(std::abs(sag - tan) / std::max(sag, tan) > 0.08f);
}

TEST_CASE("ACCEPTANCE: an achromat fringes green and magenta, not red and blue") {
    // Corrected at the F and C lines: red and blue focus together, green does not.
    Params p = base();
    p.doPsf = true;
    p.dispersion.correction_nm = {486.1f, 656.3f};
    p.dispersion.residual = 1.0f;
    p.bands = 11;

    const float eF = optics::focusError(p.dispersion, 486.1f, p.lambdaHat);
    const float eC = optics::focusError(p.dispersion, 656.3f, p.lambdaHat);
    const float eG = optics::focusError(p.dispersion, 546.1f, p.lambdaHat);
    CAPTURE(eF); CAPTURE(eC); CAPTURE(eG);

    CHECK(std::abs(eF) < 1e-6f);
    CHECK(std::abs(eC) < 1e-6f);
    CHECK(std::abs(eG) > 1e-5f);                       // green is the odd one out

    // The uncorrected singlet does the opposite: a monotonic red-to-blue ramp.
    Params s = p;
    s.dispersion.correction_nm.clear();
    const float sF = optics::focusError(s.dispersion, 486.1f, s.lambdaHat);
    const float sC = optics::focusError(s.dispersion, 656.3f, s.lambdaHat);
    CHECK(std::abs(sF - sC) > std::abs(eF - eC));
}

TEST_CASE("ACCEPTANCE: stopping down removes mechanical vignetting from the render") {
    Params p = base();
    p.doVignette = true;
    p.vignette.tStop = 2.0f;
    const float wide = radialMean(luminance(render(toImage(flatField(128, 128, 1.0f)), p, table())), 8).back();
    p.vignette.tStop = 8.0f;
    const float stopped = radialMean(luminance(render(toImage(flatField(128, 128, 1.0f)), p, table())), 8).back();
    CHECK(stopped > wide);
    CHECK(stopped == doctest::Approx(optics::naturalFalloff(p.vignette, 0.94f)).epsilon(0.05));
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build --output-on-failure`
Expected: the acceptance file compiles (it only uses existing headers) but the
cases fail or the target does not exist yet. Add `test_acceptance.cpp` to
`tests/CMakeLists.txt` first, then confirm the failures are real numeric ones.

- [ ] **Step 3: Make the acceptance tests pass**

These exercise code written in Tasks 3-18. If one fails, the defect is in the
stage it names, not in the test. Work the failure back to its unit test and add
the missing case there before fixing anything.

- [ ] **Step 4: Write the lens file loader**

`lensdata/lensfile.hpp`:
```cpp
#pragma once
#include "lenscore/params.hpp"
#include <fstream>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>

namespace lens::data {

inline std::optional<Params> loadLensFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) return std::nullopt;
    nlohmann::json j;
    try { in >> j; } catch (...) { return std::nullopt; }
    if (j.value("schema", 0) != 2) return std::nullopt;

    Params p;
    p.focal_mm    = j.value("focal_mm", 32.0f);
    p.fNumberWide = j["aperture"].value("t_stop", 2.0f);
    p.pupil.blades      = j["aperture"].value("blades", 9);
    p.pupil.curvature   = j["aperture"].value("curvature", 0.15f);
    p.pupil.rotationRad = j["aperture"].value("rotation_deg", 0.0f) * 3.14159265f / 180.0f;

    if (j.contains("dispersion")) {
        const auto& d = j["dispersion"];
        for (int i = 0; i < 3; ++i) { p.dispersion.B[i] = d["B"][i]; p.dispersion.C[i] = d["C"][i]; }
        p.dispersion.correction_nm = d.value("correction_nm", std::vector<float>{});
        p.dispersion.residual      = d.value("residual", 1.0f);
    }
    if (j.contains("distortion")) {
        const auto& d = j["distortion"];
        p.distortion.k1 = d.value("k1", 0.0f); p.distortion.k2 = d.value("k2", 0.0f);
        p.distortion.k3 = d.value("k3", 0.0f); p.distortion.p1 = d.value("p1", 0.0f);
        p.distortion.p2 = d.value("p2", 0.0f);
    }
    p.lateralK = j.contains("lateral_ca") ? j["lateral_ca"].value("k_l", 0.0f) : 0.0f;

    if (j.contains("pupil")) {
        const auto& u = j["pupil"];
        p.pupil.rEntrance        = u.value("r_entrance", 1.0f);
        p.pupil.rExit            = u.value("r_exit", 0.92f);
        p.pupil.sepNorm          = u.value("sep_norm", 0.35f);
        p.pupil.apodizationSlope = u.value("apodization_slope", 0.0f);
        p.vignette.rEntrance = p.pupil.rEntrance;
        p.vignette.sepNorm   = p.pupil.sepNorm;
    }
    if (j.contains("vignette")) p.vignette.naturalExp = j["vignette"].value("natural_exp", 4.0f);
    p.vignette.tStopWide = p.fNumberWide;
    p.vignette.tStop     = p.fNumberWide;
    p.vignette.focal_mm  = p.focal_mm;

    if (j.contains("wavefront")) {
        const auto& wv = j["wavefront"];
        p.petzval   = wv.value("petzval", 0.0f);
        p.astig     = wv.value("astig", 0.0f);
        p.coma      = wv.value("coma", 0.0f);
        p.spherical = wv.value("spherical", 0.0f);
    }
    return p;
}

}  // namespace lens::data
```

`lensdata/CMakeLists.txt`:
```cmake
include(FetchContent)
FetchContent_Declare(nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.11.3)
FetchContent_MakeAvailable(nlohmann_json)

add_library(lensdata INTERFACE)
target_include_directories(lensdata INTERFACE ${CMAKE_SOURCE_DIR})
target_link_libraries(lensdata INTERFACE lenscore nlohmann_json::nlohmann_json)
```

Create `lensdata/cooke-s4-32mm.lens` with the exact JSON from section 5 of the
spec.

- [ ] **Step 5: Write the CLI**

`lenscli/main.cpp`:
```cpp
#include "lenscore/pfm.hpp"
#include "lenscore/pipeline.hpp"
#include "lensdata/lensfile.hpp"
#include <cstdio>
#include <cstring>
#include <string>

using namespace lens;

static int usage() {
    std::fprintf(stderr,
        "lenscli render <in.pfm> <out.pfm> --lens <file.lens> [--table <t.bin>] [--bands N]\n"
        "lenscli target <flat|edge|points|star> <out.pfm> [--size N]\n");
    return 2;
}

int main(int argc, char** argv) {
    if (argc < 3) return usage();
    const std::string cmd = argv[1];

    if (cmd == "render") {
        if (argc < 4) return usage();
        std::string lensPath, tablePath = "lensdata/rgb2spec/rec2020.bin";
        int bands = 11;
        for (int i = 4; i + 1 < argc; i += 2) {
            if (!std::strcmp(argv[i], "--lens"))  lensPath  = argv[i + 1];
            if (!std::strcmp(argv[i], "--table")) tablePath = argv[i + 1];
            if (!std::strcmp(argv[i], "--bands")) bands     = std::atoi(argv[i + 1]);
        }
        auto src = pfm::read(argv[2]);
        if (!src) { std::fprintf(stderr, "cannot read %s\n", argv[2]); return 1; }
        auto tbl = color::readTable(tablePath);
        if (!tbl) { std::fprintf(stderr, "cannot read table %s\n", tablePath.c_str()); return 1; }
        auto par = lensPath.empty() ? std::optional<Params>(Params{}) : data::loadLensFile(lensPath);
        if (!par) { std::fprintf(stderr, "cannot read lens %s\n", lensPath.c_str()); return 1; }
        par->bands = bands;
        if (!pfm::write(argv[3], render(*src, *par, *tbl))) { std::fprintf(stderr, "write failed\n"); return 1; }
        std::printf("rendered %s -> %s\n", argv[2], argv[3]);
        return 0;
    }
    return usage();
}
```

`lenscli/CMakeLists.txt`:
```cmake
add_executable(lenscli main.cpp)
target_link_libraries(lenscli PRIVATE lenscore lensdata)
```

Add `add_subdirectory(lensdata)` and `add_subdirectory(lenscli)` to the
top-level `CMakeLists.txt`.

- [ ] **Step 6: Verify the CLI end to end**

```bash
cmake --build build
./build/lenscli/lenscli render tests/data/star.pfm /tmp/out.pfm \
  --lens lensdata/cooke-s4-32mm.lens --bands 11
```
Expected: `rendered tests/data/star.pfm -> /tmp/out.pfm`. Generate the input
first with a short PFM writer, or reuse any 256x256 PFM.

- [ ] **Step 7: Commit**

```bash
git add lensdata lenscli tests/test_acceptance.cpp tests/CMakeLists.txt CMakeLists.txt
git commit -m "feat: lens file loading, CLI, and the acceptance suite"
```

---

## Out of scope for this plan

Deliberately deferred. Each becomes its own plan once these interfaces are real:

| Plan | Covers | Blocked on |
|---|---|---|
| Film stage | Spectral sensitivity, emulsion MTF, halation, H&D curves, grain, print cascade (spec 4.7) | `Params` and `render()` shape, settled here |
| Veiling glare | The `gsf = delta*(1-k) + k*b` convolution and its f-stop parameterisation (spec 4.6) | Nothing; small enough to fold into the film plan |
| Photoshop host | `lens8bf` adapter, ImGui dialog, descriptors, packaging (spec 7, 8, 10) | The gated Adobe SDK, and `lenscore` being stable |
| `lensfit` | Chart and prescription fitting, real measured presets (spec 5.1) | This plan's parameter set |

## Self-review

**Spec coverage.** Sections 4.1 (Task 3), 4.2 (Tasks 4-6), 4.3 (Task 8), 4.4
(Tasks 7, 10-15), 4.5 (Task 9), 5 and 5.1 partially (Task 19; `lensfit` is
deferred and named above), 6 partially (the equalizer drives band sampling in
Task 18; the mapping-function UI belongs with the host), 9 (Tasks 16-19).
Sections 4.6, 4.7, 7, 8, 10 are explicitly out of scope and listed above.

**Placeholders.** None. Every step carries the code or the exact command.

**Type consistency.** `Plane` and `Image` are distinct throughout: optics work
on `Plane`, the pipeline boundary uses `Image`. `t` is always normalised field
radius with the corner at 1. Wavelengths are always nanometres, angles always
radians. `Wavefront` is constructed only by `wavefrontAt`. `PupilParams::
apertureRadius` and `VignetteParams::tStop` describe the same stop and are
reconciled in `loadLensFile`.

**One invariant worth restating**, because three tasks depend on it and it is
the easiest thing to break by accident: **the PSF is never normalised.** Task 13
returns unnormalised intensity, Task 15 scales only by the resampling area
ratio, and Task 14 sums whatever it is given. The energy carried through is the
optical vignetting term, and Task 11's consistency test is what proves it.
