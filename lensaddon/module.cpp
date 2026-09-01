// lensaddon -- the native half of the UXP Hybrid plugin.
//
// It does one job: take a pixel buffer and a parameter object from the panel,
// run it through lenscore, and hand the result back. It holds no UI and no
// Photoshop state.
//
// Every entry point catches everything. An exception escaping into Photoshop's
// process would take the application down with it -- see the error policy in
// docs/superpowers/specs.

#include <exception>
#include <memory>
#include <mutex>
#include <optional>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

#include "UxpAddon.h"
#include "UxpTask.h"
#include "UxpValue.h"

#include "lenscore/pipeline.hpp"
#include "lenscore/color/spectable.hpp"

// The 9.4MB spectral table is loaded once and reused. Guarded because UXP may
// call in from more than one thread over the plugin's life.
std::mutex gTableMutex;
std::optional<lens::color::SpecTable> gTable;

const lens::color::SpecTable& tableFrom(const std::string& path) {
    std::lock_guard<std::mutex> lock(gTableMutex);
    if (!gTable) {
        auto t = lens::color::readTable(path);
        if (!t) throw std::runtime_error("could not read spectral table at " + path);
        gTable = std::move(*t);
    }
    return *gTable;
}

double numberProp(addon_env env, addon_value obj, const char* name, double fallback) {
    addon_value v = nullptr;
    if (UxpAddonApis.uxp_addon_get_named_property(env, obj, name, &v) != addon_ok || !v) return fallback;
    double out = fallback;
    if (UxpAddonApis.uxp_addon_get_value_double(env, v, &out) != addon_ok) return fallback;
    return out;
}

std::string stringProp(addon_env env, addon_value obj, const char* name) {
    addon_value v = nullptr;
    if (UxpAddonApis.uxp_addon_get_named_property(env, obj, name, &v) != addon_ok || !v) return {};
    size_t len = 0;
    if (UxpAddonApis.uxp_addon_get_value_string_utf8(env, v, nullptr, 0, &len) != addon_ok) return {};
    std::string s(len + 1, '\0');
    size_t written = 0;
    if (UxpAddonApis.uxp_addon_get_value_string_utf8(env, v, s.data(), len + 1, &written) != addon_ok) return {};
    s.resize(written);
    return s;
}

// Maps the panel's controls onto lens::Params. The panel speaks in plain
// amounts; the physical parameters live here so the UI never has to know what
// a Petzval coefficient is.
lens::Params paramsFrom(addon_env env, addon_value obj) {
    lens::Params p;

    const double ca       = numberProp(env, obj, "chromaticAberration", 0.0);  // 0..100
    const double vig      = numberProp(env, obj, "vignette", 0.0);             // 0..100
    const double edge     = numberProp(env, obj, "edgeBlur", 0.0);             // 0..100
    const double astig    = numberProp(env, obj, "astigmatism", 0.0);          // 0..100
    const double coma     = numberProp(env, obj, "coma", 0.0);                 // 0..100
    const double tStop    = numberProp(env, obj, "tStop", 2.0);
    const int    bands    = int(numberProp(env, obj, "bands", 3));

    p.bands = bands < 3 ? 3 : (bands > 15 ? 15 : bands);
    p.highlightRecovery = numberProp(env, obj, "highlightRecovery", 1.0) > 0.5;

    // Amounts scale onto measured-magnitude coefficients. 100% is roughly what a
    // real lens of this class does; the panel can exceed it deliberately.
    p.lateralK  = float(ca   * 1.0e-6);
    p.petzval   = float(edge  * 0.02);
    p.astig     = float(astig * 0.02);
    p.coma      = float(coma  * 0.02);

    p.vignette.tStop     = float(tStop);
    p.vignette.tStopWide = float(numberProp(env, obj, "tStopWide", 2.0));
    p.vignette.naturalExp = float(4.0 * (vig / 100.0));

    p.doLateralCa = p.lateralK != 0.0f;
    p.doPsf       = (p.petzval != 0.0f) || (p.astig != 0.0f) || (p.coma != 0.0f);
    p.doVignette  = vig > 0.0;
    return p;
}

// renderPixels(Float32Array rgb, width, height, params, tablePath) -> Float32Array
addon_value RenderPixels(addon_env env, addon_callback_info info) {
    try {
        size_t argc = 5;
        addon_value argv[5] = {};
        Check(UxpAddonApis.uxp_addon_get_cb_info(env, info, &argc, argv, nullptr, nullptr));
        if (argc < 5) throw std::runtime_error("renderPixels expects 5 arguments");

        // Incoming pixels. This SDK exposes get_arraybuffer_info but no
        // get_typedarray_info, so the panel hands us the underlying ArrayBuffer
        // and we treat it as tightly packed float32 RGB.
        void* data = nullptr;
        size_t byteLength = 0;
        Check(UxpAddonApis.uxp_addon_get_arraybuffer_info(env, argv[0], &data, &byteLength));
        if (!data || byteLength % sizeof(float) != 0)
            throw std::runtime_error("pixels must be a float32 ArrayBuffer");
        const size_t length = byteLength / sizeof(float);

        int32_t w = 0, h = 0;
        Check(UxpAddonApis.uxp_addon_get_value_int32(env, argv[1], &w));
        Check(UxpAddonApis.uxp_addon_get_value_int32(env, argv[2], &h));
        if (w <= 0 || h <= 0) throw std::runtime_error("bad image dimensions");
        if (length < size_t(w) * size_t(h) * 3) throw std::runtime_error("pixel buffer too small for dimensions");

        lens::Params params = paramsFrom(env, argv[3]);
        const std::string tablePath = stringProp(env, argv[4], "path");

        lens::Image src(w, h);
        const float* in = static_cast<const float*>(data);
        std::copy(in, in + src.px.size(), src.px.begin());

        const lens::Image out = lens::render(src, params, tableFrom(tablePath));

        void* outData = nullptr;
        addon_value outBuffer = nullptr;
        Check(UxpAddonApis.uxp_addon_create_arraybuffer(env, out.px.size() * sizeof(float), &outData, &outBuffer));
        std::copy(out.px.begin(), out.px.end(), static_cast<float*>(outData));

        addon_value result = nullptr;
        Check(UxpAddonApis.uxp_addon_create_typedarray(
            env, addon_float32_array, out.px.size(), outBuffer, 0, &result));
        return result;
    } catch (...) {
        return CreateErrorFromException(env);
    }
}

// version() -> string. Lets the panel prove the native half actually loaded.
addon_value Version(addon_env env, addon_callback_info info) {
    try {
        const char* v = "lensaddon 0.2.0 / lenscore";
        addon_value s = nullptr;
        Check(UxpAddonApis.uxp_addon_create_string_utf8(env, v, strlen(v), &s));
        return s;
    } catch (...) {
        return CreateErrorFromException(env);
    }
}

addon_value Init(addon_env env, addon_value exports, const addon_apis& api) {
    addon_value fn = nullptr;
    if (api.uxp_addon_create_function(env, NULL, 0, RenderPixels, NULL, &fn) != addon_ok ||
        api.uxp_addon_set_named_property(env, exports, "renderPixels", fn) != addon_ok) {
        api.uxp_addon_throw_error(env, NULL, "could not export renderPixels");
    }
    if (api.uxp_addon_create_function(env, NULL, 0, Version, NULL, &fn) != addon_ok ||
        api.uxp_addon_set_named_property(env, exports, "version", fn) != addon_ok) {
        api.uxp_addon_throw_error(env, NULL, "could not export version");
    }
    return exports;
}

UXP_ADDON_INIT(Init)

// UXP calls this on unload. Registering it matters: without the macro the
// _uxp_addon_terminate symbol is never exported and the addon is only half
// wired up. Drop the cached spectral table so a reload starts clean.
void Terminate(addon_env) {
    try {
        std::lock_guard<std::mutex> lock(gTableMutex);
        gTable.reset();
    } catch (...) {
    }
}

UXP_ADDON_TERMINATE(Terminate)
