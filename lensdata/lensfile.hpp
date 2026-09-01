#pragma once
#include "lenscore/constants.hpp"
#include "lenscore/params.hpp"
#include <fstream>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>
#include <vector>

namespace lens::data {

namespace detail {
// j[key] the safe way. Plain `j[key]` is the NON-const operator[]: on a
// missing key it silently INSERTS a null node, and .value() called on that
// null then throws nlohmann::json::type_error -- straight out of a function
// whose whole contract is std::optional, never an exception. A key that is
// present but holds something other than an object (a non-object
// "dispersion" or "distortion", say) throws the same way once .value() is
// called on it. This covers both: missing key and wrong-typed key both fall
// back to an empty object, so every .value() call below sees an object and
// applies its own default instead of throwing.
inline nlohmann::json objectOrEmpty(const nlohmann::json& j, const char* key) {
    return (j.contains(key) && j.at(key).is_object()) ? j.at(key) : nlohmann::json::object();
}
}  // namespace detail

inline std::optional<Params> loadLensFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) return std::nullopt;
    nlohmann::json j;
    try { in >> j; } catch (...) { return std::nullopt; }

    // Everything from here on reads a parsed-but-untrusted document: a
    // malformed field the guards above don't anticipate (dispersion.B not a
    // 3-element numeric array, say) can still throw out of .value() or
    // operator[]. The function's contract is std::optional on any bad
    // input, never an exception, so the rest of the parse is one try block
    // rather than guarding every individual field access.
    try {
        if (!j.is_object() || j.value("schema", 0) != 2) return std::nullopt;

        // "aperture" is the one section this format treats as required (it
        // is where fNumberWide -- the working stop everything else in
        // Params gets reconciled against -- comes from), unlike dispersion/
        // distortion/pupil/vignette/wavefront below, which are genuinely
        // optional and fall back to physical defaults when absent. Missing
        // or wrong-typed, that is a malformed file, not a lens with no
        // aperture: return nullopt rather than silently defaulting.
        if (!j.contains("aperture") || !j.at("aperture").is_object()) return std::nullopt;
        const auto& aperture = j.at("aperture");

        Params p;
        p.focal_mm = j.value("focal_mm", 32.0f);
        p.fNumberWide       = aperture.value("t_stop", 2.0f);
        p.pupil.blades      = aperture.value("blades", 9);
        p.pupil.curvature   = aperture.value("curvature", 0.15f);
        p.pupil.rotationRad = aperture.value("rotation_deg", 0.0f) * kPi / 180.0f;

        if (j.contains("dispersion") && j.at("dispersion").is_object()) {
            const auto& d = j.at("dispersion");
            for (int i = 0; i < 3; ++i) {
                p.dispersion.B[i] = d.at("B").at(i);
                p.dispersion.C[i] = d.at("C").at(i);
            }
            p.dispersion.correction_nm = d.value("correction_nm", std::vector<float>{});
            p.dispersion.residual      = d.value("residual", 1.0f);
        }
        if (j.contains("distortion") && j.at("distortion").is_object()) {
            const auto& d = j.at("distortion");
            p.distortion.k1 = d.value("k1", 0.0f); p.distortion.k2 = d.value("k2", 0.0f);
            p.distortion.k3 = d.value("k3", 0.0f); p.distortion.p1 = d.value("p1", 0.0f);
            p.distortion.p2 = d.value("p2", 0.0f);
        }
        p.lateralK = detail::objectOrEmpty(j, "lateral_ca").value("k_l", 0.0f);

        if (j.contains("pupil") && j.at("pupil").is_object()) {
            const auto& u = j.at("pupil");
            p.pupil.rEntrance        = u.value("r_entrance", 1.0f);
            p.pupil.rExit            = u.value("r_exit", 0.92f);
            p.pupil.sepNorm          = u.value("sep_norm", 0.35f);
            p.pupil.apodizationSlope = u.value("apodization_slope", 0.0f);
            p.vignette.rEntrance = p.pupil.rEntrance;
            p.vignette.sepNorm   = p.pupil.sepNorm;
        }
        p.vignette.naturalExp = detail::objectOrEmpty(j, "vignette").value("natural_exp", 4.0f);
        p.vignette.tStopWide = p.fNumberWide;
        p.vignette.tStop     = p.fNumberWide;
        p.vignette.focal_mm  = p.focal_mm;

        if (j.contains("wavefront") && j.at("wavefront").is_object()) {
            const auto& wv = j.at("wavefront");
            p.petzval   = wv.value("petzval", 0.0f);
            p.astig     = wv.value("astig", 0.0f);
            p.coma      = wv.value("coma", 0.0f);
            p.spherical = wv.value("spherical", 0.0f);
        }
        return p;
    } catch (...) {
        return std::nullopt;
    }
}

}  // namespace lens::data
