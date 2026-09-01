#pragma once
#include "lenscore/constants.hpp"
#include "lenscore/params.hpp"
#include <fstream>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>
#include <vector>

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
    p.pupil.rotationRad = j["aperture"].value("rotation_deg", 0.0f) * kPi / 180.0f;

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
