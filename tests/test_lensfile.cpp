#include <doctest/doctest.h>
#include "lensdata/lensfile.hpp"
#include <cstdio>
#include <fstream>
#include <functional>
#include <string>

using namespace lens;

// Writes a scratch file, runs the assertion, and removes it -- mirrors
// test_spectable.cpp's pattern for temporary on-disk fixtures.
static void withFile(const std::string& path, const std::string& content,
                      const std::function<void()>& body) {
    { std::ofstream out(path); out << content; }
    body();
    std::remove(path.c_str());
}

TEST_CASE("loadLensFile returns nullopt, not a throw, when aperture is missing") {
    // The bug this guards: aperture used to be read via the non-const
    // operator[], which inserts a null node for a missing key, and .value()
    // on that null then threw nlohmann::json::type_error out of a function
    // whose signature promises std::optional.
    withFile("missing_aperture.lens", R"({
        "schema": 2,
        "focal_mm": 32.0
    })", [] {
        std::optional<Params> p;
        CHECK_NOTHROW(p = data::loadLensFile("missing_aperture.lens"));
        CHECK_FALSE(p.has_value());
    });
}

TEST_CASE("loadLensFile returns nullopt, not a throw, when aperture is not an object") {
    withFile("bad_aperture.lens", R"({
        "schema": 2,
        "focal_mm": 32.0,
        "aperture": "t/2.0"
    })", [] {
        std::optional<Params> p;
        CHECK_NOTHROW(p = data::loadLensFile("bad_aperture.lens"));
        CHECK_FALSE(p.has_value());
    });
}

TEST_CASE("loadLensFile tolerates a non-object dispersion without throwing") {
    // dispersion (unlike aperture) is optional: a malformed one should fall
    // back to Params' default glass, not fail the whole file.
    withFile("bad_dispersion.lens", R"({
        "schema": 2,
        "focal_mm": 32.0,
        "aperture": { "t_stop": 2.0 },
        "dispersion": 42
    })", [] {
        std::optional<Params> p;
        CHECK_NOTHROW(p = data::loadLensFile("bad_dispersion.lens"));
        REQUIRE(p.has_value());
        CHECK(p->fNumberWide == doctest::Approx(2.0f));
        CHECK(p->dispersion.residual == doctest::Approx(1.0f));
    });
}

TEST_CASE("loadLensFile returns nullopt for a non-object document") {
    withFile("not_an_object.lens", "[1, 2, 3]", [] {
        std::optional<Params> p;
        CHECK_NOTHROW(p = data::loadLensFile("not_an_object.lens"));
        CHECK_FALSE(p.has_value());
    });
}

TEST_CASE("loadLensFile returns nullopt for an unreadable path") {
    CHECK_FALSE(data::loadLensFile("/does/not/exist.lens").has_value());
}

TEST_CASE("loadLensFile returns nullopt for invalid json") {
    withFile("garbage.lens", "{not json at all", [] {
        std::optional<Params> p;
        CHECK_NOTHROW(p = data::loadLensFile("garbage.lens"));
        CHECK_FALSE(p.has_value());
    });
}

TEST_CASE("loadLensFile reads the shipped Cooke S4/i preset") {
    // LENS_SOURCE_DIR (tests/CMakeLists.txt) is the repo root: ctest runs
    // with the build directory as its working directory, not the source
    // tree, so a path relative to CMAKE_SOURCE_DIR is used instead of a
    // bare relative one that would resolve against the wrong directory.
    const auto p = data::loadLensFile(std::string(LENS_SOURCE_DIR) + "/lensdata/cooke-s4-32mm.lens");
    REQUIRE(p.has_value());
    CHECK(p->focal_mm == doctest::Approx(32.0f));
    CHECK(p->fNumberWide == doctest::Approx(2.0f));
    CHECK(p->pupil.blades == 9);
    CHECK(p->dispersion.correction_nm.size() == 2u);
    CHECK(p->petzval == doctest::Approx(0.42f));
}
