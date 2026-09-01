#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "lenscore/version.hpp"

TEST_CASE("version string is present") {
    CHECK(std::string(lens::kVersion) == "0.1.0");
}
