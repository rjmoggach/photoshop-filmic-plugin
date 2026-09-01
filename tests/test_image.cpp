#include <doctest/doctest.h>
#include "lenscore/image.hpp"
#include "lenscore/pfm.hpp"
#include <cstdio>
#include <cstddef>
#include <string>

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

TEST_CASE("pfm read rejects malformed header") {
    const std::string p = "malformed.pfm";
    {
        std::FILE* f = std::fopen(p.c_str(), "wb");
        std::fprintf(f, "XX\n");  // wrong magic
        std::fclose(f);
    }
    auto result = lens::pfm::read(p);
    CHECK(!result.has_value());
    std::remove(p.c_str());
}

TEST_CASE("pfm read rejects truncated payload") {
    const std::string p = "truncated.pfm";
    {
        std::FILE* f = std::fopen(p.c_str(), "wb");
        std::fprintf(f, "PF\n2 2\n-1.0\n");
        // Write only 1 float instead of 2*2*3 = 12 floats
        float val = 1.5f;
        std::fwrite(&val, sizeof(float), 1, f);
        std::fclose(f);
    }
    auto result = lens::pfm::read(p);
    CHECK(!result.has_value());
    std::remove(p.c_str());
}

TEST_CASE("pfm read rejects big-endian (positive scale)") {
    const std::string p = "bigendian.pfm";
    {
        std::FILE* f = std::fopen(p.c_str(), "wb");
        std::fprintf(f, "PF\n2 2\n1.0\n");  // positive scale = big-endian
        // Write the data anyway (should be rejected before reading)
        float data[12] = {0.0f};
        std::fwrite(data, sizeof(float), 12, f);
        std::fclose(f);
    }
    auto result = lens::pfm::read(p);
    CHECK(!result.has_value());
    std::remove(p.c_str());
}
