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
