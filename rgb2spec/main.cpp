#include "lenscore/color/spectable.hpp"
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

int main(int argc, char** argv) {
    const int res = (argc > 1) ? std::atoi(argv[1]) : 64;
    const char* out = (argc > 2) ? argv[2] : "lensdata/rgb2spec/rec2020.bin";
    std::printf("building %d^3 table...\n", res);
    lens::color::SpecTable t;
    try {
        t = lens::color::buildTable(res);
    } catch (const std::invalid_argument& e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
    if (!lens::color::writeTable(out, t)) { std::fprintf(stderr, "write failed: %s\n", out); return 1; }
    std::printf("wrote %s (%zu floats)\n", out, t.data.size());
    return 0;
}
