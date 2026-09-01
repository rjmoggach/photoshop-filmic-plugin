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
