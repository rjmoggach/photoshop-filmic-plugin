#include "lenscore/pfm.hpp"
#include "lenscore/pipeline.hpp"
#include "lensdata/lensfile.hpp"
#include <cstdio>
#include <cstdlib>
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
