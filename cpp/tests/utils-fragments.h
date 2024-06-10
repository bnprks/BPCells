// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <memory>
#include <random>

#include <arrayIO/vector.h>
#include <fragmentIterators/FragmentIterator.h>
#include <fragmentIterators/StoredFragments.h>

namespace Testing {

class Frag {
  public:
    uint32_t chr, start, end, cell;
};

// Create a Fragment object from a vector of fragments (sorting not required)
std::unique_ptr<BPCells::VecReaderWriterBuilder> writeFragmentTuple(
    std::vector<Frag> frag_vec,
    uint32_t min_cell_count = 0,
    bool skip_sort = false,
    uint32_t chunk_size = 1024
) {
    using namespace BPCells;

    VecReaderWriterBuilder v;
    UIntWriter w_cell = v.createUIntWriter("cell");
    UIntWriter w_start = v.createUIntWriter("start");
    UIntWriter w_end = v.createUIntWriter("end");
    UIntWriter w_end_max = v.createUIntWriter("end_max");
    UIntWriter w_chr_ptr = v.createUIntWriter("chr_ptr");
    v.writeVersion("unpacked-fragments-v1");

    std::unique_ptr<StringWriter> w_chr_names = v.createStringWriter("chr_names");
    std::unique_ptr<StringWriter> w_cell_names = v.createStringWriter("cell_names");

    // Sort by (chr, start, end, cell) to ensure consistent ordering
    if (!skip_sort) {
        std::sort(frag_vec.begin(), frag_vec.end(), [](const Frag &a, const Frag &b) {
            if (a.chr != b.chr) return a.chr < b.chr;
            if (a.start != b.start) return a.start < b.start;
            if (a.end != b.end) return a.end < b.end;
            return a.cell < b.cell;
        });
    }

    uint32_t current_chr = 0;
    w_chr_ptr.write_one(0);
    uint32_t count = 0;
    for (auto &f : frag_vec) {
        while (f.chr > current_chr) {
            w_chr_ptr.write_one(count);
            w_chr_ptr.write_one(count);
            current_chr += 1;
        }
        w_cell.write_one(f.cell);
        w_start.write_one(f.start);
        w_end.write_one(f.end);
        count++;
    }
    w_chr_ptr.write_one(count);

    for (int i = 0; i < count; i += 128) {
        w_end_max.write_one(UINT32_MAX);
    }

    std::vector<std::string> chr_names;
    for (int i = 0; i <= current_chr; i++) {
        chr_names.push_back("chr" + std::to_string(i));
    }

    for (auto &f : frag_vec)
        min_cell_count = std::max(min_cell_count, f.cell + 1);
    std::vector<std::string> cell_names;
    for (int i = 0; i < min_cell_count; i++) {
        cell_names.push_back("c" + std::to_string(i));
    }

    w_chr_names->write(VecStringReader(chr_names));
    w_cell_names->write(VecStringReader(cell_names));
    w_cell.finalize();
    w_start.finalize();
    w_end.finalize();
    w_end_max.finalize();
    w_chr_ptr.finalize();
    StoredFragments manual_frags = StoredFragments::openUnpacked(v);

    auto ret = std::make_unique<VecReaderWriterBuilder>(chunk_size);
    StoredFragmentsWriter::createUnpacked(*ret).write(manual_frags);

    return ret;
}

std::vector<Frag> generateFrags(
    uint32_t n,
    uint32_t max_chr,
    uint32_t max_coord,
    uint32_t max_cell,
    uint32_t max_width,
    uint32_t seed = 12548
) {
    std::minstd_rand gen(seed
    ); // Linear congruential engine (not super random, but probably enough for this)
    std::uniform_int_distribution<> chr(0, max_chr);
    std::uniform_int_distribution<> cell(0, max_cell);
    std::uniform_int_distribution<> width(1, max_width);
    std::uniform_int_distribution<> start(0, max_coord); // 1/5 chance of being non-zero

    std::vector<Frag> ret;
    for (uint32_t i = 0; i < n; i++) {
        uint32_t s = start(gen);
        ret.push_back(Frag{(uint32_t)chr(gen), s, s + width(gen), (uint32_t)cell(gen)});
    }
    return ret;
}

bool fragments_identical(BPCells::FragmentLoader &fragments1, BPCells::FragmentLoader &fragments2) {
    using namespace BPCells;
    fragments1.restart();
    fragments2.restart();
    FragmentIterator i1((std::unique_ptr<FragmentLoader>(&fragments1)));
    FragmentIterator i2((std::unique_ptr<FragmentLoader>(&fragments2)));
    i1.preserve_input_loader();
    i2.preserve_input_loader();

    while (true) {
        bool res1 = i1.nextChr();
        bool res2 = i2.nextChr();
        if (res1 != res2) {
            std::cout << "Different number of remaining chromosomes." << std::endl;
            return false;
        }
        if (!res1) break;
        if (i1.currentChr() != i2.currentChr()) {
            std::cout << "Different chromosome ID loaded" << std::endl;
            return false;
        }
        while (true) {
            bool res1 = i1.nextFrag();
            bool res2 = i2.nextFrag();
            if (res1 != res2) {
                std::cout << "Different number of fragments in chromosome." << std::endl;
                return false;
            }
            if (!res1) break;
            if (i1.cell() != i2.cell() || i1.start() != i2.start() || i1.end() != i2.end()) {
                printf(
                    "Mismatched fragments: %s(id=%d):%d-%d vs. %s(id=%d):%d-%d\n",
                    i1.cellNames(i1.cell()),
                    i1.cell(),
                    i1.start(),
                    i1.end(),
                    i2.cellNames(i2.cell()),
                    i2.cell(),
                    i2.start(),
                    i2.end()
                );
                return false;
            }
        }
    }
    for (uint32_t i = 0;; i++) {
        const char *n1 = i1.cellNames(i);
        const char *n2 = i2.cellNames(i);
        if ((n1 == NULL) != (n2 == NULL)) {
            std::cout << "Mismatched number of cell names" << std::endl;
            return false;
        }
        if (n1 == NULL) break;
        if (strcmp(n1, n2) != 0) {
            std::cout << "Mismatched cell names" << std::endl;
            return false;
        }
    }
    for (uint32_t i = 0;; i++) {
        const char *n1 = i1.chrNames(i);
        const char *n2 = i2.chrNames(i);
        if ((n1 == NULL) != (n2 == NULL)) {
            std::cout << "Mismatched number of chr names" << std::endl;
            return false;
        }
        if (n1 == NULL) break;
        if (strcmp(n1, n2) != 0) {
            std::cout << "Mismatched  chr names" << std::endl;
            return false;
        }
    }
    return true;
}

} // namespace Testing