// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <atomic>

#include "FootprintMatrix.h"
#include "../simd/math.h"

namespace BPCells {

Eigen::MatrixXd footprintMatrix(
    FragmentLoader &frags,
    const std::vector<uint32_t> &chr,
    const std::vector<uint32_t> &center,
    const std::vector<int32_t> &strand,
    uint32_t flank_width,
    std::unique_ptr<StringReader> &&chr_levels,
    const std::vector<uint32_t> &cell_groups,
    const std::vector<double> &cell_weights,
    std::atomic<bool> *user_interrupt
) {
    class Region {
      public:
        uint32_t chr, start, end;
        bool pos_strand;
    };

    std::vector<Region> sorted_regions;
    std::vector<Region> active_regions;
    uint32_t next_active_region = 0;

    // **********************
    // Check input arguments
    // **********************
    if (frags.cellCount() < 0)
        throw std::invalid_argument(
            "frags must have a known cell count. Consider using a cell selection to define the "
            "number of cells."
        );
    if (frags.cellCount() != (int64_t) cell_groups.size() || frags.cellCount() != (int64_t) cell_weights.size()) {
        throw std::invalid_argument(
            "frags must have same cell count as cell_groups and cell_weights"
        );
    }

    if (chr.size() != center.size() || chr.size() != strand.size())
        throw std::invalid_argument("chr, center, and strand must all be same length");

    // Check that chr name matches for all the available chrNames in frags
    for (uint32_t i = 0; i < chr_levels->size(); i++) {
        const char *chr_name_frag = frags.chrNames(i);
        const char *chr_name_args = chr_levels->get(i);
        if (chr_name_frag != NULL &&
            (chr_name_args == NULL || strcmp(chr_name_frag, chr_name_args) != 0)) {
            throw std::runtime_error(
                std::string("FootprintMatrix encountered fragment with incorrect chrLevel: ") +
                std::string(chr_name_frag) + std::string(" expected: ") + std::string(chr_name_args)
            );
        }
    }

    for (size_t i = 0; i < chr.size(); i++) {
        if (chr[i] >= chr_levels->size())
            throw std::invalid_argument(
                "FootprintMatrix: chr has values higher than length of chr_levels"
            );
        Region r;
        if (center[i] < flank_width)
            throw std::invalid_argument("FootprintMatrix: flank_width expands to negative bases");
        r.start = center[i] - flank_width;
        r.end = center[i] + flank_width + 1;
        r.chr = chr[i];
        r.pos_strand = strand[i] == 1;
        if (strand[i] != 1 && strand[i] != -1)
            throw std::invalid_argument("strand must have values of only +/- 1");
        sorted_regions.push_back(r);
    }

    uint32_t max_group = 0;
    for (auto g : cell_groups)
        max_group = std::max(max_group, g);

    // Sentinel value at end of sorted_peaks
    sorted_regions.push_back({UINT32_MAX, UINT32_MAX, UINT32_MAX, 0});

    std::sort(sorted_regions.begin(), sorted_regions.end(), [](Region a, Region b) {
        if (a.chr != b.chr) return a.chr < b.chr;
        else if (a.start != b.start) return a.start < b.start;
        else return a.end < b.end;
    });

    Eigen::MatrixXd out(max_group + 1, 2 * flank_width + 1);
    out.setZero();

    // **********************
    // Loop through fragments
    // **********************
    uint32_t prev_chr_id = 0;
    frags.restart();
    while (frags.nextChr()) {
        if (frags.currentChr() < prev_chr_id) {
            throw std::runtime_error(
                "FootprintMatrix encountered fragments with out of order chromosome IDs. Please "
                "save + load fragments before passing to FootprintMatrix to fix this issue."
            );
        }
        prev_chr_id = frags.currentChr();
        // Check that chr name matches
        const char *chr_name_frag = frags.chrNames(frags.currentChr());
        const char *chr_name_args = chr_levels->get(frags.currentChr());
        if (chr_name_frag == NULL || chr_name_args == NULL ||
            strcmp(chr_name_frag, chr_name_args) != 0) {
            throw std::runtime_error(
                std::string("FootprintMatrix encountered fragment with incorrect chrLevel: ") +
                std::string(chr_name_frag) + std::string(" expected: ") + std::string(chr_name_args)
            );
        }
        while (sorted_regions[next_active_region].chr < frags.currentChr()) {
            next_active_region++;
        }
        active_regions.clear();
        while (frags.load()) {
            uint32_t capacity = frags.capacity();
            const uint32_t *start_data = frags.startData();
            const uint32_t *end_data = frags.endData();
            const uint32_t *cell_data = frags.cellData();
            uint32_t i = 0;
            uint32_t end_max = 0;

            if (user_interrupt != nullptr && *user_interrupt) return out;

            // Loop through reads in blocks of 128 at a time
            while (i < capacity) {
                uint32_t items = std::min(128U, capacity - i);
                end_max = std::max(end_max, simd::max(end_data, items));

                // Check for new peaks to activate
                while (sorted_regions[next_active_region].chr == frags.currentChr() &&
                       sorted_regions[next_active_region].start < end_max) {

                    active_regions.push_back(sorted_regions[next_active_region]);
                    next_active_region += 1;
                }

                // For each active peak, iterate through the fragments & tally overlaps
                for (uint32_t j = 0; j < active_regions.size(); j++) {
                    Region r = active_regions[j];

                    uint32_t k = 0;

                    for (; k < items && start_data[i + k] < r.end; k++) {
                        if (start_data[i + k] >= r.start && start_data[i + k] < r.end) {
                            uint32_t cell_group = cell_groups[cell_data[i + k]];
                            uint32_t base;
                            if (r.pos_strand) base = start_data[i + k] - r.start;
                            else base = r.end - 1 - start_data[i + k];
                            out(cell_group, base) += cell_weights[cell_data[i + k]];
                        }
                        if (end_data[i + k] > r.start && end_data[i + k] <= r.end) {
                            uint32_t cell_group = cell_groups[cell_data[i + k]];
                            uint32_t base;
                            if (r.pos_strand) base = end_data[i + k] - 1 - r.start;
                            else base = r.end - end_data[i + k];
                            out(cell_group, base) += cell_weights[cell_data[i + k]];
                        }
                    }

                    // Remove the peak from active_regions if we're done, and mark the
                    // next completed peak on our list
                    if (k < items) {
                        std::swap(active_regions.back(), active_regions[j]);
                        active_regions.pop_back();
                        j -= 1;
                    }
                }
                i += items;
            }
        }
    }
    return out;
}

} // end namespace BPCells
