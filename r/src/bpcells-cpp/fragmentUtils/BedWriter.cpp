// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "BedWriter.h"
#include "../fragmentIterators/FragmentIterator.h"
#include "InsertionIterator.h"
#include "../utils/filesystem_compat.h"
#include "../utils/gzfile_wrapper.h"
#include <zlib.h>
#include <functional>
#include <cmath>
namespace BPCells {

// Write bedgraph files including chr, start, end with only 1bp insertions,
// only considering a subset of cells given by the cells vector. Leave duplicates in the bedfile.
// Args:
// - fragments: source of fragments to convert to insertions
// - output_path: The file path to save the bedgraph
// - mode: StartOnly = include only start coords, EndOnly = include only end coords, Both = include start + end coords
void writeInsertionBed(
    FragmentLoader &fragments,
    const std::string &output_path,
    const BedgraphInsertionMode &mode,
    std::atomic<bool> *user_interrupt
) {
    InsertionIterator it(fragments);
    const uint32_t buffer_size = 1 << 20;
    gzFileWrapper file(output_path, "w", buffer_size);
    while (it.nextChr()) {
        if (fragments.chrNames(it.chr()) == NULL) {
            throw std::runtime_error("writeInsertionBed: No chromosome name found for ID: " + std::to_string(it.chr()));
        }
        std::string chr_name(fragments.chrNames(it.chr()));
        uint32_t count = 0;
        while (it.nextInsertion()) {
            if (mode == BedgraphInsertionMode::StartOnly && !it.isStart()) continue;
            if (mode == BedgraphInsertionMode::EndOnly && it.isStart()) continue;
            uint32_t bytes_written = gzprintf(
                *file,
                "%s\t%d\t%d\n",
                chr_name.c_str(),
                it.coord(),
                it.coord() + 1
            );
            if (bytes_written <= 0) {
                throw std::runtime_error("writeInsertionBed: Failed to write data");
            }
            if (user_interrupt != nullptr && count++ % 65536 == 0 && *user_interrupt) return;
        }
    }
}

// Helper function to determine whether an insertion should be accepted based on mode
// ie accept start if mode is StartOnly or Both, accept end if mode is EndOnly or Both
template<BedgraphInsertionMode mode>
static inline bool keepInsertion(bool isStart) noexcept {
    if constexpr (mode == BedgraphInsertionMode::Both) return true;
    else if constexpr (mode == BedgraphInsertionMode::StartOnly) return isStart;
    else return !isStart; // EndOnly
}

template<BedgraphInsertionMode mode>
static inline void count_cpm_per_group(
    InsertionIterator &it,
    const std::vector<uint32_t> &cell_groups,
    std::vector<uint64_t> &total_insertions,
    std::atomic<bool> *user_interrupt = nullptr
) {
    uint32_t count = 0;
    while (it.nextChr()) {
        while (it.nextInsertion()) {
            if (!keepInsertion<mode>(it.isStart())) continue;
            total_insertions[cell_groups[it.cell()]] += 1;
            if (user_interrupt != nullptr && count++ % 65536 == 0 && *user_interrupt) return;
        }
    }
    it.restart();
}

// Compute normalization factors for each group based on the normalization method
// Args:
// - fragments: source of fragments to convert to insertions & calculate coverage
// - cell_groups: For each cell in fragments, the index of the pseudobulk to assign it to
// - mode: 0 = include start + end coords, 1 = include just start coords, 2 = include just end coords
// - normalization_method:  Normalization method for coverage values.
template<BedgraphInsertionMode mode>
std::vector<double> compute_group_normalization_factors(
    InsertionIterator &it,
    const std::vector<uint32_t> &cell_groups,
    const PseudobulkNormalizationMethod normalization_method,
    size_t n_groups
) {
    std::vector<double> group_norm_factors(n_groups, 1.0);
    std::vector<uint32_t> cell_counter(n_groups, 0);
    if (normalization_method == PseudobulkNormalizationMethod::CPM) {
        // Total insertions per group
        std::vector<uint64_t> total_insertions(n_groups, 0);
        // Count insertion starts if mode is start/both, ends if mode is end/both
        count_cpm_per_group<mode>(it, cell_groups, total_insertions);
        it.restart();
        for (size_t i = 0; i < n_groups; i++) {
            if (total_insertions[i] > 0) {
                group_norm_factors[i] =  1e6 / (double)total_insertions[i];
            }
        }
    } else if (normalization_method == PseudobulkNormalizationMethod::NCells) {
        for (const auto &x: cell_groups) cell_counter[x]++;
        for (size_t i = 0; i < n_groups; i++) {
            if (cell_counter[i] > 0) {
                group_norm_factors[i] = 1 / (double)cell_counter[i];
            }
        }
    }
    return group_norm_factors;
}

// Template implementation of writeInsertionBedgraph for a specific mode
template<BedgraphInsertionMode mode>
void writeInsertionBedgraph_impl(
    FragmentLoader &fragments,
    const std::vector<uint32_t> &cell_groups,
    const uint32_t& tile_width,
    const std::vector<std::string> &output_paths,
    const PseudobulkNormalizationMethod normalization_method,
    const std::vector<uint32_t> *chrom_sizes,
    std::atomic<bool> *user_interrupt
) {
    for (const auto &x : cell_groups) {
        if (x >= output_paths.size()) {
            throw std::runtime_error("writeInsertionBedgraph: found cell_group larger than number of output paths: " + 
                std::to_string(x));
        }
    }
    InsertionIterator it(fragments);
    std::vector<gzFileWrapper> files;
    const uint32_t buffer_size = 1 << 20;
    for (const auto &path : output_paths) {
        files.push_back(gzFileWrapper(path, "w", buffer_size));
    }
    std::vector<double> group_norm_factors = compute_group_normalization_factors<mode>(
        it, cell_groups, normalization_method, output_paths.size()
    );
    // Last start of the tile used from each group
    std::vector<uint32_t> last_tile_start(output_paths.size(), UINT32_MAX);
    // Running insertion tally for each start
    std::vector<uint32_t> tally(output_paths.size(), 0);

    // Write insertions to tiled bedgraph
    while (it.nextChr()) {
        if (fragments.chrNames(it.chr()) == NULL) {
            throw std::runtime_error("writeInsertionBedgraph: No chromosome name found for ID: " + 
                std::to_string(it.chr()));
        }
        uint32_t chrom_end = chrom_sizes != nullptr ? (*chrom_sizes)[it.chr()] : UINT32_MAX;
        std::string chr_name(fragments.chrNames(it.chr()));
        uint32_t count = 0;

        // reset for every chromosone
        for (size_t i = 0; i < output_paths.size(); i++) {
            last_tile_start[i] = UINT32_MAX;
            tally[i] = 0;
        }

        while (it.nextInsertion()) {
            if (!keepInsertion<mode>(it.isStart())) continue;
            // Determine the start of the tile this insertion belongs to
            // this should be integer division, but static cast just in case
            uint32_t curr_tile_start = static_cast<uint32_t>(it.coord() / tile_width) * tile_width;
            int cell_group = cell_groups[it.cell()];

            if (last_tile_start[cell_group] == curr_tile_start) {
                tally[cell_group] += 1;
                continue;
            }
        
            // Handle the case where this is not the same tile we saw last
            if (last_tile_start[cell_group] != UINT32_MAX) {
                double val = tally[cell_group] * group_norm_factors[cell_group];
                uint32_t bytes_written = gzprintf(
                    *files[cell_group],
                    "%s\t%d\t%d\t%.4f\n",
                    chr_name.c_str(),
                    last_tile_start[cell_group],
                    std::min(last_tile_start[cell_group] + tile_width, chrom_end),
                    val
                );
                if (bytes_written <= 0) {
                    throw std::runtime_error("writeInsertionBedgraph: Failed to write data");
                }
            }

            tally[cell_group] = 1;
            // last base is now the current tile
            last_tile_start[cell_group] = curr_tile_start;
            if (user_interrupt != nullptr && count++ % 65536 == 0 && *user_interrupt) return;
        }
        // Cleanup at the end of the chromosome
        for (size_t i = 0; i < files.size(); i++) {
            if (last_tile_start[i] == UINT32_MAX) continue;
            double val = tally[i] * group_norm_factors[i];
            uint32_t bytes_written = gzprintf(
                *files[i],
                "%s\t%d\t%d\t%.4f\n",
                chr_name.c_str(),
                last_tile_start[i],
                std::min(last_tile_start[i] + tile_width, chrom_end),
                val
            );
            if (bytes_written <= 0) {
                throw std::runtime_error("writeInsertionBedgraph: Failed to write data");
            }
        }
    }
};

// Write bedgraph coverage files for insertions computed from fragment pseudobulks with possible tiling and normalization.
// Args:
// - fragments: source of fragments to convert to insertions & calculate coverage
// - cell_groups: For each cell in fragments, the index of the pseudobulk to assign it to
// - tile_width: Width of each tile in the bedgraph
// - output_paths: The file path to save the bedgraph for each pseudobulk
// - mode: 0 = include start + end coords, 1 = include just start coords, 2 = include just end coords
// - normalization_method:  Normalization method for coverage values.
// - chrom_sizes: Total size of each chromosome, used to determine when to the final size of the last tile.
void writeInsertionBedgraph(
    FragmentLoader &fragments,
    const std::vector<uint32_t> &cell_groups,
    const uint32_t& tile_width,
    const std::vector<std::string> &output_paths,
    const BedgraphInsertionMode mode,
    const PseudobulkNormalizationMethod normalization_method,
    const std::vector<uint32_t> *chrom_sizes,
    std::atomic<bool> *user_interrupt
) {
    if (tile_width == 0) {
        throw std::runtime_error("writeInsertionBedgraph: tile_width must be > 0");
    }
    if (mode == BedgraphInsertionMode::Both) {
        writeInsertionBedgraph_impl<BedgraphInsertionMode::Both>(
            fragments, cell_groups, tile_width, output_paths, normalization_method, chrom_sizes, user_interrupt
        );
    } else if (mode == BedgraphInsertionMode::StartOnly) {
        writeInsertionBedgraph_impl<BedgraphInsertionMode::StartOnly>(
            fragments, cell_groups, tile_width, output_paths, normalization_method, chrom_sizes, user_interrupt
        );
    } else { // EndOnly
        writeInsertionBedgraph_impl<BedgraphInsertionMode::EndOnly>(
            fragments, cell_groups, tile_width, output_paths, normalization_method, chrom_sizes, user_interrupt
        );
    }
};

} // end namespace BPCells