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


// Write bedgraph coverage files for insertions computed from fragment pseudobulks with possible tiling.
// Args:
// - fragments: source of fragments to convert to insertions & calculate coverage
// - cell_groups: For each cell in fragments, the index of the pseudobulk to assign it to
// - tile_width: Width of each tile in the bedgraph
// - output_paths: The file path to save the bedgraph for each pseudobulk
// - mode: 0 = include start + end coords, 1 = include just start coords, 2 = include just end coords
// - normalization_method:  Normalization method for coverage values.   One of "None", "NFrags" (CPM-like), or "NCells".
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
    for (const auto &x : cell_groups) {
        if (x >= output_paths.size()) {
            throw std::runtime_error("writeBedgraph: found cell_group larger than number of output paths: " + std::to_string(x));
        }
    }

    InsertionIterator it(fragments);
    std::vector<gzFileWrapper> files;
    const uint32_t buffer_size = 1 << 20;
    for (const auto &path : output_paths) {
        files.push_back(gzFileWrapper(path, "w", buffer_size));
    }
        
    // Last start of the tile used from each group
    std::vector<uint32_t> last_tile_start(output_paths.size(), UINT32_MAX);
    // Running insertion tally for each start
    std::vector<uint32_t> tally(output_paths.size(), 0);
    // Total insertions per group (for normalization)
    std::vector<uint64_t> total_insertions(output_paths.size(), 0);
    // Count insertion starts if mode is start/both, ends if mode is end/both
    auto accept = [&](bool isStart) -> int {
        switch (mode) {
            case BedgraphInsertionMode::StartOnly: return isStart;
            case BedgraphInsertionMode::EndOnly:   return !isStart;
            case BedgraphInsertionMode::Both:      return 1;
        }
        return 0;
    };
    if (normalization_method == PseudobulkNormalizationMethod::CPM) {
        // First pass to count total insertions per group
        while (it.nextChr()) {
            while (it.nextInsertion()) {
                total_insertions[cell_groups[it.cell()]] += accept(it.isStart());
            }
        }
        it.restart();
    }
    // Count number of cells per group
    std::vector<uint32_t> cell_counter(output_paths.size(), 0);
    if (normalization_method == PseudobulkNormalizationMethod::NCells) {
        for (const auto &x: cell_groups) cell_counter[x]++;
    }
    // Write insertions to tiled bedgraph
    while (it.nextChr()) {
        if (fragments.chrNames(it.chr()) == NULL) {
            throw std::runtime_error("writeInsertionBedgraph: No chromosome name found for ID: " + std::to_string(it.chr()));
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
            if (!accept(it.isStart())) continue;
            // Determine the start of the tile this insertion belongs to
            uint32_t curr_tile_start = std::floor(it.coord() / tile_width) * tile_width;
            int cell_group = cell_groups[it.cell()];

            if (last_tile_start[cell_group] == curr_tile_start) {
                tally[cell_group] += 1;
                continue;
            }
        
            // Handle the case where this is not the same tile we saw last
            if (last_tile_start[cell_group] != UINT32_MAX) {
                // not sure what type for val
                double val = tally[cell_group];
                if (normalization_method == PseudobulkNormalizationMethod::CPM) {
                    val = (double)val * 1e6 / (double)total_insertions[cell_group];
                } else if (normalization_method == PseudobulkNormalizationMethod::NCells) {
                    val = (double)val / (double)cell_counter[cell_group];
                }
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
            double val = tally[i];
            if (normalization_method == PseudobulkNormalizationMethod::CPM) {
                val = (double)val * 1e6 / (double)total_insertions[i];
            } else if (normalization_method == PseudobulkNormalizationMethod::NCells) {
                val = (double)val / (double)cell_counter[i];
            } 
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
}

} // end namespace BPCells