// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <atomic>
#include <vector>

#include "../fragmentIterators/FragmentIterator.h"

namespace BPCells {


enum class BedgraphInsertionMode {
    Both,
    StartOnly,
    EndOnly
};

// Normalization method on insertions for bedgraph coverage values
// NFrags = Total # insertions in a pseudobulk, scaled by 1e6
// NCells = Total # cells in a pseudobulk
// None = No normalization
enum class PseudobulkNormalizationMethod {
  CPM,
  NCells,
  None
};

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
);

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
);

} // end namespace BPCells