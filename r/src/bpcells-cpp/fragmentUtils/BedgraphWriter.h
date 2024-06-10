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

// Write bedgraph coverage files for insertions computed from fragment pseudobulks.
// Note that
// Args:
// - fragments: source of fragments to convert to insertions & calculate coverage
// - cell_groups: For each cell in fragments, the index of the pseudobulk to assign it to
// - output_paths: The file path to save the bedgraph for each pseudobulk
// - mode: Which combination of start/end coordinates to include

void writeInsertionBedgraph(
    FragmentLoader &fragments,
    const std::vector<uint32_t> &cell_groups,
    const std::vector<std::string> &output_paths,
    const BedgraphInsertionMode mode = BedgraphInsertionMode::Both,
    std::atomic<bool> *user_interrupt = NULL
);

} // end namespace BPCells