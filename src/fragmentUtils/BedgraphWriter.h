#pragma once

#include <vector>

#include "../fragmentIterators/FragmentIterator.h"

namespace BPCells {

// Write bedgraph coverage files for insertions computed from fragment pseudobulks.
// Note that
// Args:
// - fragments: source of fragments to convert to insertions & calculate coverage
// - cell_groups: For each cell in fragments, the index of the pseudobulk to assign it to,
//     or UINT32_MAX if the cell should be dropped from coverage calculations
// - output_paths: The file path to save the bedgraph for each pseudobulk
// - smooth_bp: How many bp to smooth the coverage by

void writeBedgraph(
    FragmentLoader &fragments,
    std::vector<uint32_t> cell_groups,
    std::vector<std::string> output_paths,
    uint32_t smooth_bp
) {}

} // end namespace BPCells