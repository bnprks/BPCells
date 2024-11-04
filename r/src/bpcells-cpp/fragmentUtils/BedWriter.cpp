// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "BedWriter.h"
#include <zlib.h>
#include "../fragmentIterators/FragmentIterator.h"
#include "InsertionIterator.h"
#include "../utils/filesystem_compat.h"
#include "../utils/gzfile_wrapper.h"
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

// Write bedgraph coverage files for insertions computed from fragment pseudobulks.
// Args:
// - fragments: source of fragments to convert to insertions & calculate coverage
// - cell_groups: For each cell in fragments, the index of the pseudobulk to assign it to
// - output_paths: The file path to save the bedgraph for each pseudobulk
// - mode: 0 = include start + end coords, 1 = include just start coords, 2 = include just end coords
void writeInsertionBedgraph(
    FragmentLoader &fragments,
    const std::vector<uint32_t> &cell_groups,
    const std::vector<std::string> &output_paths,
    const BedgraphInsertionMode mode,
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
        
    // Last bp observed from each group
    std::vector<uint32_t> last_base(output_paths.size(), UINT32_MAX);
    
    // Running insertion tally for each group
    std::vector<uint32_t> tally(output_paths.size(), 0);
    
    while (it.nextChr()) {
        if (fragments.chrNames(it.chr()) == NULL) {
            throw std::runtime_error("writeBedgraph: No chromosome name found for ID: " + std::to_string(it.chr()));
        }
        std::string chr_name(fragments.chrNames(it.chr()));
        uint32_t count = 0;

        for (size_t i = 0; i < output_paths.size(); i++) {
            last_base[i] = UINT32_MAX;
            tally[i] = 0;
        }

        while (it.nextInsertion()) {
            if (mode == BedgraphInsertionMode::StartOnly && !it.isStart()) continue;
            if (mode == BedgraphInsertionMode::EndOnly && it.isStart()) continue;

            int cell_group = cell_groups[it.cell()];
            
            if (last_base[cell_group] == it.coord()) {
                tally[cell_group] += 1;
                continue;
            }
        
            // Handle the case where this is not the same basepair we saw last
            if (last_base[cell_group] != UINT32_MAX) {
                uint32_t bytes_written = gzprintf(
                    *files[cell_group],
                    "%s\t%d\t%d\t%d\n",
                    chr_name.c_str(),
                    last_base[cell_group],
                    last_base[cell_group] + 1,
                    tally[cell_group]
                );

                if (bytes_written <= 0) {
                    throw std::runtime_error("writeBedgraph: Failed to write data");
                }
            }

            tally[cell_group] = 1;
            last_base[cell_group] = it.coord();
            if (user_interrupt != nullptr && count++ % 65536 == 0 && *user_interrupt) return;
        }
        // Cleanup at the end of the chromosome
        for (size_t i = 0; i < files.size(); i++) {
            if (last_base[i] == UINT32_MAX) continue;
            uint32_t bytes_written = gzprintf(
                *files[i],
                "%s\t%d\t%d\t%d\n",
                chr_name.c_str(),
                last_base[i],
                last_base[i] + 1,
                tally[i]
            );
            if (bytes_written <= 0) {
                throw std::runtime_error("writeBedgraph: Failed to write data");
            }

            last_base[i] = 0;
            tally[i] = 0;
        }
    }
}

} // end namespace BPCells