// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include <vector>

#include "FragmentIterator.h"

namespace BPCells {

// Merge several fragment sources into a single sorted loader.
// Cell IDs will be adjusted to be sequential between fragments
// Chromosomes will be merged in the order of the provided chromosome ordering. It
// is an error if any input has chromosmes not in the output ordering.
// All inputs must have known cell counts, chr counts+names, and be seekable
class MergeFragments : public FragmentLoader {
  private:
    // Helper class to handle reading a fixed number of fragments at a time,
    // but in bulk.
    class ChunkedLoader {
        std::unique_ptr<FragmentLoader> frags;
        uint32_t loaded = 0; // How many fragments are currently loaded and available to copy
        const uint32_t cell_offset, chunk_size; 

      public:
        ChunkedLoader(std::unique_ptr<FragmentLoader> &&loader, uint32_t cell_offset, uint32_t chunk_size);
        // Copy fragments to the destinations (must have space for at least chunk_size elements)
        // Return how many fragments were actually copied
        uint32_t load_chunk(uint32_t *start, uint32_t *end, uint32_t *cell);
        
        // Return the next start coordinate without marking it loaded, or UINT32_MAX if no more available
        uint32_t peek_start(); 

        // Wrapper methods for internal fragments
        void seek(uint32_t chr_id, uint32_t base);
        void restart();

        int chrCount() const;
        int cellCount() const;

        const char *chrNames(uint32_t chr_id);
        const char *cellNames(uint32_t cell_id);

        bool nextChr();

        uint32_t currentChr() const;
    };

    uint32_t const load_size, chunk_size, chunks_per_input;

    std::vector<ChunkedLoader> frags;
    std::vector<uint32_t> cell_id_offset;

    std::vector<std::string> chr_order; // Order of output chromosomes
    std::vector<std::vector<uint32_t>>
        source_chr; // source_chr[i][j] is the chr ID in frags[i] with name chr_order[j]

    // Size these to 3 * chunk_size * frags.size() 
    std::vector<uint32_t> start, end, cell, start_buf, end_buf, cell_buf;

    // Heap of (start_coord, loader_idx) pairs to use in load() operation
    std::vector<std::pair<uint32_t, uint32_t>> heap;

    // Loaded is how many values are saved in the data buffers
    // Available is how many of those are ready to output
    // cur_idx is the current index that we are outputting
    uint32_t loaded, available, cur_idx; 
    uint32_t current_chr = UINT32_MAX;

  public:
    MergeFragments(
        std::vector<std::unique_ptr<FragmentLoader>> &&fragments,
        const std::vector<std::string> &chr_order,
        uint32_t load_size = 1024, // Output load size
        uint32_t chunk_size = 32 // Input load chunk size
    );

    ~MergeFragments() = default;

    bool isSeekable() const override;
    void seek(uint32_t chr_id, uint32_t base) override;
    void restart() override;

    int chrCount() const override;
    int cellCount() const override;

    const char *chrNames(uint32_t chr_id) override;
    const char *cellNames(uint32_t cell_id) override;

    bool nextChr() override;

    uint32_t currentChr() const override;

    // Load algorithm:
    // 1. Add chunks_per_input * N chunks from the N input loaders into the cell, start, end vectors
    //    - Always add the chunk with the smallest start coordinate
    //    - Track what's the largest start coordinate added from any input loader.
    //      The smallest of these coordinates well be the limit of our known-sorted data.
    // 2. Sort all the fragments
    // 3. Output fragments until we hit our start coordinate limit
    //    - This is guaranteed to leave at most N chunks left over (at most 1 chunk for each input loader)
    //    - This means if we have space for (chunks_per_input + 1)*N chunks, and add (chunks_per_input)*N chunks each time we can never overflow 
    //      our buffers, and we'll never sort a fragment more than once
    bool load() override;

    uint32_t capacity() const override;

    uint32_t *cellData() override;
    uint32_t *startData() override;
    uint32_t *endData() override;
};

} // end namespace BPCells