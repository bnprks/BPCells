#pragma once

#include <algorithm>
#include <functional>

#include "FragmentsIterator.h"
#include "InsertionsIterator.h"
#include "../lib/dary_heap.hpp"

namespace BPCells {

// Class to conveniently iterate over insertion sites in sorted order.
// Input must be a *sorted* fragments iterator
class InsertionsIterator2 : public FragmentsLoaderWrapper {
private:
    uint32_t load_size = 1024;
    uint32_t current_chr;

    uint32_t next_start = UINT32_MAX; 
    uint32_t next_end = UINT32_MAX;

    uint32_t current_coord;
    uint32_t current_cell;

    struct insertion_buffer {
        std::vector<uint32_t> coord, cell;
        void resize(size_t s) {coord.resize(s); cell.resize(s);}
        size_t size() {return coord.size();}
    };
    insertion_buffer start_buf, end_buf, end_sort_scratch;
    
    bool finished_chr = false;
    
    // Use radix sort (LSD) to sort the end buffer
    inline void sortEnds();
    
    // Load new insertions into start_buf and end_buf
    inline void loadInsertions();
    
    inline void reset_internals() {
        next_start = UINT32_MAX;
        next_end = end_buf.size();
        finished_chr = false;
    }
public:
    // Construct iterator with a given internal buffer size (must be a power of 2 >= 128)
    InsertionsIterator2(FragmentsLoader &loader) : FragmentsLoaderWrapper(loader) {
        start_buf.resize(load_size);
        end_buf.resize(load_size);
        end_sort_scratch.resize(load_size);
        reset_internals();
    }

    virtual ~InsertionsIterator2() = default;

    // Return false if there isn't a nextFragment in the current chromosome
    inline bool nextInsertion() {
        if (next_start >= start_buf.size()) {
            if (!finished_chr) loadInsertions();
            if (finished_chr) {
                if (next_end >= end_buf.size()) return false;
                current_coord = end_buf.coord[next_end] - 1;
                current_cell = end_buf.cell[next_end];
                next_end += 1;
                return true;
            }
        }
        
        // TODO: consider making start_buf padded by UINT32_MAX at the end, so I 
        // don't need to compare against finished_chr
        bool next_is_start = start_buf.coord[next_start] < end_buf.coord[next_end];
        current_coord = next_is_start ? start_buf.coord[next_start] : end_buf.coord[next_end] - 1;
        current_cell = next_is_start ? start_buf.cell[next_start] : end_buf.cell[next_end]; 
        
        next_start += next_is_start;
        next_end += !next_is_start;

        return true;
    }
    // Advance the iterator to the next chromosome. Return false if there are no more chromosomes
    inline bool nextChr() override {
        bool res = loader.nextChr();
        if (res) current_chr = loader.currentChr();
        reset_internals();
        return res;
    }
    // Access chr, start, end, cell from current fragment
    inline uint32_t chr() const {return current_chr; };
    inline uint32_t coord() const {return current_coord; };
    inline uint32_t cell() const {return current_cell; };
    
    // Move iterator to just before fragments which end after "base".
    // It's possible that fragments returned after seek will end before "base",
    // but it's guaranteed that it won't skip over any fragments ending before "base"
    inline void seek(uint32_t chr_id, uint32_t base) override {
        loader.seek(chr_id, base);
        current_chr = chr_id;
        reset_internals();
    }   

    inline void restart() override {
        loader.restart();
        reset_internals();
    }

    int32_t load(uint32_t count, FragmentArray &buffer) override {
        return loader.load(count, buffer);
    };
};


void InsertionsIterator2::sortEnds() {
    uint32_t radix_counts[4][256] = {{0}};
    bool skip_byte[4] = {false};
    // Count up how many we see of each byte combination in 1 pass of the data
    for (uint32_t coord : end_buf.coord) {
        for (int i = 0; i < 4; i++) {
            radix_counts[i][255 & ((uint32_t) coord >> (i * 8))]++;
        }
    }
    // Tally up which output index each byte combination should start at
    // So radix_counts[i][j] will turn into the sum of all radix_counts[i][k] where k < j
    for (int i = 0; i < 4; i++) {
        uint32_t running_sum = 0;
        for (int j = 0; j < 256; j++) {
            if (radix_counts[i][j] == end_buf.size()) skip_byte[i] = true;
            running_sum += radix_counts[i][j];
            radix_counts[i][j] = running_sum - radix_counts[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        // Skip if this byte is all identical
        if (skip_byte[i]) { continue; }
        for (int j = 0; j < end_buf.size(); j++) {
            uint32_t bucket = 255 & ((uint32_t) end_buf.coord[j] >> (i * 8));
            end_sort_scratch.coord[radix_counts[i][bucket]] = end_buf.coord[j];
            end_sort_scratch.cell[radix_counts[i][bucket]] = end_buf.cell[j];
            radix_counts[i][bucket]++;
        }
        std::swap(end_sort_scratch, end_buf);
    }
}

void InsertionsIterator2::loadInsertions() {
    // Increase buffer size if ends remaining is a large fraction of buffer size
    // This avoids having to re-sort data too frequently
    uint32_t ends_remaining = end_buf.size() - next_end;
    if (ends_remaining * 2 > load_size) {
        load_size = load_size * 2;
        printf("Expanding load size: %d\n", load_size);
        end_buf.resize(load_size);
        end_sort_scratch.resize(load_size);
        start_buf.resize(load_size);
    }
    
    // Read in data to fill remaining buffer
    uint32_t fragments_to_load = (load_size - ends_remaining);
    
    start_buf.resize(fragments_to_load);
    end_sort_scratch.resize(ends_remaining + fragments_to_load);
    
    FragmentArray buf = {
        &start_buf.coord[0], &end_sort_scratch.coord[ends_remaining], &start_buf.cell[0], fragments_to_load
    };
    uint32_t load_count = loader.load(fragments_to_load, buf);

    if (load_count == 0) {
        finished_chr = true;
        next_start = UINT32_MAX;
    } else {
        end_buf.resize(ends_remaining + load_count);
        end_sort_scratch.resize(ends_remaining + load_count);
        start_buf.resize(load_count);

        // Copy remaing ends into buffer
        std::memmove(&end_sort_scratch.cell[0], &end_buf.cell[next_end], ends_remaining*sizeof(uint32_t));
        std::memmove(&end_sort_scratch.coord[0], &end_buf.coord[next_end], ends_remaining*sizeof(uint32_t));
        // Copy loaded cells into buffer
        std::memmove(&end_sort_scratch.cell[ends_remaining], &start_buf.cell[0], load_count * sizeof(uint32_t));

        std::swap(end_sort_scratch, end_buf);

        sortEnds();

        next_end = 0;
        next_start = 0;
    }
}

} // end namespace BPCells