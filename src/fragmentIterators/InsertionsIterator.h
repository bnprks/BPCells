#pragma once

#include <algorithm>
#include <functional>

#include "FragmentsIterator.h"
namespace BPCells {

// Class to conveniently iterate over insertion sites in sorted order.
// Input must be a *sorted* fragments iterator
class InsertionsIterator : public FragmentsLoaderWrapper {
private:
    const uint32_t chunk_capacity;
    int32_t chunk_size;
    uint32_t last_start;
    uint32_t current_chr;
    std::vector<uint32_t> start_buf, end_buf, cell_buf;
    FragmentArray fragments_buf;

    struct insertion {
        uint32_t coord;
        uint32_t cell_id;
        friend bool operator>(const insertion &i1, const insertion &i2) {
            return i1.coord > i2.coord;
        }
    };

    // State of the heap:
    // - At start of chromosome, heap is empty
    // - While reading, last element of heap is the next insertion in order,
    //   as it was just put there by heap_pop
    std::vector<insertion> heap;

    // Load more fragments into the heap
    inline bool loadFragments() {
        chunk_size = loader.load(chunk_capacity, fragments_buf);
        if (chunk_size == 0) return false;
        for (int i = 0; i < chunk_size; i++) {
            if (fragments_buf.start[i] < last_start)
                throw std::runtime_error("Input fragments not in sorted order by start");
            if (fragments_buf.start[i] >= fragments_buf.end[i])
                throw std::runtime_error("Input fragments have end <= start");
            last_start = fragments_buf.start[i];
            heap.push_back({fragments_buf.start[i], fragments_buf.cell[i]});
            std::push_heap(heap.begin(), heap.end(), std::greater<insertion>());

            heap.push_back({fragments_buf.end[i] - 1, fragments_buf.cell[i]});
            std::push_heap(heap.begin(), heap.end(), std::greater<insertion>());
        }
        return true;
    }
public:
    // Construct iterator with a given internal buffer size (must be a power of 2 >= 128)
    InsertionsIterator(FragmentsLoader &loader, uint32_t buffer_size = 1024);

    virtual ~InsertionsIterator() = default;
    
    // Return false if there isn't a nextFragment in the current chromosome
    inline bool nextInsertion() {
        if (!heap.empty()) heap.pop_back();
        while (heap.empty() || heap[0].coord < last_start) {
            if (!loadFragments()) {
                // This ensures we will just pop all the remaining items
                last_start = 0; 
                if (heap.empty()) return false;
            }
        }
        std::pop_heap(heap.begin(), heap.end(), std::greater<insertion>());
        return true;
    }
    // Advance the iterator to the next chromosome. Return false if there are no more chromosomes
    inline bool nextChr() override {
        bool res = loader.nextChr();
        if (res) current_chr = loader.currentChr();
        heap.clear();
        last_start = 0;
        return res;
    }
    // Access chr, start, end, cell from current fragment
    inline uint32_t chr() const {return current_chr; };
    inline uint32_t coord() const {return heap.back().coord; };
    inline uint32_t cell() const {return heap.back().cell_id; };   
    
    inline uint32_t get_chunk_capacity() {return chunk_capacity;}; 

    // Move iterator to just before fragments which end after "base".
    // It's possible that fragments returned after seek will end before "base",
    // but it's guaranteed that it won't skip over any fragments ending before "base"
    inline void seek(uint32_t chr_id, uint32_t base) override {
        loader.seek(chr_id, base);
        heap.clear();
        last_start = 0;
    }   

    int32_t load(uint32_t count, FragmentArray &buffer) override;
};

} // end namespace BPCells