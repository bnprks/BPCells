#pragma once

#include <algorithm>
#include <functional>

#include "FragmentsIterator.h"
#include "../lib/dary_heap.hpp"
// This is a record of insertionsIterator where the heap is starts+ends, and 
// gets added to incrementally during nextInsertion rather than at load time
namespace BPCells {

// Class to conveniently iterate over insertion sites in sorted order.
// Input must be a *sorted* fragments iterator
class InsertionsIterator : public FragmentsLoaderWrapper {
private:
    const uint32_t chunk_capacity;
    int32_t chunk_size;
    uint32_t last_start;
    uint32_t current_chr;
    uint32_t idx;
    std::vector<uint32_t> start_buf, end_buf, cell_buf;
    FragmentArray fragments_buf;
    uint32_t last_output = 0;

    class insertion {
    private:
        uint32_t coord_;
        uint32_t cell_id_;
    public:
        insertion(uint32_t coord, uint32_t cell_id):
            coord_(coord), cell_id_(cell_id) {}
        inline uint32_t coord() const {return coord_;}
        inline uint32_t cell_id() const {return cell_id_;}
        friend bool operator>(const insertion &i1, const insertion &i2) {
            return i1.coord_ > i2.coord_;
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
        //if (heap.size() > 0) printf("Heap min start: %d", heap.front().coord);
        for (int i = 0; i < chunk_size; i++) {
            if (fragments_buf.start[i] < last_start)
                throw std::runtime_error("Input fragments not in sorted order by start");
            if (fragments_buf.start[i] >= fragments_buf.end[i])
                throw std::runtime_error("Input fragments have end <= start");
            last_start = fragments_buf.start[i];
            heap.push_back({fragments_buf.start[i], fragments_buf.cell[i]});
            dary_heap::push_heap<2>(heap.begin(), heap.end(), std::greater<insertion>());

            heap.push_back({fragments_buf.end[i] - 1, fragments_buf.cell[i]});
            dary_heap::push_heap<2>(heap.begin(), heap.end(), std::greater<insertion>());
        }
        printf("Loaded fragments; new heap size: %ld, new heap min: %d\n", heap.size(), heap.front().coord());
        return true;
    }
public:
    // Construct iterator with a given internal buffer size (must be a power of 2 >= 128)
    InsertionsIterator(FragmentsLoader &loader, uint32_t buffer_size = 1024);

    virtual ~InsertionsIterator() = default;
    
    // Return false if there isn't a nextFragment in the current chromosome
    inline bool nextInsertion() {
        if (!heap.empty()) heap.pop_back();
        while (heap.empty() || heap[0].coord() > last_start) {
            idx += 1;
            if (idx >= chunk_size) {
                chunk_size = loader.load(chunk_capacity, fragments_buf);
                idx = 0;
                if (chunk_size == 0) {
                    if (heap.empty()) return false;
                    else break;
                }
            }
            if (fragments_buf.start[idx] < last_start)
                throw std::runtime_error("Input fragments not in sorted order by start");
            if (fragments_buf.start[idx] >= fragments_buf.end[idx])
                throw std::runtime_error("Input fragments have end <= start");
            last_start = fragments_buf.start[idx];
            heap.push_back({fragments_buf.start[idx], fragments_buf.cell[idx]});
            dary_heap::push_heap<2>(heap.begin(), heap.end(), std::greater<insertion>());

            heap.push_back({fragments_buf.end[idx] - 1, fragments_buf.cell[idx]});
            dary_heap::push_heap<2>(heap.begin(), heap.end(), std::greater<insertion>());
        }
        dary_heap::pop_heap<2>(heap.begin(), heap.end(), std::greater<insertion>());
        if (heap.back().coord() < last_output) {
            throw std::runtime_error("logic error in InsertionsIteator.h");
        }
        last_output = heap.back().coord();
        return true;
    }
    // Advance the iterator to the next chromosome. Return false if there are no more chromosomes
    inline bool nextChr() override {
        bool res = loader.nextChr();
        if (res) current_chr = loader.currentChr();
        heap.clear();
        last_start = 0;
        idx = chunk_capacity;
        last_output = 0;
        return res;
    }
    // Access chr, start, end, cell from current fragment
    inline uint32_t chr() const {return current_chr; };
    inline uint32_t coord() const {return heap.back().coord(); };
    inline uint32_t cell() const {return heap.back().cell_id(); };   
    
    inline uint32_t get_chunk_capacity() {return chunk_capacity;}; 

    // Move iterator to just before fragments which end after "base".
    // It's possible that fragments returned after seek will end before "base",
    // but it's guaranteed that it won't skip over any fragments ending before "base"
    inline void seek(uint32_t chr_id, uint32_t base) override {
        loader.seek(chr_id, base);
        heap.clear();
        last_start = 0;
        last_output = 0;
        idx = chunk_capacity;
    }   

    int32_t load(uint32_t count, FragmentArray &buffer) override;
};

} // end namespace BPCells