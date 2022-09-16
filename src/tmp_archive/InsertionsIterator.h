#pragma once

#include <algorithm>
#include <cassert>
#include <functional>

#include "../lib/dary_heap.hpp"
#include "FragmentIterator.h"

namespace BPCells {
class insertion {
  private:
    uint32_t coord_;
    uint32_t cell_id_;

  public:
    insertion() : coord_(UINT32_MAX), cell_id_(UINT32_MAX) {}
    insertion(uint32_t coord, uint32_t cell_id) : coord_(coord), cell_id_(cell_id) {}
    inline uint32_t coord() const { return coord_; }
    inline uint32_t cell_id() const { return cell_id_; }
    friend bool operator<(const insertion &i1, const insertion &i2) {
        return i1.coord_ < i2.coord_;
    }
    friend bool operator>(const insertion &i1, const insertion &i2) {
        return i1.coord_ > i2.coord_;
    }
};
// A buffer that uses a heap to sort end coordinates
class HeapBuffer {
  private:
    std::vector<insertion> heap;

  public:
    HeapBuffer() { heap.push_back({UINT32_MAX, UINT32_MAX}); }
    inline void add_coord(const uint32_t coord, const uint32_t cell) {
        heap.push_back({coord, cell});
        dary_heap::push_heap<4>(heap.begin(), heap.end(), std::greater<insertion>());
    }
    inline uint32_t minCoord() const { return heap.front().coord(); }
    inline uint32_t minCell() const { return heap.front().cell_id(); }
    inline void popMin() {
        dary_heap::pop_heap<4>(heap.begin(), heap.end(), std::greater<insertion>());
        heap.pop_back();
    }
    inline bool empty() const { return heap.front().coord() == UINT32_MAX; }
    inline void clear() {
        heap.resize(0);
        heap.push_back({UINT32_MAX, UINT32_MAX});
    }
};

// A buffer that uses a circular
class CircleBuffer {
  private:
    std::vector<insertion> buf;
    // Properties:
    // end = next location to insert an item
    // end == start -> empty buffer
    // end == start-1 -> full buffer
    size_t start, end;

    inline uint32_t nextIdx(uint32_t idx) const {
        return idx + (idx < buf.size() - 1 ? 1 : -buf.size() + 1);
    }
    inline uint32_t prevIdx(uint32_t idx) const { return idx - (idx > 0 ? 1 : -buf.size() + 1); }
    inline void resize() {
        // printf("Calling resize, size=%zu, start=%zu, end=%zu\n", buf.size(), start, end);
        assert(end == prevIdx(start));
        end = buf.size() - 1;
        buf.resize(buf.size() + buf.size() / 2 + 1);
        // printf("new buf size: %zu\n", buf.size());
        for (size_t i = 0; i < start; i++) {
            end = nextIdx(end);
            buf[end] = buf[i];
        }
        // printf("Finished resize, size=%zu, start=%zu, end=%zu\n", buf.size(), start, end);
    }

  public:
    CircleBuffer() {
        start = 0;
        end = 1;
        buf.resize(2);
        buf[0] = {UINT32_MAX, UINT32_MAX};
    }
    inline void add_coord(const uint32_t coord, const uint32_t cell) {
        // printf("insertion-pre: start=%zu, end=%zu, array=", start, end);
        // for(int i = 0; i < buf.size(); i++)
        //     printf("%d ", buf[i].coord());
        // printf("\n");
        if (end == prevIdx(start)) resize();
        uint32_t cur = end;
        end = nextIdx(end);
        uint32_t neighbor = prevIdx(cur);
        buf[cur] = {coord, cell};
        while (buf[cur] < buf[neighbor] && cur != start) {
            std::swap(buf[cur], buf[neighbor]);
            cur = neighbor;
            neighbor = prevIdx(cur);
        }
        // printf("insertion-post: start=%zu, end=%zu, array=", start, end);
        // for(int i = 0; i < buf.size(); i++)
        //     printf("%d ", buf[i].coord());
        // printf("\n");
    }
    inline uint32_t minCoord() const { return buf[start].coord(); }
    inline uint32_t minCell() const { return buf[start].cell_id(); }
    inline void popMin() {
        // printf("Removing coord %d. new_min=%d\n", minCoord(), buf[nextIdx(start)].coord());
        start = nextIdx(start);
    }
    inline bool empty() const { return buf[start].coord() == UINT32_MAX; }
    inline void clear() {
        end = start;
        // printf("calling clear\n");
        add_coord(UINT32_MAX, UINT32_MAX);
    }
};

// Class to conveniently iterate over insertion sites in sorted order.
// Input must be a *sorted* fragments iterator
class InsertionsIterator : public FragmentsLoaderWrapper {
  private:
    const uint32_t chunk_capacity;
    int32_t chunk_size;
    uint32_t current_chr;
    uint32_t idx;
    std::vector<uint32_t> start_buf, end_buf, cell_buf;
    FragmentArray fragments_buf;

    insertion next_insertion = {UINT32_MAX, UINT32_MAX};

    HeapBuffer buf;
    // CircleBuffer buf;
  public:
    // Construct iterator with a given internal buffer size (must be a power of 2 >= 128)
    InsertionsIterator(FragmentsLoader &loader, uint32_t buffer_size = 1024);

    virtual ~InsertionsIterator() = default;

    // Return false if there isn't a nextFragment in the current chromosome
    inline bool nextInsertion() {
        if (idx >= chunk_size) {
            chunk_size = loader.load(chunk_capacity, fragments_buf);
            idx = 0;
            if (chunk_size == 0) {
                if (buf.empty()) return false;
                next_insertion = {buf.minCoord(), buf.minCell()};
                buf.popMin();
                return true;
            }
        }
        // printf("buf_min = %d, start_min=%d\n", buf.minCoord(), fragments_buf.start[idx]);
        if (buf.minCoord() < fragments_buf.start[idx]) {
            next_insertion = {buf.minCoord(), buf.minCell()};
            buf.popMin();
        } else {
            next_insertion = {fragments_buf.start[idx], fragments_buf.cell[idx]};
            buf.add_coord(fragments_buf.end[idx] - 1, fragments_buf.cell[idx]);
            idx++;
        }
        return true;
    }
    // Advance the iterator to the next chromosome. Return false if there are no more chromosomes
    inline bool nextChr() override {
        bool res = loader.nextChr();
        if (res) current_chr = loader.currentChr();
        buf.clear();
        idx = chunk_capacity;
        return res;
    }
    // Access chr, start, end, cell from current fragment
    inline uint32_t chr() const { return current_chr; };
    inline uint32_t coord() const { return next_insertion.coord(); };
    inline uint32_t cell() const { return next_insertion.cell_id(); };

    inline uint32_t get_chunk_capacity() { return chunk_capacity; };

    // Move iterator to just before fragments which end after "base".
    // It's possible that fragments returned after seek will end before "base",
    // but it's guaranteed that it won't skip over any fragments ending before "base"
    inline void seek(uint32_t chr_id, uint32_t base) override {
        loader.seek(chr_id, base);
        current_chr = chr_id;
        buf.clear();
        idx = chunk_capacity;
    }

    inline void restart() override {
        loader.restart();
        idx = chunk_capacity;
        buf.clear();
    }

    int32_t load(uint32_t count, FragmentArray buffer) override;
};

} // end namespace BPCells