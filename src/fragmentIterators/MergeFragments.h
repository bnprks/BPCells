#pragma once

#include "../lib/dary_heap.hpp"
#include <cstring>
#include <vector>
#include <unordered_map>

#include "FragmentIterator.h"

namespace BPCells {

// Merge several fragment sources into a single sorted loader.
// Cell IDs will be adjusted to be sequential between fragments
// Chromosomes will be merged in the order of the provided chromosome ordering. It
// is an error if any input has chromosmes not in the output ordering.
// All inputs must have known cell counts, chr counts+names, and be seekable
class MergeFragments : public FragmentLoader {
  private:
    std::vector<FragmentIterator> frags;
    std::vector<uint32_t> heap; // Heap of indices into the fragments array
    std::vector<uint32_t> cell_id_offset;

    std::vector<std::string> chr_order; // Order of output chromosomes
    std::vector<std::vector<uint32_t>>
        source_chr; // source_chr[i][j] is the chr ID in frags[i] with name chr_order[j]

    std::vector<uint32_t> start, end, cell;
    uint32_t loaded;
    uint32_t current_chr = UINT32_MAX;

    // Order by start, breaking ties by the fragment loader index
    inline bool compare(uint32_t idx1, uint32_t idx2) {
        if (frags[idx1].start() != frags[idx2].start())
            return frags[idx1].start() > frags[idx2].start();
        return idx1 > idx2;
    }

  public:
    MergeFragments(
        std::vector<std::unique_ptr<FragmentLoader>> &&fragments, const std::vector<std::string> &chr_order, uint32_t load_size = 1024
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

    bool load() override;

    uint32_t capacity() const override;

    uint32_t *cellData() override;
    uint32_t *startData() override;
    uint32_t *endData() override;
};

} // end namespace BPCells