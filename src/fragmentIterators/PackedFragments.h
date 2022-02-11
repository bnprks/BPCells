#pragma once

#include <algorithm>

#include "../bitpacking/packing_utils.h"
#include "FragmentsIterator.h"

namespace BPCells {

class PackedFragments : public FragmentsLoader {
public:
    // Read a fragment TSV with columns chr, start, end, cell_id [optional others] in that order.
    // cell and chromosme IDs are assigned in sequential order from the order they're seen
    PackedFragments(const std::vector<packed_frags> frags, 
        const std::vector<std::string> cell_names, const std::vector<std::string> chr_names);
    
    PackedFragments() = delete;

    bool isSeekable() const override;
    
    void seek(uint32_t chr_id, uint32_t base) override;

    void restart() override;

    int chrCount() const override;
    int cellCount() const override;

    const char* chrNames(uint32_t chr_id) const override;
    const char* cellNames(uint32_t cell_id) const override;
    
    // Advance the iterator to the next chromosome. Return false if there are no more chromosomes
    bool nextChr() override;
    // Return chromosome ID of current fragment
    uint32_t currentChr() const override;
private:
    const std::vector<packed_frags> frags;
    const std::vector<std::string> chr_names, cell_names;

    int32_t current_chr;
    uint32_t current_block;

    int32_t load(uint32_t count, FragmentArray buffer) override;
};



class PackedFragmentsWriter : public FragmentsWriter {
public:
    PackedFragmentsWriter() = default;
    PackedFragmentsWriter(const PackedFragmentsWriter&) = delete;
    PackedFragmentsWriter& operator=(const PackedFragmentsWriter& other) = delete;

    bool write(FragmentsIterator &fragments, void (*checkInterrupt)(void) = NULL) override;

    std::vector<packed_frags> getPackedFrags();
    std::vector<std::string> getChrNames();
    std::vector<std::string> getCellNames();
private:
    struct vec_packed_frags {
        std::vector<uint32_t> 
            start_data, start_idx, start_start,
            end_data, end_idx, end_max,
            cell_data, cell_idx;
        size_t len;
    };
    std::vector<std::string> chr_names, cell_names;
    std::vector<struct vec_packed_frags> frags;
    std::vector<packed_frags> packed_frags_vec;

    vec_uint32_t convert_to_vec(std::vector<uint32_t> &v);
};

} // end namespace BPCells