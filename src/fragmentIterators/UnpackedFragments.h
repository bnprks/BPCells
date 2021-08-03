#pragma once

#include <algorithm>

#include "../bitpacking/packing_utils.h"
#include "FragmentsIterator.h"

namespace BPCells {

class UnpackedFragments : public FragmentsLoader {
public:
    UnpackedFragments(
        const std::vector<vec_uint32_t> start, const std::vector<vec_uint32_t> end,
        const std::vector<vec_uint32_t> cell_id,  
        const std::vector<std::string> cell_names, const std::vector<std::string> chr_names);
    
    UnpackedFragments() = delete;

    bool isSeekable() const override;
    // Move iterator to just before fragments which end after "base".
    // It's possible that fragments returned after seek will end before "base",
    // but it's guaranteed that it won't skip over any fragments ending before "base"
    void seek(uint32_t chr_id, uint32_t base) override;

    // Reset the iterator to start from the beginning
    void restart() override;

    // Return the number of cells/chromosomes, or return -1 if this number is 
    // not known ahead of time
    int chrCount() const override;
    int cellCount() const override;

    // Return name for a given chr_id or cell_id. Only valid to call
    // for chromosme or cell_ids that have been actually returned by the iterator
    const char* chrNames(uint32_t chr_id) const override;
    const char* cellNames(uint32_t cell_id) const override;
    
    // Advance the iterator to the next chromosome. Return false if there are no more chromosomes
    bool nextChr() override;

    // Return chromosome ID of current fragment
    uint32_t currentChr() const override;

    // Return number of items loaded. Should repeatedly return 0 at the end of a chromosome.
    // Return -1 for error
    int32_t load(uint32_t count, FragmentArray &buffer) override;
private:
    const std::vector<vec_uint32_t> start;
    const std::vector<vec_uint32_t> end;
    const std::vector<vec_uint32_t> cell_id;
    
    const std::vector<std::string> chr_names, cell_names;

    int32_t current_chr;
    uint32_t current_fragment;
};


class UnpackedFragmentsWriter : public FragmentsWriter {
public:
    UnpackedFragmentsWriter() = default;
    UnpackedFragmentsWriter(const UnpackedFragmentsWriter&) = delete;
    UnpackedFragmentsWriter& operator=(const UnpackedFragmentsWriter& other) = delete;

    bool write(FragmentsIterator &fragments, void (*checkInterrupt)(void) = NULL) override;

    struct unpacked_frags {
        std::vector<uint32_t> start_data, end_data, cell_data;
        size_t len;
    };
    const std::vector<UnpackedFragmentsWriter::unpacked_frags> *getFrags();
    std::vector<std::string> getChrNames();
    std::vector<std::string> getCellNames();


private:
    std::vector<std::string> chr_names, cell_names;
    std::vector<struct unpacked_frags> frags;
};



} // end namespace BPCells