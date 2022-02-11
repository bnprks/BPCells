#pragma once

#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>

#include "../bitpacking/bp128.h"
#include "FragmentsIterator.h"
#include "../arrayIO/array_interfaces.h"

namespace BPCells {

class PackedFragments3: public FragmentsLoader {
private:
    using ReaderPtr = std::unique_ptr<UIntReader>;
    // cell, start, end = concatenated fragment data across chromosomes
    // end_max[i] = max(ends[chr_start:i*128])
    // chr_ptr[2*i] = start index of chr i in data arrays
    // chr_ptr[2*i + 1] = end index of chr i in data arrays
    ReaderPtr cell_data, cell_idx, start_data, start_idx, start_starts, 
        end_data, end_idx, end_max, chr_ptr;
    std::unique_ptr<StringReader> chr_names, cell_names;

    // end_max_buf holds the end_max values for the current chromosome, but to 
    // deal with possible overhang from high end values in the previous chromosomes we 
    // ensure end_max_buf[0] <= end_max_buf[1]
    std::vector<uint32_t> end_max_buf;
    uint32_t current_chr = UINT32_MAX;
    uint32_t current_idx = UINT32_MAX;
    uint32_t chr_start_ptr, chr_end_ptr;

    uint32_t cell_buf[128], start_buf[128], end_buf[128];
    uint32_t prev_cell_idx, prev_start_idx, prev_end_idx;

    // Read end_max_buf from end_max iterator, making it equal to the values between
    // start_idx and end_idx
    void readEndMaxBuf(uint32_t start_idx, uint32_t end_idx);

    void load128(uint32_t *cell_out, uint32_t *start_out, uint32_t *end_out);
public:
    PackedFragments3(ReaderPtr &&cell_data, ReaderPtr &&cell_idx, ReaderPtr &&start_data, ReaderPtr &&start_idx, ReaderPtr &&start_starts, 
        ReaderPtr &&end_data, ReaderPtr &&end_idx, ReaderPtr &&end_max, ReaderPtr &&chr_ptr, 
        std::unique_ptr<StringReader> &&chr_names,
        std::unique_ptr<StringReader> &&cell_names);
                
    bool isSeekable() const override;
    void seek(uint32_t chr_id, uint32_t base) override;

    void restart() override;

    int chrCount() const override;
    int cellCount() const override;
    const char* chrNames(uint32_t chr_id) const override;
    const char* cellNames(uint32_t cell_id) const override;

    bool nextChr() override;
    uint32_t currentChr() const override;
    
    int32_t load(uint32_t count, FragmentArray buffer) override;
};

class PackedFragmentsWriter3: public FragmentsWriter {
private: 
    using WriterPtr = std::unique_ptr<UIntWriter>;
    WriterPtr cell_data, cell_idx, start_data, start_idx, start_starts, 
        end_data, end_idx, end_max, chr_ptr;
    std::unique_ptr<StringWriter> chr_names, cell_names;

    void pack128(const uint32_t* cell_in, const uint32_t* start_in, uint32_t* end_in,
                uint32_t &cur_cell_idx, uint32_t &cur_start_idx, uint32_t &cur_end_idx);
public:
    PackedFragmentsWriter3(WriterPtr &&cell_data, WriterPtr &&cell_idx, WriterPtr &&start_data, WriterPtr &&start_idx, WriterPtr &&start_starts, 
        WriterPtr &&end_data, WriterPtr &&end_idx, WriterPtr &&end_max, WriterPtr &&chr_ptr, 
        std::unique_ptr<StringWriter> &&chr_names,
        std::unique_ptr<StringWriter> &&cell_names);
    
    bool write(FragmentsIterator &fragments, void (*checkInterrupt)(void) = NULL) override;
};

} // end namespace BPCells