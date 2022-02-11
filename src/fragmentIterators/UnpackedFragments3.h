#pragma once

#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>

#include "FragmentsIterator.h"
#include "../arrayIO/array_interfaces.h"

namespace BPCells {

class UnpackedFragments3: public FragmentsLoader {
private:
    using ReaderPtr = std::unique_ptr<UIntReader>;
    // cell, start, end = concatenated fragment data across chromosomes
    // end_max[i] = max(ends[chr_start:i*128])
    // chr_ptr[2*i] = start index of chr i in data arrays
    // chr_ptr[2*i + 1] = end index of chr i in data arrays
    ReaderPtr cell, start, end, end_max, chr_ptr;
    std::unique_ptr<StringReader> chr_names, cell_names;

    // end_max_buf holds the end_max values for the current chromosome, but to 
    // deal with possible overhang from high end values in the previous chromosomes we 
    // ensure end_max_buf[0] <= end_max_buf[1]
    std::vector<uint32_t> end_max_buf;
    uint32_t current_chr = UINT32_MAX;
    uint32_t current_idx = UINT32_MAX;
    uint32_t chr_start_ptr, chr_end_ptr;

    // Read end_max_buf from end_max iterator, making it equal to the values between
    // start_idx and end_idx
    void readEndMaxBuf(uint32_t start_idx, uint32_t end_idx);
public:
    UnpackedFragments3(ReaderPtr &&cell, ReaderPtr &&start, ReaderPtr &&end, 
        ReaderPtr &&end_max, ReaderPtr && chr_ptr, 
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

class UnpackedFragmentsWriter3: public FragmentsWriter {
private: 
    using WriterPtr = std::unique_ptr<UIntWriter>;
    WriterPtr cell, start, end, end_max, chr_ptr;
    std::unique_ptr<StringWriter> chr_names, cell_names;
public:
    UnpackedFragmentsWriter3(WriterPtr &&cell, WriterPtr &&start, WriterPtr &&end, 
        WriterPtr &&end_max, WriterPtr && chr_ptr,
        std::unique_ptr<StringWriter> &&chr_names, std::unique_ptr<StringWriter> &&cell_names);
    
    bool write(FragmentsIterator &fragments, void (*checkInterrupt)(void) = NULL) override;
};

} // end namespace BPCells