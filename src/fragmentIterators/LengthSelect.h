#pragma once
#include "FragmentIterator.h"

namespace BPCells {

// Transform a fragments iterator by subsetting to fragments of a given length range
class LengthSelect : public FragmentLoaderWrapper {
private:
    uint32_t loaded = 0;
    uint32_t min_len, max_len;
public:
    // min_len and max_len provide inclusive limits on the size of fragments
    LengthSelect(FragmentLoader &loader, uint32_t min_len, uint32_t max_len=UINT32_MAX);

    ~LengthSelect() = default;

    bool load() override;
    uint32_t capacity() const override;
};


} // end namespace BPCells