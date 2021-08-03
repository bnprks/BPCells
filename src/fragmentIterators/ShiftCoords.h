#pragma once

#include "FragmentsIterator.h"
#include "../bitpacking/simd_vec.h"
namespace BPCells {

// Shift start & end coordinates of a FragmentsLoader, and expose result
// as a new loader
class ShiftCoords : public FragmentsLoaderWrapper{
private:
    int32_t shift_start, shift_end;
public:
    
    ShiftCoords(FragmentsLoader &loader, int32_t shift_start, int32_t shift_end);

    ~ShiftCoords() = default;
    
    int32_t load(uint32_t count, FragmentArray &buffer) override;

    // Move loader to just before fragments which end after "base".
    // It's possible that fragments returned after seek will end before "base",
    // but it's guaranteed that it won't skip over any fragments ending before "base"
    void seek(uint32_t chr_id, uint32_t base) override;
};


} // end namespace BPCells