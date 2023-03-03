#pragma once

#include "../bitpacking/simd_vec.h"
#include "FragmentIterator.h"
namespace BPCells {

// Shift start & end coordinates of a FragmentsLoader, and expose result
// as a new loader
class ShiftCoords : public FragmentLoaderWrapper {
  private:
    int32_t shift_start, shift_end;

  public:
    ShiftCoords(std::unique_ptr<FragmentLoader> &&loader, int32_t shift_start, int32_t shift_end);

    ~ShiftCoords() = default;

    // Move loader to just before fragments which end after "base".
    // It's possible that fragments returned after seek will end before "base",
    // but it's guaranteed that it won't skip over any fragments ending before "base"
    void seek(uint32_t chr_id, uint32_t base) override;
    bool load() override;
};

} // end namespace BPCells