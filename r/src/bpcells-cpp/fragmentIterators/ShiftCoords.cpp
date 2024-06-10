// Copyright 2021 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "ShiftCoords.h"
#include "../simd/math.h"

namespace BPCells {

ShiftCoords::ShiftCoords(
    std::unique_ptr<FragmentLoader> &&loader, int32_t shift_start, int32_t shift_end
)
    : FragmentLoaderWrapper(std::move(loader))
    , shift_start(shift_start)
    , shift_end(shift_end) {}

bool ShiftCoords::load() {
    if (!loader->load()) return false;
    uint32_t *start = loader->startData();
    uint32_t *end = loader->endData();
    uint32_t capacity = loader->capacity();

    simd::add(start, shift_start, capacity);
    simd::add(end, shift_end, capacity);
    return true;
}

// Move loader to just before fragments which end after "base".
// It's possible that fragments returned after seek will end before "base",
// but it's guaranteed that it won't skip over any fragments ending before "base"
void ShiftCoords::seek(uint32_t chr_id, uint32_t base) {
    int32_t m = std::min(shift_start, shift_end);
    uint32_t base2;
    // Handle shifts causing underflow
    if (m < 0) base2 = std::min(base, base + m);
    else base2 = base + m;
    loader->seek(chr_id, base2);
}

} // end namespace BPCells
