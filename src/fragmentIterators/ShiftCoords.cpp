#include "ShiftCoords.h"

namespace BPCells {

ShiftCoords::ShiftCoords(FragmentLoader &loader, int32_t shift_start, int32_t shift_end)
    : FragmentLoaderWrapper(loader)
    , shift_start(shift_start)
    , shift_end(shift_end) {}

bool ShiftCoords::load() {
    if (!loader.load()) return false;
    uint32_t *start = loader.startData();
    uint32_t *end = loader.endData();
    uint32_t capacity = loader.capacity();

    vec start_vec = splat(shift_start);
    vec end_vec = splat(shift_end);
    // Use SIMD ops to shift values 4 at a time, followed by
    // a cleanup loop
    uint32_t i;
    for (i = 0; i + 4 <= capacity; i += 4) {
        vec in_start = BPCells::load((vec *)&start[i]);
        store((vec *)&start[i], add(start_vec, in_start));

        vec in_end = BPCells::load((vec *)&end[i]);
        store((vec *)&end[i], add(end_vec, in_end));
    }
    for (; i < capacity; i++) {
        start[i] += shift_start;
        end[i] += shift_end;
    }
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
    loader.seek(chr_id, base2);
}

} // end namespace BPCells
