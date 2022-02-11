#include "ShiftCoords.h"

namespace BPCells {

ShiftCoords::ShiftCoords(FragmentsLoader &loader, int32_t shift_start, int32_t shift_end) :
    FragmentsLoaderWrapper(loader), 
    shift_start(shift_start), shift_end(shift_end)
        {} 



int32_t ShiftCoords::load(uint32_t count, FragmentArray buffer) {
    int32_t res = loader.load(count, buffer);
    if (res <= 0) return res;
    vec start = splat(shift_start);
    vec end = splat(shift_end);
    // Use SIMD ops to shift values 4 at a time, followed by
    // a cleanup loop
    int i;
    for (i = 0; i + 4 < res; i += 4) {
        vec in_start = BPCells::load((vec *) &buffer.start[i]);
        BPCells::store((vec *) &buffer.start[i], add(start, in_start));

        vec in_end = BPCells::load((vec *) &buffer.end[i]);
        BPCells::store((vec *) &buffer.end[i], add(end, in_end));
    }
    for(; i < res; i++) {
        buffer.start[i] += shift_start;
        buffer.end[i] += shift_end;
    }
    return res;
};

// Move loader to just before fragments which end after "base".
// It's possible that fragments returned after seek will end before "base",
// but it's guaranteed that it won't skip over any fragments ending before "base"
void ShiftCoords::seek(uint32_t chr_id, uint32_t base) {
    int32_t m = std::min(shift_start, shift_end);
    uint32_t base2;
    // Handle shifts causing underflow
    if (m < 0) base2 = std::min(base, base + m);
    else base2 =  base + m;
    loader.seek(chr_id, base2);
};

} // end namespace BPCells