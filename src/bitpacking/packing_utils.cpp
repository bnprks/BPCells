#include "packing_utils.h"

namespace BPCells {


void unpack_128(const bp128_array &src, uint32_t *dst, const uint32_t idx) {
    uint32_t bits = (src.idx[idx+1] - src.idx[idx]) / 4;
    simdunpack(&src.data[src.idx[idx]], dst, bits);
};

void unpack_128_d1(const bp128_d1_array &src, uint32_t *dst, const uint32_t idx) {
    uint32_t bits = (src.idx[idx+1] - src.idx[idx]) / 4;
    simdunpackd1(src.starts[idx], &src.data[src.idx[idx]], dst, bits);
};

void unpack_128_for(const bp128_for_array &src, uint32_t *dst, const uint32_t idx) {
    uint32_t bits = (src.idx[idx+1] - src.idx[idx]) / 4;
    simdunpackFOR(src.mins[idx], &src.data[src.idx[idx]], dst, bits);
};

void unpack_128_frags(const packed_frags &src, uint32_t *starts, uint32_t *ends, uint32_t* cell_ids, const uint32_t idx) {
    unpack_128_d1(src.starts, starts, idx);
    unpack_128(src.ends, ends, idx);
    simdadd(starts, ends, ends);
    // Original design considered using FOR packing for cell_ids.
    // Since we are currently just storing pos-sorted without
    // much cell chunking I'll disable the FOR encoding
    unpack_128(src.cell_ids, cell_ids, idx);
}


} // end namespace BPCells