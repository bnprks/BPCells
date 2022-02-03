#include <cstdio>
#include "simd_vec.h"
#include "bp128.h"

using namespace BPCells;

// static uint32_t a[4] = {0x01020304, 0x05060708, 0x090a0b0c, 0x0d0e0f10};
// static uint32_t b[4] = {0x11121314, 0x15161718, 0x191a1b1c, 0x1d1e1f01};

// static uint32_t scratch[4];
// void assert_mem(vec v, uint32_t a0, uint32_t a1, uint32_t a2, uint32_t a3) {
//     store((vec *) scratch, v);
//     if (
//         scratch[0] != a0 || 
//         scratch[1] != a1 ||
//         scratch[2] != a2 ||
//         scratch[3] != a3
//     ) {
//         throw new std::runtime_error("mismatch in memory!");
//     }
// }

int main() {
    printf("SIMD mode: %d\n", _SIMDBP128_MODE_);
    if (test_bitpacking()) {
        printf("all okay!\n");
    }
    return 0;
}