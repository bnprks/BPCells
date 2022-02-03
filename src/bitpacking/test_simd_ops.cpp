#include <cstdio>
#include "simd_vec.h"

using namespace BPCells;

static uint32_t a[4] = {0x01020304, 0x05060708, 0x090a0b0c, 0x0d0e0f10};

static uint32_t b[4] = {12, 16, 21, 39};
static uint32_t c[4] = {3, 4, 5, 6};

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

void print_ints(BPCells::vec i) {
    uint32_t array[4];
    BPCells::store((BPCells::vec *) &array, i);
    printf("%d\t%d\t%d\t%d\n", array[0], array[1], array[2], array[3]);
}

int main() {
    vec av = load((vec *) a);
    
    printf("SIMD_MODE: %d\n", _SIMDBP128_MODE_);
    printf("a: "); print_ints(av); printf("\n");
    printf("move_l_1: "); print_ints(move_l_1(av));
    printf("move_l_2: "); print_ints(move_l_2(av));
    printf("move_l_3: "); print_ints(move_l_3(av));
    printf("move_r_1: "); print_ints(move_r_1(av));
    printf("move_r_2: "); print_ints(move_r_2(av));
    printf("move_r_3: "); print_ints(move_r_3(av));
    
    printf("\n");
    vec bv = load((vec *) b);
    vec cv = load((vec *) c);
    printf("b: "); print_ints(bv);
    printf("c: "); print_ints(cv); printf("\n");
    printf("delta(b, c): "); print_ints(delta(bv, cv));
    printf("prefixSum(b, c): "); print_ints(prefixSum(bv, cv));
    
    return 0;
}