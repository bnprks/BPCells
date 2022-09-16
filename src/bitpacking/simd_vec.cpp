#include "simd_vec.h"

namespace BPCells {

#define BP128_UNROLL_LOOP1(start, counter, body)                                                   \
    counter = start;                                                                               \
    body

#define BP128_UNROLL_LOOP4(start, counter, body)                                                   \
    BP128_UNROLL_LOOP1((start + 0), counter, body)                                                 \
    BP128_UNROLL_LOOP1((start + 1), counter, body)                                                 \
    BP128_UNROLL_LOOP1((start + 2), counter, body)                                                 \
    BP128_UNROLL_LOOP1((start + 3), counter, body)

#define BP128_UNROLL_LOOP32(counter, body)                                                         \
    BP128_UNROLL_LOOP4(0, counter, body)                                                           \
    BP128_UNROLL_LOOP4(4, counter, body)                                                           \
    BP128_UNROLL_LOOP4(8, counter, body)                                                           \
    BP128_UNROLL_LOOP4(12, counter, body)                                                          \
    BP128_UNROLL_LOOP4(16, counter, body)                                                          \
    BP128_UNROLL_LOOP4(20, counter, body)                                                          \
    BP128_UNROLL_LOOP4(24, counter, body)                                                          \
    BP128_UNROLL_LOOP4(28, counter, body)

// Write end - start to out for 128 integers
void simdsubtract(uint32_t *end, const uint32_t *start) {
    const vec *_start = (vec *)start;
    vec *_end = (vec *)end;
    vec Reg;
    [[maybe_unused]] int i;
    BP128_UNROLL_LOOP32(i, {
        Reg = sub(load(_end), load(_start++));
        store(_end++, Reg);
    })
}

// Write end + start to end for 128 integers
void simdadd(uint32_t *end, const uint32_t *start) {
    const vec *_start = (vec *)start;
    vec *_end = (vec *)end;
    vec Reg;
    [[maybe_unused]] int i;
    BP128_UNROLL_LOOP32(i, {
        Reg = add(load(_end), load(_start++));
        store(_end++, Reg);
    })
}

// Get the maximum value of 128 integers
uint32_t simdmax(const uint32_t *in) {
    const vec *_in = (vec *)in;
    vec Reg = splat(0);
    [[maybe_unused]] int i;
    BP128_UNROLL_LOOP32(i, { Reg = max(Reg, load(_in++)); })
    // last of _tmp2 holds the max of maxs
    const vec _tmp1 = max(move_l_2(Reg), Reg);
    const vec _tmp2 = max(move_l_1(_tmp1), _tmp1);
    return extract_left(_tmp2);
}

} // end namespace BPCells