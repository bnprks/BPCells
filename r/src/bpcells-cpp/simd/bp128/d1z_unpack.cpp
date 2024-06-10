// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

// Compilation speed is a concern, hence why each bp128 function is split up (better
// parallelization). This is because 32-step unrolled loop * 32 bitwidths * 6 functions * 6-8
// architectures = a lot of code When in debug mode AND optimize is turned off, disable dynamic
// dispatch and compile for just baseline architecture
#if !defined(NDEBUG) && !defined(__OPTIMIZE__)
#define HWY_COMPILE_ONLY_STATIC 1
#endif

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "simd/bp128/d1z_unpack.cpp"
#include <hwy/foreach_target.h>

#include "bp128-helper-inl.h"
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::bp128::HWY_NAMESPACE {

void unpack_d1z(
    uint32_t initvalue,
    const uint32_t *HWY_RESTRICT in,
    uint32_t *HWY_RESTRICT out,
    const uint32_t bit
) {
    using namespace hwy::HWY_NAMESPACE;
    using D = Full128<uint32_t>;
    using Vec = Vec<D>;
    D d;
    Rebind<int32_t, D> d_signed;

    Vec prevOffset = Set(d, initvalue);
    Vec one = Set(d, 1);
    unpack(in, out, bit, [d, d_signed, &prevOffset, one](auto v) {
        // Zigzag decode: (i >> 1) ^ -(i & 1)
        v = Xor(
            ShiftRight<1>(v), BitCast(d, Neg(BitCast(d_signed, And(v, one))))
        );

        // OutReg = bitwise_xor(shift_r(OutReg, 1), sub(zero, bitwise_and(OutReg, one)));
        // Delta decode: Input: [a,b,c,d]; [?,?,?,h]; Output [a+h, a+b+h, a+b+c+h, a+b+c+d+h]
        const auto tmp1 = Add(ShiftLeftLanes<2>(d, v), v);
        const auto tmp2 = Add(ShiftLeftLanes<1>(d, tmp1), tmp1);
        const auto res = Add(tmp2, BroadcastLane<3>(prevOffset));

        // Update offset for next iteration
        prevOffset = res;

        return res;
    });
}

} // namespace BPCells::simd::bp128::HWY_NAMESPACE

HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace BPCells::simd::bp128 {

HWY_EXPORT(unpack_d1z);
void unpack_d1z(
    uint32_t initvalue,
    const uint32_t *HWY_RESTRICT in,
    uint32_t *HWY_RESTRICT out,
    const uint32_t bit
) {
    HWY_DYNAMIC_DISPATCH(unpack_d1z)(initvalue, in, out, bit);
}

} // namespace BPCells::simd::bp128
#endif