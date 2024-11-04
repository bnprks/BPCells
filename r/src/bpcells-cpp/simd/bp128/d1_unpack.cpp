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
#define HWY_TARGET_INCLUDE "simd/bp128/d1_unpack.cpp"
#include <hwy/foreach_target.h>

#include "bp128-helper-inl.h"
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::bp128::HWY_NAMESPACE {

BP128_UNPACK_DECL(
    unpack_d1,
    // Setup Code
    Vec prevOffset = Set(d, initvalue);
    ,
    // Transform Code (OutReg is the data)
    {
        // Do a prefix sum:
        // Input: [a,b,c,d]; [?,?,?,h]; Output [a+h, a+b+h, a+b+c+h, a+b+c+d+h]
        const auto tmp1 = Add(ShiftLeftLanes<2>(d, OutReg), OutReg);
        const auto tmp2 = Add(ShiftLeftLanes<1>(d, tmp1), tmp1);
        OutReg = Add(tmp2, BroadcastLane<3>(prevOffset));
        prevOffset = OutReg;
    },
    // Call args in parens
    (initvalue, in, out),
    // Arguments
    uint32_t initvalue,
    const uint32_t *HWY_RESTRICT in,
    uint32_t *HWY_RESTRICT out
)

} // namespace BPCells::simd::bp128::HWY_NAMESPACE

HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace BPCells::simd::bp128 {

HWY_EXPORT(unpack_d1);
void unpack_d1(
    uint32_t initvalue,
    const uint32_t *HWY_RESTRICT in,
    uint32_t *HWY_RESTRICT out,
    const uint32_t bit
) {
    HWY_DYNAMIC_DISPATCH(unpack_d1)(initvalue, in, out, bit);
}

} // namespace BPCells::simd::bp128
#endif