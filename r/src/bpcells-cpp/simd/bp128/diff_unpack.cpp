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
#define HWY_TARGET_INCLUDE "simd/bp128/diff_unpack.cpp"
#include <hwy/foreach_target.h>

#include "bp128-helper-inl.h"
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::bp128::HWY_NAMESPACE {

BP128_UNPACK_DECL(
    unpack_diff,
    // Setup Code
    // (none)
    ,
    // Transform Code (OutReg is the data)
    {
        OutReg = Add(OutReg, LoadU(d, ref));
        ref += 4;
    },
    // Call args in parens
    (ref, in, out),
    // Arguments
    const uint32_t *HWY_RESTRICT ref,
    const uint32_t *HWY_RESTRICT in,
    uint32_t *HWY_RESTRICT out
)

} // namespace BPCells::simd::bp128::HWY_NAMESPACE

HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace BPCells::simd::bp128 {

HWY_EXPORT(unpack_diff);
void unpack_diff(
    const uint32_t *HWY_RESTRICT ref,
    const uint32_t *HWY_RESTRICT in,
    uint32_t *HWY_RESTRICT out,
    const uint32_t bit
) {
    HWY_DYNAMIC_DISPATCH(unpack_diff)(ref, in, out, bit);
}

} // namespace BPCells::simd::bp128
#endif