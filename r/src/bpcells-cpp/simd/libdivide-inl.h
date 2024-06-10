// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

// Per-target include guard
#if defined(BPCELLS_SIMD_LIBDIVIDE_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef BPCELLS_SIMD_LIBDIVIDE_INL_H_
#undef BPCELLS_SIMD_LIBDIVIDE_INL_H_
#else
#define BPCELLS_SIMD_LIBDIVIDE_INL_H_
#endif

#include "../../vendor/libdivide/libdivide.h"
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::HWY_NAMESPACE {

// Just copy the libdivide algorithm from `libdivide_u32_do` translated to highway instructions
template <class V>
HWY_INLINE V libdivide_u32(V numers, const struct libdivide::libdivide_u32_t *denom) {
    using namespace hwy::HWY_NAMESPACE;

    DFromV<V> d;

    uint8_t more = denom->more;
    if (!denom->magic) {
        return ShiftRightSame(numers, more);
    } else {
        V magic = Set(d, denom->magic);

        // Get the high u32 from a multiplication
        auto even = BitCast(d, ShiftRight<32>(MulEven(magic, numers)));
        auto odd = BitCast(d, MulOdd(magic, numers));
        auto q = OddEven(odd, even);

        if (more & libdivide::LIBDIVIDE_ADD_MARKER) {
            auto t = Add(ShiftRight<1>(Sub(numers, q)), q);
            uint32_t shift = more & libdivide::LIBDIVIDE_32_SHIFT_MASK;
            return ShiftRightSame(t, shift);
        } else {
            // All upper bits are 0,
            // don't need to mask them off.
            return ShiftRightSame(q, more);
        }
    }
}

} // namespace BPCells::simd::HWY_NAMESPACE

HWY_AFTER_NAMESPACE();

#endif // BPCELLS_SIMD_TRANSFORM_INL_H_