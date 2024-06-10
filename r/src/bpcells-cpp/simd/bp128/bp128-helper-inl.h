// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

// Per-target include guard
#if defined(BPCELLS_SIMD_BP128_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef BPCELLS_SIMD_BP128_INL_H_
#undef BPCELLS_SIMD_BP128_INL_H_
#else
#define BPCELLS_SIMD_BP128_INL_H_
#endif

#include <cstring>
#include <hwy/highway.h>
#include <utility>
// Helpers for BP128 packing

HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::bp128::HWY_NAMESPACE {

// ##############################################################################
// #                  MANUAL LOOP UNROLL TEMPLATES                              #
// ##############################################################################

// How this works: The main "unroll" function takes a one-parameter lambda with signature void
// f(auto i). This `i` can be converted into a constexpr integer automatically.
//
// There are two requirements that add up to this solution:
//   1. The performance of the BP128 function is dependent on full loop unrolling, which is
//   otherwise hard to
//      force a compiler to do. The alternative is a mix of unrolling macros + templates, which I've
//      moved away from for clarity at the call-site.
//   2. The hwy::ShiftLeft and hwy::ShiftRight functions require a comptime-known shift as a
//   template parameter.
//
// From some testing on godbolt.org, it seems like this lambda style will be fully inlined for
// optimization levels -O1 and above. I'm fine with this given that it enables much less code
// duplication for all the BP128 variants
template <int N> struct ComptimeInteger {
    constexpr operator int() const { return N; }
};

template <int... inds, class F> void unroll_helper(F &&f, std::integer_sequence<int, inds...>) {
    // This is C++17 "pack expansion" with an expression called once per integer input
    (f(ComptimeInteger<inds>{}), ...);
}

template <class F> constexpr void unroll32(F &&f) {
    unroll_helper(std::forward<F>(f), std::make_integer_sequence<int, 32>{});
}

// ##############################################################################
// #                    BP128 CODEC CORE TEMPLATES                              #
// ##############################################################################

/**
 * @brief Generic BP128 unpacking. Unpacks 128 integers with known bitwidth.
 *
 * @tparam B Number of bits per integer
 * @tparam F Lambda function type
 * @param in Bitpacked data
 * @param out Uncompressed data
 * @param unpack_transform One-argument function that transforms the result just prior to output
 */
template <unsigned B, class F>
void unpack(const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, F &&unpack_transform) {
    using namespace hwy::HWY_NAMESPACE;
    using D = Full128<uint32_t>;
    using Vec = Vec<D>;

    D d;

    // Handle special cases here, since we can't give partial specialization for functions
    if constexpr (B == 0) {
        unroll32([&out, &unpack_transform, d](auto i) {
            StoreU(unpack_transform(Zero(d)), d, out + i*4);
        });
        return;
    }
    if constexpr (B == 32) {
        memcpy(out, in, sizeof(uint32_t) * 128);
        return;
    }

    Vec InReg, OutReg;
    Vec mask = Set(d, B < 32 ? (1U << B) - 1 : 0xffffffff);

    unroll32([&in, &out, &unpack_transform, d, mask, &InReg, &OutReg](auto i) {
        constexpr int shift = (i * B) & 31;
        if constexpr (shift == 0) {
            InReg = LoadU(d, in);
            in += 4;
            OutReg = InReg;
        } else {
            OutReg = ShiftRight<shift>(InReg);
        }

        if constexpr (shift > 32 - B) {
            InReg = LoadU(d, in);
            in += 4;
            OutReg = Or(OutReg, ShiftLeft<32 - shift>(InReg));
        }

        OutReg = And(OutReg, mask);

        OutReg = unpack_transform(OutReg);

        StoreU(OutReg, d, out);
        out += 4;
    });
}

/**
 * @brief Generic BP128 unpacking. Unpacks 128 integers with known bitwidth.
 *
 * @tparam B Number of bits per integer
 * @tparam MASK If true, then apply masking prior to packing for safety with too-large values
 * @tparam F Lambda function type
 * @param in Uncompressed data
 * @param out Bitpacked data
 * @param pack_transform One-argument function that transforms the result prior to bitpacking
 */
template <unsigned B, bool MASK, class F>
void pack(const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, F &&pack_transform) {
    // Handle special cases here, since we can't give partial specialization for functions
    if constexpr (B == 0) {
        return;
    }
    if constexpr (B == 32) {
        memcpy(out, in, sizeof(uint32_t) * 128);
        return;
    }
    using namespace hwy::HWY_NAMESPACE;
    using D = Full128<uint32_t>;
    using Vec = Vec<D>;

    D d;
    Vec InReg, OutReg;
    Vec mask = Set(d, B < 32 ? (1U << B) - 1 : 0xffffffff);

    unroll32([&in, &out, &pack_transform, d, &InReg, &OutReg, mask](auto i) {
        constexpr int shift = (i * B) & 31;
        InReg = LoadU(d, in);
        in += 4;

        InReg = pack_transform(InReg);

        if constexpr (MASK) {
            InReg = And(InReg, mask);
        }

        if constexpr (shift == 0) {
            OutReg = InReg;
        } else {
            OutReg = Or(OutReg, ShiftLeft<shift>(InReg));
        }

        if constexpr (shift >= 32 - B) {
            StoreU(OutReg, d, out);
            out += 4;
            if (shift > 32 - B) {
                OutReg = ShiftRight<32 - shift>(InReg);
            }
        }
    });
}

template <unsigned B, class F>
void pack_mask(const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, F &&pack_transform) {
    pack<B, true, F>(in, out, std::forward<F>(pack_transform));
}

// From
// https://github.com/lemire/simdcomp/blob/dd9317ff7df8ad6350d492e6a54600a9a69e7919/src/simdcomputil.c#L16-L29
//  Return number of bits needed to represent v
inline uint32_t bits(const uint32_t v) {
#ifdef _MSC_VER
    unsigned long answer;
    if (v == 0) {
        return 0;
    }
    _BitScanReverse(&answer, v);
    return answer + 1;
#else
    return v == 0 ? 0 : 32 - __builtin_clz(v); /* assume GCC-like compiler if not microsoft */
#endif
}

/**
 * @brief Calculate maximum bits needed to store a chunk of 128 numbers
 *
 * @tparam F
 * @param in
 * @param pack_transform
 * @return uint32_t
 */
template <class F> uint32_t maxbits(const uint32_t *HWY_RESTRICT in, F &&pack_transform) {
    using namespace hwy::HWY_NAMESPACE;
    using D = Full128<uint32_t>;
    using Vec = Vec<D>;

    D d;
    Vec accumulator = Zero(d);
    // Rather than "max" operations, we can do bitwise or, since we just care about the highest set
    // bit

    unroll32([&in, &pack_transform, d, &accumulator](auto i) {
        accumulator = Or(accumulator, pack_transform(LoadU(d, in)));
        in += 4;
    });

    // Reduce Or across all lanes
    const auto tmp1 = Or(ShiftRightLanes<2>(d, accumulator), accumulator);
    const auto tmp2 = Or(ShiftRightLanes<1>(d, tmp1), tmp1);
    uint32_t val = GetLane(tmp2);
    return bits(val);
}

#define BP128_SWITCH_CASE(i, function, ...)                                                        \
    case i:                                                                                        \
        function<i>(__VA_ARGS__);                                                                  \
        break;

#define BP128_SWITCH_CASE4(start, function, ...)                                                   \
    BP128_SWITCH_CASE((start + 0), function, __VA_ARGS__)                                          \
    BP128_SWITCH_CASE((start + 1), function, __VA_ARGS__)                                          \
    BP128_SWITCH_CASE((start + 2), function, __VA_ARGS__)                                          \
    BP128_SWITCH_CASE((start + 3), function, __VA_ARGS__)

#define BP128_SWITCH_CASE32(function, ...)                                                         \
    BP128_SWITCH_CASE4(1, function, __VA_ARGS__)                                                   \
    BP128_SWITCH_CASE4(5, function, __VA_ARGS__)                                                   \
    BP128_SWITCH_CASE4(9, function, __VA_ARGS__)                                                   \
    BP128_SWITCH_CASE4(13, function, __VA_ARGS__)                                                  \
    BP128_SWITCH_CASE4(17, function, __VA_ARGS__)                                                  \
    BP128_SWITCH_CASE4(21, function, __VA_ARGS__)                                                  \
    BP128_SWITCH_CASE4(25, function, __VA_ARGS__)                                                  \
    BP128_SWITCH_CASE4(29, function, __VA_ARGS__)

// This overload allows runtime selection of number of bits
template <class F>
void pack_mask(
    const uint32_t *HWY_RESTRICT in,
    uint32_t *HWY_RESTRICT out,
    const uint32_t bit,
    F &&unpack_transform
) {
    switch (bit) { BP128_SWITCH_CASE32(pack_mask, in, out, unpack_transform) }
}

// This overload allows runtime selection of number of bits
template <class F>
void unpack(
    const uint32_t *HWY_RESTRICT in,
    uint32_t *HWY_RESTRICT out,
    const uint32_t bit,
    F &&unpack_transform
) {
    switch (bit) {
    case 0:
        unpack<0>(in, out, unpack_transform);
        break;
    
    BP128_SWITCH_CASE32(unpack, in, out, unpack_transform)
    }
}

} // namespace BPCells::simd::bp128::HWY_NAMESPACE

HWY_AFTER_NAMESPACE();

#endif // BPCELLS_SIMD_BP128_INL_H_