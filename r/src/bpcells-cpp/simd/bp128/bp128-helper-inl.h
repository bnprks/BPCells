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
// #                           REFERENCE RESOURCES                              #
// ##############################################################################
// - Original BP-128 paper by Daniel Lemire: https://arxiv.org/pdf/1209.2137
// - Check out pack128 and unpack128 from the tests in `test-bp128.cpp` for an easy-to-read implementation
//   of the vanilla BP-128 logic
// - Daniel Lemire repo with some of these functions implemented: https://github.com/lemire/simdcomp
//    - A bit tough to read, but was an important reference for the original implementation


// ##############################################################################
// #                              HELPER MACROS                                 #
// ##############################################################################

/**
 * Note: A previous attempt did this all with nice template functions and closures
 * (immediately after the Genentech code incorporation).
 * That proved to not inline everything as desired, so now there's a big mess of macros to handle
 * code generation. The macros have very specific, somewhat unusual argument requirements, though
 * these are documented by each macro and should be easy to follow from examples.
 */

// BP128_UNROLL_LOOP32(counter, body) will repeat `body` 32 times, with interleaved `counter = i`
// lines
#define BP128_UNROLL_LOOP1(start, counter, body)                                                   \
    {                                                                                              \
        [[maybe_unused]] constexpr int counter = start;                                            \
        body                                                                                       \
    }

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

// BP_SWITCH_CASE31 will switch over a variable for case 1: ... case 31:
// ARGS_IN_PARENS should be a list of arguments already enclosed in parens, e.g. `(arg1, arg2,
// arg3)`
#define BP128_SWITCH_CASE(i, function, ARGS_IN_PARENS)                                             \
    case i:                                                                                        \
        function<i> ARGS_IN_PARENS;                                                                \
        break;

#define BP128_SWITCH_CASE4(start, function, ARGS_IN_PARENS)                                        \
    BP128_SWITCH_CASE((start + 0), function, ARGS_IN_PARENS)                                       \
    BP128_SWITCH_CASE((start + 1), function, ARGS_IN_PARENS)                                       \
    BP128_SWITCH_CASE((start + 2), function, ARGS_IN_PARENS)                                       \
    BP128_SWITCH_CASE((start + 3), function, ARGS_IN_PARENS)

#define BP128_SWITCH_CASE31(function, ARGS_IN_PARENS)                                              \
    BP128_SWITCH_CASE4(1, function, ARGS_IN_PARENS)                                                \
    BP128_SWITCH_CASE4(5, function, ARGS_IN_PARENS)                                                \
    BP128_SWITCH_CASE4(9, function, ARGS_IN_PARENS)                                                \
    BP128_SWITCH_CASE4(13, function, ARGS_IN_PARENS)                                               \
    BP128_SWITCH_CASE4(17, function, ARGS_IN_PARENS)                                               \
    BP128_SWITCH_CASE4(21, function, ARGS_IN_PARENS)                                               \
    BP128_SWITCH_CASE4(25, function, ARGS_IN_PARENS)                                               \
    BP128_SWITCH_CASE(29, function, ARGS_IN_PARENS)                                                \
    BP128_SWITCH_CASE(30, function, ARGS_IN_PARENS)                                                \
    BP128_SWITCH_CASE(31, function, ARGS_IN_PARENS)

// ##############################################################################
// #                     BP128 CODEC CORE MACROS                                #
// ##############################################################################

/** BP128_PACK_DECL helps to declare bitpacking functions with full unrolling.
 *
 * Declares FN_NAME_static, which is templated with a static bitwidth and
 * FN_NAME, which takes a dynamic bitwidth as its final parameter.
 *
 * No masking of input data is performed, so inputs that require higher than specified number
 * of bits will mess up adjacent values.
 *
 * @param FN_NAME (variable name) Name to use for the output function
 * @param setup (code block) Code to run setup just prior to starting the bitpacking loop
 * @param pack_transform (code block) Code block that applies transforms on `InReg` prior to
 * bitpacking, re-assinging the result to the `InReg` variable
 * @param CALL_ARGS_IN_PARENS (parameter list in parens) List of the parameter names in function
 * signature without types, e.g. `(in, out)`
 * @param ... Function signature. Must be at least: `const uint32_t *HWY_RESTRICT in, uint32_t
 * *HWY_RESTRICT out`
 *
 */
#define BP128_PACK_DECL(FN_NAME, setup, pack_transform, CALL_ARGS_IN_PARENS, ...)                  \
    template <unsigned BITS> void FN_NAME##_static(__VA_ARGS__) {                                  \
        using namespace hwy::HWY_NAMESPACE;                                                        \
        using D = Full128<uint32_t>;                                                               \
        using Vec = Vec<D>;                                                                        \
        D d;                                                                                       \
        Vec InReg, OutReg;                                                                         \
        setup;                                                                                     \
        BP128_UNROLL_LOOP32(i, {                                                                   \
            constexpr int shift = (i * BITS) & 31;                                                 \
            InReg = LoadU(d, in);                                                                  \
            in += 4;                                                                               \
            pack_transform;                                                                        \
            if constexpr (shift == 0) {                                                            \
                OutReg = InReg;                                                                    \
            } else {                                                                               \
                OutReg = Or(OutReg, ShiftLeft<shift>(InReg));                                      \
            }                                                                                      \
            if constexpr (shift >= 32 - BITS) {                                                    \
                StoreU(OutReg, d, out);                                                            \
                out += 4;                                                                          \
                if (shift > 32 - BITS) {                                                           \
                    OutReg = ShiftRight<32 - shift>(InReg);                                        \
                }                                                                                  \
            }                                                                                      \
        })                                                                                         \
    }                                                                                              \
                                                                                                   \
    void FN_NAME(__VA_ARGS__, const uint32_t bit) {                                                \
        switch (bit) {                                                                             \
        case 0:                                                                                    \
            break;                                                                                 \
            /* Cases 1 to 31 */                                                                    \
            BP128_SWITCH_CASE31(FN_NAME##_static, CALL_ARGS_IN_PARENS)                             \
        case 32:                                                                                   \
            memcpy(out, in, sizeof(uint32_t) * 128);                                               \
        }                                                                                          \
    }

/** BP128_UNPACK_DECL helps to declare bit-unpacking functions with full unrolling.
 *
 * Declares FN_NAME_static, which is templated with a static bitwidth and
 * FN_NAME, which takes a dynamic bitwidth as its final parameter.
 *
 * @param FN_NAME (variable name) Name to use for the output function
 * @param setup (code block) Code to run setup just prior to starting the bitpacking loop
 * @param unpack_transform (code block) Code block that applies transforms on `OutReg` after
 * unpacking, re-assinging the result to the `OutReg` variable
 * @param CALL_ARGS_IN_PARENS (parameter list in parens) List of the parameter names in function
 * signature without types, e.g. `(in, out)`
 * @param ... Function signature. Must be at least: `const uint32_t *HWY_RESTRICT in, uint32_t
 * *HWY_RESTRICT out`
 *
 */
#define BP128_UNPACK_DECL(FN_NAME, setup, unpack_transform, CALL_ARGS_IN_PARENS, ...)              \
    template <unsigned BITS> void FN_NAME##_static(__VA_ARGS__) {                                           \
        using namespace hwy::HWY_NAMESPACE;                                                        \
        using D = Full128<uint32_t>;                                                               \
        using Vec = Vec<D>;                                                                        \
        D d;                                                                                       \
        Vec InReg, OutReg;                                                                         \
        Vec mask = Set(d, BITS < 32 ? (1U << BITS) - 1 : 0xffffffff);                              \
        setup;                                                                                     \
        BP128_UNROLL_LOOP32(i, {                                                                   \
            constexpr int shift = (i * BITS) & 31;                                                 \
            if constexpr (shift == 0) {                                                            \
                InReg = LoadU(d, in);                                                              \
                in += 4;                                                                           \
                OutReg = InReg;                                                                    \
            } else {                                                                               \
                OutReg = ShiftRight<shift>(InReg);                                                 \
            }                                                                                      \
            if constexpr (shift > 32 - BITS) {                                                     \
                InReg = LoadU(d, in);                                                              \
                in += 4;                                                                           \
                OutReg = Or(OutReg, ShiftLeft<32 - shift>(InReg));                                 \
            }                                                                                      \
            OutReg = And(OutReg, mask);                                                            \
            unpack_transform;                                                                      \
            StoreU(OutReg, d, out);                                                                \
            out += 4;                                                                              \
        })                                                                                         \
    }                                                                                              \
                                                                                                   \
    void FN_NAME(__VA_ARGS__, uint32_t bit) {                                                      \
        switch (bit) {                                                                             \
        case 0: {                                                                                  \
            using namespace hwy::HWY_NAMESPACE;                                                    \
            using D = Full128<uint32_t>;                                                           \
            using Vec = Vec<D>;                                                                    \
            D d;                                                                                   \
            setup;                                                                                 \
            BP128_UNROLL_LOOP32(i, {                                                               \
                Vec OutReg = Zero(d);                                                              \
                unpack_transform;                                                                  \
                StoreU(OutReg, d, out + i * 4);                                                    \
            })                                                                                     \
            break;                                                                                 \
        }                                                                                          \
            BP128_SWITCH_CASE31(FN_NAME##_static, CALL_ARGS_IN_PARENS)                                      \
        case 32:                                                                                   \
            memcpy(out, in, sizeof(uint32_t) * 128);                                               \
            break;                                                                                 \
        }                                                                                          \
    }

/** BP128_MAXBITS_DECL helps to declare bitpacking functions to detect maximum bitwidth with full
 * unrolling
 *
 * @param FN_NAME (variable name) Name to use for the output function
 * @param setup (code block) Code to run setup just prior to starting the bitpacking loop
 * @param pack_transform (code block) Code block that applies transforms on `InReg` prior to
 * bitpacking, re-assinging the result to the `InReg` variable
 * @param ... Function signature. Must be at least: `const uint32_t *HWY_RESTRICT in`
 *
 */
#define BP128_MAXBITS_DECL(FN_NAME, setup, pack_transform, ...)                                    \
    uint32_t FN_NAME(__VA_ARGS__) {                                                                \
        using namespace hwy::HWY_NAMESPACE;                                                        \
        using D = Full128<uint32_t>;                                                               \
        using Vec = Vec<D>;                                                                        \
                                                                                                   \
        D d;                                                                                       \
        Vec accumulator = Zero(d);                                                                 \
        setup;                                                                                     \
        /* Rather than "max" operations, we can do bitwise or, since we just care about the        \
         * highest set bit */                                                                      \
                                                                                                   \
        BP128_UNROLL_LOOP32(i, {                                                                   \
            Vec InReg = LoadU(d, in);                                                              \
            pack_transform;                                                                        \
            accumulator = Or(accumulator, InReg);                                                  \
            in += 4;                                                                               \
        })                                                                                         \
                                                                                                   \
        /* Reduce Or across all lanes */                                                           \
        const auto tmp1 = Or(ShiftRightLanes<2>(d, accumulator), accumulator);                     \
        const auto tmp2 = Or(ShiftRightLanes<1>(d, tmp1), tmp1);                                   \
        uint32_t val = GetLane(tmp2);                                                              \
        return bits(val);                                                                          \
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

} // namespace BPCells::simd::bp128::HWY_NAMESPACE

HWY_AFTER_NAMESPACE();

#endif // BPCELLS_SIMD_BP128_INL_H_