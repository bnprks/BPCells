// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <hwy/detect_targets.h>
#include <hwy/targets.h>
#include <hwy/highway.h>

#include <simd/bp128.h>

using namespace BPCells;

// Use vector extensions for a simple, loopy reference implementation.
// On Clang + GCC, this is confirmed to give bitwise correct results for the format
typedef uint32_t uint32x4 __attribute__((vector_size(16), aligned(4)));
typedef int32_t int32x4 __attribute__((vector_size(16), aligned(4)));

void pack128(const uint32x4 *in, uint32x4 *out, uint8_t B) {
    uint32x4 InReg, OutReg;
    uint32_t mask = (B == 32) ? 0xffffffff : ((1U << B) - 1);
    uint32x4 mask_vec = {mask, mask, mask, mask};
    if (B == 0) return;
    for (int i = 0; i < 32; i++) {
        unsigned int shift = (i * B) & 31;
        InReg = *in++;
        InReg = InReg & mask_vec;

        if (shift == 0) OutReg = InReg;
        else OutReg = OutReg | (InReg << shift);

        if (shift >= 32 - B) *out++ = OutReg;
        if (shift > 32 - B) OutReg = InReg >> (32 - shift);
    }
}

void unpack128(const uint32x4 *in, uint32x4 *out, uint8_t B) {
    uint32x4 InReg, OutReg;
    uint32_t mask = (B == 32) ? 0xffffffff : ((1U << B) - 1);
    uint32x4 mask_vec = {mask, mask, mask, mask};
    for (int i = 0; i < 32; i++) {
        unsigned int shift = (i * B) & 31;
        if (shift == 0) OutReg = InReg = *in++;
        else OutReg = InReg >> shift;

        if (shift > 32 - B) {
            InReg = *in++;
            OutReg |= InReg << (32 - shift);
        }
        *out++ = OutReg & mask;
    }
}

void encode_d1(uint32_t prev_val, uint32_t *inout) {
    for (int i = 0; i < 128; i++) {
        uint32_t tmp = inout[i];
        inout[i] = inout[i] - prev_val;
        prev_val = tmp;
    }
}

void decode_d1(uint32_t prev_val, uint32_t *inout) {
    for (int i = 0; i < 128; i++) {
        inout[i] += prev_val;
        prev_val = inout[i];
    }
}

void encode_zigzag(uint32_t *inout) {
    for (int i = 0; i < 128; i++) {
        inout[i] = (int32_t(inout[i]) >> 31) ^ (inout[i] << 1);
    }
}

void decode_zigzag(uint32_t *inout) {
    for (int i = 0; i < 128; i++) {
        inout[i] = (inout[i] >> 1) ^ -(inout[i] & 1);
    }
}

void encode_for(uint32_t ref_val, uint32_t *inout) {
    for (int i = 0; i < 128; i++) {
        inout[i] -= ref_val;
    }
}

void decode_for(uint32_t ref_val, uint32_t *inout) {
    for (int i = 0; i < 128; i++) {
        inout[i] += ref_val;
    }
}

void encode_diff(const uint32_t *ref_vals, uint32_t *inout) {
    for (int i = 0; i < 128; i++) {
        inout[i] -= ref_vals[i];
    }
}

void decode_diff(const uint32_t *ref_vals, uint32_t *inout) {
    for (int i = 0; i < 128; i++) {
        inout[i] += ref_vals[i];
    }
}

// Test helpers
uint32_t random_uint32_t() {
    return (
        (uint32_t)(rand() & 255) + (uint32_t)((rand() & 255) << 8) +
        (uint32_t)((rand() & 255) << 16) + (uint32_t)((rand() & 255) << 24)
    );
}
void fill_buf(uint32_t *buf, int bits) {
    uint32_t mask = (uint32_t)((UINT64_C(1) << bits) - 1);

    for (int i = 0; i < 128; i++) {
        buf[i] = random_uint32_t() & mask;
    }
}
bool equal_buf(uint32_t *buf1, uint32_t *buf2) {
    for (int i = 0; i < 128; i++) {
        if (buf1[i] != buf2[i]) return false;
    }
    return true;
}

bool equal_buf(uint32_t *buf1, uint32_t *buf2, uint8_t bits) {
    for (int i = 0; i < 4 * bits; i++) {
        if (buf1[i] != buf2[i]) return false;
    }
    return true;
}

TEST(Bitpacking, BP128) {
    uint32_t input_buf[128];
    uint32_t packed_buf[128];
    uint32_t ref_packed_buf[128];
    uint32_t output_buf[128];

    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        for (int bits = 0; bits <= 32; bits++) {
            fill_buf(input_buf, bits);
            fill_buf(packed_buf, 32); // fill the other bufs with garbage to avoid accidental success

            EXPECT_FALSE(equal_buf(input_buf, output_buf))
                << "Input and output aren't different before test, bits=" << bits;

            int simd_bits = simd::bp128::maxbits(input_buf);
            EXPECT_EQ(simd_bits, bits) << "Wrong number of simdmaxbits";

            EXPECT_FALSE(equal_buf(packed_buf, ref_packed_buf))
                << "Packed ref and obs aren't different before test, bits=" << bits;
            simd::bp128::pack(input_buf, packed_buf, simd_bits);
            simd::bp128::unpack(packed_buf, output_buf, simd_bits);
            
            pack128((uint32x4 *)input_buf, (uint32x4 *)ref_packed_buf, bits);
            EXPECT_TRUE(equal_buf(packed_buf, ref_packed_buf, simd_bits))
                << "Packed ref and obs are't identical after test, bits=" << bits;
            
            if (bits > 0) {
                EXPECT_TRUE(equal_buf(input_buf, output_buf))
                    << "Input buffer doesn't match output buffer, bits=" << bits;
            }

        }
    }
}

TEST(Bitpacking, BP128d1) {
    uint32_t input_buf[128];
    uint32_t packed_buf[128];
    uint32_t ref_packed_buf[128];
    uint32_t output_buf[128];

    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        for (int bits = 0; bits <= 32; bits++) {
            fill_buf(input_buf, bits);
            fill_buf(packed_buf, 32); // fill the other bufs with garbage to avoid accidental success
            uint32_t start_val = rand() & ((2 << 12) - 1);
            decode_d1(start_val, input_buf);

            EXPECT_FALSE(equal_buf(input_buf, output_buf))
                << "Input and output aren't different before test, bits=" << bits;

            int simd_bits = simd::bp128::maxbits_d1(start_val, input_buf);
            EXPECT_EQ(simd_bits, bits) << "Wrong number of simdmaxbits";

            simd::bp128::pack_d1(start_val, input_buf, packed_buf, simd_bits);
            simd::bp128::unpack_d1(start_val, packed_buf, output_buf, simd_bits);
            if (bits > 0) {
                EXPECT_TRUE(equal_buf(input_buf, output_buf))
                    << "Input buffer doesn't match output buffer, bits=" << bits;
            }

            EXPECT_FALSE(equal_buf(packed_buf, ref_packed_buf))
                << "Packed ref and obs aren't different before test, bits=" << bits;
            if (bits != 32) encode_d1(start_val, input_buf);
            pack128((uint32x4 *)input_buf, (uint32x4 *)ref_packed_buf, bits);
            EXPECT_TRUE(equal_buf(packed_buf, ref_packed_buf, simd_bits))
                << "Packed ref and obs are't identical after test, bits=" << bits;
        }
    }
}

TEST(Bitpacking, BP128d1z) {
    uint32_t input_buf[128];
    uint32_t packed_buf[128];
    uint32_t ref_packed_buf[128];
    uint32_t output_buf[128];

    // Test vanilla bitpacking
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        for (int bits = 0; bits <= 32; bits++) {
            fill_buf(input_buf, bits);
            fill_buf(packed_buf, 32); // fill the other bufs with garbage to avoid accidental success
            uint32_t start_val = rand() & ((2 << 12) - 1);
            decode_zigzag(input_buf);
            decode_d1(start_val, input_buf);

            EXPECT_FALSE(equal_buf(input_buf, output_buf))
                << "Input and output aren't different before test, bits=" << bits;

            int simd_bits = simd::bp128::maxbits_d1z(start_val, input_buf);
            EXPECT_EQ(simd_bits, bits) << "Wrong number of simdmaxbits";

            simd::bp128::pack_d1z(start_val, input_buf, packed_buf, simd_bits);
            simd::bp128::unpack_d1z(start_val, packed_buf, output_buf, simd_bits);
            if (bits > 0) {
                EXPECT_TRUE(equal_buf(input_buf, output_buf))
                    << "Input buffer doesn't match output buffer, bits=" << bits;
            }

            EXPECT_FALSE(equal_buf(packed_buf, ref_packed_buf))
                << "Packed ref and obs aren't different before test, bits=" << bits;
            if (bits != 32) {
                encode_d1(start_val, input_buf);
                encode_zigzag(input_buf);
            }
            pack128((uint32x4 *)input_buf, (uint32x4 *)ref_packed_buf, bits);
            EXPECT_TRUE(equal_buf(packed_buf, ref_packed_buf, simd_bits))
                << "Packed ref and obs are't identical after test, bits=" << bits;
        }
    }
}

TEST(Bitpacking, BP128for) {
    uint32_t input_buf[128];
    uint32_t packed_buf[128];
    uint32_t ref_packed_buf[128];
    uint32_t output_buf[128];

    // Test FOR bitpacking
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        for (int bits = 0; bits <= 32; bits++) {
            fill_buf(input_buf, bits);
            fill_buf(packed_buf, 32); // fill the other bufs with garbage to avoid accidental success
            uint32_t start_val = rand() & ((2 << 12) - 1);
            decode_for(start_val, input_buf);

            EXPECT_FALSE(equal_buf(input_buf, output_buf))
                << "Input and output aren't different before test, bits=" << bits;

            EXPECT_FALSE(equal_buf(packed_buf, ref_packed_buf))
                << "Packed ref and obs aren't different before test, bits=" << bits;
            
            int simd_bits = simd::bp128::maxbits_FOR(start_val, input_buf);
            EXPECT_EQ(simd_bits, bits) << "Wrong number of simdmaxbits";

            simd::bp128::pack_FOR(start_val, input_buf, packed_buf, simd_bits);
            simd::bp128::unpack_FOR(start_val, packed_buf, output_buf, simd_bits);
            if (bits > 0) {
                EXPECT_TRUE(equal_buf(input_buf, output_buf))
                    << "Input buffer doesn't match output buffer, bits=" << bits;
            }

            if (bits != 32) {
                encode_for(start_val, input_buf);
            }
            pack128((uint32x4 *)input_buf, (uint32x4 *)ref_packed_buf, bits);
            EXPECT_TRUE(equal_buf(packed_buf, ref_packed_buf, simd_bits))
                << "Packed ref and obs are't identical after test, bits=" << bits;
        }
    }
}

TEST(Bitpacking, BP128diff) {
    uint32_t input_buf[128];
    uint32_t packed_buf[128];
    uint32_t ref_packed_buf[128];
    uint32_t output_buf[128];
    uint32_t offsets[128];

    // Test FOR bitpacking
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        hwy::SetSupportedTargetsForTest(target);
        for (int bits = 0; bits <= 32; bits++) {
            fill_buf(input_buf, bits);
            fill_buf(packed_buf, 32); // fill the other bufs with garbage to avoid accidental success
            fill_buf(offsets, 12);
            decode_diff(offsets, input_buf);

            EXPECT_FALSE(equal_buf(input_buf, output_buf))
                << "Input and output aren't different before test, bits=" << bits;

            EXPECT_FALSE(equal_buf(packed_buf, ref_packed_buf))
                << "Packed ref and obs aren't different before test, bits=" << bits;
            
            int simd_bits = simd::bp128::maxbits_diff(offsets, input_buf);
            EXPECT_EQ(simd_bits, bits) << "Wrong number of simdmaxbits";

            simd::bp128::pack_diff(offsets, input_buf, packed_buf, simd_bits);
            simd::bp128::unpack_diff(offsets, packed_buf, output_buf, simd_bits);
            if (bits > 0) {
                EXPECT_TRUE(equal_buf(input_buf, output_buf))
                    << "Input buffer doesn't match output buffer, bits=" << bits;
            }

            if (bits != 32) {
                encode_diff(offsets, input_buf);
            }
            pack128((uint32x4 *)input_buf, (uint32x4 *)ref_packed_buf, bits);
            EXPECT_TRUE(equal_buf(packed_buf, ref_packed_buf, simd_bits))
                << "Packed ref and obs are't identical after test, bits=" << bits;
        }
    }
}