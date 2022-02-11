#include <gtest/gtest.h>
#include <bitpacking/bp128.h>

using namespace BPCells;

char** my_argv;
int my_argc;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    
    my_argc = argc;
    my_argv = argv;

    return RUN_ALL_TESTS();
}

TEST(Bitpacking, ExpectedArchitechture) {
    // Check that the test is being run on the expected simd mode
    ASSERT_EQ(my_argc, 3);
    ASSERT_STREQ(my_argv[1], "--arch");
    ASSERT_EQ(my_argv[2][0] - '0', _SIMDBP128_MODE_);
}

uint32_t random_uint32_t() {
    return (
        (uint32_t) (rand() & 255) +
        (uint32_t) ((rand() & 255) << 8) +
        (uint32_t) ((rand() & 255) << 16) +
        (uint32_t) ((rand() & 255) << 24)
    );
}
void fill_buf(uint32_t *buf, int bits) {
    uint32_t mask = (uint32_t) ((UINT64_C(1) << bits) - 1);
    
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
void print_buf(uint32_t *buf) {
    for (int i = 0; i < 128; i++) {
        printf("%d,\t", buf[i]);
    }
}

TEST(Bitpacking, BP128) {
    uint32_t input_buf[128];
    uint32_t packed_buf[128];
    uint32_t output_buf[128];

    // Test vanilla bitpacking
    for (int bits = 0; bits <= 32; bits++) {
        fill_buf(input_buf, bits);
        EXPECT_FALSE(equal_buf(input_buf, output_buf)) << 
            "Input and output aren't different before test, bits=" << bits;
            
        int simd_bits = simdmaxbits(input_buf);
        EXPECT_EQ(simd_bits, bits) << "Wrong number of simdmaxbits";
        
        simdpack(input_buf, packed_buf, simd_bits);
        simdunpack(packed_buf, output_buf, simd_bits);
        EXPECT_TRUE(equal_buf(input_buf, output_buf)) <<
            "Input buffer doesn't match output buffer, bits=" << bits;
    }
}

TEST(Bitpacking, BP128d1) {
    uint32_t input_buf[128];
    uint32_t packed_buf[128];
    uint32_t output_buf[128];

    for (int bits = 0; bits <= 32; bits++) {
        fill_buf(input_buf, bits);
        uint32_t start_val = rand() & ((2<<12) - 1);
        input_buf[0] += start_val;
        for (int i = 1; i < 128; i++)
            input_buf[i] += input_buf[i-1];

        EXPECT_FALSE(equal_buf(input_buf, output_buf)) <<
            "Input and output aren't different before test, bits=" << bits;
        
        int simd_bits = simdmaxbitsd1(start_val, input_buf);
        EXPECT_EQ(simd_bits, bits) << "Wrong number simdmaxbits, start_val=" << start_val;
        
        simdpackd1(start_val, input_buf, packed_buf, simd_bits);
        simdunpackd1(start_val, packed_buf, output_buf, simd_bits);
        
        EXPECT_TRUE(equal_buf(input_buf, output_buf)) << 
            "Input buffer doesn't match output buffer, bits=" << bits << 
            ", start_val=" << start_val;
    }
}

TEST(Bitpacking, BP128d1z) {
    uint32_t input_buf[128];
    uint32_t packed_buf[128];
    uint32_t output_buf[128];

    for (int bits = 0; bits <= 32; bits++) {
        fill_buf(input_buf, bits);
        for (int i = 0; i < 128; i++) {
            // Do a zigzag decode on the random offsets
            input_buf[i] = ((input_buf[i] >> 1) ^ -(input_buf[i] & 1));
        }
        uint32_t start_val = rand() & ((2<<12) - 1);
        input_buf[0] += start_val;
        for (int i = 1; i < 128; i++) {
            input_buf[i] += input_buf[i-1];
        }

        EXPECT_FALSE(equal_buf(input_buf, output_buf)) <<
            "Input and output aren't different before test, bits=" << bits;
        
        int simd_bits = simdmaxbitsd1z(start_val, input_buf);
        EXPECT_EQ(simd_bits, bits) << "Wrong number simdmaxbits, start_val=" << start_val;
        
        simdpackd1z(start_val, input_buf, packed_buf, simd_bits);
        simdunpackd1z(start_val, packed_buf, output_buf, simd_bits);
        
        EXPECT_TRUE(equal_buf(input_buf, output_buf)) << 
            "Input buffer doesn't match output buffer, bits=" << bits << 
            ", start_val=" << start_val;
    }
}

TEST(Bitpacking, BP128FOR) {
    uint32_t input_buf[128];
    uint32_t packed_buf[128];
    uint32_t output_buf[128];

    for (int bits = 0; bits <= 32; bits++) {
        fill_buf(input_buf, bits);
        uint32_t start_val = rand() & ((2<<12) - 1);
        
        for (int i = 0; i < 128; i++)
            input_buf[i] += start_val;

        EXPECT_FALSE(equal_buf(input_buf, output_buf)) <<
            "Input and output aren't different before test, bits=" << bits;
        
        uint32_t simd_bits, min_value;
        simdmaxbitsFORwithmin(input_buf, simd_bits, min_value);
        uint32_t real_min = INT32_MAX;
        for (int i = 0; i < 128; i++) {
            real_min = std::min(real_min, input_buf[i]);
        }
        EXPECT_EQ(real_min, min_value) << "Wrong minimum value, bits:" << bits;
        EXPECT_EQ(simd_bits, bits) << "Wrong number simdmaxbits, bits=" << bits <<
            ", start_val=" << start_val;
        
        EXPECT_EQ(simdmaxbitsFOR(min_value, input_buf), bits) <<
            "Wrong number from simdmaxbitsFOR, bits=" << bits <<
            ", start_val=" << start_val;
        
        simdpackFOR(start_val, input_buf, packed_buf, simd_bits);
        simdunpackFOR(start_val, packed_buf, output_buf, simd_bits);
        EXPECT_TRUE(equal_buf(input_buf, output_buf)) <<
            "Input buffer doesn't match output buffer, bits=" << bits <<
            ", start_val=" << start_val;
    }
}