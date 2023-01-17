#include <algorithm>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <lib/sleef_wrapper.h>

using namespace BPCells;
using ::testing::DoubleEq;
using ::testing::FloatEq;

char **my_argv;
int my_argc;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    my_argc = argc;
    my_argv = argv;

    return RUN_ALL_TESTS();
}

TEST(SimdSleef, ExpectedArchitechture) {
    // Check that the test is being run on the expected simd mode
    ASSERT_EQ(my_argc, 3);
    ASSERT_STREQ(my_argv[1], "--arch");
    ASSERT_EQ(my_argv[2][0] - '0', BPCELLS_SLEEF_MODE);
}

// Testing helper functions
template <typename VecOp, typename Op> void test_unary_op_float(VecOp vec_op, Op op) {
    float in[BPCELLS_VEC_FLOAT_SIZE];
    float out[BPCELLS_VEC_FLOAT_SIZE];
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        in[i] = 1.1 * (i + 1);
    }
    store_float(out, vec_op(load_float(in)));

    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        EXPECT_THAT(out[i], FloatEq(op(in[i])));
    }
}

template <typename VecOp, typename Op> void test_unary_op_double(VecOp vec_op, Op op) {
    double in[BPCELLS_VEC_DOUBLE_SIZE];
    double out[BPCELLS_VEC_DOUBLE_SIZE];
    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        in[i] = 1.1 * (i + 1);
    }
    store_double(out, vec_op(load_double(in)));

    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        EXPECT_THAT(out[i], DoubleEq(op(in[i])));
    }
}

template <typename VecOp, typename Op> void test_binary_op_float(VecOp vec_op, Op op) {
    float in1[BPCELLS_VEC_FLOAT_SIZE];
    float in2[BPCELLS_VEC_FLOAT_SIZE];
    float out[BPCELLS_VEC_FLOAT_SIZE];
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        in1[i] = 1.1 * (i + 1);
        in2[i] = 2.2 * (i + 1);
    }
    store_float(out, vec_op(load_float(in1), load_float(in2)));

    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        EXPECT_THAT(out[i], FloatEq(op(in1[i], in2[i])));
    }
}

template <typename VecOp, typename Op> void test_binary_op_double(VecOp vec_op, Op op) {
    double in1[BPCELLS_VEC_DOUBLE_SIZE];
    double in2[BPCELLS_VEC_DOUBLE_SIZE];
    double out[BPCELLS_VEC_DOUBLE_SIZE];
    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        in1[i] = 1.1 * (i + 1);
        in2[i] = 2.2 * (i + 1);
    }
    store_double(out, vec_op(load_double(in1), load_double(in2)));

    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        EXPECT_THAT(out[i], DoubleEq(op(in1[i], in2[i])));
    }
}

template <typename VecOp, typename Op> void test_ternary_op_float(VecOp vec_op, Op op) {
    float in1[BPCELLS_VEC_FLOAT_SIZE];
    float in2[BPCELLS_VEC_FLOAT_SIZE];
    float in3[BPCELLS_VEC_FLOAT_SIZE];
    float out[BPCELLS_VEC_FLOAT_SIZE];
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        in1[i] = 1.1 * (i + 1);
        in2[i] = 2.2 * (i + 1);
        in3[i] = 3.3 * (i + 1);
    }
    store_float(out, vec_op(load_float(in1), load_float(in2), load_float(in3)));

    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        EXPECT_THAT(out[i], FloatEq(op(in1[i], in2[i], in3[i])));
    }
}

template <typename VecOp, typename Op> void test_ternary_op_double(VecOp vec_op, Op op) {
    double in1[BPCELLS_VEC_DOUBLE_SIZE];
    double in2[BPCELLS_VEC_DOUBLE_SIZE];
    double in3[BPCELLS_VEC_DOUBLE_SIZE];
    double out[BPCELLS_VEC_DOUBLE_SIZE];
    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        in1[i] = 1.1 * (i + 1);
        in2[i] = 2.2 * (i + 1);
        in3[i] = 3.3 * (i + 1);
    }
    store_double(out, vec_op(load_double(in1), load_double(in2), load_double(in3)));

    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        EXPECT_THAT(out[i], DoubleEq(op(in1[i], in2[i], in3[i])));
    }
}

TEST(SimdSleef, LoadStore) {
    float a_f[BPCELLS_VEC_FLOAT_SIZE + 1];
    float out_f[BPCELLS_VEC_FLOAT_SIZE + 2];

    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE + 1; i++) {
        a_f[i] = 1.0 * (i + 1);
    }
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE + 2; i++) {
        out_f[i] = 0.0;
    }

    // Do one load + store
    vec_float vf = load_float(a_f);
    store_float(out_f + 1, vf);
    EXPECT_EQ(out_f[0], 0);
    EXPECT_EQ(out_f[BPCELLS_VEC_FLOAT_SIZE + 1], 0);
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        EXPECT_EQ(out_f[i + 1], 1.0 * (i + 1));
    }
    // Do load + store at different alignment
    vf = load_float(a_f + 1);
    store_float(out_f + 2, vf);
    EXPECT_EQ(out_f[0], 0);
    EXPECT_EQ(out_f[BPCELLS_VEC_FLOAT_SIZE + 1], 1.0 * (BPCELLS_VEC_FLOAT_SIZE + 1));
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        EXPECT_EQ(out_f[i + 2], 1.0 * (i + 2));
    }

    double a_d[BPCELLS_VEC_DOUBLE_SIZE + 1];
    double out_d[BPCELLS_VEC_DOUBLE_SIZE + 2];

    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE + 1; i++) {
        a_d[i] = 1.0 * (i + 1);
    }
    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE + 2; i++) {
        out_d[i] = 0.0;
    }

    // Do one load + store
    vec_double vd = load_double(a_d);
    store_double(out_d + 1, vd);
    EXPECT_EQ(out_d[0], 0);
    EXPECT_EQ(out_d[BPCELLS_VEC_DOUBLE_SIZE + 1], 0);
    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        EXPECT_EQ(out_d[i + 1], 1.0 * (i + 1));
    }

    // Do load + store at different alignment
    vd = load_double(a_d + 1);
    store_double(out_d + 2, vd);
    EXPECT_EQ(out_d[0], 0);
    EXPECT_EQ(out_d[BPCELLS_VEC_DOUBLE_SIZE + 1], 1.0 * (BPCELLS_VEC_DOUBLE_SIZE + 1));
    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        EXPECT_EQ(out_d[i + 2], 1.0 * (i + 2));
    }
}

TEST(SimdSleef, Splat) {
    float out_f[BPCELLS_VEC_FLOAT_SIZE];

    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        out_f[i] = 0.0;
    }
    store_float(out_f, splat_float(1.0));
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        EXPECT_EQ(out_f[i], 1.0);
    }

    double out_d[BPCELLS_VEC_DOUBLE_SIZE + 2];

    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        out_d[i] = 0.0;
    }
    store_double(out_d, splat_double(1.0));
    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        EXPECT_EQ(out_d[i], 1.0);
    }
}

TEST(SimdSleef, StoreConvert) {
    float in1[BPCELLS_VEC_FLOAT_SIZE];
    double out1[BPCELLS_VEC_FLOAT_SIZE + 2];
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE + 2; i++) {
        out1[i] = 0.0;
    }
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        in1[i] = 1.0 * (i + 1);
    }
    store_float_to_double(out1 + 1, load_float(in1));
    EXPECT_EQ(out1[0], 0.0);
    EXPECT_EQ(out1[BPCELLS_VEC_FLOAT_SIZE + 1], 0.0);
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        EXPECT_EQ(out1[i + 1], 1.0 * (i + 1));
    }
    // Try again at a different output offset
    EXPECT_EQ(out1[0], 0.0);
    EXPECT_EQ(out1[1], 1.0);
    store_float_to_double(out1 + 2, load_float(in1));
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        EXPECT_EQ(out1[i + 2], 1.0 * (i + 1));
    }

    // Test on doubles
    double in2[BPCELLS_VEC_DOUBLE_SIZE];
    float out2[BPCELLS_VEC_DOUBLE_SIZE + 2];
    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE + 2; i++) {
        out2[i] = 0.0;
    }
    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        in2[i] = 1.0 * (i + 1);
    }
    store_double_to_float(out2 + 1, load_double(in2));
    EXPECT_EQ(out2[0], 0.0);
    EXPECT_EQ(out2[BPCELLS_VEC_DOUBLE_SIZE + 1], 0.0);
    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        EXPECT_EQ(out2[i + 1], 1.0 * (i + 1));
    }
    // Try again at a different output offset
    EXPECT_EQ(out2[0], 0.0);
    EXPECT_EQ(out2[1], 1.0);
    store_double_to_float(out2 + 2, load_double(in2));
    for (int i = 0; i < BPCELLS_VEC_DOUBLE_SIZE; i++) {
        EXPECT_EQ(out2[i + 2], 1.0 * (i + 1));
    }
}

TEST(SimdSleef, LoadConvert) {
    double in[BPCELLS_VEC_FLOAT_SIZE];
    float out[BPCELLS_VEC_FLOAT_SIZE];
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        in[i] = 1.0 * (i + 1);
        out[i] = -1.0;
    }
    store_float(out, load_double_to_float(in));
    for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
        EXPECT_EQ(out[i], 1.0 * (i + 1));
    }
}

TEST(SimdSleef, Log1p) {
    test_unary_op_float(log1p_f, log1pf);
    test_unary_op_double(log1p_d, log1p);
}

TEST(SimdSleef, Expm1) {
    test_unary_op_float(expm1_f, expm1f);
    test_unary_op_double(expm1_d, expm1);
}

TEST(SimdSleef, Pow) {
    test_binary_op_float(pow_f, pow);
    test_binary_op_double(pow_d, pow);
}

TEST(SimdSleef, FusedMultiplyAdd) {
    test_ternary_op_float(fma_f, [](float a, float b, float c) { return a * b + c; });
    test_ternary_op_double(fma_d, [](double a, double b, double c) { return a * b + c; });
}

TEST(SimdSleef, Add) {
    test_binary_op_float(add_f, [](float a, float b) { return a + b; });
    test_binary_op_double(add_d, [](double a, double b) { return a + b; });
}

TEST(SimdSleef, Mul) {
    test_binary_op_float(mul_f, [](float a, float b) { return a * b; });
    test_binary_op_double(mul_d, [](double a, double b) { return a * b; });
}

TEST(SimdSleef, Min) {
    test_binary_op_float(min_f, [](float a, float b) { return std::min(a, b); });
    test_binary_op_double(min_d, [](double a, double b) { return std::min(a, b); });
}

TEST(SimdSleef, Max) {
    test_binary_op_float(max_f, [](float a, float b) { return std::max(a, b); });
    test_binary_op_double(max_d, [](double a, double b) { return std::max(a, b); });
}

TEST(SimdSleef, Neg) {
    test_unary_op_float(neg_f, [](float a) { return -a; });
    test_unary_op_double(neg_d, [](double a) { return -a; });
}

TEST(SimdSleef, Rsqrt) {
    test_unary_op_float(rsqrt_f, [](float a) { return 1 / sqrtf(a); });
    test_unary_op_double(rsqrt_d, [](double a) { return 1 / sqrt(a); });
}
