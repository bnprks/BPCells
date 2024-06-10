// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <random>
#include <sstream>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <arrayIO/vector.h>
#include <matrixIterators/CSparseMatrix.h>
#include <matrixIterators/ConcatenateMatrix.h>
#include <matrixIterators/MatrixIndexSelect.h>
#include <matrixIterators/MatrixIterator.h>
#include <matrixIterators/StoredMatrix.h>

#include <matrixTransforms/Log1p.h>
#include <matrixTransforms/Pow.h>
#include <matrixTransforms/Min.h>
#include <matrixTransforms/Scale.h>
#include <matrixTransforms/Shift.h>

#include <Eigen/Core>

using namespace BPCells;
using namespace ::testing;
using namespace Eigen;


MatrixXd generate_dense_mat(uint32_t n_row, uint32_t n_col, uint32_t seed = 125124) {
    std::mt19937 gen(seed); // Standard mersenne_twister_engine
    std::uniform_int_distribution<> distrib(1, 20);

    MatrixXd mat(n_row, n_col);
    for (int i = 0; i < n_row; i++) {
        for (int j = 0; j < n_col; j++) {
            mat(i, j) = (double) distrib(gen);
        }
    }
    return mat;
}

MatrixXd generate_dense_vec(uint32_t n, uint32_t seed = 125124) {
    std::mt19937 gen(seed); // Standard mersenne_twister_engine
    std::uniform_int_distribution<> distrib(1, 20);

    VectorXd vec(n);
    for (int i = 0; i < n; i++) {
        vec(i) = (double) distrib(gen);
    }
    return vec;
}

SparseMatrix<double> generate_mat(uint32_t n_row, uint32_t n_col, uint32_t seed = 125124) {
    std::mt19937 gen(seed); // Standard mersenne_twister_engine
    std::uniform_int_distribution<> distrib(1, 20);
    std::uniform_int_distribution<> nonzero(0, 4); // 1/5 chance of being non-zero

    std::vector<Triplet<double>> triplets;

    for (int i = 0; i < n_row; i++) {
        for (int j = 0; j < n_col; j++) {
            if (0 == nonzero(gen)) {
                triplets.push_back({i, j, (double) distrib(gen)});
            }
        }
    }

    SparseMatrix<double> mat(n_row, n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

template <typename T> Map<T> get_map(T &mat) { return Map<T>(mat.data(), mat.rows(), mat.cols()); }

template <> Map<SparseMatrix<double>> get_map<SparseMatrix<double>>(SparseMatrix<double> &mat) {
    return Map<SparseMatrix<double>>(
        mat.rows(),
        mat.cols(),
        mat.nonZeros(),
        (int *)mat.outerIndexPtr(),
        (int *)mat.innerIndexPtr(),
        (double *)mat.valuePtr()
    );
}

TEST(MatrixMath, Multiply) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    MatrixXd b_right = generate_dense_mat(50, 4);
    MatrixXd b_left = generate_dense_mat(4, 100);

    CSparseMatrix mat_1(get_map(m1));

    MatrixXd r1 = mat_1.denseMultiplyLeft(get_map<MatrixXd>(b_left));
    MatrixXd r2 = mat_1.denseMultiplyRight(get_map<MatrixXd>(b_right));

    EXPECT_EQ(r1, b_left * m1);
    EXPECT_EQ(r2, m1 * b_right);

    VectorXd v_right = b_right.col(0);
    VectorXd v_left = b_left.row(0);

    VectorXd r3 = mat_1.vecMultiplyLeft(get_map<VectorXd>(v_left));
    VectorXd r4 = mat_1.vecMultiplyRight(get_map<VectorXd>(v_right));

    VectorXd ans3 = b_left({0}, all) * m1;
    EXPECT_EQ(r3, ans3);

    VectorXd ans4 = m1 * b_right.col(0);
    EXPECT_EQ(r4, ans4);
}

TEST(MatrixMath, Stats) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    // SparseMatrix<double> m1 = generate_mat(5, 3, 125123);

    // Getting useful represetntations of input
    MatrixXd d1(m1);
    ArrayXXd a1(d1);

    CSparseMatrix mat_1(get_map(m1));

    ArrayXXd row_stats(3, m1.rows());
    ArrayXXd col_stats(3, m1.cols());
    row_stats.row(0) = (a1 > 0).rowwise().count().cast<double>();
    row_stats.row(1) = a1.rowwise().mean();
    row_stats.row(2) =
        (a1.colwise() - row_stats.row(1).transpose()).square().rowwise().sum() / (a1.cols() - 1);

    col_stats.row(0) = (a1 > 0).colwise().count().cast<double>();
    col_stats.row(1) = a1.colwise().mean();
    col_stats.row(2) = (a1.rowwise() - col_stats.row(1)).square().colwise().sum() / (a1.rows() - 1);

    for (int i = 1; i <= 3; i++) {
        StatsResult r = mat_1.computeMatrixStats((Stats)i, (Stats)i);

        EXPECT_EQ(r.row_stats.rows(), i);
        EXPECT_EQ(r.row_stats.cols(), m1.rows());

        EXPECT_EQ(r.col_stats.rows(), i);
        EXPECT_EQ(r.col_stats.cols(), m1.cols());

        EXPECT_TRUE(r.row_stats.isApprox(row_stats.topRows(i)));
        EXPECT_TRUE(r.col_stats.isApprox(col_stats.topRows(i)));
    }
}

TEST(MatrixMath, Log1p) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    MatrixXd ans = MatrixXd(m1).array().log1p();

    Log1p mat_1_trans(std::make_unique<CSparseMatrix<double>>(get_map(m1)));

    CSparseMatrixWriter<double> res;
    res.write(mat_1_trans);

    EXPECT_TRUE(MatrixXd(res.getMat()).isApprox(ans));
}

TEST(MatrixMath, Log1pSIMD) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    MatrixXd ans = MatrixXd(m1).array().log1p();

    Log1pSIMD mat_1_trans(std::make_unique<CSparseMatrix<double>>(get_map(m1)));

    CSparseMatrixWriter<double> res;
    res.write(mat_1_trans);

    EXPECT_TRUE(MatrixXd(res.getMat()).isApprox(ans, Eigen::NumTraits<float>::dummy_precision()));
}

TEST(MatrixMath, Expm1) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    MatrixXd ans = MatrixXd(m1).array().expm1();

    Expm1 mat_1_trans(std::make_unique<CSparseMatrix<double>>(get_map(m1)));

    CSparseMatrixWriter<double> res;
    res.write(mat_1_trans);

    EXPECT_TRUE(MatrixXd(res.getMat()).isApprox(ans));
}

TEST(MatrixMath, Expm1SIMD) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    MatrixXd ans = MatrixXd(m1).array().expm1();

    Expm1SIMD mat_1_trans(std::make_unique<CSparseMatrix<double>>(get_map(m1)));

    CSparseMatrixWriter<double> res;
    res.write(mat_1_trans);

    EXPECT_TRUE(MatrixXd(res.getMat()).isApprox(ans, Eigen::NumTraits<float>::dummy_precision()));
}

TEST(MatrixMath, Pow) {
    double exp = 3.0;
    
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    MatrixXd ans = MatrixXd(m1).array().pow(exp);

    ArrayXd global_params(1);
    global_params = exp;

    Pow mat_1_trans(std::make_unique<CSparseMatrix<double>>(get_map(m1)), {{}, {}, global_params});

    CSparseMatrixWriter<double> res;
    res.write(mat_1_trans);

    EXPECT_TRUE(MatrixXd(res.getMat()).isApprox(ans));
}

TEST(MatrixMath, SquareSIMD) {
    double exp = 2.0;
    
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    MatrixXd ans = MatrixXd(m1).array().pow(exp);

    ArrayXd global_params(1);
    global_params = exp;

    SquareSIMD mat_1_trans(std::make_unique<CSparseMatrix<double>>(get_map(m1)), {{}, {}, global_params});

    CSparseMatrixWriter<double> res;
    res.write(mat_1_trans);

    EXPECT_TRUE(MatrixXd(res.getMat()).isApprox(ans));
}

TEST(MatrixMath, Square) {
    double exp = 2.0;
    
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    MatrixXd ans = MatrixXd(m1).array().pow(exp);

    Square mat_1_trans(std::make_unique<CSparseMatrix<double>>(get_map(m1)), {});

    CSparseMatrixWriter<double> res;
    res.write(mat_1_trans);

    EXPECT_TRUE(MatrixXd(res.getMat()).isApprox(ans));
}

TEST(MatrixMath, Min) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);

    double c = 4.0;

    ArrayXXd cmp(m1.rows(), m1.cols());
    cmp = c;
    MatrixXd ans = MatrixXd(m1).array().min(cmp);

    ArrayXd global_params(1);
    global_params = c;

    Min mat_1_trans(std::make_unique<CSparseMatrix<double>>(get_map(m1)), {{}, {}, global_params});

    CSparseMatrixWriter<double> res;
    res.write(mat_1_trans);

    EXPECT_TRUE(MatrixXd(res.getMat()).isApprox(ans));
}

void checkMultiplyOps(MatrixLoader<double> &mat, MatrixXd res) {
    MatrixXd mat_right = generate_dense_mat(mat.cols(), 4, 125412);
    MatrixXd mat_left = generate_dense_mat(4, mat.rows(), 8432);

    VectorXd vec_right = mat_right.col(0);
    VectorXd vec_left = mat_left.row(0).transpose();

    EXPECT_TRUE((res * mat_right).isApprox(mat.denseMultiplyRight(get_map<MatrixXd>(mat_right))));
    EXPECT_TRUE((res * vec_right).isApprox(mat.vecMultiplyRight(get_map<VectorXd>(vec_right))));

    EXPECT_TRUE((mat_left * res).isApprox(mat.denseMultiplyLeft(get_map<MatrixXd>(mat_left))));
    EXPECT_TRUE((vec_left.transpose() * res)
                    .transpose()
                    .isApprox(mat.vecMultiplyLeft(get_map<VectorXd>(vec_left))));
}

TEST(MatrixMath, Scale) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);

    ArrayXXd scale_row = generate_dense_vec(m1.rows(), 1513).transpose();
    ArrayXXd scale_col = generate_dense_vec(m1.cols(), 14582).transpose();

    // Scale row + col
    MatrixXd ans1 = (MatrixXd(m1).array().rowwise() * scale_col.row(0)).colwise() *
                    scale_row.row(0).transpose();
    Scale s1(std::make_unique<CSparseMatrix<double>>(get_map(m1)), TransformFit{scale_row, scale_col});
    checkMultiplyOps(s1, ans1);

    CSparseMatrixWriter<double> r1;
    s1.restart();
    r1.write(s1);
    EXPECT_TRUE(MatrixXd(r1.getMat()).isApprox(ans1));

    // Scale just row
    MatrixXd ans2 = MatrixXd(m1).array().colwise() * scale_row.row(0).transpose();
    Scale s2(std::make_unique<CSparseMatrix<double>>(get_map(m1)), TransformFit{scale_row, {}});
    checkMultiplyOps(s2, ans2);

    CSparseMatrixWriter<double> r2;
    s2.restart();
    r2.write(s2);
    EXPECT_TRUE(MatrixXd(r2.getMat()).isApprox(ans2));

    // Scale just col
    MatrixXd ans3 = MatrixXd(m1).array().rowwise() * scale_col.row(0);
    Scale s3(std::make_unique<CSparseMatrix<double>>(get_map(m1)), TransformFit{{}, get_map<ArrayXXd>(scale_col)});
    checkMultiplyOps(s3, ans3);

    CSparseMatrixWriter<double> r3;
    s3.restart();
    r3.write(s3);
    EXPECT_TRUE(MatrixXd(r3.getMat()).isApprox(ans3));
}

class SimpleRowShift : public MatrixTransformDense {
  public:
    SimpleRowShift(std::unique_ptr<MatrixLoader<double>> &&mat, TransformFit fit) : MatrixTransformDense(std::move(mat), fit) {}

    bool loadZeroSubtracted(MatrixLoader<double> &loader) override { return loader.load(); }
    void loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) override {
        for (uint32_t i = 0; i < count; i++) {
            values[i] = fit.row_params(0, start_row + i);
        }
    }
};

TEST(MatrixMath, TransformDense) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);

    ArrayXXd shift_row = generate_dense_vec(m1.rows(), 1513).transpose();

    // Shift row
    MatrixXd ans1 = MatrixXd(m1).array().colwise() + shift_row.row(0).transpose();
    SimpleRowShift s1(std::make_unique<CSparseMatrix<double>>(get_map(m1)), TransformFit{get_map<ArrayXXd>(shift_row), {}});

    checkMultiplyOps(s1, ans1);

    CSparseMatrixWriter<double> r1;
    s1.restart();
    r1.write(s1);

    EXPECT_TRUE(MatrixXd(r1.getMat()).isApprox(ans1));
}

TEST(MatrixMath, TransformDenseReorder) {
    // Regression test for issue #4: segfault on transforms of reordered matrices
    
    // Make enough rows so that we'll have more than one read chunk per col
    SparseMatrix<double> m1 = generate_mat(10000, 2, 125123);

    // Put in reverse order
    std::vector<uint32_t> row_select;
    for (int i = 0; i < m1.rows(); i++) {
        row_select.push_back(m1.rows() - 1 - i);
    }

    std::unique_ptr<MatrixLoader<double>> mat = std::make_unique<CSparseMatrix<double>>(get_map(m1));
    mat = std::make_unique<MatrixRowSelect<double>>(std::move(mat), row_select);

    ArrayXXd shift_row = generate_dense_vec(m1.rows(), 1513).transpose();
    // Shift row 
    mat = std::make_unique<SimpleRowShift>(std::move(mat), TransformFit{get_map<ArrayXXd>(shift_row), {}});

    MatrixXd ans1 = MatrixXd(m1).array().colwise().reverse().colwise() + shift_row.row(0).transpose();

    checkMultiplyOps(*mat, ans1);

    CSparseMatrixWriter<double> r1;
    mat->restart();
    r1.write(*mat);

    EXPECT_TRUE(MatrixXd(r1.getMat()).isApprox(ans1));
}

TEST(MatrixMath, Shift) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);

    ArrayXXd shift_row = generate_dense_vec(m1.rows(), 1513).transpose();
    ArrayXXd shift_col = generate_dense_vec(m1.cols(), 1242).transpose();

    // Shift row
    MatrixXd ans1 = MatrixXd(m1).array().colwise() + shift_row.row(0).transpose();
    std::unique_ptr<MatrixLoader<double>> mat = std::make_unique<CSparseMatrix<double>>(get_map(m1));
    mat = std::make_unique<ShiftRows>(std::move(mat), TransformFit{shift_row, {}});

    checkMultiplyOps(*mat, ans1);

    CSparseMatrixWriter<double> r1;
    mat->restart();
    r1.write(*mat);

    EXPECT_TRUE(MatrixXd(r1.getMat()).isApprox(ans1));

    // Shift col
    MatrixXd ans2 = MatrixXd(m1).array().rowwise() + shift_col.row(0);
    mat = std::make_unique<CSparseMatrix<double>>(get_map(m1));
    mat = std::make_unique<ShiftCols>(std::move(mat), TransformFit{{}, shift_col});

    checkMultiplyOps(*mat, ans2);

    CSparseMatrixWriter<double> r2;
    mat->restart();
    r2.write(*mat);

    EXPECT_TRUE(MatrixXd(r2.getMat()).isApprox(ans2));
}