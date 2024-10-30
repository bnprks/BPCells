// Copyright 2024 BPCells contributors
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

#include <matrixIterators/CSparseMatrix.h>
#include <matrixIterators/FilterZeros.h>

#include <Eigen/Core>

using namespace BPCells;
using namespace ::testing;
using namespace Eigen;

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

Map<SparseMatrix<double>> get_map(SparseMatrix<double> &mat) {
    return Map<SparseMatrix<double>>(
        mat.rows(),
        mat.cols(),
        mat.nonZeros(),
        (int *)mat.outerIndexPtr(),
        (int *)mat.innerIndexPtr(),
        (double *)mat.valuePtr()
    );
}

bool matrices_are_identical(SparseMatrix<double> a, SparseMatrix<double> b) {
    if (a.rows() != b.rows()) return false;
    if (a.cols() != b.cols()) return false;
    EXPECT_TRUE(a.isCompressed());
    EXPECT_TRUE(b.isCompressed());
    if (a.nonZeros() != b.nonZeros()) return false;
    for (size_t i = 0; i < a.nonZeros(); i++) {
        if (a.innerIndexPtr()[i] != b.innerIndexPtr()[i]) return false;
        if (a.valuePtr()[i] != b.valuePtr()[i]) return false;
    }
    if (a.outerSize() != b.outerSize()) return false;
    for (size_t i = 0; i < a.outerSize(); i++) {
        if (a.outerIndexPtr()[i] != b.outerIndexPtr()[i]) return false;
    } 
    return true;
}

TEST(MatrixIterators, FilterZeros) {
    // A matrix without any explicit zeros should not be transformed
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    FilterZeros<double> m1_filt(std::make_unique<CSparseMatrix<double>>(get_map(m1), std::unique_ptr<StringReader>(), std::unique_ptr<StringReader>(), 5));

    CSparseMatrixWriter<double> res_m1;
    res_m1.write(m1_filt);
    EXPECT_TRUE(matrices_are_identical(m1, res_m1.getMat()));

    // Introduce explicit zeros into m2
    SparseMatrix<double> m2 = m1;
    for (auto &x : m2.coeffs()) {
        if (x < 10) x = 0;
    }
    SparseMatrix<double> m2_pruned = m2.pruned();
    
    // The explicit zeros should be filtered out
    FilterZeros<double> m2_filt(std::make_unique<CSparseMatrix<double>>(get_map(m2), std::unique_ptr<StringReader>(), std::unique_ptr<StringReader>(), 5));
    CSparseMatrixWriter<double> res_m2;
    res_m2.write(m2_filt);
    EXPECT_FALSE(matrices_are_identical(m2, res_m2.getMat()));
    EXPECT_TRUE(matrices_are_identical(m2_pruned, res_m2.getMat()));
}
