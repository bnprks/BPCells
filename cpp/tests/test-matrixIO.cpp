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
#include <matrixIterators/StoredMatrixWriter.h>

#include <matrixIterators/ImportMatrixHDF5.h>

#include <Eigen/Core>

using namespace BPCells;
using namespace ::testing;
using namespace Eigen;

template <typename T> bool matrix_identical(MatrixLoader<T> &mat1, MatrixLoader<T> &mat2) {
    mat1.restart();
    mat2.restart();
    MatrixIterator<T> i1((std::unique_ptr<MatrixLoader<T>>(&mat1)));
    MatrixIterator<T> i2((std::unique_ptr<MatrixLoader<T>>(&mat2)));
    i1.preserve_input_loader();
    i2.preserve_input_loader();

    while (true) {
        bool res1 = i1.nextCol();
        bool res2 = i2.nextCol();
        if (res1 != res2) {
            return false;
        }
        if (!res1) break;
        if (i1.currentCol() != i2.currentCol()) {
            return false;
        }
        while (true) {
            bool res1 = i1.nextValue();
            bool res2 = i2.nextValue();
            if (res1 != res2) {
                return false;
            }
            if (!res1) break;
            if (i1.row() != i2.row() || i1.col() != i2.col() || i1.val() != i2.val()) {
                return false;
            }
        }
    }
    return true;
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

Map<SparseMatrix<double>> get_map(const SparseMatrix<double> &mat) {
    return Map<SparseMatrix<double>>(
        mat.rows(),
        mat.cols(),
        mat.nonZeros(),
        (int *)mat.outerIndexPtr(),
        (int *)mat.innerIndexPtr(),
        (double *)mat.valuePtr()
    );
}

TEST(MatrixIO, UnpackedVec) {
    const Eigen::SparseMatrix<double> orig_mat = generate_mat(10, 10);

    MatrixConverterLoader<double, uint32_t> mat_i(std::make_unique<CSparseMatrix<double>>(get_map(orig_mat))
    );

    VecReaderWriterBuilder vb1(1024);
    auto w1 = StoredMatrixWriter<uint32_t>::createUnpacked(vb1);
    w1.write(mat_i);

    auto loader_double = MatrixConverterLoader<uint32_t, double>(
        std::make_unique<StoredMatrix<uint32_t>>(StoredMatrix<uint32_t>::openUnpacked(vb1))
    );

    CSparseMatrixWriter<double> w2;
    w2.write(loader_double);

    EXPECT_TRUE(w2.getMat().isApprox(orig_mat));
}

TEST(MatrixIO, PackedVec) {
    const Eigen::SparseMatrix<double> orig_mat = generate_mat(10, 10);

    MatrixConverterLoader<double, uint32_t> mat_i(std::make_unique<CSparseMatrix<double>>(get_map(orig_mat))
    );

    VecReaderWriterBuilder vb1(1024);
    auto w1 = StoredMatrixWriter<uint32_t>::createPacked(vb1);
    w1.write(mat_i);

    auto loader_double = MatrixConverterLoader<uint32_t, double>(
        std::make_unique<StoredMatrix<uint32_t>>(StoredMatrix<uint32_t>::openPacked(vb1))
    );

    CSparseMatrixWriter<double> w2;
    w2.write(loader_double);
    EXPECT_TRUE(w2.getMat().isApprox(orig_mat));
}

TEST(MatrixIO, SeekCSparse) {
    std::vector<Triplet<double>> triplets;
    const uint32_t n_row = 6;
    const uint32_t n_col = n_row - 1;

    // Matrix where each column j starts with an entry of j at row j+1
    for (int j = 0; j < n_col; j++) {
        for (int i = j + 1; i < n_row; i++) {
            triplets.push_back({i, j, (double) j});
        }
    }
    SparseMatrix<double> mat(n_row, n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());

    MatrixIterator<double> it(std::make_unique<CSparseMatrix<double>>(get_map(mat)));
    for (auto j : {4, 1, 3, 0, 2}) {
        it.seekCol(j);
        EXPECT_TRUE(it.nextValue());
        EXPECT_EQ(it.row(), j + 1);
        EXPECT_EQ(it.val(), j);
    }
}

TEST(MatrixIO, SeekStoredVec) {
    std::vector<Triplet<double>> triplets;
    const uint32_t n_row = 6;
    const uint32_t n_col = n_row - 1;

    // Matrix where each column j starts with an entry of j at row j+1
    for (int j = 0; j < n_col; j++) {
        for (int i = j + 1; i < n_row; i++) {
            triplets.push_back({i, j, (double) j});
        }
    }
    SparseMatrix<double> mat(n_row, n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());

    MatrixConverterLoader<double, uint32_t> mat_int(std::make_unique<CSparseMatrix<double>>(get_map(mat)));

    VecReaderWriterBuilder vb(1024);
    auto w = StoredMatrixWriter<uint32_t>::createUnpacked(vb);
    w.write(mat_int);

    std::unique_ptr<MatrixLoader<uint32_t>> mat_int2 =
        std::make_unique<StoredMatrix<uint32_t>>(StoredMatrix<uint32_t>::openUnpacked(vb));
    std::unique_ptr<MatrixLoader<double>> mat_double =
        std::make_unique<MatrixConverterLoader<uint32_t, double>>(std::move(mat_int2));

    MatrixIterator<double> it(std::move(mat_double));
    for (auto j : {4, 1, 3, 0, 2}) {
        it.seekCol(j);
        EXPECT_TRUE(it.nextValue());
        EXPECT_EQ(it.row(), j + 1);
        EXPECT_EQ(it.val(), j);
    }
}

TEST(MatrixIO, ColSelectCSparse) {
    std::vector<Triplet<double>> triplets;
    const uint32_t n_row = 6;
    const uint32_t n_col = n_row - 1;

    for (int j = 0; j < n_col; j++) {
        for (int i = 0; i < n_row; i++) {
            triplets.push_back({i, j, (double) j + n_col * i});
        }
    }

    SparseMatrix<double> mat(n_row, n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());

    std::unique_ptr<MatrixLoader<double>> loader = std::make_unique<CSparseMatrix<double>>(get_map(mat));
    loader = std::make_unique<MatrixColSelect<double>>(std::move(loader), std::vector<uint32_t>{0,4,2});

    CSparseMatrixWriter<double> writer;
    writer.write(*loader);

    // Test that the matrix is the same
    EXPECT_EQ(MatrixXd(writer.getMat()), MatrixXd(mat)(all, {0, 4, 2}));

    MatrixIterator<double> it(std::move(loader));
    // Check that seeking columns works
    for (auto j : {2, 0, 1}) {
        it.seekCol(j);
        EXPECT_TRUE(it.nextValue());
        EXPECT_EQ(it.row(), 0);
        EXPECT_EQ(it.val(), std::vector<uint32_t>({0, 4, 2})[j]);
    }
}

TEST(MatrixIO, RowSelectCSparse) {
    std::vector<Triplet<double>> triplets;
    const uint32_t n_row = 6;
    const uint32_t n_col = n_row - 1;

    for (int j = 0; j < n_col; j++) {
        for (int i = 0; i < n_row; i++) {
            triplets.push_back({i, j, (double) j + n_col * i});
        }
    }

    SparseMatrix<double> mat(n_row, n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());

    MatrixRowSelect<double> select_1(std::make_unique<CSparseMatrix<double>>(get_map(mat)), {0, 4, 2});

    CSparseMatrixWriter<double> writer1;
    writer1.write(select_1);

    EXPECT_EQ(MatrixXd(writer1.getMat()), MatrixXd(mat)({0, 4, 2}, all));

    MatrixRowSelect<double> select_2(std::make_unique<CSparseMatrix<double>>(get_map(mat)), {0, 2, 4});

    CSparseMatrixWriter<double> writer2;
    writer2.write(select_2);

    EXPECT_EQ(MatrixXd(writer2.getMat()), MatrixXd(mat)({0, 2, 4}, all));
}

TEST(MatrixIO, ConcatRows) {
    SparseMatrix<double> m1 = generate_mat(3000, 10, 12512);
    SparseMatrix<double> m2 =
        generate_mat(1, 10, 7345); // Very few rows to try getting 0-entry columns
    SparseMatrix<double> m3 = generate_mat(256, 10, 3864);
    SparseMatrix<double> mx = generate_mat(8, 5, 92568);

    MatrixXd concat_dense(m1.rows() + m2.rows() + m3.rows(), m1.cols());
    concat_dense << MatrixXd(m1), MatrixXd(m2), MatrixXd(m3);
    SparseMatrix<double> concat = concat_dense.sparseView();

    std::vector<std::unique_ptr<MatrixLoader<double>>> error_mat_vec;
    error_mat_vec.push_back(std::make_unique<CSparseMatrix<double>>(get_map(m1)));
    error_mat_vec.push_back(std::make_unique<CSparseMatrix<double>>(get_map(mx)));
    EXPECT_ANY_THROW(ConcatRows<double>(std::move(error_mat_vec), 0));

    CSparseMatrixWriter<double> res;
    std::vector<std::unique_ptr<MatrixLoader<double>>> mat_vec;
    mat_vec.push_back(std::make_unique<CSparseMatrix<double>>(get_map(m1)));
    mat_vec.push_back(std::make_unique<CSparseMatrix<double>>(get_map(m2)));
    mat_vec.push_back(std::make_unique<CSparseMatrix<double>>(get_map(m3)));
    ConcatRows<double> my_concat(std::move(mat_vec), 0);
    res.write(my_concat);

    EXPECT_TRUE(res.getMat().isApprox(concat));
}

TEST(MatrixIO, ConcatCols) {
    SparseMatrix<double> m1 = generate_mat(10, 3000, 12512);
    SparseMatrix<double> m2 =
        generate_mat(10, 1, 7345); // Very few rows to try getting 0-entry columns
    SparseMatrix<double> m3 = generate_mat(10, 256, 3864);
    SparseMatrix<double> mx = generate_mat(5, 8, 92568);

    MatrixXd concat_dense(m1.rows(), m1.cols() + m2.cols() + m3.cols());
    concat_dense << MatrixXd(m1), MatrixXd(m2), MatrixXd(m3);
    SparseMatrix<double> concat = concat_dense.sparseView();

    std::vector<std::unique_ptr<MatrixLoader<double>>> error_mat_vec;
    error_mat_vec.push_back(std::make_unique<CSparseMatrix<double>>(get_map(m1)));
    error_mat_vec.push_back(std::make_unique<CSparseMatrix<double>>(get_map(mx)));
    EXPECT_ANY_THROW(ConcatCols<double>(std::move(error_mat_vec), 0));

    CSparseMatrixWriter<double> res;
    std::vector<std::unique_ptr<MatrixLoader<double>>> mat_vec;
    mat_vec.push_back(std::make_unique<CSparseMatrix<double>>(get_map(m1)));
    mat_vec.push_back(std::make_unique<CSparseMatrix<double>>(get_map(m2)));
    mat_vec.push_back(std::make_unique<CSparseMatrix<double>>(get_map(m3)));
    ConcatCols<double> my_concat(std::move(mat_vec), 0);

    res.write(my_concat);

    EXPECT_TRUE(res.getMat().isApprox(concat));
}

void test_order_rows(SparseMatrix<double> m, uint32_t load_size) {
    std::vector<uint64_t> col(m.outerIndexPtr(), m.outerIndexPtr() + m.cols() + 1);
    std::vector<uint32_t> row(m.innerIndexPtr(), m.innerIndexPtr() + m.nonZeros());
    std::vector<uint32_t> row_orig(row);
    std::vector<double> val(m.valuePtr(), m.valuePtr() + m.nonZeros());
    std::vector<double> val_orig(val);

    // Within each column, shuffle the rows
    std::mt19937 gen;
    gen.seed(125123);
    for (uint32_t i = 0; i < m.cols(); i++)
        std::shuffle(row.data() + col[i], row.data() + col[i + 1], gen);
    gen.seed(125123);
    for (uint32_t i = 0; i < m.cols(); i++)
        std::shuffle(val.data() + col[i], val.data() + col[i + 1], gen);

    uint32_t row_count = m.rows();

    StoredMatrix<double> orig(
        UIntReader(std::make_unique<VecUIntReader>(row_orig.data(), row_orig.size()), 16),
        DoubleReader(std::make_unique<VecNumReader<double>>(val_orig.data(), val_orig.size()), 16),
        ULongReader(std::make_unique<VecNumReader<uint64_t>>(col.data(), col.size()), 16),
        row_count,
        std::make_unique<VecStringReader>(std::vector<std::string>()),
        std::make_unique<VecStringReader>(std::vector<std::string>())
    );

    StoredMatrix<double> unordered(
        UIntReader(std::make_unique<VecUIntReader>(row.data(), row.size()), 16),
        DoubleReader(std::make_unique<VecNumReader<double>>(val.data(), val.size()), 16),
        ULongReader(std::make_unique<VecNumReader<uint64_t>>(col.data(), col.size()), 16),
        row_count,
        std::make_unique<VecStringReader>(std::vector<std::string>()),
        std::make_unique<VecStringReader>(std::vector<std::string>())
    );
    ASSERT_FALSE(matrix_identical(unordered, orig));

    OrderRows<double> ordered(std::unique_ptr<StoredMatrix<double>>(&unordered), load_size);
    ordered.preserve_input_loader();

    ASSERT_TRUE(matrix_identical(ordered, orig));
}

TEST(MatrixIO, OrderRows) {
    test_order_rows(generate_mat(5, 5), 100);
    test_order_rows(generate_mat(1000, 100), 16);
    test_order_rows(generate_mat(1000, 100), 1024);
}