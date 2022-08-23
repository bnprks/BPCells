#include <filesystem>
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

#include <Eigen/Core>

namespace fs = std::filesystem;
using namespace BPCells;
using namespace ::testing;
using namespace Eigen;

template<typename T>
bool matrix_identical(MatrixLoader<T> &mat1, MatrixLoader<T> &mat2) {
    mat1.restart(); mat2.restart();
    MatrixIterator<T> i1(mat1);
    MatrixIterator<T> i2(mat2);

    while(true) {
        bool res1 = i1.nextCol();
        bool res2 = i2.nextCol();
        if(res1 != res2) {
            return false;
        }
        if (!res1) break;
        if(i1.currentCol() != i2.currentCol()) {
            return false;
        }
        while(true) {
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
                triplets.push_back({i, j, distrib(gen)});
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
        (double *)mat.valuePtr());
}

TEST(MatrixIO, UnpackedVec) {
    const Eigen::SparseMatrix<double> orig_mat = generate_mat(10, 10);
    CSparseMatrix mat_d(get_map(orig_mat));

    MatrixConverterLoader<double, uint32_t> mat_i(mat_d);

    VecReaderWriterBuilder vb1(1024);
    auto w1 = StoredMatrixWriter<uint32_t>::createUnpacked(vb1);
    w1.write(mat_i);

    auto loader = StoredMatrix<uint32_t>::openUnpacked(vb1);
    auto loader_double = MatrixConverterLoader<uint32_t, double>(loader);

    CSparseMatrixWriter w2;
    w2.write(loader_double);

    EXPECT_TRUE(w2.getMat().isApprox(orig_mat));
}

TEST(MatrixIO, PackedVec) {
    const Eigen::SparseMatrix<double> orig_mat = generate_mat(10, 10);
    CSparseMatrix mat_d(get_map(orig_mat));

    MatrixConverterLoader<double, uint32_t> mat_i(mat_d);

    VecReaderWriterBuilder vb1(1024);
    auto w1 = StoredMatrixWriter<uint32_t>::createPacked(vb1);
    w1.write(mat_i);

    auto loader = StoredMatrix<uint32_t>::openPacked(vb1, 1024);
    auto loader_double = MatrixConverterLoader<uint32_t, double>(loader);

    CSparseMatrixWriter w2;
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
            triplets.push_back({i, j, j});
        }
    }
    SparseMatrix<double> mat(n_row, n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());

    CSparseMatrix mat_l(get_map(mat));
    MatrixIterator it(mat_l);
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
            triplets.push_back({i, j, j});
        }
    }
    SparseMatrix<double> mat(n_row, n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());

    CSparseMatrix mat_double(get_map(mat));
    MatrixConverterLoader<double, uint32_t> mat_int(mat_double);

    VecReaderWriterBuilder vb(1024);
    auto w = StoredMatrixWriter<uint32_t>::createUnpacked(vb);
    w.write(mat_int);

    auto loader = StoredMatrix<uint32_t>::openUnpacked(vb);
    auto loader_double = MatrixConverterLoader<uint32_t, double>(loader);

    MatrixIterator it(loader_double);
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
            triplets.push_back({i, j, j + n_col * i});
        }
    }

    SparseMatrix<double> mat(n_row, n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());

    CSparseMatrix mat_l(get_map(mat));
    MatrixColSelect<double> mat_col_select(mat_l, {0, 4, 2});

    CSparseMatrixWriter writer;
    writer.write(mat_col_select);

    // Test that the matrix is the same
    EXPECT_EQ(MatrixXd(writer.getMat()), MatrixXd(mat)(all, {0, 4, 2}));

    MatrixIterator it(mat_col_select);
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
            triplets.push_back({i, j, j + n_col * i});
        }
    }

    SparseMatrix<double> mat(n_row, n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());

    CSparseMatrix mat_l(get_map(mat));
    MatrixRowSelect<double> select_1(mat_l, {0, 4, 2});

    CSparseMatrixWriter writer1;
    writer1.write(select_1);

    EXPECT_EQ(MatrixXd(writer1.getMat()), MatrixXd(mat)({0, 4, 2}, all));

    CSparseMatrix mat_l_2(get_map(mat));
    MatrixRowSelect<double> select_2(mat_l_2, {0, 2, 4});

    CSparseMatrixWriter writer2;
    writer2.write(select_2);

    EXPECT_EQ(MatrixXd(writer2.getMat()), MatrixXd(mat)({0, 2, 4}, all));
}

TEST(MatrixIO, ConcatRows) {
    SparseMatrix<double> m1 = generate_mat(3000, 10, 12512);
    SparseMatrix<double> m2 = generate_mat(1, 10, 7345); // Very few rows to try getting 0-entry columns
    SparseMatrix<double> m3 = generate_mat(256, 10, 3864);
    SparseMatrix<double> mx = generate_mat(8, 5, 92568);

    MatrixXd concat_dense(m1.rows() + m2.rows() + m3.rows(), m1.cols());
    concat_dense << MatrixXd(m1), MatrixXd(m2), MatrixXd(m3);
    SparseMatrix<double> concat = concat_dense.sparseView();

    CSparseMatrix mat_1(get_map(m1));
    CSparseMatrix mat_2(get_map(m2));
    CSparseMatrix mat_3(get_map(m3));
    CSparseMatrix mat_x(get_map(mx));

    EXPECT_ANY_THROW(ConcatRows<double>({&mat_1, &mat_x}));

    CSparseMatrixWriter res;
    ConcatRows<double> my_concat({&mat_1, &mat_2, &mat_3});
    res.write(my_concat);

    EXPECT_TRUE(res.getMat().isApprox(concat));
}

TEST(MatrixIO, ConcatCols) {
    SparseMatrix<double> m1 = generate_mat(10, 3000, 12512);
    SparseMatrix<double> m2 = generate_mat(10, 1, 7345); // Very few rows to try getting 0-entry columns
    SparseMatrix<double> m3 = generate_mat(10, 256, 3864);
    SparseMatrix<double> mx = generate_mat(5, 8, 92568);

    MatrixXd concat_dense(m1.rows(), m1.cols() + m2.cols() + m3.cols());
    concat_dense << MatrixXd(m1), MatrixXd(m2), MatrixXd(m3);
    SparseMatrix<double> concat = concat_dense.sparseView();

    CSparseMatrix mat_1(get_map(m1));
    CSparseMatrix mat_2(get_map(m2));
    CSparseMatrix mat_3(get_map(m3));
    CSparseMatrix mat_x(get_map(mx));

    EXPECT_ANY_THROW(ConcatCols<double>({&mat_1, &mat_x}));

    CSparseMatrixWriter res;
    ConcatCols<double> my_concat({&mat_1, &mat_2, &mat_3});
    res.write(my_concat);

    EXPECT_TRUE(res.getMat().isApprox(concat));
}

void test_order_rows(SparseMatrix<double> m, uint32_t load_size) {
    std::vector<uint32_t> col(m.outerIndexPtr(), m.outerIndexPtr() + m.cols() + 1);
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
    StoredMatrix<double> unordered(
        UIntReader(std::make_unique<VecUIntReader>(row.data(), row.size()), 16),
        DoubleReader(std::make_unique<VecNumReader<double>>(val.data(), val.size()), 16),
        UIntReader(std::make_unique<VecUIntReader>(col.data(), col.size()), 16),
        row_count,
        std::make_unique<VecStringReader>(std::vector<std::string>()),
        std::make_unique<VecStringReader>(std::vector<std::string>())
    );

    OrderRows<double> ordered(unordered, load_size);

    StoredMatrix<double> orig(
        UIntReader(std::make_unique<VecUIntReader>(row_orig.data(), row_orig.size()), 16),
        DoubleReader(std::make_unique<VecNumReader<double>>(val_orig.data(), val_orig.size()), 16),
        UIntReader(std::make_unique<VecUIntReader>(col.data(), col.size()), 16),
        row_count,
        std::make_unique<VecStringReader>(std::vector<std::string>()),
        std::make_unique<VecStringReader>(std::vector<std::string>())
    );
    
    ASSERT_FALSE(matrix_identical(unordered, orig));
    ASSERT_TRUE(matrix_identical(ordered, orig));
}

TEST(MatrixIO, OrderRows) {
    test_order_rows(generate_mat(5, 5), 100);
    test_order_rows(generate_mat(1000, 100), 16);
    test_order_rows(generate_mat(1000, 100), 1024);
    
}