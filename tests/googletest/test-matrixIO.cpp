#include <filesystem>
#include <random>
#include <sstream>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <matrixIterators/MatrixIterator.h>
#include <matrixIterators/CSparseMatrix.h>
#include <matrixIterators/StoredMatrix.h>
#include <matrixIterators/MatrixIndexSelect.h>
#include <arrayIO/vector.h>

#include <Eigen/Core>

namespace fs = std::filesystem;
using namespace BPCells;
using namespace ::testing;
using namespace Eigen;

SparseMatrix<double> generate_mat(uint32_t n_row, uint32_t n_col) {
    std::random_device rd;
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
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
    
    SparseMatrix<double> mat(n_row,n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

Map<SparseMatrix<double>> get_map(const SparseMatrix<double> &mat) {
   return Map<SparseMatrix<double>>(
        mat.rows(), 
        mat.cols(), 
        mat.nonZeros(), 
        (int *) mat.outerIndexPtr(), 
        (int *) mat.innerIndexPtr(), 
        (double *) mat.valuePtr()
    );
}

TEST(MatrixIO, UnpackedVec) {
    const Eigen::SparseMatrix<double> orig_mat = generate_mat(10, 10);
    CSparseMatrix mat_d (get_map(orig_mat));

    MatrixConverterLoader<double, uint32_t> mat_i(mat_d);
    MatrixIterator it1(mat_i);

    VecReaderWriterBuilder vb1(1024);
    auto w1 = StoredMatrixWriter::createUnpacked(vb1);
    w1.write(it1);

    auto loader = StoredMatrix<uint32_t>::openUnpacked(vb1);
    auto loader_double = MatrixConverterLoader<uint32_t, double>(loader);
    
    CSparseMatrixWriter w2;
    MatrixIterator it2(loader_double);
    w2.write(it2);
    EXPECT_TRUE(w2.getMat().isApprox(orig_mat));
}

TEST(MatrixIO, PackedVec) {
    const Eigen::SparseMatrix<double> orig_mat = generate_mat(10, 10);
    CSparseMatrix mat_d (get_map(orig_mat));

    MatrixConverterLoader<double, uint32_t> mat_i(mat_d);
    MatrixIterator it1(mat_i);

    VecReaderWriterBuilder vb1(1024);
    auto w1 = StoredMatrixWriter::createPacked(vb1);
    w1.write(it1);

    auto loader = StoredMatrix<uint32_t>::openPacked(vb1, 1024);
    auto loader_double = MatrixConverterLoader<uint32_t, double>(loader);
    
    CSparseMatrixWriter w2;
    MatrixIterator it2(loader_double);
    w2.write(it2);
    EXPECT_TRUE(w2.getMat().isApprox(orig_mat));
}

TEST(MatrixIO, SeekCSparse) {
    std::vector<Triplet<double>> triplets;
    const uint32_t n_row = 6;
    const uint32_t n_col = n_row - 1;

    // Matrix where each column j starts with an entry of j at row j+1
    for (int j = 0; j < n_col; j++) {
        for (int i = j+1; i < n_row; i++) {
            triplets.push_back({i, j, j});
        }
    }
    SparseMatrix<double> mat(n_row,n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());

    CSparseMatrix mat_l(get_map(mat));
    MatrixIterator it(mat_l);
    for (auto j : {4,1,3,0,2}) {
        it.seekCol(j);
        EXPECT_TRUE(it.nextValue());
        EXPECT_EQ(it.row(), j+1);
        EXPECT_EQ(it.val(), j);
    }
}

TEST(MatrixIO, SeekStoredVec) {
    std::vector<Triplet<double>> triplets;
    const uint32_t n_row = 6;
    const uint32_t n_col = n_row - 1;

    // Matrix where each column j starts with an entry of j at row j+1
    for (int j = 0; j < n_col; j++) {
        for (int i = j+1; i < n_row; i++) {
            triplets.push_back({i, j, j});
        }
    }
    SparseMatrix<double> mat(n_row,n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());

    CSparseMatrix mat_double(get_map(mat));
    MatrixConverterLoader<double, uint32_t> mat_int(mat_double);
    MatrixIterator it1(mat_int);

    VecReaderWriterBuilder vb(1024);
    auto w = StoredMatrixWriter::createUnpacked(vb);
    w.write(it1);

    auto loader = StoredMatrix<uint32_t>::openUnpacked(vb);
    auto loader_double = MatrixConverterLoader<uint32_t, double>(loader);

    
    MatrixIterator it(loader_double);
    for (auto j : {4,1,3,0,2}) {
        it.seekCol(j);
        EXPECT_TRUE(it.nextValue());
        EXPECT_EQ(it.row(), j+1);
        EXPECT_EQ(it.val(), j);
    }
}


TEST(MatrixIO, ColSelectCSparse) {
    std::vector<Triplet<double>> triplets;
    const uint32_t n_row = 6;
    const uint32_t n_col = n_row - 1;

    for (int j = 0; j < n_col; j++) {
        for (int i = 0; i < n_row; i++) {
            triplets.push_back({i, j, j+n_col*i});
        }
    }
    
    SparseMatrix<double> mat(n_row,n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    
    CSparseMatrix mat_l(get_map(mat));
    MatrixColSelect<double> mat_col_select(mat_l, {0,4,2});
    MatrixIterator it(mat_col_select);

    CSparseMatrixWriter writer;
    writer.write(it);

    // Test that the matrix is the same
    EXPECT_EQ(MatrixXd(writer.getMat()), MatrixXd(mat)(all, {0,4,2}));

    // Check that seeking columns works
    for (auto j : {2,0,1}) {
        it.seekCol(j);
        EXPECT_TRUE(it.nextValue());
        EXPECT_EQ(it.row(), 0);
        EXPECT_EQ(it.val(), std::vector<uint32_t>({0,4,2})[j]);
    }
}

TEST(MatrixIO, RowSelectCSparse) {
    std::vector<Triplet<double>> triplets;
    const uint32_t n_row = 6;
    const uint32_t n_col = n_row - 1;

    for (int j = 0; j < n_col; j++) {
        for (int i = 0; i < n_row; i++) {
            triplets.push_back({i, j, j+n_col*i});
        }
    }
    
    SparseMatrix<double> mat(n_row,n_col);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    
    CSparseMatrix mat_l(get_map(mat));
    MatrixRowSelect<double> select_1(mat_l, {0,4,2});
    MatrixIterator it(select_1);

    CSparseMatrixWriter writer1;
    writer1.write(it);

    EXPECT_EQ(MatrixXd(writer1.getMat()), MatrixXd(mat)({0,4,2}, all));

    CSparseMatrix mat_l_2(get_map(mat));
    MatrixRowSelect<double> select_2(mat_l_2, {0,2,4});
    MatrixIterator it2(select_2);

    CSparseMatrixWriter writer2;
    writer2.write(it2);

    EXPECT_EQ(MatrixXd(writer2.getMat()), MatrixXd(mat)({0,2,4}, all));
}