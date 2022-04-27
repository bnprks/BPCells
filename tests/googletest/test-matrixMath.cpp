#include <filesystem>
#include <random>
#include <sstream>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <matrixIterators/MatrixIterator.h>
#include <matrixIterators/ConcatenateMatrix.h>
#include <matrixIterators/CSparseMatrix.h>
#include <matrixIterators/StoredMatrix.h>
#include <matrixIterators/MatrixIndexSelect.h>
#include <arrayIO/vector.h>

#include <Eigen/Core>

namespace fs = std::filesystem;
using namespace BPCells;
using namespace ::testing;
using namespace Eigen;

MatrixXd generate_dense_mat(uint32_t n_row, uint32_t n_col, uint32_t seed=125124) {
    std::mt19937 gen(seed); //Standard mersenne_twister_engine
    std::uniform_int_distribution<> distrib(1, 20);
    
    MatrixXd mat(n_row, n_col);
    for (int i = 0; i < n_row; i++) {
        for (int j = 0; j < n_col; j++) {
            mat(i, j) = distrib(gen);
        }
    }
    return mat;
}

MatrixXd generate_dense_vec(uint32_t n, uint32_t seed=125124) {
    std::mt19937 gen(seed); //Standard mersenne_twister_engine
    std::uniform_int_distribution<> distrib(1, 20);
    
    VectorXd vec(n);
    for (int i = 0; i < n; i++) {
            vec(i) = distrib(gen);
    }
    return vec;
}

SparseMatrix<double> generate_mat(uint32_t n_row, uint32_t n_col, uint32_t seed=125124) {
    std::mt19937 gen(seed); //Standard mersenne_twister_engine
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

Map<MatrixXd> get_map_dense(MatrixXd &mat) {
    return Map<MatrixXd>(mat.data(), mat.rows(), mat.cols());
}

Map<VectorXd> get_map_vec(VectorXd &mat) {
    return Map<VectorXd>(mat.data(), mat.rows(), mat.cols());
}

TEST(MatrixMath, Multiply) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    MatrixXd b_right = generate_dense_mat(50, 4);
    MatrixXd b_left = generate_dense_mat(4, 100);
    
    CSparseMatrix mat_1(get_map(m1));

    MatrixXd r1 = mat_1.denseMultiplyLeft(get_map_dense(b_left));
    MatrixXd r2 = mat_1.denseMultiplyRight(get_map_dense(b_right));

    EXPECT_EQ(r1, b_left * m1);
    EXPECT_EQ(r2, m1 * b_right);

    VectorXd v_right = b_right.col(0);
    VectorXd v_left = b_left.row(0);

    VectorXd r3 = mat_1.vecMultiplyLeft(get_map_vec(v_left));
    VectorXd r4 = mat_1.vecMultiplyRight(get_map_vec(v_right));
    
    VectorXd ans3 = b_left({0}, all) * m1;
    EXPECT_EQ(r3, ans3);

    VectorXd ans4 = m1 * b_right.col(0);
    EXPECT_EQ(r4, ans4);

}


TEST(MatrixMath, Stats) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    //SparseMatrix<double> m1 = generate_mat(5, 3, 125123);
    
    // Getting useful represetntations of input
    MatrixXd d1(m1);
    ArrayXXd a1(d1);

    CSparseMatrix mat_1(get_map(m1));

    ArrayXXd row_stats(3, m1.rows());
    ArrayXXd col_stats(3, m1.cols());
    row_stats.row(0) = (a1 > 0).rowwise().count().cast<double>();
    row_stats.row(1) = a1.rowwise().mean();
    row_stats.row(2) = (a1.colwise() - row_stats.row(1).transpose()).square().rowwise().sum() / (a1.cols() - 1);

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