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

#include <matrixTransforms/Log1p.h>
#include <matrixTransforms/Scale.h>
#include <matrixTransforms/Shift.h>

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


TEST(MatrixMath, Log1p) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    MatrixXd ans = MatrixXd(m1).array().log1p();

    CSparseMatrix mat_1(get_map(m1));
    Log1p mat_1_trans(mat_1);

    CSparseMatrixWriter res;
    res.write(mat_1_trans);

    EXPECT_TRUE(MatrixXd(res.getMat()).isApprox(ans));
}

void checkMultiplyOps(MatrixLoader<double> &mat, MatrixXd res) {
    MatrixXd mat_right = generate_dense_mat(mat.cols(), 4, 125412);
    MatrixXd mat_left = generate_dense_mat(4, mat.rows(), 8432);

    VectorXd vec_right = mat_right.col(0);
    VectorXd vec_left = mat_left.row(0).transpose();

    EXPECT_TRUE((res * mat_right).isApprox(mat.denseMultiplyRight(get_map_dense(mat_right))));
    EXPECT_TRUE((res * vec_right).isApprox(mat.vecMultiplyRight(get_map_vec(vec_right))));

    EXPECT_TRUE((mat_left * res).isApprox(mat.denseMultiplyLeft(get_map_dense(mat_left))));
    EXPECT_TRUE((vec_left.transpose() * res).transpose().isApprox(mat.vecMultiplyLeft(get_map_vec(vec_left))));
}

TEST(MatrixMath, Scale) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    CSparseMatrix mat_1(get_map(m1));
    
    VectorXd scale_row = generate_dense_vec(m1.rows(), 1513);
    VectorXd scale_col = generate_dense_vec(m1.cols(), 14582);
    
    // Scale row + col
    MatrixXd ans1 = (MatrixXd(m1).array().rowwise() * scale_col.array().transpose()).colwise() * scale_row.array();
    Scale s1(mat_1, TransformFit{scale_row.transpose(), scale_col.transpose(), {}});
    checkMultiplyOps(s1, ans1);

    CSparseMatrixWriter r1;
    s1.restart();
    r1.write(s1);
    EXPECT_TRUE(MatrixXd(r1.getMat()).isApprox(ans1));

    // Scale just row
    MatrixXd ans2 = MatrixXd(m1).array().colwise() * scale_row.array();
    Scale s2(mat_1, TransformFit{scale_row.transpose(), {}, {}});
    checkMultiplyOps(s2, ans2);

    CSparseMatrixWriter r2;
    s2.restart();
    r2.write(s2);
    EXPECT_TRUE(MatrixXd(r2.getMat()).isApprox(ans2));

    // Scale just col
    MatrixXd ans3 = MatrixXd(m1).array().rowwise() * scale_col.array().transpose();
    Scale s3(mat_1, TransformFit{{}, scale_col.transpose(), {}});
    checkMultiplyOps(s3, ans3);

    CSparseMatrixWriter r3;
    s3.restart();
    r3.write(s3);
    EXPECT_TRUE(MatrixXd(r3.getMat()).isApprox(ans3));
}

class SimpleRowShift : public MatrixTransformDense {
public:
    SimpleRowShift(MatrixLoader<double> &mat, TransformFit fit) : MatrixTransformDense(mat, fit) {}

    bool loadZeroSubtracted() override {return loader.load();}
    void loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) override {
        for (uint32_t i = 0; i < count; i++) {
            values[i] = fit.row_params(0, start_row + i);
        }
    }
};

TEST(MatrixMath, TransformDense) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    CSparseMatrix mat_1(get_map(m1));
    
    VectorXd shift_row = generate_dense_vec(m1.rows(), 1513);
    VectorXd shift_col = generate_dense_vec(m1.cols(), 1242);
    
    // Shift row
    MatrixXd ans1 = MatrixXd(m1).array().colwise() + shift_row.array();
    SimpleRowShift s1(mat_1, TransformFit{shift_row.transpose(), {}, {}});
    
    checkMultiplyOps(s1, ans1);
    
    CSparseMatrixWriter r1;
    s1.restart();
    r1.write(s1);
    
    EXPECT_TRUE(MatrixXd(r1.getMat()).isApprox(ans1));
}

TEST(MatrixMath, Shift) {
    SparseMatrix<double> m1 = generate_mat(100, 50, 125123);
    CSparseMatrix mat_1(get_map(m1));
    
    VectorXd shift_row = generate_dense_vec(m1.rows(), 1513);
    VectorXd shift_col = generate_dense_vec(m1.cols(), 1242);
    
    // Shift row
    MatrixXd ans1 = MatrixXd(m1).array().colwise() + shift_row.array();
    ShiftRows s1(mat_1, TransformFit{shift_row.transpose(), {}, {}});
    
    checkMultiplyOps(s1, ans1);
    
    CSparseMatrixWriter r1;
    s1.restart();
    r1.write(s1);
    
    EXPECT_TRUE(MatrixXd(r1.getMat()).isApprox(ans1));

    // Shift col
    MatrixXd ans2 = MatrixXd(m1).array().rowwise() + shift_col.transpose().array();
    ShiftCols s2(mat_1, TransformFit{{}, shift_col.transpose(), {}});
    
    checkMultiplyOps(s2, ans2);
    
    CSparseMatrixWriter r2;
    s2.restart();
    r2.write(s2);
    
    EXPECT_TRUE(MatrixXd(r2.getMat()).isApprox(ans2));
}