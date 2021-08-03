#include "MatrixOps.h"

namespace BPCells {

// Calculate matrix multiply A*B where B is a dense matrix
Eigen::MatrixXd denseMultiplyRight(MatrixIterator<double> &A, const Eigen::Map<Eigen::MatrixXd> B,
                                   bool transpose_A) {
    // Use transposed output so that results write contiguously in memory.
    // Note that left multiply will have better performance properties due
    // to better cache locality on the input matrix.
    // So if the performance is key and the size of the input is small,
    // it may be worth paying the price of additional transpose operations
    Eigen::Index output_rows = transpose_A ? A.cols() : A.rows();
    Eigen::MatrixXd res(B.cols(), output_rows);
    res.setZero();
    
    if (!transpose_A) {
        while (A.nextCol()) {
            while (A.nextValue()) {
                res.col(A.row()) += A.val() * B.row(A.col());
            }
        }
    } else { // Output transposed version
        while (A.nextCol()) {
            while (A.nextValue()) {
                res.col(A.col()) += A.val() * B.row(A.row());
            }
        }
        
    }
    return res.transpose();
}

// Calculate matrix multiply B*A where B is a dense matrix
Eigen::MatrixXd denseMultiplyLeft(MatrixIterator<double> &A, const Eigen::Map<Eigen::MatrixXd> B,
                                  bool transpose_A) {
    Eigen::Index output_cols = transpose_A ? A.rows() : A.cols();
    Eigen::MatrixXd res(B.rows(), output_cols);
    res.setZero();
    
    if (!transpose_A) {
        while (A.nextCol()) {
            while (A.nextValue()) {
                res.col(A.col()) += A.val() * B.col(A.row());
            }
        }
    } else {
        while (A.nextCol()) {
            while (A.nextValue()) {
                res.col(A.row()) += A.val() * B.col(A.col());
            }
        }
    }
    return res;
}

// Calculate vector multiply Av where v is a dense (column) vector
Eigen::VectorXd vecMultiplyRight(MatrixIterator<double> &A, const Eigen::Map<Eigen::VectorXd> v,
                                bool transpose_A) {                                
    if (transpose_A) return vecMultiplyLeft(A, v, false);
    
    Eigen::VectorXd res(A.rows());
    res.setZero();
    while (A.nextCol()) {
        while (A.nextValue()) {
            res(A.row()) += A.val() * v(A.col());
        }
    }
    return res;
}

// Calculate vector multiply vA where v is a dense (row) vector
Eigen::VectorXd vecMultiplyLeft(MatrixIterator<double> &A, const Eigen::Map<Eigen::VectorXd> v,
                                bool transpose_A) {
    if (transpose_A) return vecMultiplyRight(A, v, false);
    Eigen::VectorXd res(A.cols());
    res.setZero();
    while (A.nextCol()) {
        while (A.nextValue()) {
            res(A.col()) += A.val() * v(A.row());
        }
    }
    return res;
}

} // end namespace BPCells