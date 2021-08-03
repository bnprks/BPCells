#pragma once
#include <vector>

#include "MatrixIterator.h"
#include <RcppEigen.h>

namespace BPCells {

// Calculate matrix multiply A*B where B is a dense matrix
Eigen::MatrixXd denseMultiplyRight(MatrixIterator<double> &A, const Eigen::Map<Eigen::MatrixXd> B, 
                                   bool transpose_A = false);

// Calculate matrix multiply B*A where B is a dense matrix
Eigen::MatrixXd denseMultiplyLeft(MatrixIterator<double> &A, const Eigen::Map<Eigen::MatrixXd> B,
                                  bool transpose_A = false);

// Calculate vector multiply Av where v is a dense (column) vector
Eigen::VectorXd vecMultiplyRight(MatrixIterator<double> &A, const Eigen::Map<Eigen::VectorXd> v,
                                 bool transpose_A = false);

// Calculate vector multiply vA where v is a dense (row) vector
Eigen::VectorXd vecMultiplyLeft(MatrixIterator<double> &A, const Eigen::Map<Eigen::VectorXd> v,
                                 bool transpose_A = false);

template<typename T>
std::vector<T> colSums(MatrixIterator<T> &m) {
    std::vector<T> sums(m.cols());
    while (m.nextCol()) {
        while(m.nextValue()) {
            sums[m.col()] += m.val();
        }
    }
    return sums;
}

template<typename T>
std::vector<T> rowSums(MatrixIterator<T> &m) {
    std::vector<T> sums(m.rows());

    while (m.nextCol()) {
        while(m.nextValue()) {
            sums[m.row()] += m.val();
        }
    }
    return sums;
}

} // end namespace BPCells