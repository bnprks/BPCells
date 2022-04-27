#pragma once


#ifndef RCPP_EIGEN
#include <Eigen/SparseCore>
#else
#include <RcppEigen.h>
#endif
// [[Rcpp::depends(RcppEigen)]]

namespace BPCells {

// Enum to request stats to calculate in MatrixLoader::computeMatrixStats
enum class Stats {
    None = 0,
    NonZeroCount = 1,
    Mean = 2,
    Variance = 3
};

// Result class for row + column stats.
// Each stat is a row in the matrix
class StatsResult {
public:
    Eigen::ArrayXXd row_stats;
    Eigen::ArrayXXd col_stats;
    
    Eigen::ArrayXXd rowNonzeros();
    Eigen::ArrayXXd rowMean();
    Eigen::ArrayXXd rowVariance();

    Eigen::ArrayXXd colNonzeros();
    Eigen::ArrayXXd colMean();
    Eigen::ArrayXXd colVariance();

    StatsResult transpose();
};


} // end namespace BPCells