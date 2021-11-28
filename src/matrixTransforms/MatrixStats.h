#pragma once

#include <stdexcept>
#include <vector>

#include <RcppEigen.h>

#include "../matrixIterators/MatrixIterator.h"

namespace BPCells {

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
};


// Calculate stats on the rows or columns of a matrix in a single pass.
// For each of rows and columns, the user can choose to calculate no stats,
// just means, or mean + variance. 
// Outputs results to matrices row_output and col_output, which are col-major matrices,
// with one column per # rows or # columns as appropriate, and one row per output statistic
// Available statistics:
// 1. None
// 2. NonZeroCount
// 3. Mean
// 4. Variance
// Each of the later statistics is always calculated simultaneously to the earlier statistics 
StatsResult computeMatrixStats(MatrixLoader<double> &mat, Stats row_stats, Stats col_stats,
                        bool transpose = false, uint32_t buffer_size = 1024, void (*checkInterrupt)(void) = NULL);

} // end namespace BPCells