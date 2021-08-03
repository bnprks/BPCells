#include "MatrixStats.h"

namespace BPCells {

inline void update_mean(double val, double &count, double &running_mean) {
    count += 1;
    running_mean += (val - running_mean) / count;
}

inline void update_mean_and_variance(double val, double &count, double &running_mean, double &running_variance) {
    count += 1;
    double delta = val - running_mean;
    running_mean += delta / count;
    double delta2 = val - running_mean;
    running_variance += delta * delta2;
}

Eigen::ArrayXXd StatsResult::rowNonzeros() {
    if (row_stats.rows() < 1) throw std::runtime_error("Nonzeros not calculated in this StatsResult");
    return row_stats.row(0);
};
Eigen::ArrayXXd StatsResult::rowMean() {
    if (row_stats.rows() < 2) throw std::runtime_error("Nonzeros not calculated in this StatsResult");
    return row_stats.row(1);
};
Eigen::ArrayXXd StatsResult::rowVariance() {
    if (row_stats.rows() < 3) throw std::runtime_error("Nonzeros not calculated in this StatsResult");
    return row_stats.row(2);
};

Eigen::ArrayXXd StatsResult::colNonzeros() {
    if (col_stats.rows() < 1) throw std::runtime_error("Nonzeros not calculated in this StatsResult");
    return col_stats.row(0);
};
Eigen::ArrayXXd StatsResult::colMean() {
    if (col_stats.rows() < 2) throw std::runtime_error("Nonzeros not calculated in this StatsResult");
    return col_stats.row(1);
};
Eigen::ArrayXXd StatsResult::colVariance() {
    if (col_stats.rows() < 3) throw std::runtime_error("Nonzeros not calculated in this StatsResult");
    return col_stats.row(2);
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
                        bool transpose, uint32_t buffer_size) {
    if (transpose) {
        StatsResult res = computeMatrixStats(mat, col_stats, row_stats, false, buffer_size);
        return StatsResult {
            res.col_stats,
            res.row_stats
        };
    }

    StatsResult res {
        Eigen::ArrayXXd((int) row_stats, mat.rows()),
        Eigen::ArrayXXd((int) col_stats, mat.cols())
    };

    res.row_stats.setZero();
    res.col_stats.setZero();

    std::vector<uint32_t> row_buf(buffer_size);
    std::vector<double> val_buf(buffer_size);

    SparseVector<double> buf;
    buf.idx = &row_buf[0];
    buf.val = &val_buf[0];
    buf.capacity = buffer_size;

    uint32_t current_col;
    uint32_t num_loaded;
    while (mat.nextCol()) {
        current_col = mat.currentCol();

        while ( (num_loaded = mat.load(buffer_size, buf)) ) {
            // Update row stats if needed
            if (row_stats == Stats::NonZeroCount) {
                for (size_t i = 0; i < num_loaded; i++) {
                    res.row_stats(0, row_buf[i]) += 1;
                }
            } else if (row_stats == Stats::Mean) {
                for (size_t i = 0; i < num_loaded; i++) {
                    update_mean(
                        val_buf[i], 
                        res.row_stats(0, row_buf[i]), 
                        res.row_stats(1, row_buf[i])
                    );
                }
            } else if (row_stats == Stats::Variance) {
                for (size_t i = 0; i < num_loaded; i++) {
                    update_mean_and_variance(
                        val_buf[i], 
                        res.row_stats(0, row_buf[i]), 
                        res.row_stats(1, row_buf[i]), 
                        res.row_stats(2, row_buf[i])
                    );
                }
            }
            // Update col stats if needed
            if (col_stats == Stats::NonZeroCount) {
                for (size_t i = 0; i < num_loaded; i++) {
                    res.col_stats(0, current_col) += 1;
                }
            } else if(col_stats == Stats::Mean) {
                for (size_t i = 0; i < num_loaded; i++) {
                    update_mean(
                        val_buf[i], 
                        res.col_stats(0, current_col),
                        res.col_stats(1, current_col)
                    );
                }
            } else if (col_stats == Stats::Variance) {
                for (size_t i = 0; i < num_loaded; i++) {
                    update_mean_and_variance(
                        val_buf[i], 
                        res.col_stats(0, current_col), 
                        res.col_stats(1, current_col),
                        res.col_stats(2, current_col)
                    );
                }
            }
        }
    }

    // fix up mean + variance accounting for the zero entries
    if (col_stats == Stats::Mean) {
        auto nnz = res.col_stats.row(0);
        auto mean = res.col_stats.row(1);
        res.col_stats.row(1) = mean * nnz / mat.rows();
    } else if (col_stats == Stats::Variance) {
        // m2 += mean^2 * nnz*(nrow - nnz)/nrow;
        auto nnz = res.col_stats.row(0);
        auto mean = res.col_stats.row(1);
        auto m2 = res.col_stats.row(2);
        res.col_stats.row(2) = (m2 + mean.square() * nnz * ((double) mat.rows() - nnz) / mat.rows()) /
            (mat.rows() - 1);
        res.col_stats.row(1) = mean * nnz / mat.rows();
    }

    if (row_stats == Stats::Mean) {
        auto nnz = res.row_stats.row(0);
        auto mean = res.row_stats.row(1);
        res.row_stats.row(1) = mean * nnz / mat.cols();
    } else if (row_stats == Stats::Variance) {
        // m2 += mean^2 * nnz*(ncol - nnz)/ncol;
        auto nnz = res.row_stats.row(0);
        auto mean = res.row_stats.row(1);
        auto m2 = res.row_stats.row(2);
        res.row_stats.row(2) = (m2 + mean.square() * nnz * ((double) mat.cols() - nnz) / mat.cols()) /
            (mat.cols() - 1);
        res.row_stats.row(1) = mean * nnz / mat.cols();
    }

    return res;
}


} // end namespace BPCells