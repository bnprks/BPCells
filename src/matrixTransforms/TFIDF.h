#pragma once

#include "MatrixTransform.h"

namespace BPCells {

class TFIDF : public SparseTransform<TFIDF> {
public:
    using SparseTransform::SparseTransform;

    static inline TransformFit fit(MatrixLoader<double> &mat, const Eigen::ArrayXd &global_params, bool transpose) {
        auto stats = computeMatrixStats(mat, Stats::Mean, Stats::Mean, transpose);
        auto rows = transpose ? mat.cols() : mat.rows();
        TransformFit f {
            stats.rowMean(),
            rows * stats.colMean(),
            global_params
        };
        return f;
    }

    inline double transform(double val, const double *row_params, const double *col_params, const double *global_params) {
        return std::log1p(global_params[0] * val / row_params[0] / col_params[0]);
    }
};

} // end namespace BPCells