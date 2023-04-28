#pragma once
#include <atomic>
#include "MatrixTransform.h"

namespace BPCells {

// Scale rows and/or columns of a matrix
// out[i,j] = in[i,j] * row_params[0,i] * col_params[0,j]
// If row_params or col_params have 0 rows, then skip scaling along that dimension
class Scale : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;

    Eigen::MatrixXd denseMultiplyRight(
        const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt = NULL
    ) override;
    Eigen::MatrixXd denseMultiplyLeft(
        const Eigen::Map<Eigen::MatrixXd> B, std::atomic<bool> *user_interrupt = NULL
    ) override;
    // Calculate matrix-vector product A*v where A (this) is sparse and B is a dense matrix.
    Eigen::VectorXd vecMultiplyRight(
        const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt = NULL
    ) override;
    Eigen::VectorXd vecMultiplyLeft(
        const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt = NULL
    ) override;
};

} // end namespace BPCells