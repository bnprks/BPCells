#pragma once

#include "MatrixTransform.h"

namespace BPCells {

// Shift rows a matrix
// out[i,j] = in[i,j] + row_params[0,i]
class ShiftRows : public MatrixTransformDense {
  public:
    ShiftRows(MatrixLoader<double> &mat, TransformFit fit);

    bool loadZeroSubtracted() override;
    void loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) override;

    Eigen::MatrixXd denseMultiplyRight(
        const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void) = NULL
    ) override;
    Eigen::MatrixXd denseMultiplyLeft(
        const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void) = NULL
    ) override;
    // Calculate matrix-vector product A*v where A=this and B is a dense matrix.
    Eigen::VectorXd vecMultiplyRight(
        const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void) = NULL
    ) override;
    Eigen::VectorXd vecMultiplyLeft(
        const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void) = NULL
    ) override;
};

class ShiftCols : public MatrixTransformDense {
  public:
    ShiftCols(MatrixLoader<double> &mat, TransformFit fit);

    bool loadZeroSubtracted() override;
    void loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) override;

    Eigen::MatrixXd denseMultiplyRight(
        const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void) = NULL
    ) override;
    Eigen::MatrixXd denseMultiplyLeft(
        const Eigen::Map<Eigen::MatrixXd> B, void (*checkInterrupt)(void) = NULL
    ) override;
    // Calculate matrix-vector product A*v where A=this and B is a dense matrix.
    Eigen::VectorXd vecMultiplyRight(
        const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void) = NULL
    ) override;
    Eigen::VectorXd vecMultiplyLeft(
        const Eigen::Map<Eigen::VectorXd> v, void (*checkInterrupt)(void) = NULL
    ) override;
};

} // end namespace BPCells