#pragma once

#include "MatrixTransform.h"

namespace BPCells {

// Pearson residuals of an SCTransform model
// First row of row params is 1/theta
// mu = exp(t(row_params[1:,:]) * col_params[:,:])
// Transform = (X - mu) / sqrt(mu + mu^2/theta)
class SCTransformPearson : public MatrixTransformDense {
  private:
    //uint32_t cached_col = UINT32_MAX;
    //Eigen::ArrayXf col_mu; // Cached values of mu for current column
    Eigen::ArrayXf theta_inv;
    Eigen::ArrayXXf col_mat;
    Eigen::ArrayXXf row_mat;
    Eigen::Array<float, 2048, 1> mu_tmp;

    //void ensure_cached_mu(uint32_t col);
  public:
    SCTransformPearson(MatrixLoader<double> &loader, TransformFit fit);

    // On seekCol or nextCol, update col_mu in a single operation
    // bool nextCol() override;
    // void seekCol(uint32_t col) override;

    bool loadZeroSubtracted(MatrixLoader<double> &loader) override;
    void loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) override;
};

} // end namespace BPCells