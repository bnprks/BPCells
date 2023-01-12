#pragma once

#include "MatrixTransform.h"

namespace BPCells {

// Pearson residuals of an SCTransform model
// First row of row params is 1/theta
// mu = exp(t(row_params[1:,:]) * col_params[:,:])
// Transform = (X - mu) / sqrt(mu + mu^2/theta)
class SCTransformPearson : public MatrixTransformDense {
  private:
    //Eigen::VectorXd col_mu; // Cached values of mu for current column
  public:
    using MatrixTransformDense::MatrixTransformDense;

    // On seekCol or nextCol, update col_mu in a single operation
    //bool nextCol() override;
    //void seekCol(uint32_t col) override;

    bool loadZeroSubtracted(MatrixLoader<double> &loader) override;
    void loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) override;
};


} // end namespace BPCells