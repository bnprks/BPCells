// Copyright 2021 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <array>
#include <atomic>
#include <vector>

#include "../matrixIterators/OrderRows.h"
#include "../matrixIterators/MatrixIterator.h"

namespace BPCells {

class TransformFit {
  public:
    Eigen::ArrayXXd row_params;
    Eigen::ArrayXXd col_params;
    Eigen::ArrayXd global_params;
};

class MatrixTransform : public MatrixLoaderWrapper<double> {
  protected:
    TransformFit fit;

  public:
    enum class RecalculateFit { None, Rows, Cols };

    MatrixTransform(std::unique_ptr<MatrixLoader<double>> &&loader);
    MatrixTransform(std::unique_ptr<MatrixLoader<double>> &&loader, TransformFit fit);
    // Constructor argument conventions:
    // MatrixTransform(MatrixLoader<double> &mat): Fit a transform, then iterate over mat
    // MatrixTransform(MatrixLoader<double> &mat, TransformFit fit, RecalculateFit recalculate =
    // RecalculateFit::None)
    //     Use the given transform parameters on the matrix, while optionally re-calculating
    //     transform parameters for either none, rows, or columns. (rows or columns could
    //     be useful for re-projecting matrices to match a given normalization)
    virtual ~MatrixTransform() = default;

    // Return the fit object calculated at construction time
    TransformFit getFit();
};

// Framework for implementing dense matrix transformations. (i.e. at least some
// zero values are transformed into non-zero values)
//
// This class implements the basic loading primitives while inserting the newly
// non-zero values into the data stream. Note that these primitives will likely have much
// worse performance on dense matrices, but are provided for ease-of-use.
//
// To provide efficient implementations of both loading and matrix/vector products,
// Child classes must implement both loadZero and loadZeroSubtracted. These functions
// involve calculating transformed values as if the underlying sparse matrix was all zeros.
class MatrixTransformDense : public MatrixTransform {
  private:
    static inline const uint32_t buf_size = 1024;
    std::array<double, buf_size> val_data;
    std::array<uint32_t, buf_size> row_data;
    uint32_t loader_idx = UINT32_MAX; // Index of loader data for this->load() output
    uint32_t loader_capacity =
        0; // Capacity of loader data. loader_capacity==0 signals no more data for this column
    uint32_t loader_col = UINT32_MAX;
    uint32_t current_row = UINT32_MAX;
    uint32_t current_col = UINT32_MAX;

  public:
    MatrixTransformDense(std::unique_ptr<MatrixLoader<double>> &&mat, TransformFit fit);

    // Reset the iterator to start from the beginning
    void restart() override;

    // Seek to a specific column without reading data
    // Next call should be to load(). col must be < cols()
    void seekCol(uint32_t col) override;

    // Advance to the next column, or return false if there
    // are no more columns
    bool nextCol() override;

    // Return the index of the current column
    uint32_t currentCol() const override;

    bool load() override;
    // Number of loaded entries available
    uint32_t capacity() const override;

    // Pointers to the loaded entries
    uint32_t *rowData() override;
    double *valData() override;

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

    // Calculate row/column sums of the matrix
    std::vector<double> colSums(std::atomic<bool> *user_interrupt = NULL) override;
    std::vector<double> rowSums(std::atomic<bool> *user_interrupt = NULL) override;

  protected:
    // Perform a normal load from the underlying matrix, then subtract transform(0)
    // from each entry and return false if there are no more non-zero values to load
    // from the underlying matrix
    virtual bool loadZeroSubtracted(MatrixLoader<double> &loader) = 0;
    // Load a range of transform(0) values into an output vector.
    // The output values should represent rows `start_row` to `start_row + count`,
    // and come from column `col`
    virtual void loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) = 0;

    // Perform matrix-matrix or matrix-vector products as if all
    // entries in the underlying matrix were 0, adding the results into
    // an existing allocated dense matrix
    virtual void denseMultiplyRightZero(
        Eigen::MatrixXd &out,
        const Eigen::Map<Eigen::MatrixXd> B,
        std::atomic<bool> *user_interrupt = NULL
    );
    virtual void denseMultiplyLeftZero(
        Eigen::MatrixXd &out,
        const Eigen::Map<Eigen::MatrixXd> B,
        std::atomic<bool> *user_interrupt = NULL
    );

    virtual void vecMultiplyRightZero(
        Eigen::VectorXd &out,
        const Eigen::Map<Eigen::VectorXd> v,
        std::atomic<bool> *user_interrupt = NULL
    );
    virtual void vecMultiplyLeftZero(
        Eigen::VectorXd &out,
        const Eigen::Map<Eigen::VectorXd> v,
        std::atomic<bool> *user_interrupt = NULL
    );
};

} // end namespace BPCells