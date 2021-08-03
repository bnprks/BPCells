#pragma once

#include <RcppEigen.h>

#include "../matrixIterators/MatrixIterator.h"
#include "../matrixIterators/MatrixOps.h"

#include "MatrixStats.h"

namespace BPCells {
    
class TransformFit {
public:
    Eigen::ArrayXXd row_params;
    Eigen::ArrayXXd col_params;
    Eigen::ArrayXd global_params;
};


class MatrixTransform {
public:    
    enum class RecalculateFit {
        None,
        Rows,
        Cols
    };
    // Constructor argument conventions:
    // MatrixTransform(MatrixLoader<double> &mat): Fit a transform, then iterate over mat
    // MatrixTransform(MatrixLoader<double> &mat, TransformFit fit, RecalculateFit recalculate = RecalculateFit::None)
    //     Use the given transform parameters on the matrix, while optionally re-calculating
    //     transform parameters for either none, rows, or columns. (rows or columns could
    //     be useful for re-projecting matrices to match a given normalization)
    virtual ~MatrixTransform() = default;

    // Multiply the transformed matrix by a dense matrix
    virtual Eigen::MatrixXd denseMultiplyRight(Eigen::Map<Eigen::MatrixXd> B) = 0;
    virtual Eigen::MatrixXd denseMultiplyLeft(Eigen::Map<Eigen::MatrixXd> B) = 0;

    // Multiply the transformed matrix by a dense vector
    virtual Eigen::VectorXd vecMultiplyRight(const Eigen::Map<Eigen::VectorXd> v) = 0;
    virtual Eigen::VectorXd vecMultiplyLeft(const Eigen::Map<Eigen::VectorXd> v) = 0;

    // Return the fit object calculated at construction time
    virtual TransformFit getFit() = 0;

    // Return mean + variance of either rows or columns.
    // This is used to make PCA easier to implement.
    // The return value will be a 2-row matrix with 1st containing means and
    // second containing variance
    using StatsResultMatrix = Eigen::Matrix<double, 2, Eigen::Dynamic>;
    virtual StatsResultMatrix rowStats(uint32_t buffer_size) = 0;
    virtual StatsResultMatrix colStats(uint32_t buffer_size) = 0;
};

// A CRTP base class designed to make for very succinct implementations
// of sparse matrix transformations.
//
// Child classes should implement the following methods:
// - static TransformFit fit(MatrixLoader<double> &mat, const Eigen::ArrayXd &global_params)
//       Fit all row+col parameters on a given matrix, returning the result.
//       Global parameters are listed in the input if required, and will be copied
//       to the TransformFit object automatically
//
// - double transform(double val, const double *row_params, const double *col_params, const double *global_params)
//       Transform a single value, given the relevant parameters for its row and column
//
// There are just two minor caveats:
// - Some work will be wasted in fitting transforms for re-projection
// - The global parameters must not be calculated and are passed in by the user
//   then saved exactly
//
// See the LSI imlementation for an example of how to implement a sparse transform
template<typename ConcreteTransform>
class SparseTransform : public MatrixTransform, public MatrixLoaderWrapper<double> {
private:
    TransformFit fit;
    bool transpose;
    // Transform a value given the row and column indices
    inline double transform(double val, uint32_t row, uint32_t col) {
        double *row_params = &fit.row_params(0, row);
        double *col_params = &fit.col_params(0, col);
        return static_cast<ConcreteTransform*>(this)->transform(
            val, row_params, col_params, &fit.global_params(0)
        );
    }
public:    
    // Must provide at least a matrix to construct
    SparseTransform() = delete;
    
    // Fit a new transform to the matrix
    SparseTransform(MatrixLoader<double> &mat, const Eigen::ArrayXd &params, bool transpose=false) : 
        MatrixLoaderWrapper(mat), fit(ConcreteTransform::fit(mat, params, transpose)), transpose(transpose) {};

    SparseTransform(MatrixLoader<double> &mat,
                    TransformFit fit, MatrixTransform::RecalculateFit recalculate,
                    bool transpose) :
        MatrixLoaderWrapper(mat), transpose(transpose) {
        this->fit = fit;
        if (recalculate != MatrixTransform::RecalculateFit::None) {
            TransformFit new_fit = ConcreteTransform::fit(mat, fit.global_params, transpose);
            if(recalculate == MatrixTransform::RecalculateFit::Rows) {
                this->fit.row_params = new_fit.row_params;
            }
            if (recalculate == MatrixTransform::RecalculateFit::Cols) {
                this->fit.col_params = new_fit.col_params;
            }
        }
        if (!transpose) {
            assert(this->fit.row_params.cols() == loader.rows());
            assert(this->fit.col_params.cols() == loader.cols());
        } else {
            assert(this->fit.row_params.cols() == loader.cols());
            assert(this->fit.col_params.cols() == loader.rows());
        }
    }

    // Multiply the transformed matrix by a dense matrix
    Eigen::MatrixXd denseMultiplyRight(Eigen::Map<Eigen::MatrixXd> B) override {
        loader.restart();
        MatrixIterator<double> A(*this);
        return BPCells::denseMultiplyRight(A, B, transpose);
    };
    Eigen::MatrixXd denseMultiplyLeft(Eigen::Map<Eigen::MatrixXd> B) override {
        loader.restart();
        MatrixIterator<double> A(*this);
        return BPCells::denseMultiplyLeft(A, B, transpose);
    };

    // Multiply the transformed matrix by a dense vector
    Eigen::VectorXd vecMultiplyRight(const Eigen::Map<Eigen::VectorXd> v) override {
        loader.restart();
        MatrixIterator<double> A(*this);
        if (transpose) 
            return BPCells::vecMultiplyLeft(A, v);
        else
            return BPCells::vecMultiplyRight(A, v);
    };
    Eigen::VectorXd vecMultiplyLeft(const Eigen::Map<Eigen::VectorXd> v) override {
        loader.restart();
        MatrixIterator<double> A(*this);
        if (transpose) 
            return BPCells::vecMultiplyRight(A, v);
        else
            return BPCells::vecMultiplyLeft(A, v);
    };

    // Return the fit object calculated at construction time
    TransformFit getFit() override {return fit;};

    int32_t load(uint32_t count, SparseVector<double> buffer) override {
        int32_t ret = loader.load(count, buffer);
        if (transpose) {
            for (int i = 0; i < ret; i++) {
                buffer.val[i] = transform(buffer.val[i], current_col, buffer.idx[i]);
            }
        } else {
            for (int i = 0; i < ret; i++) {
                buffer.val[i] = transform(buffer.val[i], buffer.idx[i], current_col);
            }
        }
        return ret;
    };

    // Return mean + variance of either rows or columns.
    // This is used to make PCA easier to implement
    MatrixTransform::StatsResultMatrix rowStats(uint32_t buffer_size = 1024) override {
        loader.restart();
        
        StatsResult stats = computeMatrixStats(*this, Stats::Variance, Stats::None, transpose, buffer_size);
        
        return stats.row_stats.bottomRows<2>();
    };
    MatrixTransform::StatsResultMatrix colStats(uint32_t buffer_size = 1024) override {
        loader.restart();

        StatsResult stats = computeMatrixStats(*this, Stats::None, Stats::Variance, transpose, buffer_size);

        return stats.col_stats.bottomRows<2>();
    };
};



} // end namespace BPCells