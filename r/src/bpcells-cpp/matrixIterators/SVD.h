#pragma once

#include <Eigen/Core>

#include "MatrixIterator.h"

namespace BPCells {

namespace SVD {
class UserInterruptException : std::exception {
  public:
    char const *what() const noexcept override { return "Terminated at user request"; }
};
} // namespace SVD

class SpectraMatOp {
  private:
    MatrixLoader<double> *mat;
    std::atomic<bool> *user_interrupt;
    const bool tall; // If true, calculate t(A)*A; else calculate A*t(A)

  public:
    SpectraMatOp(MatrixLoader<double> *mat, std::atomic<bool> *user_interrupt);
    using Scalar = double;
    Eigen::Index rows() const;
    Eigen::Index cols() const;
    // y_out = M * x_in
    void perform_op(const double *x_in, double *y_out) const;
};

class SVDResult {
  public:
    Eigen::VectorXd d; 
    Eigen::MatrixXd u, v;
    int num_iterations, num_operations, num_converged; // Number of iters, multiply ops, and converged singular values
    bool success = false;
};

// Parameters: 
// mat - matrix to load 
// k - number of components to calculate
// n_cv - convergence speed parameter. Higher values use more memory, and more
//        operations per iteration, but converge faster
// maxit - maximum iterations
// transpose - if true, treat mat as transposed
// tol - precision for eigenvalues
SVDResult
svd(MatrixLoader<double> *mat,
    int k,
    int n_cv,
    int maxit,
    double tol,
    bool transpose,
    std::atomic<bool> *user_interrupt);

} // end namespace BPCells