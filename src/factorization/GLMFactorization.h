#include <algorithm>
#include <memory>
#include <vector>

#ifndef RCPP_EIGEN
#include <Eigen/Core>
#else
#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <RcppEigen.h>
#endif
// [[Rcpp::depends(RcppEigen)]]

#include "../lib/sleef_wrapper.h"
#include "../matrixIterators/MatrixIterator.h"

using namespace Eigen;

namespace BPCells::poisson2 {

// Statistics for fitting a single Poisson GLM vector
class FitStats {
  public:
    int iterations = 0;          // Number of iterations
    int alpha_lower = 0;         // Number of times alpha is lowered
    int alpha_raise = 0;         // Number of times alpha is raised
    double min_alpha = 1;        // Minimum alpha seen during backtracking
    int compute_derivatives = 0; // Times computing loss + gradients
    int compute_loss = 0;        // Times computing just loss
};

// Parameters to control fitting procedure
class FitParams {
  public:
    int max_it;            // Maximum iterations
    double abstol, reltol; // absolute and relative tolerances for convergence criteria
    double c1 = 1e-3;      // Fraction of first-order improvement required for Armijo condition
    double armijo_tol =
        1e-6;             // Relative error tolerance when computing loss change across iterations
    double ridge_penalty; // Ridge regression penalty
};

// Base-class solver for a GLM.
// The constructor specifies a matrix of GLMs to solve,
// then the methods can solve one particular vector at a time.
class GlmSolver {
  public:
    GlmSolver(int K, FitParams fit_params);

    // Compute the GLM best fit for problem `i`.
    // Best-fit vector returned in `out`, and fitting stats returned
    // in `stats
    virtual void fit_vector(int i, VectorXd &out, FitStats &stats) final;

  protected:
    // Holds intermediate gradient results
    class Derivatives {
      public:
        VectorXd gradient; // Dim K x 1
        MatrixXd hessian;  // Dim K x K
        Derivatives(int K) : gradient(K), hessian(K, K) {}
        void setZero() {
            gradient.setZero();
            hessian.setZero();
        }
        void setZero(int k) {
            gradient.setZero(k);
            hessian.setZero(k, k);
        }
    };

    // Set active problem (e.g. which GLM loss and loss_and_derivatives affect)
    // This should modify the value of `beta` along with other internal variables of the subclass
    virtual void set_problem(int i) = 0;
    // Set `out` to the value of beta
    virtual void get_beta(VectorXd &out) = 0;
    // Return loss for the current problem and `beta` value
    virtual double loss() = 0;
    // Return loss while adding derivatives into `out` (for the current problem and `beta` value)
    virtual double loss_and_derivatives(Derivatives &out) = 0;


    // Current value of beta for calculating loss/derivatives
    // Should not include any offset variables
    VectorXd beta; // Dim K x 1.

  private:
    const FitParams fit_params;

    // Working memory re-used across fit_vector calls
    Derivatives derivatives;
    VectorXd update;    // Dim K x 1
    VectorXd prev_beta; // Dim K x 1
    LLT<MatrixXd> cholesky;

    // Get ridge penalty loss and derivatives.
    // Penalty = ridge_weight / 2 * sum(beta^2)
    // Derivatives are added into `out`
    double regularized_loss_and_derivatives(Derivatives &out);
    double regularized_loss();
};

class PoissonSolverSimd : public GlmSolver {
  public:
    PoissonSolverSimd(
        const MatrixXd &X,                         // Dim (K + J) x N
        std::shared_ptr<const MatrixXd> XY,        // Dim (K + J) x M
        std::shared_ptr<const MatrixXd> beta_init, // Dim (K + J) x M
        std::vector<int> fixed_dims,               // Length J
        FitParams fit_params
    );

    // In a rearranged matrix, We take a K X N matrix and turn it into a 
    // BPCELLS_VEC_FLOAT_SIZE x (K*N/BPCELLS_VEC_FLOAT_SIZE) matrix, padded out
    // with zeros until N is a multiple of BPCELLS_VEC_FLOAT_SIZE.
    // The end result has BPCELLS_VEC_FLOAT_SIZE chunks of each row stored contiguously
    // in memory, and each set of BPCELLS_VEC_FLOAT_SIZE columns across all rows is stored
    // in a contiguous chunk. This allows vectorized loading and maintains memory locality
    // during processing
    using MatRearranged = Eigen::Matrix<float, BPCELLS_VEC_FLOAT_SIZE, Dynamic>;

  private:
    static std::shared_ptr<MatRearranged> rearrange_matrix(const MatrixXd &x);
    static std::vector<int> sorted_vector(std::vector<int> v);

    const int J;
    const int K;
    const std::vector<int> fixed_dims;
    int Xpadding;
    VectorXd offset_beta; // Dim J x 1
    VectorXd XtY;         // Dim K x 1
    VectorXd offset_XtY;  // Dim J x 1

    std::shared_ptr<const MatRearranged> X, offset_X;
    std::shared_ptr<const MatrixXd> XY, beta_init;

    Eigen::Matrix<float, BPCELLS_VEC_FLOAT_SIZE, 1> loss_scratch; // Dim BPCELLS_VEC_FLOAT_SIZE x 1
    MatRearranged gradient_scratch;                               // Dim BPCELLS_VEC_FLOAT_SIZE x K
    MatRearranged hessian_scratch;                                // Dim K * (K + 1) / 2 x 1
    VectorXd hessian_accumulator; // Dim K * (K + 1) / 2 x 1
  
    template <bool CalcDerivatives>
    double loss_and_derivatives_impl(Derivatives *out);
  protected:
    // Set active problem (e.g. which GLM loss and loss_and_derivatives affect)
    // This should modify the value of `beta` along with other internal variables of the subclass
    void set_problem(int i) override;
    // Set `out` to the value of beta
    void get_beta(VectorXd &out) override;
    // Return loss for the current problem and `beta` value
    double loss() override {
      return loss_and_derivatives_impl<false>(NULL);
    }
    // Return loss while adding derivatives into `out` (for the current problem and `beta` value)
    double loss_and_derivatives(Derivatives &out) override {
      return loss_and_derivatives_impl<true>(&out);
    }
};

class PoissonSolverEigen : public GlmSolver {
  public:
    PoissonSolverEigen(
        const MatrixXd &X,                         // Dim (K + J) x N
        std::shared_ptr<const MatrixXd> XY,        // Dim (K + J) x M
        std::shared_ptr<const MatrixXd> beta_init, // Dim (K + J) x M
        std::vector<int> fixed_dims,               // Length J
        FitParams fit_params
    );
  private:
    static std::vector<int> sorted_vector(std::vector<int> v);

    const int J;
    const int K;
    const std::vector<int> fixed_dims;

    VectorXd offset_beta; // Dim J x 1
    VectorXd XtY;         // Dim K x 1
    VectorXd offset_XtY;  // Dim J x 1

    std::shared_ptr<const MatrixXf> X, offset_X;
    std::shared_ptr<const MatrixXd> XY, beta_init;

    VectorXf mu; // Dim N x 1

  protected:
    // Set active problem (e.g. which GLM loss and loss_and_derivatives affect)
    // This should modify the value of `beta` along with other internal variables of the subclass
    void set_problem(int i) override;
    // Set `out` to the value of beta
    void get_beta(VectorXd &out) override;
    // Return loss for the current problem and `beta` value
    double loss() override;
    // Return loss while adding derivatives into `out` (for the current problem and `beta` value)
    double loss_and_derivatives(Derivatives &out) override;
};

} // namespace BPCells