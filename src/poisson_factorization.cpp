#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>
#include <RcppEigen.h>

#include "lib/sleef_wrapper.h"
#include "matrixIterators/MatrixIterator.h"

using namespace Rcpp;
using namespace Eigen;
using namespace BPCells;

namespace BPCells::poisson {

class FitParams {
  public:
    int max_it;            // Maximum iterations
    double abstol, reltol; // absolute and relative tolerances for convergence criteria
    double c1 = 1e-3;      // Fraction of first-order improvement required for Armijo condition
    double armijo_tol =
        1e-6;             // Relative error tolerance when computing loss change across iterations
    double ridge_penalty; // Ridge regression penalty
};

enum class compute_target {
    loss = 1 << 0,      // Calculate log likelihood loss
    gradients = 1 << 1, // Calculate gradient and hessian
    all = (1 << 0) | (1 << 1)
};

class StepInputs {
  public:
    VectorXd XtY;  // Dim p x 1
    MatrixXf X;    // Dim n x p, but pad rows with 0 to a multiple of BPCELLS_VEC_FLOAT_SIZE
    VectorXd beta; // Dim p x 1
    int Xpadding;  // Number of padding rows added to X (to correct the loss calculation)
    static StepInputs init(int n, int p) {
        StepInputs x;
        x.XtY = VectorXd::Zero(p);
        int padded_rows =
            ((n + BPCELLS_VEC_FLOAT_SIZE - 1) / BPCELLS_VEC_FLOAT_SIZE) * BPCELLS_VEC_FLOAT_SIZE;
        x.X = MatrixXf::Zero(padded_rows, p);
        x.beta = VectorXd::Zero(p);
        x.Xpadding = padded_rows - n;
        return x;
    }
};

class StepOutputs {
  public:
    double loss;
    VectorXd gradient;                                            // Dim p x 1
    MatrixXd hessian;                                             // Dim p x p
    Eigen::Matrix<float, BPCELLS_VEC_FLOAT_SIZE, 1> loss_scratch; // Dim BPCELLS_VEC_FLOAT_SIZE x 1
    Eigen::Matrix<float, BPCELLS_VEC_FLOAT_SIZE, Dynamic>
        gradient_scratch; // Dim BPCELLS_VEC_FLOAT_SIZE x p
    Eigen::Matrix<float, BPCELLS_VEC_FLOAT_SIZE, Dynamic>
        hessian_scratch;          // Dim BPCELLS_VEC_FLOAT_SIZE x p * (p + 1) / 2
    VectorXd hessian_accumulator; // Dim p * (p + 1) / 2 x 1
    static StepOutputs init(int p) {
        StepOutputs x;
        x.gradient = VectorXd::Zero(p);
        x.hessian = MatrixXd::Zero(p, p);
        x.gradient_scratch =
            Eigen::Matrix<float, BPCELLS_VEC_FLOAT_SIZE, Dynamic>::Zero(BPCELLS_VEC_FLOAT_SIZE, p);
        x.hessian_scratch = Eigen::Matrix<float, BPCELLS_VEC_FLOAT_SIZE, Dynamic>::Zero(
            BPCELLS_VEC_FLOAT_SIZE, (p * p + p) / 2
        );
        x.hessian_accumulator = VectorXd::Zero((p * p + p) / 2);
        return x;
    }
};

class FitScratch {
  public:
    VectorXd update;    // Dim p x 1
    VectorXd prev_beta; // Dim p x 1, output parameter for previous beta
    LLT<MatrixXd> cholesky;
    static FitScratch init(int p) {
        FitScratch x;
        x.update = VectorXd::Zero(p);
        x.prev_beta = VectorXd::Zero(p);
        x.cholesky = LLT<MatrixXd>(p);
        return x;
    }
};

// Statistics for fitting a single Poisson GLM vector
class FitStats {
  public:
    int iterations = 0;       // Number of iterations
    int alpha_lower = 0;      // Number of times alpha is lowered
    int alpha_raise = 0;      // Number of times alpha is raised
    double min_alpha = 1;     // Minimum alpha seen during backtracking
    int compute_gradient = 0; // Times computing just gradients
    int compute_loss = 0;     // Times computing just loss
    int compute_both = 0;     // Times computing loss + gradients simultaneously
};

// Calculate loss and/or gradient+hessian for a poisson optimization step
template <compute_target Target>
inline void loss_and_gradients(
    const StepInputs &in, StepOutputs &out, const FitParams params, FitStats &stats
) {
    // Update stats
    if constexpr (Target == compute_target::gradients) stats.compute_gradient += 1;
    else if constexpr (Target == compute_target::loss) stats.compute_loss += 1;
    else if constexpr (Target == compute_target::all) stats.compute_both += 1;

    if constexpr (int(Target) & int(compute_target::gradients)) {
        // gradient = - t(X) * y
        out.gradient = -in.XtY;
        out.hessian.setZero();
        out.gradient_scratch.setZero();
        out.hessian_scratch.setZero();
        out.hessian_accumulator.setZero();
    }
    if constexpr (int(Target) & int(compute_target::loss)) {
        out.loss = -in.beta.dot(in.XtY) + params.ridge_penalty / 2 * in.beta.array().square().sum();
        out.loss_scratch.setZero();
    }
    for (int i = 0; i + BPCELLS_VEC_FLOAT_SIZE <= in.X.rows(); i += BPCELLS_VEC_FLOAT_SIZE) {
        // Calculate lambda = exp(X * beta)
        vec_float lambda = splat_float(0.0);
        for (int p = 0; p < in.beta.rows(); p++) {
            vec_float Xp = load_float(in.X.col(p).data() + i);
            lambda = fma_f(Xp, splat_float(in.beta(p)), lambda);
        }
        lambda = exp_f(lambda);

        // Loss calculation
        if constexpr (int(Target) & int(compute_target::loss)) {
            store_float(
                out.loss_scratch.data(), add_f(lambda, load_float(out.loss_scratch.data()))
            );
        }

        // Gradient + Hessian calculation
        if constexpr (int(Target) & int(compute_target::gradients)) {
            // gradient += t(lambda) * X
            // hessian += Xt * diag(lambda) * X
            float *hessian_tmp = out.hessian_scratch.data();
            for (int p = 0; p < in.beta.rows(); p++) {
                float *gradient_tmp = out.gradient_scratch.col(p).data();
                vec_float XpLambda = mul_f(lambda, load_float(in.X.col(p).data() + i));
                store_float(gradient_tmp, add_f(XpLambda, load_float(gradient_tmp)));
                for (int q = 0; q <= p; q++) {
                    vec_float Xq = load_float(in.X.col(q).data() + i);
                    store_float(hessian_tmp, fma_f(XpLambda, Xq, load_float(hessian_tmp)));
                    hessian_tmp += BPCELLS_VEC_FLOAT_SIZE;
                }
            }
        }

        // Store to the main accumulators to mitigate some of the precision loss during summation
        if (i % (64 * BPCELLS_VEC_FLOAT_SIZE) == 64 * BPCELLS_VEC_FLOAT_SIZE - 1) {
            // if (true) {
            if constexpr (int(Target) & int(compute_target::loss)) {
                out.loss += out.loss_scratch.sum();
                out.loss_scratch.setZero();
            }
            if constexpr (int(Target) & int(compute_target::gradients)) {
                out.gradient += out.gradient_scratch.cast<double>().colwise().sum();
                out.hessian_accumulator.array() +=
                    out.hessian_scratch.cast<double>().colwise().sum().array();
                out.gradient_scratch.setZero();
                out.hessian_scratch.setZero();
            }
        }
    }
    // Final cleanup
    if constexpr (int(Target) & int(compute_target::loss)) {
        out.loss += out.loss_scratch.sum();
        out.loss -= in.Xpadding;
    }
    if constexpr (int(Target) & int(compute_target::gradients)) {
        out.gradient += out.gradient_scratch.cast<double>().colwise().sum();
        out.hessian_accumulator.array() +=
            out.hessian_scratch.cast<double>().colwise().sum().array();
        // Copy hessian accumulator to the actual lower triangle of the matrix
        int i = 0;
        for (int p = 0; p < in.beta.rows(); p++) {
            for (int q = 0; q <= p; q++) {
                out.hessian(p, q) = out.hessian_accumulator(i);
                i += 1;
            }
        }
        // Add ridge penalty
        out.gradient += in.beta * params.ridge_penalty;
        out.hessian.diagonal().array() += params.ridge_penalty;
    }
}

// Fit a GLM vector
// Preconditions:
//   - step_in should be fully initialized
// Postconditions:
//   - step_out contains a fully fit beta vector
void glm_fit_vector(
    StepInputs &step_in,
    StepOutputs &step_out,
    FitScratch &fit_scratch,
    const FitParams &fit_params,
    FitStats &stats
) {
    double alpha = 1;
    double prev_loss = INFINITY;
    double minimum_loss_change;

    poisson::loss_and_gradients<compute_target::all>(step_in, step_out, fit_params, stats);

    int it;
    for (it = 0; it < fit_params.max_it; it++) {
        // Convergence criteria: abs(beta_new - beta_old) <= (abstol + reltol * absolute(beta_old))
        if (it > 0) {
            auto abs_diff = (fit_scratch.prev_beta - step_in.beta).array().abs();
            auto max_abs_diff =
                fit_params.abstol + fit_params.reltol * fit_scratch.prev_beta.array().abs();
            if ((abs_diff < max_abs_diff).all()) {
                it += 1;
                break;
            }
        }

        fit_scratch.cholesky.compute(step_out.hessian);
        fit_scratch.update = fit_scratch.cholesky.solve(-step_out.gradient);
        minimum_loss_change = fit_scratch.update.dot(step_out.gradient) * fit_params.c1;
        prev_loss = step_out.loss;
        fit_scratch.prev_beta = step_in.beta;
        step_in.beta += alpha * fit_scratch.update;

        if (it + 1 < fit_params.max_it && alpha == 1) {
            // Optimistically compute the gradients assuming we will not need backtracking
            poisson::loss_and_gradients<compute_target::all>(step_in, step_out, fit_params, stats);
        } else {
            // We just need loss to begin backtracking
            poisson::loss_and_gradients<compute_target::loss>(step_in, step_out, fit_params, stats);
        }

        bool armijo = step_out.loss <= prev_loss + alpha * minimum_loss_change ||
                      abs(step_out.loss - prev_loss) / prev_loss < fit_params.armijo_tol;
        if (it > 0 && (alpha != 1 || !armijo)) {
            // Perform line search with Armijo condition
            // See: https://en.wikipedia.org/wiki/Backtracking_line_search and Armijo conditions
            // In most cases we will should skip this (alpha == 1 and Armijo condition is
            // satisfied), but for unruly fits we may hit these conditions

            if (!armijo) {
                // Lower our alpha until Armijo condition is satisfied
                while (alpha > 1e-9 && !armijo) {
                    alpha *= 0.5;
                    stats.alpha_lower += 1;
                    step_in.beta = fit_scratch.prev_beta + alpha * fit_scratch.update;
                    poisson::loss_and_gradients<compute_target::loss>(
                        step_in, step_out, fit_params, stats
                    );
                    armijo = step_out.loss <= prev_loss + alpha * minimum_loss_change ||
                             abs(step_out.loss - prev_loss) / prev_loss < fit_params.armijo_tol;
                }
                if (it + 1 < fit_params.max_it) {
                    poisson::loss_and_gradients<compute_target::gradients>(
                        step_in, step_out, fit_params, stats
                    );
                }
            } else {
                // Raise our alpha until alpha == 1 or Armijo condition is not satisified
                while (alpha < 1.0 && armijo) {
                    alpha = std::min(1.0, alpha * 2);
                    stats.alpha_raise += 1;
                    step_in.beta = fit_scratch.prev_beta + alpha * fit_scratch.update;
                    poisson::loss_and_gradients<compute_target::loss>(
                        step_in, step_out, fit_params, stats
                    );
                    armijo = step_out.loss <= prev_loss + alpha * minimum_loss_change ||
                             abs(step_out.loss - prev_loss) / prev_loss < fit_params.armijo_tol;
                }
                if (!armijo) {
                    alpha *= 0.5;
                    stats.alpha_raise -= 1; // Roll back our tested alpha that went too far
                    step_in.beta = fit_scratch.prev_beta + alpha * fit_scratch.update;
                    if (it + 1 < fit_params.max_it) {
                        poisson::loss_and_gradients<compute_target::all>(
                            step_in, step_out, fit_params, stats
                        );
                    }
                }
            }
            stats.min_alpha = std::min(alpha, stats.min_alpha);
        }
    }
    stats.iterations = it;
}

} // namespace BPCells::poisson

// Return the loss, gradient, and hessian from the given beta
// [[Rcpp::export]]
SEXP glm_check_gradient_cpp(
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::Map<Eigen::MatrixXd> XtY,
    const Eigen::Map<Eigen::MatrixXd> beta_init,
    double ridge_penalty
) {
    using namespace BPCells::poisson;
    int p = X.cols();
    int n = X.rows();

    FitParams fit_params;
    fit_params.ridge_penalty = ridge_penalty;

    StepInputs step_in = StepInputs::init(n, p);
    step_in.X.topRows(n) = X.cast<float>();
    step_in.XtY = XtY.col(0);
    step_in.beta = beta_init.col(0);

    StepOutputs step_out = StepOutputs::init(p);

    FitStats stats;

    poisson::loss_and_gradients<compute_target::all>(step_in, step_out, fit_params, stats);

    return List::create(
        Named("loss") = step_out.loss,
        Named("gradient") = step_out.gradient,
        Named("hessian") = step_out.hessian
    );
}

// Fit a series of poisson GLM problems with a shared model matrix
// Parameters:
//   - X: Shared model matrix (dimension n x p)
//   - YtX: XPtr<IterableMatrix> (dimension p x m), result of t(Y) * X
//   - beta_init: Initial guess for beta (dimension p x m)
//   - ridge_penalty: multiplier for ridge penalty in loss function
//   - max_it: Maximum number of Newton-Raphson iterations
//   - abstol, reltol: Convergence criteria: test abs(beta_new - beta_old) <= (abstol + reltol *
//   absolute(beta_old))
//                     (this is like numpy.allclose)
// Return list of:
//   - "beta": Matrix of fit betas (dimension m x p)
//   - "num_iters": Number of Netwon-Raphson iterations performed
// [[Rcpp::export]]
SEXP glm_fit_matrix_cpp(
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::Map<Eigen::MatrixXd> XtY,
    const Eigen::Map<Eigen::MatrixXd> beta_init,
    double ridge_penalty,
    int max_it,
    double abstol,
    double reltol
) {
    using namespace BPCells::poisson;

    int p = X.cols();
    int n = X.rows();
    int m = beta_init.cols();

    FitParams fit_params;
    fit_params.max_it = max_it;
    fit_params.abstol = abstol;
    fit_params.reltol = reltol;
    fit_params.ridge_penalty = ridge_penalty;

    StepInputs step_in = StepInputs::init(n, p);
    step_in.X.topRows(n) = X.cast<float>();

    StepOutputs step_out = StepOutputs::init(p);
    FitScratch fit_scratch = FitScratch::init(p);

    MatrixXd beta_transpose(p, m);

    // We'll return a vector for each fitting stat
    std::vector<int> iterations;
    std::vector<int> alpha_lower;
    std::vector<int> alpha_raise;
    std::vector<double> min_alpha;
    std::vector<int> compute_gradient;
    std::vector<int> compute_loss;
    std::vector<int> compute_both;

    for (int i = 0; i < m; i++) {
        FitStats stats;
        step_in.XtY = XtY.col(i);
        step_in.beta = beta_init.col(i);
        glm_fit_vector(step_in, step_out, fit_scratch, fit_params, stats);
        beta_transpose.col(i) = step_in.beta;
        // Track all the stats
        iterations.push_back(stats.iterations);
        alpha_lower.push_back(stats.alpha_lower);
        alpha_raise.push_back(stats.alpha_raise);
        min_alpha.push_back(stats.min_alpha);
        compute_gradient.push_back(stats.compute_gradient);
        compute_loss.push_back(stats.compute_loss);
        compute_both.push_back(stats.compute_both);
    }

    // Return fit with stats on fit computations
    return List::create(
        Named("beta") = beta_transpose.transpose(),
        Named("fit_stats") = DataFrame::create(
            Named("iterations") = iterations,
            Named("alpha_lower") = alpha_lower,
            Named("alpha_raise") = alpha_raise,
            Named("min_alpha") = min_alpha,
            Named("compute_gradient") = compute_gradient,
            Named("compute_loss") = compute_loss,
            Named("compute_both") = compute_both
        )
    );
}

// Fit a series of poisson GLM problems with a shared model matrix
// Parameters:
//   - X: Shared model matrix (dimension n x p)
//   - Y: XPtr<IterableMatrix> (dimension n x m)
//   - beta_init: Initial guess for beta (dimension m x p)
//   - ridge_penalty: multiplier for ridge penalty in loss function
//   - max_it: Maximum number of Newton-Raphson iterations
//   - abstol, reltol: Convergence criteria: test abs(beta_new - beta_old) <= (abstol + reltol *
//   absolute(beta_old))
//                     (this is like numpy.allclose)
// Return list of:
//   - "beta": Matrix of fit betas (dimension m x p)
//   - "num_iters": Number of Netwon-Raphson iterations performed
// [[Rcpp::export]]
SEXP poisson_glm_matrix_cpp(
    const Eigen::Map<Eigen::MatrixXd> X,
    SEXP Y,
    const Eigen::Map<Eigen::MatrixXd> beta_init,
    double ridge_penalty,
    int max_it,
    double abstol,
    double reltol
) {
    int p = X.cols();
    int n = X.rows();
    int m = beta_init.rows();

    // Note: Remember that Eigen defaults to col-major storage, hence our choice of
    // orientations for gradient_from_data and X1 when we expect column-based slicing to be needed
    XPtr<BPCells::MatrixLoader<double>> Y_loader(Y);
    MatrixXd gradient_from_data; // dimension p x m
    {
        MatrixXd Xt_tmp(X.transpose());
        gradient_from_data = -1 * Y_loader->denseMultiplyLeft(
                                      Map<MatrixXd>(Xt_tmp.data(), Xt_tmp.rows(), Xt_tmp.cols()),
                                      &Rcpp::checkUserInterrupt
                                  );
    }

    // Appending a row of all 1 to X helps reduce our Newton-Raphson math to a single matrix
    // multiply
    MatrixXd X1(p + 1, n); // p+1 x n
    X1.topRows(p) = X.transpose();
    X1.bottomRows(1).setOnes();

    LLT<MatrixXd> cholesky; // Share memory for cholesky decomposition among all regressions
    MatrixXd hessian_gradient(p, p + 1); // Hessian for first p columns, followed by gradient vector
    VectorXd update(p);                  // Update vector

    MatrixXd beta = beta_init.transpose(); // p x n
    int num_iters = 0;
    int num_backtracks = 0;
    bool backtrack_error = false;
    bool iters_error = false;

    VectorXd y(n);
    Y_loader->restart();
    for (int i = 0; i < m; i++) {
        // Calculate log_likelihood of starting solution
        double last_neg_log_likelihood =
            (X * beta.col(i)).array().exp().sum() +
            (gradient_from_data.transpose() * beta.col(i)).array().sum();
        int it;
        for (it = 0; it < max_it; it++) {
            num_iters += 1;
            // Calculate hessian and gradient
            // lambda = exp(X * beta)
            // hessian = Xt * lambda * X
            // gradient = Xt * lambda - Xt * y
            auto lambda = (X * beta.col(i)).array().exp().matrix().asDiagonal();
            hessian_gradient.noalias() = X1.topRows(p) * lambda * X1.transpose();
            hessian_gradient.rightCols<1>() += gradient_from_data.col(i);
            // Add ridge penalty
            hessian_gradient.rightCols<1>() += ridge_penalty * beta.col(i);   // gradient
            hessian_gradient.leftCols(p).diagonal().array() += ridge_penalty; // hessian

            // Compute Newton-Raphson update direction
            // beta_n+1 = beta_n - Hessian^-1 * gradient
            // Hessian * (beta_n+1 - beta_n) = -gradient
            cholesky.compute(hessian_gradient.leftCols(p));
            update = cholesky.solve(hessian_gradient.rightCols<1>());

            beta.col(i) -= update;
            // Check for convergence
            if ((update.array().abs() < abstol + reltol * beta.col(i).array().abs()).all()) break;

            // Armijo backtracking line search
            int backtracks;
            int max_bactracks = 30;
            double c1 =
                1e-3; // As suggestend in
                      // https://optimization.cbe.cornell.edu/index.php?title=Line_search_methods#Armijo_.28Sufficient_Decrease.29_Condition
            double minimum_improvement = c1 * update.dot(hessian_gradient.rightCols<1>());
            for (backtracks = 0; backtracks < max_bactracks; backtracks++) {
                double neg_log_likelihood = (X * beta.col(i)).array().exp().sum() -
                                            (y.transpose() * X * beta.col(i)).array().sum();
                if (neg_log_likelihood <= last_neg_log_likelihood + minimum_improvement) {
                    last_neg_log_likelihood = neg_log_likelihood;
                    break;
                } else {
                    update /= 2;
                    minimum_improvement /= 2;
                    beta.col(i) += update; // This might have some numerical stability concerns
                }
            }
            num_backtracks += backtracks;
            if (backtracks == max_bactracks && !backtrack_error) {
                backtrack_error = true;
                printf("Warning: exceeded maximum backtracks on problem %d\n", i);
            }
        }
        num_iters += it;
        if (it == max_it && !iters_error) {
            iters_error = true;
            printf("Warning: exceeded maximum iterations on problem %d\n", i);
        }
    }

    return List::create(
        Named("beta") = beta.transpose(),
        Named("num_iters") = num_iters,
        Named("num_backtracks") = num_backtracks
    );
}
