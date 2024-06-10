// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#if 0

#include "GLMFactorization.h"

using namespace Eigen;

namespace BPCells::poisson2 {

GlmSolver::GlmSolver(int K, FitParams fit_params)
    : beta(K)
    , fit_params(fit_params)
    , derivatives(K)
    , update(K)
    , prev_beta(K)
    , cholesky(K) {}

void GlmSolver::fit_vector(int i, VectorXd &out, FitStats &stats) {
    set_problem(i);

    double alpha = 1;
    double prev_loss = INFINITY;
    double minimum_loss_change;

    derivatives.setZero();
    double curr_loss = regularized_loss_and_derivatives(derivatives);
    stats.compute_derivatives += 1;

    int it;
    for (it = 0; it < fit_params.max_it; it++) {
        // Convergence criteria: abs(beta_new - beta_old) <= (abstol + reltol * absolute(beta_old))
        if (it > 0) {
            auto abs_diff = (prev_beta - beta).array().abs();
            auto max_abs_diff = fit_params.abstol + fit_params.reltol * prev_beta.array().abs();
            if ((abs_diff < max_abs_diff).all()) {
                it += 1;
                break;
            }
        }

        cholesky.compute(derivatives.hessian);
        update = cholesky.solve(-derivatives.gradient);
        minimum_loss_change = update.dot(derivatives.gradient) * fit_params.c1;
        prev_loss = curr_loss;
        prev_beta = beta;
        beta += alpha * update;

        if (it + 1 < fit_params.max_it && alpha == 1) {
            // Optimistically compute the gradients assuming we will not need backtracking
            derivatives.setZero();
            curr_loss = regularized_loss_and_derivatives(derivatives);
            stats.compute_derivatives += 1;
        } else {
            // We just need loss to begin backtracking
            curr_loss = regularized_loss();
            stats.compute_loss += 1;
        }

        bool armijo = curr_loss <= prev_loss + alpha * minimum_loss_change ||
                      abs(curr_loss - prev_loss) / prev_loss < fit_params.armijo_tol;
        if (alpha != 1 || !armijo) {
            // Perform line search with Armijo condition
            // See: https://en.wikipedia.org/wiki/Backtracking_line_search and Armijo conditions
            // In most cases we will should skip this (alpha == 1 and Armijo condition is
            // satisfied), but for unruly fits we may hit these conditions

            if (!armijo) {
                // Lower our alpha until Armijo condition is satisfied
                while (alpha > 1e-9 && !armijo) {
                    alpha *= 0.5;
                    stats.alpha_lower += 1;
                    beta = prev_beta + alpha * update;
                    curr_loss = regularized_loss();
                    stats.compute_loss += 1;
                    armijo = curr_loss <= prev_loss + alpha * minimum_loss_change ||
                             abs(curr_loss - prev_loss) / prev_loss < fit_params.armijo_tol;
                }
                if (it + 1 < fit_params.max_it) {
                    derivatives.setZero();
                    curr_loss = regularized_loss_and_derivatives(derivatives);
                    stats.compute_derivatives += 1;
                }
            } else {
                // Raise our alpha until alpha == 1 or Armijo condition is not satisified
                while (alpha < 1.0 && armijo) {
                    alpha = std::min(1.0, alpha * 2);
                    stats.alpha_raise += 1;
                    beta = prev_beta + alpha * update;
                    curr_loss = regularized_loss();
                    stats.compute_loss += 1;
                    armijo = curr_loss <= prev_loss + alpha * minimum_loss_change ||
                             abs(curr_loss - prev_loss) / prev_loss < fit_params.armijo_tol;
                }
                if (!armijo) {
                    alpha *= 0.5;
                    stats.alpha_raise -= 1; // Roll back our tested alpha that went too far
                    beta = prev_beta + alpha * update;
                    if (it + 1 < fit_params.max_it) {
                        derivatives.setZero();
                        curr_loss = regularized_loss_and_derivatives(derivatives);
                        stats.compute_derivatives += 1;
                    }
                }
            }
            stats.min_alpha = std::min(alpha, stats.min_alpha);
        }
    }
    stats.iterations = it;
    get_beta(out);
}

double GlmSolver::loss_and_derivatives(int i, VectorXd &gradient_out, MatrixXd &hessian_out) {
    set_problem(i);
    derivatives.setZero();
    double l = regularized_loss_and_derivatives(derivatives);
    gradient_out = derivatives.gradient;
    hessian_out = derivatives.hessian;
    return l;
}

double GlmSolver::regularized_loss_and_derivatives(Derivatives &out) {
    double loss = loss_and_derivatives(out);
    loss += fit_params.ridge_penalty / 2 * beta.array().square().sum();

    out.gradient += beta * fit_params.ridge_penalty;
    out.hessian.diagonal().array() += fit_params.ridge_penalty;

    return loss;
}

double GlmSolver::regularized_loss() {
    return loss() + fit_params.ridge_penalty / 2 * beta.array().square().sum();
}

PoissonSolverSimd::PoissonSolverSimd(
    const MatrixXd &X,                         // Dim (K + J) x N
    std::shared_ptr<const MatrixXd> XY,        // Dim (K + J) x M
    std::shared_ptr<const MatrixXd> beta_init, // Dim (K + J) x M
    std::vector<int> fixed_dims,               // Length J
    FitParams fit_params
)
    : GlmSolver(X.rows() - fixed_dims.size(), fit_params)
    , J(fixed_dims.size())
    , K(X.rows() - fixed_dims.size())
    , fixed_dims(sorted_vector(fixed_dims))
    , offset_beta(J)
    , XtY(K)
    , offset_XtY(J)
    , XY(XY)
    , beta_init(beta_init)
    , loss_scratch(BPCELLS_VEC_FLOAT_SIZE)
    , gradient_scratch(BPCELLS_VEC_FLOAT_SIZE, K)
    , hessian_scratch(BPCELLS_VEC_FLOAT_SIZE, (K + K * K) / 2)
    , hessian_accumulator((K * K + K) / 2) {

    // Check there aren't any duplicates in fixed_dims
    for (size_t i = 1; i < fixed_dims.size(); i++) {
        if (fixed_dims[i - 1] == fixed_dims[i]) {
            throw std::runtime_error("Error in PoissonSolverSimd: repeated entries in fixed_dims");
        }
    }

    // Split and re-order the X matrix
    MatrixXd X_var(K, X.cols());
    MatrixXd X_off(J, X.cols());
    size_t k = 0;
    size_t j = 0;
    for (int i = 0; i < K + J; i++) {
        if (j < this->fixed_dims.size() && this->fixed_dims[j] == i) {
            X_off.row(j) = X.row(i);
            j++;
        } else {
            X_var.row(k) = X.row(i);
            k++;
        }
    }
    this->X = rearrange_matrix(X_var);
    this->offset_X = rearrange_matrix(X_off);
    Xpadding = this->X->cols() / K * BPCELLS_VEC_FLOAT_SIZE - X.cols();
}

std::shared_ptr<PoissonSolverSimd::MatRearranged>
PoissonSolverSimd::rearrange_matrix(const MatrixXd &X) {
    int n = X.cols();
    int k = X.rows();
    int padded_rows =
        ((n + BPCELLS_VEC_FLOAT_SIZE - 1) / BPCELLS_VEC_FLOAT_SIZE) * BPCELLS_VEC_FLOAT_SIZE;

    auto X_ptr = std::make_shared<PoissonSolverSimd::MatRearranged>(
        BPCELLS_VEC_FLOAT_SIZE, k * padded_rows / BPCELLS_VEC_FLOAT_SIZE
    );
    X_ptr->setZero();
    for (int i = 0; i < padded_rows / BPCELLS_VEC_FLOAT_SIZE; i++) {
        int copy_size = std::min(BPCELLS_VEC_FLOAT_SIZE, n - i * BPCELLS_VEC_FLOAT_SIZE);
        X_ptr->block(0, i * k, copy_size, k) =
            X.block(0, i * BPCELLS_VEC_FLOAT_SIZE, k, copy_size).transpose().cast<float>();
    }
    return X_ptr;
}

std::vector<int> PoissonSolverSimd::sorted_vector(std::vector<int> v) {
    std::sort(v.begin(), v.end());
    return v;
}

// Set active problem (e.g. which GLM loss and loss_and_derivatives affect)
// This should modify the value of `beta` along with other internal variables of the subclass
void PoissonSolverSimd::set_problem(int i) {
    size_t k = 0;
    size_t j = 0;
    for (int d = 0; d < K + J; d++) {
        if (j < fixed_dims.size() && fixed_dims[j] == d) {
            offset_beta(j) = (*beta_init)(d, i);
            offset_XtY(j) = (*XY)(d, i);
            j++;
        } else {
            beta(k) = (*beta_init)(d, i);
            XtY(k) = (*XY)(d, i);
            k++;
        }
    }
}
// Set `out` to the value of beta
void PoissonSolverSimd::get_beta(VectorXd &out) {
    out.resize(K + J);
    size_t k = 0;
    size_t j = 0;
    for (int d = 0; d < K + J; d++) {
        if (j < fixed_dims.size() && fixed_dims[j] == d) {
            out(d) = offset_beta(j);
            j++;
        } else {
            out(d) = beta(k);
            k++;
        }
    }
}

// Return loss while adding derivatives into `out` (for the current problem and `beta` value)
template <bool CalcDerivatives>
double PoissonSolverSimd::loss_and_derivatives_impl(Derivatives *out) {
    double loss = -beta.dot(XtY) + -offset_beta.dot(offset_XtY);

    loss_scratch.setZero();
    if constexpr (CalcDerivatives) {
        gradient_scratch.setZero();
        hessian_scratch.setZero();
        hessian_accumulator.setZero();

        // gradient = - t(X) * y
        out->gradient -= XtY;
    }

    int N = X->cols() / K;
    for (int i = 0; i < N; i++) {
        const float *x = X->data() + i * K * BPCELLS_VEC_FLOAT_SIZE;
        const float *x_off = offset_X->data() + i * J * BPCELLS_VEC_FLOAT_SIZE;
        // Calculate mu = exp(X * beta + O_x * O_beta)
        vec_float mu = splat_float(0.0);
        for (int k = 0; k < K; k++) {
            vec_float Xk = load_float(x + k * BPCELLS_VEC_FLOAT_SIZE);
            mu = fma_f(Xk, splat_float(beta(k)), mu);
        }
        for (int j = 0; j < J; j++) {
            vec_float Oj = load_float(x_off + j * BPCELLS_VEC_FLOAT_SIZE);
            mu = fma_f(Oj, splat_float(offset_beta(j)), mu);
        }
        mu = exp_f(mu);

        // Loss calculation
        store_float(loss_scratch.data(), add_f(mu, load_float(loss_scratch.data())));

        if constexpr (CalcDerivatives) {
            // Gradient + Hessian calculation

            // gradient += t(mu) * X
            // hessian += Xt * diag(mu) * X
            float *hessian_tmp = hessian_scratch.data();
            for (int k = 0; k < K; k++) {
                float *gradient_tmp = gradient_scratch.col(k).data();
                vec_float XkMu = mul_f(mu, load_float(x + k * BPCELLS_VEC_FLOAT_SIZE));
                store_float(gradient_tmp, add_f(XkMu, load_float(gradient_tmp)));
                for (int j = 0; j <= k; j++) {
                    vec_float Xj = load_float(x + j * BPCELLS_VEC_FLOAT_SIZE);
                    store_float(hessian_tmp, fma_f(XkMu, Xj, load_float(hessian_tmp)));
                    hessian_tmp += BPCELLS_VEC_FLOAT_SIZE;
                }
            }
        }

        // Store to the main accumulators to mitigate some of the precision loss during summation
        if (i % 64 == 64 - 1) {
            loss += loss_scratch.sum();
            loss_scratch.setZero();

            if constexpr (CalcDerivatives) {
                out->gradient += gradient_scratch.cast<double>().colwise().sum();
                hessian_accumulator.array() +=
                    hessian_scratch.cast<double>().colwise().sum().array();
                gradient_scratch.setZero();
                hessian_scratch.setZero();
            }
        }
    }
    // Final cleanup

    loss += loss_scratch.sum();
    loss -= Xpadding;

    if (CalcDerivatives) {
        out->gradient += gradient_scratch.cast<double>().colwise().sum();
        hessian_accumulator.array() += hessian_scratch.cast<double>().colwise().sum().array();
        // Copy hessian accumulator to the actual lower triangle of the matrix
        int i = 0;
        for (int k = 0; k < K; k++) {
            for (int j = 0; j <= k; j++) {
                out->hessian(k, j) = hessian_accumulator(i);
                i += 1;
            }
        }
    }
    return loss;
}

} // namespace BPCells::poisson2

#endif