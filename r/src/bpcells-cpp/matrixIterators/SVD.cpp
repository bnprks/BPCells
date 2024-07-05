#include "SVD.h"

#include <algorithm>
#include <Spectra/SymEigsSolver.h>
#include <Eigen/Core>
namespace BPCells {

SpectraMatOp::SpectraMatOp(MatrixLoader<double> *mat, std::atomic<bool> *user_interrupt)
    : mat(mat)
    , user_interrupt(user_interrupt) 
    , tall(mat->rows() > mat->cols()) {}

Eigen::Index SpectraMatOp::rows() const {
    if (tall) return (Eigen::Index)mat->cols();
    else return (Eigen::Index)mat->rows();
}

Eigen::Index SpectraMatOp::cols() const { return rows(); }
// y_out = M * x_in
void SpectraMatOp::perform_op(const double *x_in, double *y_out) const {
    Eigen::Map<Eigen::VectorXd> x_map((double *) x_in, cols());

    Eigen::VectorXd x2;
    if (tall) x2 = mat->vecMultiplyRight(x_map, user_interrupt);
    else x2 = mat->vecMultiplyLeft(x_map, user_interrupt);
    Eigen::Map<Eigen::VectorXd> x2_map(x2.data(), x2.size());

    Eigen::Map<Eigen::VectorXd> out(y_out, rows());
    if (tall) out = mat->vecMultiplyLeft(x2_map, user_interrupt);
    else out = mat->vecMultiplyRight(x2_map, user_interrupt);

    if (user_interrupt != NULL && *user_interrupt) {
        throw SVD::UserInterruptException();
    }
}

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
    std::atomic<bool> *user_interrupt) {

    SVDResult res;

    SpectraMatOp op(mat, user_interrupt);
    Spectra::SymEigsSolver<SpectraMatOp> eigs(op, k, n_cv);
    try {
        eigs.init();
        res.num_converged = eigs.compute(
            Spectra::SortRule::LargestMagn, maxit, tol, Spectra::SortRule::LargestMagn
        );
    } catch (const SVD::UserInterruptException &e) {
        return res;
    }

    if (eigs.info() != Spectra::CompInfo::Successful) return res;

    res.success = res.num_converged == k;
    res.d = eigs.eigenvalues().array().sqrt().matrix();
    res.num_iterations = eigs.num_iterations();
    res.num_operations = eigs.num_operations();

    
    auto d_inv = res.d.array().inverse().matrix().asDiagonal();

    // Calculate u from v or vice versa
    if (op.rows() == mat->rows()) {
        res.u = eigs.eigenvectors();
        // Calculate D^-1 * Ut * M = Vt
        Eigen::MatrixXd tmp = d_inv * res.u.transpose();
        Eigen::Map<Eigen::MatrixXd> tmp_map(tmp.data(), tmp.rows(), tmp.cols());
        res.v = mat->denseMultiplyLeft(tmp_map, user_interrupt).transpose();
    } else {
        res.v = eigs.eigenvectors();
        // Calculate M * V * D^-1 = U
        Eigen::MatrixXd tmp = res.v * d_inv;
        Eigen::Map<Eigen::MatrixXd> tmp_map(tmp.data(), tmp.rows(), tmp.cols());
        res.u = mat->denseMultiplyRight(tmp_map, user_interrupt);
    }
    res.num_operations += k;

    if (transpose) {
        std::swap(res.u, res.v);
    }

    return res;
}

} // namespace BPCells