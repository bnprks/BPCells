// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "SCTransform.h"
#include <atomic>

#include "../simd/sctransform.h"

namespace BPCells {

SCTransformPearsonSIMD::SCTransformPearsonSIMD(
    std::unique_ptr<MatrixLoader<double>> &&loader, TransformFit fit
)
    : MatrixTransformDense(std::move(loader), fit)
    , theta_inv(fit.row_params.row(0).cast<float>())
    , cell_read_counts(fit.col_params.row(0).cast<float>())
    , gene_beta(fit.row_params.row(1).cast<float>())
    , sd_inv_max(fit.global_params(0))
    , clip_min(fit.global_params(1))
    , clip_max(fit.global_params(2)) {}

bool SCTransformPearsonSIMD::loadZeroSubtracted(MatrixLoader<double> &loader) {
    if (!loader.load()) return false;

    simd::sctransform_zero_subtracted(
        loader.valData(),
        cell_read_counts(currentCol()),
        loader.rowData(),
        gene_beta.data(),
        theta_inv.data(),
        {this->sd_inv_max, this->clip_min, this->clip_max},
        loader.capacity()
    );

    return true;
}

void SCTransformPearsonSIMD::loadZero(
    double *values, uint32_t count, uint32_t start_row, uint32_t col
) {
    simd::sctransform_load_zero(
        values,
        cell_read_counts(col),
        gene_beta.data() + start_row,
        theta_inv.data() + start_row,
        {this->sd_inv_max, this->clip_min, this->clip_max},
        count
    );
}

inline void sctransform_vec_multiply_right_zero_impl(
    Eigen::VectorXd &out,
    const Eigen::Map<Eigen::VectorXd> v,
    std::atomic<bool> *user_interrupt,
    const Eigen::ArrayXf &cell_read_counts,
    const Eigen::ArrayXf &theta_inv,
    const Eigen::ArrayXf &gene_beta,
    const simd::SCTransformClipParam &clip_params
) {
    Eigen::VectorXf out_float(out.rows());
    out_float.setZero();

    // Logically work in the orientation where rows are genes and cells are columns
    uint32_t nrows = gene_beta.size();
    uint32_t ncols = cell_read_counts.size();

    for (uint32_t col = 0; col < ncols; col++) {
        if (user_interrupt != NULL && *user_interrupt) return;
        // Periodically flush our single-precision accumulator to avoid
        // excessive loss of precision during summation
        if (col % 64 == 0) {
            out += out_float.cast<double>();
            out_float.setZero();
        }

        simd::sctransform_multiply_right_zero(
            out_float.data(),
            v(col),
            cell_read_counts(col),
            gene_beta.data(),
            theta_inv.data(),
            clip_params,
            nrows
        );
    }
    // Output final results
    out += out_float.cast<double>();
}

inline void sctransform_vec_multiply_left_zero_impl(
    Eigen::VectorXd &out,
    const Eigen::Map<Eigen::VectorXd> v,
    std::atomic<bool> *user_interrupt,
    const Eigen::ArrayXf &cell_read_counts,
    const Eigen::ArrayXf &theta_inv,
    const Eigen::ArrayXf &gene_beta,
    const simd::SCTransformClipParam &clip_params
) {
    Eigen::VectorXf out_float(out.rows());
    out_float.setZero();

    // Logically work in the orientation where rows are genes and cells are columns
    uint32_t nrows = gene_beta.size();
    uint32_t ncols = cell_read_counts.size();

    for (uint32_t row = 0; row < nrows; row++) {
        if (user_interrupt != NULL && *user_interrupt) return;
        // Periodically flush our single-precision accumulator to avoid
        // excessive loss of precision during summation
        if (row % 64 == 0) {
            out += out_float.cast<double>();
            out_float.setZero();
        }

        simd::sctransform_multiply_left_zero(
            out_float.data(),
            v(row),
            cell_read_counts.data(),
            gene_beta(row),
            theta_inv(row),
            clip_params,
            ncols
        );
    }
    // Output final results
    out += out_float.cast<double>();
}

// Calculate matrix-vector product A*v where A (this) is sparse and B is a dense matrix.
void SCTransformPearsonSIMD::vecMultiplyRightZero(
    Eigen::VectorXd &out, const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    sctransform_vec_multiply_right_zero_impl(
        out,
        v,
        user_interrupt,
        cell_read_counts,
        theta_inv,
        gene_beta,
        {this->sd_inv_max, this->clip_min, this->clip_max}
    );
}

void SCTransformPearsonSIMD::vecMultiplyLeftZero(
    Eigen::VectorXd &out, const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    sctransform_vec_multiply_left_zero_impl(
        out,
        v,
        user_interrupt,
        cell_read_counts,
        theta_inv,
        gene_beta,
        {this->sd_inv_max, this->clip_min, this->clip_max}
    );
}

// ###################################
//  Transposed SIMD implementation
// ###################################
// Basically copy-paste from non-SIMD version with some small edits to account
// for the transposition


SCTransformPearsonTransposeSIMD::SCTransformPearsonTransposeSIMD(
    std::unique_ptr<MatrixLoader<double>> &&loader, TransformFit fit
)

    : MatrixTransformDense(std::move(loader), fit)
    , theta_inv(fit.col_params.row(0).cast<float>())
    , cell_read_counts(fit.row_params.row(0).cast<float>())
    , gene_beta(fit.col_params.row(1).cast<float>())
    , sd_inv_max(fit.global_params(0))
    , clip_min(fit.global_params(1))
    , clip_max(fit.global_params(2)) {}

bool SCTransformPearsonTransposeSIMD::loadZeroSubtracted(MatrixLoader<double> &loader) {
    if (!loader.load()) return false;

    simd::sctransform_zero_subtracted_transpose(
        loader.valData(),
        cell_read_counts.data(),
        loader.rowData(),
        gene_beta(currentCol()),
        theta_inv(currentCol()),
        {this->sd_inv_max, this->clip_min, this->clip_max},
        loader.capacity()
    );
    return true;
}

void SCTransformPearsonTransposeSIMD::loadZero(
    double *values, uint32_t count, uint32_t start_row, uint32_t col
) {
    simd::sctransform_load_zero_transpose(
        values,
        cell_read_counts.data() + start_row,
        gene_beta(col),
        theta_inv(col),
        {this->sd_inv_max, this->clip_min, this->clip_max},
        count
    );
}

// Calculate matrix-vector product A*v where A (this) is sparse and B is a dense matrix.
void SCTransformPearsonTransposeSIMD::vecMultiplyRightZero(
    Eigen::VectorXd &out, const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    sctransform_vec_multiply_left_zero_impl(
        out,
        v,
        user_interrupt,
        cell_read_counts,
        theta_inv,
        gene_beta,
        {this->sd_inv_max, this->clip_min, this->clip_max}
    );
}

void SCTransformPearsonTransposeSIMD::vecMultiplyLeftZero(
    Eigen::VectorXd &out, const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    sctransform_vec_multiply_right_zero_impl(
        out,
        v,
        user_interrupt,
        cell_read_counts,
        theta_inv,
        gene_beta,
        {this->sd_inv_max, this->clip_min, this->clip_max}
    );
}


// ###################################
//  Non-SIMD variants
// ###################################

bool SCTransformPearson::loadZeroSubtracted(MatrixLoader<double> &loader) {
    if (!loader.load()) return false;

    uint32_t *row_data = loader.rowData();
    double *val_data = loader.valData();
    uint32_t capacity = loader.capacity();

    double cell_reads = fit.col_params(0, currentCol());
    double sd_inv_max = fit.global_params(0);
    double clip_max = fit.global_params(2);
    double clip_min = fit.global_params(1);

    for (uint32_t i = 0; i < capacity; i++) {
        double mu = cell_reads * fit.row_params(1, row_data[i]);
        double theta_inv = fit.row_params(0, row_data[i]);
        double sd_i = std::min(sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));

        double zero_val = -mu * sd_i;
        // val = (X-mu)/sd
        double val = val_data[i] * sd_i + zero_val;
        // val = clamp(val, min, max)
        val = std::max(std::min(val, clip_max), clip_min);
        // val = val - max(-mu/sd, clip_min)
        val_data[i] = val - std::max(zero_val, clip_min);
    }
    return true;
}

void SCTransformPearson::loadZero(
    double *values, uint32_t count, uint32_t start_row, uint32_t col
) {
    double cell_reads = fit.col_params(0, col);
    double sd_inv_max = fit.global_params(0);
    double clip_min = fit.global_params(1);

    for (uint32_t i = 0; i < count; i++) {
        double mu = cell_reads * fit.row_params(1, start_row + i);
        double theta_inv = fit.row_params(0, start_row + i);
        double sd_i = std::min(sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));
        values[i] = std::max(clip_min, -mu * sd_i);
    }
}

bool SCTransformPearsonTranspose::loadZeroSubtracted(MatrixLoader<double> &loader) {
    if (!loader.load()) return false;

    uint32_t *row_data = loader.rowData();
    double *val_data = loader.valData();
    uint32_t capacity = loader.capacity();

    double gene_beta = fit.col_params(1, currentCol());
    double theta_inv = fit.col_params(0, currentCol());

    double sd_inv_max = fit.global_params(0);
    double clip_max = fit.global_params(2);
    double clip_min = fit.global_params(1);

    for (uint32_t i = 0; i < capacity; i++) {
        double mu = fit.row_params(0, row_data[i]) * gene_beta;
        double sd_i = std::min(sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));

        double zero_val = -mu * sd_i;
        // val = (X-mu)/sd
        double val = val_data[i] * sd_i + zero_val;
        // val = clamp(val, min, max)
        val = std::max(std::min(val, clip_max), clip_min);
        // val = val - max(-mu/sd, clip_min)
        val_data[i] = val - std::max(zero_val, clip_min);
    }
    return true;
}

void SCTransformPearsonTranspose::loadZero(
    double *values, uint32_t count, uint32_t start_row, uint32_t col
) {
    double gene_beta = fit.col_params(1, col);
    double theta_inv = fit.col_params(0, col);
    double sd_inv_max = fit.global_params(0);
    double clip_min = fit.global_params(1);

    for (uint32_t i = 0; i < count; i++) {
        double mu = fit.row_params(0, start_row + i) * gene_beta;
        double sd_i = std::min(sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));
        values[i] = std::max(clip_min, -mu * sd_i);
    }
}

} // end namespace BPCells