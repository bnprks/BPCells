#include <atomic>
#include "SCTransform.h"

namespace BPCells {

inline vec_float sd_inv(const vec_float &mu, const vec_float &theta_inv, const vec_float &max_val) {
    // min(max_val, 1/sqrt(mu + mu*mu*theta_inv))
    return min_f(rsqrt_f(fma_f(mul_f(mu, theta_inv), mu, mu)), max_val);
}

SCTransformPearsonSIMD::SCTransformPearsonSIMD(std::unique_ptr<MatrixLoader<double>> &&loader, TransformFit fit)
    : MatrixTransformDense(std::move(loader), fit)
    , theta_inv(fit.row_params.row(0).cast<float>())
    , cell_read_counts(fit.col_params.row(0).cast<float>())
    , gene_beta(fit.row_params.row(1).cast<float>())
    , sd_inv_max(fit.global_params(0))
    , clip_min(fit.global_params(1))
    , clip_max(fit.global_params(2)) {}

bool SCTransformPearsonSIMD::loadZeroSubtracted(MatrixLoader<double> &loader) {
    if (!loader.load()) return false;

    uint32_t *row_data = loader.rowData();
    double *val_data = loader.valData();
    uint32_t capacity = loader.capacity();

    vec_float col_factor = splat_float(cell_read_counts(currentCol()));
    vec_float sd_inv_max = splat_float(this->sd_inv_max);
    vec_float clip_max = splat_float(this->clip_max);
    vec_float clip_min = splat_float(this->clip_min);

    float mu_buf[BPCELLS_VEC_FLOAT_SIZE];
    float theta_inv_buf[BPCELLS_VEC_FLOAT_SIZE];

    uint32_t i;
    for (i = 0; i + BPCELLS_VEC_FLOAT_SIZE <= capacity; i += BPCELLS_VEC_FLOAT_SIZE) {
        for (uint32_t j = 0; j < BPCELLS_VEC_FLOAT_SIZE; j++) {
            mu_buf[j] = gene_beta(row_data[i + j]);
            theta_inv_buf[j] = this->theta_inv(row_data[i + j]);
        }
        vec_float mu = mul_f(col_factor, load_float(mu_buf));
        vec_float theta_inv = load_float(theta_inv_buf);
        vec_float val = load_double_to_float(val_data + i);
        vec_float sd_i = sd_inv(mu, theta_inv, sd_inv_max);
        vec_float zero_val = mul_f(neg_f(mu), sd_i);
        // val = (X-mu)/sd 
        val = add_f(mul_f(val, sd_i), zero_val);
        // val = clamp(val, min, max)
        val = max_f(min_f(val, clip_max), clip_min);
        // val = val - max(-mu/sd, clip_min)
        val = sub_f(val, max_f(zero_val, clip_min));
        
        store_float_to_double(val_data + i, val);
    }
    for (; i < capacity; i++) {
        double mu = cell_read_counts(currentCol()) * gene_beta(row_data[i]);
        double theta_inv = this->theta_inv(row_data[i]);
        double sd_i = std::min(this->sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));

        double zero_val = -mu*sd_i;
        // val = (X-mu)/sd 
        double val = val_data[i] * sd_i + zero_val;
        // val = clamp(val, min, max)
        val = std::max(std::min(val, this->clip_max), this->clip_min);
        // val = val - max(-mu/sd, clip_min)
        val_data[i] = val - std::max(zero_val, this->clip_min);
    }
    return true;
}

void SCTransformPearsonSIMD::loadZero(
    double *values, uint32_t count, uint32_t start_row, uint32_t col
) {
    uint32_t i;
    vec_float col_factor = splat_float(cell_read_counts(col));
    vec_float sd_inv_max = splat_float(this->sd_inv_max);
    vec_float clip_min = splat_float(this->clip_min);

    for (i = 0; i + BPCELLS_VEC_FLOAT_SIZE <= count; i += BPCELLS_VEC_FLOAT_SIZE) {
        vec_float mu = mul_f(col_factor, load_float(gene_beta.data() + start_row + i));
        vec_float theta_inv = load_float(this->theta_inv.data() + start_row + i);
        // val = max(clip_min, -mu / sd)
        vec_float val = mul_f(neg_f(mu), sd_inv(mu, theta_inv, sd_inv_max));
        store_float_to_double(values + i, max_f(val, clip_min));
    }
    for (; i < count; i++) {
        double mu = cell_read_counts(col) * gene_beta(start_row + i);
        double theta_inv = this->theta_inv(start_row + i);
        double sd_i = std::min(this->sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));
        values[i] = std::max(this->clip_min, -mu * sd_i);
    }
}

// Calculate matrix-vector product A*v where A (this) is sparse and B is a dense matrix.
void SCTransformPearsonSIMD::vecMultiplyRightZero(
    Eigen::VectorXd &out, const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    Eigen::VectorXf out_float(out.rows());
    out_float.setZero();

    uint32_t nrows = rows();
    uint32_t ncols = cols();

    vec_float sd_inv_max = splat_float(this->sd_inv_max);
    vec_float clip_min = splat_float(this->clip_min);

    for (uint32_t col = 0; col < ncols; col++) {
        if (user_interrupt != NULL && *user_interrupt) return;
        // Periodically flush our single-precision accumulator to avoid
        // excessive loss of precision during summation
        if (col % 64 == 0) {
            out += out_float.cast<double>();
            out_float.setZero();
        }
        uint32_t row;
        vec_float col_factor = splat_float(cell_read_counts(col));
        vec_float v_vec = splat_float(v(col));
        for (row = 0; row + BPCELLS_VEC_FLOAT_SIZE <= nrows; row += BPCELLS_VEC_FLOAT_SIZE) {
            vec_float mu = mul_f(col_factor, load_float(gene_beta.data() + row));
            vec_float theta_inv = load_float(this->theta_inv.data() + row);
            // val = max(clip_min, -mu / sd)
            vec_float val = mul_f(neg_f(mu), sd_inv(mu, theta_inv, sd_inv_max));
            val = max_f(val, clip_min);

            vec_float out_vec = fma_f(v_vec, val, load_float(out_float.data() + row));
            store_float(out_float.data() + row, out_vec);
        }
        for (; row < nrows; row++) {
            double mu = cell_read_counts(col) * gene_beta(row);
            double theta_inv = this->theta_inv(row);
            double sd_i = std::min(this->sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));
            out_float(row) += v(col) * std::max(this->clip_min, -mu * sd_i);
        }
    }
    // Output final results
    out += out_float.cast<double>();
}

void SCTransformPearsonSIMD::vecMultiplyLeftZero(
    Eigen::VectorXd &out, const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    Eigen::VectorXf v_float(v.cast<float>());

    uint32_t nrows = rows();
    uint32_t ncols = cols();

    vec_float sd_inv_max = splat_float(this->sd_inv_max);
    vec_float clip_min = splat_float(this->clip_min);

    float out_buf[BPCELLS_VEC_FLOAT_SIZE];
    for (uint32_t col = 0; col < ncols; col++) {
        if (user_interrupt != NULL && *user_interrupt) return;
        uint32_t row;
        vec_float col_factor = splat_float(cell_read_counts(col));
        vec_float out_vec = splat_float(0.0);
        for (row = 0; row + BPCELLS_VEC_FLOAT_SIZE <= nrows; row += BPCELLS_VEC_FLOAT_SIZE) {
            vec_float mu = mul_f(col_factor, load_float(gene_beta.data() + row));
            vec_float theta_inv = load_float(this->theta_inv.data() + row);
            // val = max(clip_min, -mu / sd)
            vec_float val = mul_f(neg_f(mu), sd_inv(mu, theta_inv, sd_inv_max));
            val = max_f(val, clip_min);

            vec_float v_vec = load_float(v_float.data() + row);
            out_vec = fma_f(v_vec, val, out_vec);

            // Periodically flush our single-precision accumulator to avoid
            // excessive loss of precision during summation
            if (row % (64 * BPCELLS_VEC_FLOAT_SIZE) == 0) {
                store_float(out_buf, out_vec);
                for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
                    out(col) += out_buf[i];
                }
                out_vec = splat_float(0.0);
            }
        }
        for (; row < nrows; row++) {
            double mu = cell_read_counts(col) * gene_beta(row);
            double theta_inv = this->theta_inv(row);
            double sd_i = std::min(this->sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));
            out(col) += v(row) * std::max(this->clip_min, -mu * sd_i);
        }
        store_float(out_buf, out_vec);
        for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
            out(col) += out_buf[i];
        }
    }
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

    uint32_t *row_data = loader.rowData();
    double *val_data = loader.valData();
    uint32_t capacity = loader.capacity();

    vec_float gene_beta = splat_float(this->gene_beta(currentCol()));
    vec_float theta_inv = splat_float(this->theta_inv(currentCol()));

    vec_float sd_inv_max = splat_float(this->sd_inv_max);
    vec_float clip_max = splat_float(this->clip_max);
    vec_float clip_min = splat_float(this->clip_min);

    float cell_count_buf[BPCELLS_VEC_FLOAT_SIZE];

    uint32_t i;
    for (i = 0; i + BPCELLS_VEC_FLOAT_SIZE <= capacity; i += BPCELLS_VEC_FLOAT_SIZE) {
        for (uint32_t j = 0; j < BPCELLS_VEC_FLOAT_SIZE; j++) {
            cell_count_buf[j] = cell_read_counts(row_data[i + j]);
        }
        vec_float mu = mul_f(gene_beta, load_float(cell_count_buf));
        vec_float val = load_double_to_float(val_data + i);
        vec_float sd_i = sd_inv(mu, theta_inv, sd_inv_max);
        vec_float zero_val = mul_f(neg_f(mu), sd_i);
        // val = (X-mu)/sd 
        val = add_f(mul_f(val, sd_i), zero_val);
        // val = clamp(val, min, max)
        val = max_f(min_f(val, clip_max), clip_min);
        // val = val - max(-mu/sd, clip_min)
        val = sub_f(val, max_f(zero_val, clip_min));
        store_float_to_double(val_data + i, val);
    }
    for (; i < capacity; i++) {
        double mu = cell_read_counts(row_data[i]) * this->gene_beta(currentCol());
        double theta_inv = this->theta_inv(currentCol());
        double sd_i = std::min(this->sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));
        double zero_val = -mu*sd_i;
        // val = (X-mu)/sd 
        double val = val_data[i] * sd_i + zero_val;
        // val = clamp(val, min, max)
        val = std::max(std::min(val, this->clip_max), this->clip_min);
        // val = val - max(-mu/sd, clip_min)
        val_data[i] = val - std::max(zero_val, this->clip_min);
    }
    return true;
}

void SCTransformPearsonTransposeSIMD::loadZero(
    double *values, uint32_t count, uint32_t start_row, uint32_t col
) {
    uint32_t i;
    vec_float gene_beta = splat_float(this->gene_beta(currentCol()));
    vec_float theta_inv = splat_float(this->theta_inv(currentCol()));
    vec_float sd_inv_max = splat_float(this->sd_inv_max);
    vec_float clip_min = splat_float(this->clip_min);

    for (i = 0; i + BPCELLS_VEC_FLOAT_SIZE <= count; i += BPCELLS_VEC_FLOAT_SIZE) {
        vec_float mu = mul_f(gene_beta, load_float(cell_read_counts.data() + start_row + i));
        // val = max(clip_min, -mu / sd)
        vec_float val = mul_f(neg_f(mu), sd_inv(mu, theta_inv, sd_inv_max));
        store_float_to_double(values + i, max_f(val, clip_min));
    }
    for (; i < count; i++) {
        double mu = cell_read_counts(start_row + i) * this->gene_beta(col);
        double theta_inv = this->theta_inv(col);
        double sd_i = std::min(this->sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));
        values[i] = std::max(this->clip_min, -mu * sd_i);
    }
}

void SCTransformPearsonTransposeSIMD::vecMultiplyLeftZero(
    Eigen::VectorXd &out, const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    // To convert for transpose, all we need to do is flip the `row` and `col` variables
    // and swap vecMultiplyLeftZero with vecMultiplyRightZero
    Eigen::VectorXf out_float(out.rows());
    out_float.setZero();

    uint32_t nrows = rows();
    uint32_t ncols = cols();

    vec_float sd_inv_max = splat_float(this->sd_inv_max);
    vec_float clip_min = splat_float(this->clip_min);

    for (uint32_t row = 0; row < nrows; row++) {
        if (user_interrupt != NULL && row % 128 == 0 && *user_interrupt) return;
        // Periodically flush our single-precision accumulator to avoid
        // excessive loss of precision during summation
        if (row % 64 == 0) {
            out += out_float.cast<double>();
            out_float.setZero();
        }
        uint32_t col;
        vec_float cell_reads = splat_float(cell_read_counts(row));
        vec_float v_vec = splat_float(v(row));
        for (col = 0; col + BPCELLS_VEC_FLOAT_SIZE <= ncols; col += BPCELLS_VEC_FLOAT_SIZE) {
            vec_float mu = mul_f(cell_reads, load_float(gene_beta.data() + col));
            vec_float theta_inv = load_float(this->theta_inv.data() + col);
            // val = max(clip_min, -mu / sd)
            vec_float val = mul_f(neg_f(mu), sd_inv(mu, theta_inv, sd_inv_max));
            val = max_f(val, clip_min);

            vec_float out_vec = fma_f(v_vec, val, load_float(out_float.data() + col));
            store_float(out_float.data() + col, out_vec);
        }
        for (; col < ncols; col++) {
            double mu = cell_read_counts(row) * gene_beta(col);
            double theta_inv = this->theta_inv(col);
            double sd_i = std::min(this->sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));
            out_float(col) += v(row) * std::max(this->clip_min, -mu * sd_i);
        }
    }
    // Output final results
    out += out_float.cast<double>();
}

void SCTransformPearsonTransposeSIMD::vecMultiplyRightZero(
    Eigen::VectorXd &out, const Eigen::Map<Eigen::VectorXd> v, std::atomic<bool> *user_interrupt
) {
    // To convert for transpose, all we need to do is flip the `row` and `col` variables
    // and swap vecMultiplyLeftZero with vecMultiplyRightZero
    Eigen::VectorXf v_float(v.cast<float>());

    uint32_t nrows = rows();
    uint32_t ncols = cols();

    vec_float sd_inv_max = splat_float(this->sd_inv_max);
    vec_float clip_min = splat_float(this->clip_min);

    float out_buf[BPCELLS_VEC_FLOAT_SIZE];
    for (uint32_t row = 0; row < nrows; row++) {
        if (user_interrupt != NULL && row % 128 == 0 && *user_interrupt) return;
        uint32_t col;
        vec_float cell_reads = splat_float(cell_read_counts(row));
        vec_float out_vec = splat_float(0.0);
        for (col = 0; col + BPCELLS_VEC_FLOAT_SIZE <= ncols; col += BPCELLS_VEC_FLOAT_SIZE) {
            vec_float mu = mul_f(cell_reads, load_float(gene_beta.data() + col));
            vec_float theta_inv = load_float(this->theta_inv.data() + col);
            // val = max(clip_min, -mu / sd)
            vec_float val = mul_f(neg_f(mu), sd_inv(mu, theta_inv, sd_inv_max));
            val = max_f(val, clip_min);

            vec_float v_vec = load_float(v_float.data() + col);
            out_vec = fma_f(v_vec, val, out_vec);

            // Periodically flush our single-precision accumulator to avoid
            // excessive loss of precision during summation
            if (col % (64 * BPCELLS_VEC_FLOAT_SIZE) == 0) {
                store_float(out_buf, out_vec);
                for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
                    out(row) += out_buf[i];
                }
                out_vec = splat_float(0.0);
            }
        }
        for (; col < ncols; col++) {
            double mu = cell_read_counts(row) * gene_beta(col);
            double theta_inv = this->theta_inv(col);
            double sd_i = std::min(this->sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));
            out(row) += v(col) * std::max(this->clip_min, -mu * sd_i);
        }
        store_float(out_buf, out_vec);
        for (int i = 0; i < BPCELLS_VEC_FLOAT_SIZE; i++) {
            out(row) += out_buf[i];
        }
    }
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

    for (uint32_t i = 0; i < capacity; i ++) {
        double mu = cell_reads * fit.row_params(1, row_data[i]);
        double theta_inv = fit.row_params(0,row_data[i]);
        double sd_i = std::min(sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));
        
        double zero_val = -mu*sd_i;
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

    for (uint32_t i = 0; i < capacity; i ++) {
        double mu = fit.row_params(0, row_data[i]) * gene_beta;
        double sd_i = std::min(sd_inv_max, 1 / sqrt(mu + mu * mu * theta_inv));

        double zero_val = -mu*sd_i;
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