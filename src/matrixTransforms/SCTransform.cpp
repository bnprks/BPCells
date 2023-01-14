#include "SCTransform.h"

namespace BPCells {
// bool SCTransformPearson::nextCol() {
//     if (!MatrixTransformDense::nextCol()) return false;

//     col_mu = fit.row_params.block(1, 0, fit.row_params.rows()-1,
//     fit.row_params.cols()).matrix().transpose() *
//         fit.col_params.matrix().col(currentCol());
// }

// void SCTransformPearson::seekCol(uint32_t col) {
//     MatrixTransformDense::seekCol(col);
//     if (currentCol() < cols()) {
//         col_mu = fit.row_params.block(1, 0, fit.row_params.rows()-1,
//         fit.row_params.cols()).matrix().transpose() *
//             fit.col_params.matrix().col(currentCol());
//     }

// }

inline vec_float sd_inv(const vec_float &mu, const vec_float &theta_inv) {
    return rsqrt_f(fma_f(mul_f(mu, theta_inv), mu, mu));
}

SCTransformPearson::SCTransformPearson(MatrixLoader<double> &loader, TransformFit fit)
    : MatrixTransformDense(loader, fit)
    , theta_inv(fit.row_params.row(0).cast<float>())
    , col_vec((fit.row_params(2,0) * fit.col_params.row(1)).exp().cast<float>())
    , row_vec(fit.row_params.row(1).exp().cast<float>()) {}

void SCTransformPearson::ensure_cached_mu(uint32_t col) {
    if (cached_col == col || col >= col_vec.cols()) {
        return;
    }
    col_mu = (row_vec.matrix() * col_vec.col(col).matrix()).array().exp();
    cached_col = col;
}

bool SCTransformPearson::loadZeroSubtracted(MatrixLoader<double> &loader) {
    if (!loader.load()) return false;
    //ensure_cached_mu(currentCol());

    uint32_t *row_data = loader.rowData();
    double *val_data = loader.valData();
    uint32_t capacity = loader.capacity();

    vec_float col_factor = splat_float(col_vec(currentCol()));
    float mu_buf[BPCELLS_VEC_FLOAT_SIZE];
    float theta_inv_buf[BPCELLS_VEC_FLOAT_SIZE];
    
    uint32_t i;
    for (i = 0; i + BPCELLS_VEC_FLOAT_SIZE <= capacity; i += BPCELLS_VEC_FLOAT_SIZE) {
        for (uint32_t j = 0; j < BPCELLS_VEC_FLOAT_SIZE; j++) {
            mu_buf[j] = row_vec(row_data[i+j]);
            theta_inv_buf[j] = this->theta_inv(row_data[i+j]);
        }
        vec_float mu = mul_f(col_factor, load_float(mu_buf));
        vec_float theta_inv = load_float(theta_inv_buf);
        vec_float val = load_double_to_float(val_data + i);
        store_float_to_double(val_data + i, mul_f(val, sd_inv(mu, theta_inv)));
    }
    for (; i < capacity; i++) {
        double mu = col_vec(currentCol()) * row_vec(row_data[i]);
        double theta_inv = this->theta_inv(row_data[i]);
        val_data[i] /= sqrt(mu + mu * mu * theta_inv);
    }
    return true;
}

void SCTransformPearson::loadZero(
    double *values, uint32_t count, uint32_t start_row, uint32_t col
) {
    // ensure_cached_mu(col);

    uint32_t i;
    vec_float col_factor = splat_float(col_vec(col));
    for (i = 0; i + BPCELLS_VEC_FLOAT_SIZE <= count; i += BPCELLS_VEC_FLOAT_SIZE) {
        vec_float mu = mul_f(col_factor, load_float(row_vec.data() + start_row + i));
        vec_float theta_inv = load_float(this->theta_inv.data() + start_row + i);
        store_float_to_double(values + i, mul_f(neg_f(mu), sd_inv(mu, theta_inv)));
    }
    for (; i < count; i++) {
        double mu = col_vec(col) * row_vec(start_row + i);
        double theta_inv = this->theta_inv(start_row + i);
        values[i] = -mu / sqrt(mu + mu * mu * theta_inv);
    }
    // Eigen::internal::set_is_malloc_allowed(true);
}

} // end namespace BPCells