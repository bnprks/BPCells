#include "SCTransform.h"

namespace BPCells { 
// bool SCTransformPearson::nextCol() {
//     if (!MatrixTransformDense::nextCol()) return false;

//     col_mu = fit.row_params.block(1, 0, fit.row_params.rows()-1, fit.row_params.cols()).matrix().transpose() *
//         fit.col_params.matrix().col(currentCol());
// }

// void SCTransformPearson::seekCol(uint32_t col) {
//     MatrixTransformDense::seekCol(col);
//     if (currentCol() < cols()) {
//         col_mu = fit.row_params.block(1, 0, fit.row_params.rows()-1, fit.row_params.cols()).matrix().transpose() *
//             fit.col_params.matrix().col(currentCol());
//     }

// }

bool SCTransformPearson::loadZeroSubtracted(MatrixLoader<double> &loader) {
    if (!loader.load()) return false;
    
    uint32_t *row_data = loader.rowData();
    double *val_data = loader.valData();
    uint32_t capacity = loader.capacity();
    for (uint32_t i = 0; i < capacity; i++) {
        double mu = (fit.row_params.col(row_data[i]).tail(fit.row_params.rows()-1)).matrix()\
            .dot(fit.col_params.matrix().col(currentCol()));
        mu = exp(mu);
        double theta_inv = fit.row_params(0, row_data[i]);
        val_data[i] /= sqrt(mu + mu * mu * theta_inv);
    }
    return true;
}

void SCTransformPearson::loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) {
    Eigen::ArrayXd mu = (fit.row_params.block(1, start_row, fit.row_params.rows()-1, count).matrix().transpose() *
         fit.col_params.matrix().col(col)).array().exp().array();

    
    
    for (uint32_t i = 0; i < count; i++) {
        double m = mu(i);
        double theta_inv = fit.row_params(0, start_row+i);
        values[i] = -m / sqrt(m + m*m*theta_inv);
    }
}

} // end namespace BPCells