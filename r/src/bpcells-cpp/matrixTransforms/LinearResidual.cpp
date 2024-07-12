// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "LinearResidual.h"
#include <atomic>

namespace BPCells {

bool LinearResidual::loadZeroSubtracted(MatrixLoader<double> &loader) {
    return loader.load();
}

void LinearResidual::loadZero(
    double *values, uint32_t count, uint32_t start_row, uint32_t col
) {
    for (uint32_t i = 0; i < count; i++) {
        double zero_val = 0;
        // col_params = t(Q)
        // row_params = Qty
        // val = X - t(Qty) %*% t(Q)
        // val = X - t(row_params) %*% col_params
        Eigen::ArrayXXd prod_res = fit.row_params.col(start_row + i) * fit.col_params.col(col);
        zero_val = -prod_res.sum();
        values[i] = zero_val;
    }
}

} // namespace BPCells
