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
    if (!loader.load()) return false;

    uint32_t *row_data = loader.rowData();
    double *val_data = loader.valData();
    uint32_t capacity = loader.capacity();

    for (uint32_t i = 0; i < capacity; i++) {
        double y0 = 0;
        for (uint32_t j = 0; j < fit.col_params.rows(); j++) {
            y0 += fit.col_params(j, currentCol()) * fit.row_params(j, row_data[i]);
        }
        val_data[i] -= y0;
    }
    return true;
}

void LinearResidual::loadZero(
    double *values, uint32_t count, uint32_t start_row, uint32_t col
) {

    for (uint32_t i = 0; i < count; i++) {
        double y = 0;
        for (uint32_t j = 0; j < fit.col_params.rows(); j++) {
            y -= fit.col_params(j, currentCol()) * fit.row_params(j, start_row + i);
        }
        values[i] = y;
    }
}

} // namespace BPCells


