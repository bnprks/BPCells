// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "Min.h"

namespace BPCells {

bool Min::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = std::min(val_data[i], fit.global_params(0));
    }
    return true;
}

bool MinByRow::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    uint32_t *row_data = rowData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = std::min(val_data[i], fit.row_params(0, row_data[i]));
    }
    return true;
}

bool MinByCol::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();

    double col_min = fit.col_params(0, currentCol());

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = std::min(val_data[i], col_min);
    }
    return true;
}

} // end namespace BPCells
