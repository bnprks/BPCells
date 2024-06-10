// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "Pow.h"

#include "../simd/math.h"

namespace BPCells {

bool Pow::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();
    const double exponent = fit.global_params(0);

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = std::pow(val_data[i], exponent);
    }
    return true;
}

bool Square::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] *= val_data[i];
    }
    return true;
}


bool SquareSIMD::load() {
    if (!loader->load()) return false;
    
    BPCells::simd::square(valData(), capacity());

    return true;
}

} // end namespace BPCells
