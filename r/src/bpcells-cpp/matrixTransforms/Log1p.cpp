// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "Log1p.h"
#include "../simd/math.h"

namespace BPCells {

bool Log1p::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = log1p(val_data[i]);
    }
    return true;
}

bool Log1pSIMD::load() {
    if (!loader->load()) return false;
    BPCells::simd::log1p_downcast(valData(), capacity());
    return true;
}

bool Expm1::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = expm1(val_data[i]);
    }
    return true;
}

bool Expm1SIMD::load() {
    if (!loader->load()) return false;
    BPCells::simd::expm1_downcast(valData(), capacity());
    return true;
}

} // end namespace BPCells
