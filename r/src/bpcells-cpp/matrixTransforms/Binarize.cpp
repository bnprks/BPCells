// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "Binarize.h"
#include <cstdio>

namespace BPCells {

bool Binarize::load() {
    if (!loader->load()) return false;

    /*
    ** Set value to one if it's greater than or equal to
    ** threshold; otherwise set it to zero.
    */
    double *val_data = valData();
    const uint32_t cap = capacity();
    const double threshold = fit.global_params(0);
    const uint32_t strict_inequality = (uint32_t)fit.global_params(1);

    if(strict_inequality == 1) {
      for (uint32_t i = 0; i < cap; i++) {
          val_data[i] = val_data[i] > threshold ? 1 : 0;
      }
    }
    else {
      for (uint32_t i = 0; i < cap; i++) {
          val_data[i] = val_data[i] >= threshold ? 1 : 0;
      }
    }
    return true;
}

} // end namespace BPCells
