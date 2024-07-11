// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <atomic>
#include "MatrixTransform.h"

namespace BPCells {

// QR decomposition
class LinearResidual : public MatrixTransformDense {
public:
    using MatrixTransformDense::MatrixTransformDense;

    bool loadZeroSubtracted(MatrixLoader<double> &loader) override;
    void loadZero(double *values, uint32_t count, uint32_t start_row, uint32_t col) override;
};

} // namespace BPCells
