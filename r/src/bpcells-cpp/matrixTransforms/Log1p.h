// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include "MatrixTransform.h"

namespace BPCells {

class Log1p : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};

// This class uses Sleef SIMD math function to
// speed up Log1p calculation ~4x on machines with AVX2 (x86) or NEON (ARM)
// The one caveat is that this is done in single precision rather than double
class Log1pSIMD : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};

class Expm1 : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};

// This class uses Sleef SIMD math function to speed up calculations
// The one caveat is that this is done in single precision rather than double
class Expm1SIMD : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};

} // end namespace BPCells