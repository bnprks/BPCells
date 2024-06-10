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

class Min : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};

class MinByRow : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};

class MinByCol : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};


} // end namespace BPCells