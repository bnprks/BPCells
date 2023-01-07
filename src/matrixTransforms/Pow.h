#pragma once

#include <cmath>
#include "MatrixTransform.h"

namespace BPCells {

class Pow : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};

class Square : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};

} // end namespace BPCells