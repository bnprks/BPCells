#pragma once

#include "MatrixTransform.h"

namespace BPCells {

class Min : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};

} // end namespace BPCells