#pragma once

#include <cmath>
#include "../lib/sleef_wrapper.h"
#include "MatrixTransform.h"

namespace BPCells {

class Round : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};

} // end namespace BPCells
