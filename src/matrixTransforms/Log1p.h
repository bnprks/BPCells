#pragma once

#include "../lib/sleef/sleef_wrapper.h"
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

} // end namespace BPCells