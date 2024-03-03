#pragma once

#include "../lib/sleef_wrapper.h"
#include "MatrixTransform.h"

namespace BPCells {

class Pow : public MatrixTransform {
  public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};

// This class uses Sleef SIMD math function to speed up calculations
// The one caveat is that this is done in single precision rather than double
class PowSIMD : public MatrixTransform {
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