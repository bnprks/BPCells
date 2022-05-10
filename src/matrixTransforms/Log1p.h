#pragma once

#include "MatrixTransform.h"
#include "../lib/sleef/sleef_wrapper.h"

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

// This class keeps a small cache of the smallest input values in each column,
// and caches the results of log1p on them in a length-4 array.
// In common normalization cases, this cache 
class Log1pCache : public MatrixTransform {
private:
    double in[4], out[4];

public:
    using MatrixTransform::MatrixTransform;

    bool nextCol() override;
    bool load() override;
};

} // end namespace BPCells