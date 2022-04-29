#pragma once

#include "MatrixTransform.h"


namespace BPCells {

class Log1p : public MatrixTransform {
public:
    using MatrixTransform::MatrixTransform;

    bool load() override;
};

} // end namespace BPCells