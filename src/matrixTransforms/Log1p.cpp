#include "Log1p.h"


namespace BPCells {

bool Log1p::load() {
    if (!loader.load()) return false;
    
    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = std::log1p(val_data[i]);
    }
    return true;
}

} // end namespace BPCells