#include "Min.h"


namespace BPCells {

bool Min::load() {
    if (!loader.load()) return false;
    
    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = std::min(val_data[i], fit.global_params(0));
    }
    return true;
}


} // end namespace BPCells