#include "Round.h"
#include <cstdio>

namespace BPCells {

bool Round::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();
    // const uint32_t digits = fit.global_params(0);  digits is unsupported at this time.

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = nearbyint(val_data[i]);
    }
    return true;
}

} // end namespace BPCells
