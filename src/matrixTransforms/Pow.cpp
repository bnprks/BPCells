#include "Pow.h"

namespace BPCells {

bool Pow::load() {
    if (!loader.load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();
    const double exponent = fit.global_params(0);

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = std::pow(val_data[i], exponent);
    }
    return true;
}

bool Square::load() {
    if (!loader.load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] *= val_data[i];
    }
    return true;
}

} // end namespace BPCells
