#include "Binarize.h"

namespace BPCells {

bool Binarize::load() {
    if (!loader->load()) return false;

    /*
    ** Set value to one if it's greater than or equal to
    ** threshold; otherwise set it to zero.
    */
    double *val_data = valData();
    const uint32_t cap = capacity();
    const double threshold = fit.global_params(0);

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = val_data[i] >= threshold ? 1 : 0;
    }
    return true;
}

} // end namespace BPCells
