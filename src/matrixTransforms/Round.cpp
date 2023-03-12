#include "Round.h"
#include <cfenv>
#include <cmath>

namespace BPCells {

bool Round::load() {
    if (!loader->load()) return false;

    // Set rounding mode
    if (std::fegetround() != FE_TONEAREST) std::fesetround(FE_TONEAREST);

    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = std::nearbyint(val_data[i]);
    }
    return true;
}

} // end namespace BPCells
