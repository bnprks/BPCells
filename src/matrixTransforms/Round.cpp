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
    const uint32_t digits = (uint32_t)fit.global_params(0);

    if(digits == 0) {
      for (uint32_t i = 0; i < cap; i++) {
          val_data[i] = std::nearbyint(val_data[i]);
      }
    }
    else {
      const double factor = pow(10.0, (double)digits);
      const double inv_factor = 1.0 / factor;
      for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = std::nearbyint(val_data[i] * factor) * inv_factor;
      }
    }

    return true;
}

} // end namespace BPCells
