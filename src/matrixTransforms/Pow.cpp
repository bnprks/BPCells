#include "Pow.h"

namespace BPCells {

bool Pow::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();
    const double exponent = fit.global_params(0);

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = std::pow(val_data[i], exponent);
    }
    return true;
}

bool PowSIMD::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();

    float buf[BPCELLS_VEC_FLOAT_SIZE];
    vec_float exp = splat_float(fit.global_params(0));
    uint32_t i;
    for (i = 0; i + BPCELLS_VEC_FLOAT_SIZE <= cap; i += BPCELLS_VEC_FLOAT_SIZE) {
        for (uint32_t j = 0; j < BPCELLS_VEC_FLOAT_SIZE; j++) {
            buf[j] = val_data[i + j];
        }
        store_float(buf, pow_f(load_float(buf), exp));
        for (uint32_t j = 0; j < BPCELLS_VEC_FLOAT_SIZE; j++) {
            val_data[i + j] = buf[j];
        }
    }
    for (; i < cap; i++) {
        val_data[i] = powf(val_data[i], fit.global_params(0));
    }

    return true;
}

bool Square::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] *= val_data[i];
    }
    return true;
}

} // end namespace BPCells
