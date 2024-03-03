#include "Log1p.h"

namespace BPCells {

bool Log1p::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = log1p(val_data[i]);
    }
    return true;
}

bool Log1pSIMD::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();

    float buf[BPCELLS_VEC_FLOAT_SIZE];
    uint32_t i;
    for (i = 0; i + BPCELLS_VEC_FLOAT_SIZE <= cap; i += BPCELLS_VEC_FLOAT_SIZE) {
        for (uint32_t j = 0; j < BPCELLS_VEC_FLOAT_SIZE; j++) {
            buf[j] = val_data[i + j];
        }
        store_float(buf, log1p_f(load_float(buf)));
        for (uint32_t j = 0; j < BPCELLS_VEC_FLOAT_SIZE; j++) {
            val_data[i + j] = buf[j];
        }
    }
    for (; i < cap; i++) {
        val_data[i] = log1p(val_data[i]);
    }

    return true;
}


bool Expm1::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = expm1(val_data[i]);
    }
    return true;
}

bool Expm1SIMD::load() {
    if (!loader->load()) return false;

    double *val_data = valData();
    const uint32_t cap = capacity();

    float buf[BPCELLS_VEC_FLOAT_SIZE];
    uint32_t i;
    for (i = 0; i + BPCELLS_VEC_FLOAT_SIZE <= cap; i += BPCELLS_VEC_FLOAT_SIZE) {
        for (uint32_t j = 0; j < BPCELLS_VEC_FLOAT_SIZE; j++) {
            buf[j] = val_data[i + j];
        }
        store_float(buf, expm1_f(load_float(buf)));
        for (uint32_t j = 0; j < BPCELLS_VEC_FLOAT_SIZE; j++) {
            val_data[i + j] = buf[j];
        }
    }
    for (; i < cap; i++) {
        val_data[i] = expm1(val_data[i]);
    }

    return true;
}

} // end namespace BPCells
