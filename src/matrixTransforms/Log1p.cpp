#include "Log1p.h"


namespace BPCells {

bool Log1p::load() {
    if (!loader.load()) return false;
    
    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        val_data[i] = (double) log1pf((float) val_data[i]);
    }
    return true;
}

bool Log1pSIMD::load() {
    if (!loader.load()) return false;
    
    double *val_data = valData();
    const uint32_t cap = capacity();

    float buf[BPCELLS_F32_VEC_SIZE];
    uint32_t i;
    for (i = 0; i + BPCELLS_F32_VEC_SIZE <= cap; i += BPCELLS_F32_VEC_SIZE) {
        for (uint32_t j = 0; j < BPCELLS_F32_VEC_SIZE; j++) {
            buf[j] = val_data[i+j];
        }
        bpcells_log1pf_vec(buf, buf);
        for (uint32_t j = 0; j < BPCELLS_F32_VEC_SIZE; j++) {
            val_data[i+j] = buf[j];
        }
    }
    for (;i < cap; i++) {
        val_data[i] = log1p(val_data[i]);
    }

    return true;
}

bool Log1pCache::nextCol() {
    in[0] = INFINITY;
    in[1] = INFINITY;
    in[2] = INFINITY;
    in[3] = INFINITY;
    return loader.nextCol();
}

bool Log1pCache::load() {
    if (!loader.load()) return false;
    
    double *val_data = valData();
    const uint32_t cap = capacity();

    for (uint32_t i = 0; i < cap; i++) {
        // Look up if the value is in cache
        bool had_hit = val_data[i] == in[0] | val_data[i] == in[1] | val_data[i] == in[2] | val_data[i] == in[3];
        double res = 0;
        
        if (had_hit) {
            res = val_data[i] == in[0] ? out[0] : res;
            res = val_data[i] == in[1] ? out[1] : res;
            res = val_data[i] == in[2] ? out[2] : res;
            res = val_data[i] == in[3] ? out[3] : res;
            
            val_data[i] = res; 
            continue;
        } 
        
        // No hit, check if we need to add to cache
        res = log1p(val_data[i]);
        if (val_data[i] < in[3]) {
            uint32_t j;
            for (j = 3; j > 0 && val_data[i] < in[j-1]; j--) {
                in[j] = in[j-1];
                out[j] = out[j-1];
            }
            in[j] = val_data[i];
            out[j] = res;
        }
        val_data[i] = res;
    }
    return true;
}


} // end namespace BPCells