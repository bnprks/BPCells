// Copyright 2024 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

// Compilation speed is a concern, hence why each bp128 function is split up (better
// parallelization). This is because 32-step unrolled loop * 32 bitwidths * 6 functions * 6-8
// architectures = a lot of code When in debug mode AND optimize is turned off, disable dynamic
// dispatch and compile for just baseline architecture
#if !defined(NDEBUG) && !defined(__OPTIMIZE__)
#define HWY_COMPILE_ONLY_STATIC 1
#endif

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "simd/bp128/bp128-Nx128.cpp"
#include <hwy/foreach_target.h>

#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::bp128::HWY_NAMESPACE {

void unpack_bp128(const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
void pack_bp128(const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
uint32_t maxbits_bp128(const uint32_t *HWY_RESTRICT in);

extern void unpack_d1(uint32_t initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
extern void pack_d1(uint32_t initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
extern uint32_t maxbits_d1(uint32_t initvalue, const uint32_t *HWY_RESTRICT in);

extern void unpack_d1z(uint32_t initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
extern void pack_d1z(uint32_t initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
extern uint32_t maxbits_d1z(uint32_t initvalue, const uint32_t *HWY_RESTRICT in);

extern void unpack_FOR(uint32_t offset, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
extern void pack_FOR(uint32_t offset, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
extern uint32_t maxbits_FOR(uint32_t offset, const uint32_t *HWY_RESTRICT in);

extern void unpack_diff(const uint32_t *HWY_RESTRICT ref, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
extern void pack_diff(const uint32_t *HWY_RESTRICT ref, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t bit);
extern uint32_t maxbits_diff(const uint32_t *HWY_RESTRICT ref, const uint32_t *HWY_RESTRICT in);

void pack_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit) {
    for (uint32_t i = 0; i < n; i++) {
        uint32_t b = maxbits_bp128(in);
        pack_bp128(in, out, b);
        
        bit[i] = b;

        in += 128,
        out += b*4;
    }
}

void unpack_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit) {
    for (uint32_t i = 0; i < n; i++) {
        unpack_bp128(in, out, bit[i]);

        in += bit[i]*4;
        out += 128;
    }
}

void pack_d1_Nx128(uint32_t n, uint32_t *HWY_RESTRICT initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit) {
    for (uint32_t i = 0; i < n; i++) {
        uint32_t b = maxbits_d1(in[0], in);
        pack_d1(in[0], in, out, b);
        
        bit[i] = b;
        initvalue[i] = in[0];

        in += 128,
        out += b*4;
    }
}

void unpack_d1_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit) {
    for (uint32_t i = 0; i < n; i++) {
        unpack_d1(initvalue[i], in, out, bit[i]);

        in += bit[i]*4;
        out += 128;
    }
}

void pack_d1z_Nx128(uint32_t n, uint32_t *HWY_RESTRICT initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit) {
    for (uint32_t i = 0; i < n; i++) {
        uint32_t b = maxbits_d1z(in[0], in);
        pack_d1z(in[0], in, out, b);
        
        bit[i] = b;
        initvalue[i] = in[0];

        in += 128,
        out += b*4;
    }
}

void unpack_d1z_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit) {
    for (uint32_t i = 0; i < n; i++) {
        unpack_d1z(initvalue[i], in, out, bit[i]);

        in += bit[i]*4;
        out += 128;
    }
}

void pack_FOR_Nx128(uint32_t n, uint32_t offset, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit) {
    for (uint32_t i = 0; i < n; i++) {
        uint32_t b = maxbits_FOR(offset, in);
        pack_FOR(offset, in, out, b);
        
        bit[i] = b;

        in += 128,
        out += b*4;
    }
}

void unpack_FOR_Nx128(uint32_t n, uint32_t offset, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit) {
    for (uint32_t i = 0; i < n; i++) {
        unpack_FOR(offset, in, out, bit[i]);

        in += bit[i]*4;
        out += 128;
    }
}

void pack_diff_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT ref, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit) {
    for (uint32_t i = 0; i < n; i++) {
        uint32_t b = maxbits_diff(ref + 128*i, in);
        pack_diff(ref + 128*i, in, out, b);
        
        bit[i] = b;

        in += 128,
        out += b*4;
    }
}


void unpack_diff_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT ref, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit) {
    for (uint32_t i = 0; i < n; i++) {
        unpack_diff(ref + 128*i, in, out, bit[i]);

        in += bit[i]*4;
        out += 128;
    }
}


} // namespace BPCells::simd::bp128::HWY_NAMESPACE

HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace BPCells::simd::bp128 {

HWY_EXPORT(pack_Nx128);
HWY_EXPORT(unpack_Nx128);
HWY_EXPORT(pack_d1_Nx128);
HWY_EXPORT(unpack_d1_Nx128);
HWY_EXPORT(pack_d1z_Nx128);
HWY_EXPORT(unpack_d1z_Nx128);
HWY_EXPORT(pack_FOR_Nx128);
HWY_EXPORT(unpack_FOR_Nx128);
HWY_EXPORT(pack_diff_Nx128);
HWY_EXPORT(unpack_diff_Nx128);

void pack_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit) {
    HWY_DYNAMIC_DISPATCH(pack_Nx128)(n, in, out, bit);
}

void unpack_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit) {
    HWY_DYNAMIC_DISPATCH(unpack_Nx128)(n, in, out, bit);
}

void pack_d1_Nx128(uint32_t n, uint32_t *HWY_RESTRICT initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit) {
    HWY_DYNAMIC_DISPATCH(pack_d1_Nx128)(n, initvalue, in, out, bit);
}

void unpack_d1_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit) {
    HWY_DYNAMIC_DISPATCH(unpack_d1_Nx128)(n, initvalue, in, out, bit);
}

void pack_d1z_Nx128(uint32_t n, uint32_t *HWY_RESTRICT initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit) {
    HWY_DYNAMIC_DISPATCH(pack_d1z_Nx128)(n, initvalue, in, out, bit);
}

void unpack_d1z_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT initvalue, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit) {
    HWY_DYNAMIC_DISPATCH(unpack_d1z_Nx128)(n, initvalue, in, out, bit);
}

void pack_FOR_Nx128(uint32_t n, uint32_t offset, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit) {
    HWY_DYNAMIC_DISPATCH(pack_FOR_Nx128)(n, offset, in, out, bit);
}

void unpack_FOR_Nx128(uint32_t n, uint32_t offset, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit) {
    HWY_DYNAMIC_DISPATCH(unpack_FOR_Nx128)(n, offset, in, out, bit);
}

void pack_diff_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT ref, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, uint32_t *HWY_RESTRICT bit) {
    HWY_DYNAMIC_DISPATCH(pack_diff_Nx128)(n, ref, in, out, bit);
}

void unpack_diff_Nx128(uint32_t n, const uint32_t *HWY_RESTRICT ref, const uint32_t *HWY_RESTRICT in, uint32_t *HWY_RESTRICT out, const uint32_t *HWY_RESTRICT bit) {
    HWY_DYNAMIC_DISPATCH(unpack_diff_Nx128)(n, ref, in, out, bit);
}

} // namespace BPCells::simd::bp128
#endif

