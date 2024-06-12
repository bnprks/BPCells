// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "simd/dense-multiply.cpp"
#include <hwy/foreach_target.h>

#include <hwy/contrib/math/math-inl.h>
#include <hwy/highway.h>


HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::HWY_NAMESPACE {


/// @brief Helper function for denseMultiplyRight inner loop for sparse matrix A * dense matrix B
/// @param row_data Row index for each entry (length count)
/// @param val_data Value for each entry (length count)
/// @param count Number of matrix entries
/// @param res Dense matrix in memory for the result (row-major, entries index [c, c+dim) cover row c)
/// @param B_row Values in the relevant row of B (length dim)
/// @param dim Number of columns in the matrix B
void denseMultiplyRightHelper(const uint32_t *HWY_RESTRICT row_data, const double *HWY_RESTRICT val_data, uint32_t count, double *HWY_RESTRICT res, const double *HWY_RESTRICT B_row, uint32_t dim) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<double> t;
    const size_t N = Lanes(t);
    // Loop ordering strategy:
    // Principles:
    //   - Expect up to 10 outstanding memory requests at a time: https://lemire.me/blog/2018/11/05/measuring-the-memory-level-parallelism-of-a-system-using-a-small-c-program/
    //   - L1 cache of 32KB (512 lines) is probably fewer than count
    //   - Compiler doesn't know row_data has no repeats, so stores across entries might have dependencies
    //   - Therefore, try to allow for multiple outstanding memory requests by working across `dim` in the inner loop
    //   - Assume at least 16 vector registers available
    // Loop layout: 
    //   - individual entries in block of 32 (to keep everything easily within cache while amortizing loads from B_row)
    //   - Dims in blocks (to increase memory parallelism)
    // Results:
    //   - On AVX2, I see ~1.2x faster than non-vectorized code. Less than I was hoping for, but enough to not throw away this work
    //   - In disk-backed benchmark I see just ~1.1x speedup
    const size_t ENTRY_BLOCK = 32;
    const size_t DIM_BLOCK = 4;
    for (uint32_t ee = 0; ee < count; ee += ENTRY_BLOCK) {
        uint32_t dd = 0;
        // Try a block of vectors at a time (Up to 3 memory ops in flight per loop)
        for (;dd + DIM_BLOCK*N <= dim; dd += DIM_BLOCK*N) {
            auto b0 = LoadU(t, B_row + dd + 0*N);
            auto b1 = LoadU(t, B_row + dd + 1*N);
            auto b2 = LoadU(t, B_row + dd + 2*N);
            auto b3 = LoadU(t, B_row + dd + 3*N);
            for (uint32_t e = ee; e < ee + ENTRY_BLOCK && e < count; e++) {
                auto val = Set(t, val_data[e]);
                auto r0 = LoadU(t, res + dim*row_data[e] + dd + 0*N);
                auto r1 = LoadU(t, res + dim*row_data[e] + dd + 1*N);
                auto r2 = LoadU(t, res + dim*row_data[e] + dd + 2*N);
                auto r3 = LoadU(t, res + dim*row_data[e] + dd + 3*N);
                
                StoreU(MulAdd(val, b0, r0), t, res + dim*row_data[e] + dd + 0*N);
                StoreU(MulAdd(val, b1, r1), t, res + dim*row_data[e] + dd + 1*N);
                StoreU(MulAdd(val, b2, r2), t, res + dim*row_data[e] + dd + 2*N);
                StoreU(MulAdd(val, b3, r3), t, res + dim*row_data[e] + dd + 3*N);       
            }
        }
        // Finish out 1 vector at a time
        for (;dd + N <= dim; dd += N) {
            const uint32_t d = dd;
            for (uint32_t e = ee; e < ee + ENTRY_BLOCK && e < count; e++) {
                auto val = Set(t, val_data[e]);
                auto b = LoadU(t, B_row + d);
                auto r = LoadU(t, res + dim*row_data[e] + d);
                StoreU(MulAdd(val, b, r), t, res + dim*row_data[e] + d);
            }
        }
        // Finish out 1 element at a time
        for (;dd < dim; dd += 1) {
            const uint32_t d = dd;
            for (uint32_t e = ee; e < ee + ENTRY_BLOCK && e < count; e++) {
                res[dim*row_data[e] + d] += val_data[e] * B_row[d];
            }
        }
    }
}


/// @brief Helper function for denseMultiplyRight inner loop for sparse matrix A * dense matrix B
/// @param row_data Row index for each entry (length count)
/// @param val_data Value for each entry (length count)
/// @param count Number of matrix entries
/// @param res_col Values in the relevant col of res (length dim)
/// @param B Dense matrix in memory holding B (column-major; entries index [c, c+dim) cover column c)
/// @param dim Number of rows in the matrix B
void denseMultiplyLeftHelper(const uint32_t *row_data, const double *val_data, uint32_t count, double *res_col, const double *B, uint32_t dim) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<double> t;
    const size_t N = Lanes(t);

    // See looping logic from denseMultiplyRightHelper
    // Only difference is that here B requires a load while res can stay in registers,
    // so our writes are decreased but otherwise the blocking structure should still be favorable
    // Results:
    //   - On AVX2, I see ~1.5x faster than non-vectorized code. Less than I was hoping for, but enough to not throw away this work
    //   - In disk-backed benchmark I see only ~1.2x speedup
    const size_t ENTRY_BLOCK = 32;
    const size_t DIM_BLOCK = 4;
    
    for (uint32_t ee = 0; ee < count; ee += ENTRY_BLOCK) {
        uint32_t dd = 0;
        // Try a block of vectors at a time (Up to 3 memory ops in flight per loop)
        for (;dd + DIM_BLOCK*N <= dim; dd += DIM_BLOCK*N) {
            // Load from memory into (hopefully) registers
            auto r0 = LoadU(t, res_col + dd + 0*N);
            auto r1 = LoadU(t, res_col + dd + 1*N);
            auto r2 = LoadU(t, res_col + dd + 2*N);
            auto r3 = LoadU(t, res_col + dd + 3*N);
            
            // Perform updates for block of entries
            for (uint32_t e = ee; e < ee + ENTRY_BLOCK && e < count; e++) {
                auto val = Set(t, val_data[e]);
                auto b0 = LoadU(t, B + dim*row_data[e] + dd + 0*N);
                auto b1 = LoadU(t, B + dim*row_data[e] + dd + 1*N);
                auto b2 = LoadU(t, B + dim*row_data[e] + dd + 2*N);
                auto b3 = LoadU(t, B + dim*row_data[e] + dd + 3*N);
                r0 = MulAdd(val, b0, r0);
                r1 = MulAdd(val, b1, r1);
                r2 = MulAdd(val, b2, r2);
                r3 = MulAdd(val, b3, r3);
            }
            // Store back to memory
            StoreU(r0, t, res_col + dd + 0*N);
            StoreU(r1, t, res_col + dd + 1*N);
            StoreU(r2, t, res_col + dd + 2*N);
            StoreU(r3, t, res_col + dd + 3*N);
        }

        // Finish out 1 vector at a time
        for (;dd + N <= dim; dd += N) {
            const uint32_t d = dd;
            auto res = LoadU(t, res_col + d);
            for (uint32_t e = ee; e < ee + ENTRY_BLOCK && e < count; e++) {
                auto val = Set(t, val_data[e]);
                auto b = LoadU(t, B + dim*row_data[e] + d);
                res = MulAdd(val, b, res);
            }
            StoreU(res, t, res_col + d);
        }
        // Finish out 1 element at a time
        for (;dd < dim; dd += 1) {
            const uint32_t d = dd;
            auto res = res_col[d];
            for (uint32_t e = ee; e < ee + ENTRY_BLOCK && e < count; e++) {
                res += val_data[e] * B[dim*row_data[e] + d];
            }
            res_col[d] = res;
        }
    }
}


} // namespace BPCells::simd::HWY_NAMESPACE
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace BPCells::simd {

HWY_EXPORT(denseMultiplyRightHelper);
HWY_EXPORT(denseMultiplyLeftHelper);

void denseMultiplyRightHelper(const uint32_t *row_data, const double *val_data, uint32_t count, double *res, const double *B_row, uint32_t dim) {
    return HWY_DYNAMIC_DISPATCH(denseMultiplyRightHelper)(row_data, val_data, count, res, B_row, dim);
}
void denseMultiplyLeftHelper(const uint32_t *row_data, const double *val_data, uint32_t count, double *res_col, const double *B, uint32_t dim) {
    return HWY_DYNAMIC_DISPATCH(denseMultiplyLeftHelper)(row_data, val_data, count, res_col, B, dim);
}

} // namespace BPCells::simd
#endif