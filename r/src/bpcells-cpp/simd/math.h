// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once


namespace BPCells::simd {

void log1p(float *inout, size_t n);
void log1p(double *inout, size_t n);
void log1p_downcast(double *inout, size_t n);

void expm1(float *inout, size_t n);
void expm1(double *inout, size_t n);
void expm1_downcast(double *inout, size_t n);

void square(float *inout, size_t n);
void square(double *inout, size_t n);
void square_downcast(double *inout, size_t n);

// Return the max value of `in`
uint32_t max(const uint32_t *in, size_t n);

// Write inout[i] = inout[i] + a[i]
void add(uint32_t *inout, const uint32_t *a, size_t n);

// Write inout[i] = inout[i] + a
void add(uint32_t *inout, const int32_t a, size_t n);

// Write inout[i] = inout[i] - a[i]
void sub(uint32_t *inout, const uint32_t *a, size_t n);

}