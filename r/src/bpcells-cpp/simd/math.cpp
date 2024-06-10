// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "simd/math.cpp"
#include <hwy/foreach_target.h>

#include <hwy/contrib/math/math-inl.h>
#include <hwy/highway.h>

#include "transform-inl.h"

HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::HWY_NAMESPACE {

// Log1p
void log1p_float(float *inout, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    transform1<float>(inout, n, [](auto d, const auto v) HWY_ATTR { return Log1p(d, v); });
}

void log1p_double(double *inout, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    transform1<double>(inout, n, [](auto d, const auto v) HWY_ATTR { return Log1p(d, v); });
}

void log1p_downcast(double *inout, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    transform1_downcast(inout, n, [](auto d, const auto v) HWY_ATTR { return Log1p(d, v); });
}

// Expm1
void expm1_float(float *inout, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    transform1<float>(inout, n, [](auto d, const auto v) HWY_ATTR { return Expm1(d, v); });
}

void expm1_double(double *inout, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    transform1<double>(inout, n, [](auto d, const auto v) HWY_ATTR { return Expm1(d, v); });
}

void expm1_downcast(double *inout, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    transform1_downcast(inout, n, [](auto d, const auto v) HWY_ATTR { return Expm1(d, v); });
}

// Square
void square_float(float *inout, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    transform1<float>(inout, n, [](auto d, const auto v) HWY_ATTR { return Mul(v, v); });
}

void square_double(double *inout, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    transform1<double>(inout, n, [](auto d, const auto v) HWY_ATTR { return Mul(v, v); });
}

void square_downcast(double *inout, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    transform1_downcast(inout, n, [](auto d, const auto v) HWY_ATTR { return Mul(v, v); });
}

// Return the max value of `in`
uint32_t max(const uint32_t *in, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<uint32_t> d;
    size_t i = 0;
    const size_t N = Lanes(d);
    auto m = Zero(d);
    for (; i + N <= n; i += N) {
        m = Max(m, LoadU(d, in + i));
    }
    uint32_t ret = GetLane(MaxOfLanes(d, m));
    for (; i < n; i++) {
        ret = std::max(ret, in[i]);
    }
    return ret;
}

// Write inout[i] = inout[i] + a[i]
void add_vec(uint32_t *inout, const uint32_t *a, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    transform2(inout, n, a, [](auto d, const auto v1, const auto v2) { return Add(v1, v2); });
}

// Write inout[i] = inout[i] + a
void add_const(uint32_t *inout, const int32_t a, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    transform1(inout, n, [a](auto d, const auto v) {
        RebindToSigned<decltype(d)> d_signed;
        auto v2 = BitCast(d, Set(d_signed, a));
        return Add(v, v2);
    });
}

// Write inout[i] = inout[i] - a[i]
void sub(uint32_t *inout, const uint32_t *a, size_t n) {
    using namespace hwy::HWY_NAMESPACE;
    transform2(inout, n, a, [](auto d, const auto v1, const auto v2) { return Sub(v1, v2); });
}

} // namespace BPCells::simd::HWY_NAMESPACE
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace BPCells::simd {

HWY_EXPORT(log1p_float);
HWY_EXPORT(log1p_double);
HWY_EXPORT(log1p_downcast);

void log1p(float *inout, size_t n) { return HWY_DYNAMIC_DISPATCH(log1p_float)(inout, n); }

void log1p(double *inout, size_t n) { return HWY_DYNAMIC_DISPATCH(log1p_double)(inout, n); }

void log1p_downcast(double *inout, size_t n) {
    return HWY_DYNAMIC_DISPATCH(log1p_downcast)(inout, n);
}

HWY_EXPORT(expm1_float);
HWY_EXPORT(expm1_double);
HWY_EXPORT(expm1_downcast);

void expm1(float *inout, size_t n) { return HWY_DYNAMIC_DISPATCH(expm1_float)(inout, n); }

void expm1(double *inout, size_t n) { return HWY_DYNAMIC_DISPATCH(expm1_double)(inout, n); }

void expm1_downcast(double *inout, size_t n) {
    return HWY_DYNAMIC_DISPATCH(expm1_downcast)(inout, n);
}

HWY_EXPORT(square_float);
HWY_EXPORT(square_double);
HWY_EXPORT(square_downcast);

void square(float *inout, size_t n) { return HWY_DYNAMIC_DISPATCH(square_float)(inout, n); }

void square(double *inout, size_t n) { return HWY_DYNAMIC_DISPATCH(square_double)(inout, n); }

void square_downcast(double *inout, size_t n) {
    return HWY_DYNAMIC_DISPATCH(square_downcast)(inout, n);
}

HWY_EXPORT(max);
HWY_EXPORT(add_vec);
HWY_EXPORT(add_const);
HWY_EXPORT(sub);

uint32_t max(const uint32_t *in, size_t n) { return HWY_DYNAMIC_DISPATCH(max)(in, n); }

// Write inout[i] = inout[i] + a[i]
void add(uint32_t *inout, const uint32_t *a, size_t n) {
    HWY_DYNAMIC_DISPATCH(add_vec)(inout, a, n);
}

// Write inout[i] = inout[i] + a
void add(uint32_t *inout, const int32_t a, size_t n) {
    HWY_DYNAMIC_DISPATCH(add_const)(inout, a, n);
}

// Write inout[i] = inout[i] - a[i]
void sub(uint32_t *inout, const uint32_t *a, size_t n) { HWY_DYNAMIC_DISPATCH(sub)(inout, a, n); }

} // namespace BPCells::simd
#endif