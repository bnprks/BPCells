// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

// Per-target include guard
#if defined(BPCELLS_SIMD_TRANSFORM_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef BPCELLS_SIMD_TRANSFORM_INL_H_
#undef BPCELLS_SIMD_TRANSFORM_INL_H_
#else
#define BPCELLS_SIMD_TRANSFORM_INL_H_
#endif

#include <hwy/highway.h>

// Helpers for applying transformation functions across spase/dense vectors

HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::HWY_NAMESPACE {

namespace detail {

// Load index offsets for a `GatherIndex` function call.
template <class D>
inline hwy::HWY_NAMESPACE::Vec<hwy::HWY_NAMESPACE::RebindToSigned<D>> LoadIndexU(D d, const uint32_t *HWY_RESTRICT p) {
    using namespace hwy::HWY_NAMESPACE;

    if constexpr (std::is_same_v<Rebind<int32_t, D>, RebindToSigned<D>>) {
        const RebindToSigned<D> d;
        return LoadU(d, (const int32_t*) p);
    } else {
        static_assert(std::is_same_v<Rebind<int64_t, D>, RebindToSigned<D>>);
        const RebindToSigned<D> d;
        const Half<Rebind<int32_t, D>> d_half;
        return PromoteTo(d, LoadU(d_half, (const int32_t*) p));
    }
}

} // namespace detail

// Set inout[idx] = op(d, inout[idx])
template <typename T, class Op>
inline void transform1(T *HWY_RESTRICT inout, size_t n, const Op &&op) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<T> d;
    size_t i = 0;
    const size_t N = Lanes(d);
    for (; i + N <= n; i += N) {
        auto data = LoadU(d, inout + i);
        StoreU(op(d, data), d, inout + i);
    }

    CappedTag<T, 1> d1;
    for (; i < n; i++) {
        auto data = LoadU(d1, inout + i);
        StoreU(op(d1, data), d1, inout + i);
    }
}

// Set inout[idx] = (double) op(d, (float) inout[idx])
// Loses some precision, but increases throughput
template <typename Op>
inline void transform1_downcast(double *HWY_RESTRICT inout, size_t n, const Op &&op) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<float> df;
    ScalableTag<double> dd;
    const size_t Nd = Lanes(dd);
    const size_t Nf = Lanes(df);

    size_t i = 0;
    Half<decltype(df)> d_half;
    for (; i + Nf <= n; i += Nf) {
        auto data1 = LoadU(dd, inout + i);
        auto data2 = LoadU(dd, inout + i + Nd);

        auto data = Combine(df, DemoteTo(d_half, data2), DemoteTo(d_half, data1));
        data = op(df, data);

        StoreU(PromoteTo(dd, LowerHalf(d_half, data)), dd, inout + i);
        StoreU(PromoteTo(dd, UpperHalf(d_half, data)), dd, inout + i + Nd);
    }

    CappedTag<double, 1> d1;
    for (; i < n; i++) {
        auto data = LoadU(d1, inout + i);
        StoreU(op(d1, data), d1, inout + i);
    }
}

// Set inout[idx] = op(d, inout[idx], in1[idx])
template <typename T, typename Op>
inline void transform2(T *HWY_RESTRICT inout, size_t n, const T *HWY_RESTRICT in1, const Op &&op) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<T> d;
    size_t i = 0;
    const size_t N = Lanes(d);
    for (; i + N <= n; i += N) {
        auto data = LoadU(d, inout + i);
        auto arg1 = LoadU(d, in1 + i);
        StoreU(op(d, data, arg1), d, inout + i);
    }

    CappedTag<T, 1> d1;
    for (; i < n; i++) {
        auto data = LoadU(d1, inout + i);
        auto arg1 = LoadU(d1, in1 + i);
        StoreU(op(d1, data, arg1), d1, inout + i);
    }
}

// Set inout[idx] = (double) op(d, (float) inout[idx], in1[idx])
template <typename Op>
inline void transform2_downcast(
    double *HWY_RESTRICT inout, size_t n, const float *HWY_RESTRICT in1, const Op &&op
) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<float> df;
    ScalableTag<double> dd;
    const size_t Nd = Lanes(dd);
    const size_t Nf = Lanes(df);

    size_t i = 0;
    Half<decltype(df)> d_half;
    for (; i + Nf <= n; i += Nf) {
        auto data1 = LoadU(dd, inout + i);
        auto data2 = LoadU(dd, inout + i + Nd);
        auto arg1 = LoadU(df, in1 + i);

        auto data = Combine(df, DemoteTo(d_half, data2), DemoteTo(d_half, data1));
        data = op(df, data, arg1);

        StoreU(PromoteTo(dd, LowerHalf(d_half, data)), dd, inout + i);
        StoreU(PromoteTo(dd, UpperHalf(d_half, data)), dd, inout + i + Nd);
    }

    CappedTag<double, 1> d1d;
    CappedTag<float,  1> d1f;
    for (; i < n; i++) {
        auto data = LoadU(d1d, inout + i);
        auto arg1 = PromoteTo(d1d, LoadU(d1f, in1 + i));
        StoreU(op(d1d, data, arg1), d1d, inout + i);
    }
}


// Set inout[i] = op(d, inout[i], in1[idx[i]])
template <typename T, typename Op>
inline void transform2_sparse(
    T *HWY_RESTRICT inout,
    const uint32_t *HWY_RESTRICT idx,
    size_t n,
    const T *HWY_RESTRICT in1,
    const Op &&op
) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<T> d;

    size_t i = 0;
    const size_t N = Lanes(d);
    for (; i + N <= n; i += N) {
        auto indices = detail::LoadIndexU(d, idx + i);
        auto data = LoadU(d, inout + i);
        auto arg1 = GatherIndex(d, in1, indices);
        StoreU(op(d, data, arg1), d, inout + i);
    }

    CappedTag<T, 1> d1;
    for (; i < n; i++) {
        auto indices = detail::LoadIndexU(d1, idx + i);
        auto data = LoadU(d1, inout + i);
        auto arg1 = GatherIndex(d1, in1, indices);
        StoreU(op(d1, data, arg1), d1, inout + i);
    }
}

// Set inout[i] = op(d, inout[i], in1[idx[i]])
template <typename Op>
inline void transform2_sparse_downcast(
    double *HWY_RESTRICT inout,
    const uint32_t *HWY_RESTRICT idx,
    size_t n,
    const float *HWY_RESTRICT in1,
    const Op &&op
) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<float> df;
    ScalableTag<double> dd;
    const size_t Nd = Lanes(dd);
    const size_t Nf = Lanes(df);

    size_t i = 0;
    Half<decltype(df)> d_half;
    for (; i + Nf <= n; i += Nf) {
        auto indices = detail::LoadIndexU(df, idx + i);

        auto data1 = LoadU(dd, inout + i);
        auto data2 = LoadU(dd, inout + i + Nd);

        auto arg1 = GatherIndex(df, in1, indices);

        auto data = Combine(df, DemoteTo(d_half, data2), DemoteTo(d_half, data1));
        data = op(df, data, arg1);

        StoreU(PromoteTo(dd, LowerHalf(d_half, data)), dd, inout + i);
        StoreU(PromoteTo(dd, UpperHalf(d_half, data)), dd, inout + i + Nd);
    }

    CappedTag<double, 1> d1d;
    CappedTag<float,  1> d1f;
    for (; i < n; i++) {
        auto indices = detail::LoadIndexU(d1f, idx + i);

        auto data = LoadU(d1d, inout + i);

        auto arg1 = PromoteTo(d1d, GatherIndex(d1f, in1, indices));

        StoreU(op(d1d, data, arg1), d1d, inout + i);
    }
}

// Set inout[idx] = op(d, inout[idx], in1[idx], in2[idx])
template <typename T, typename Op>
inline void transform3(
    T *HWY_RESTRICT inout,
    size_t n,
    const T *HWY_RESTRICT in1,
    const T *HWY_RESTRICT in2,
    const Op &&op
) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<T> d;
    size_t i = 0;
    const size_t N = Lanes(d);
    for (; i + N <= n; i += N) {
        auto data = LoadU(d, inout + i);
        auto arg1 = LoadU(d, in1 + i);
        auto arg2 = LoadU(d, in2 + i);
        StoreU(op(d, data, arg1, arg2), d, inout + i);
    }

    CappedTag<T, 1> d1;
    for (; i < n; i++) {
        auto data = LoadU(d1, inout + i);
        auto arg1 = LoadU(d1, in1 + i);
        auto arg2 = LoadU(d1, in2 + i);
        StoreU(op(d1, data, arg1, arg2), d1, inout + i);
    }
}

// Set inout[idx] = (double) op(d, (float) inout[idx], in1[idx], in2[idx])
template <typename Op>
inline void transform3_downcast(
    double *HWY_RESTRICT inout,
    size_t n,
    const float *HWY_RESTRICT in1,
    const float *HWY_RESTRICT in2,
    const Op &&op
) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<float> df;
    ScalableTag<double> dd;
    const size_t Nd = Lanes(dd);
    const size_t Nf = Lanes(df);

    size_t i = 0;
    Half<decltype(df)> d_half;
    for (; i + Nf <= n; i += Nf) {
        auto data1 = LoadU(dd, inout + i);
        auto data2 = LoadU(dd, inout + i + Nd);
        auto arg1 = LoadU(df, in1 + i);
        auto arg2 = LoadU(df, in2 + i);

        auto data = Combine(df, DemoteTo(d_half, data2), DemoteTo(d_half, data1));
        data = op(df, data, arg1, arg2);

        StoreU(PromoteTo(dd, LowerHalf(d_half, data)), dd, inout + i);
        StoreU(PromoteTo(dd, UpperHalf(d_half, data)), dd, inout + i + Nd);
    }

    CappedTag<double, 1> d1d;
    CappedTag<float,  1> d1f;
    for (; i < n; i++) {
        auto data = LoadU(d1d, inout + i);
        auto arg1 = PromoteTo(d1d, LoadU(d1f, in1 + i));
        auto arg2 = PromoteTo(d1d, LoadU(d1f, in2 + i));
        StoreU(op(d1d, data, arg1, arg2), d1d, inout + i);
    }
}

// Set inout[i] = op(d, inout[i], in1[idx[i]], in2[idx[i]])
template <typename T, typename Op>
inline void transform3_sparse(
    T *HWY_RESTRICT inout,
    const uint32_t *HWY_RESTRICT idx,
    size_t n,
    const T *HWY_RESTRICT in1,
    const T *HWY_RESTRICT in2,
    const Op &&op
) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<T> d;

    size_t i = 0;
    const size_t N = Lanes(d);
    for (; i + N <= n; i += N) {
        auto indices = detail::LoadIndexU(d, idx + i);
        auto data = LoadU(d, inout + i);
        auto arg1 = GatherIndex(d, in1, indices);
        auto arg2 = GatherIndex(d, in2, indices);
        StoreU(op(d, data, arg1, arg2), d, inout + i);
    }

    CappedTag<T, 1> d1;
    for (; i < n; i++) {
        auto indices = detail::LoadIndexU(d1, idx + i);
        auto data = LoadU(d1, inout + i);
        auto arg1 = GatherIndex(d1, in1, indices);
        auto arg2 = GatherIndex(d1, in2, indices);
        StoreU(op(d1, data, arg1, arg2), d1, inout + i);
    }
}

// Set inout[i] = op(d, inout[i], in1[idx[i]], in2[idx[i]])
template <typename Op>
inline void transform3_sparse_downcast(
    double *HWY_RESTRICT inout,
    const uint32_t *HWY_RESTRICT idx,
    size_t n,
    const float *HWY_RESTRICT in1,
    const float *HWY_RESTRICT in2,
    const Op &&op
) {
    using namespace hwy::HWY_NAMESPACE;
    ScalableTag<float> df;
    ScalableTag<double> dd;
    const size_t Nd = Lanes(dd);
    const size_t Nf = Lanes(df);

    size_t i = 0;
    Half<decltype(df)> d_half;
    for (; i + Nf <= n; i += Nf) {
        auto indices = detail::LoadIndexU(df, idx + i);

        auto data1 = LoadU(dd, inout + i);
        auto data2 = LoadU(dd, inout + i + Nd);

        auto arg1 = GatherIndex(df, in1, indices);
        auto arg2 = GatherIndex(df, in2, indices);

        auto data = Combine(df, DemoteTo(d_half, data2), DemoteTo(d_half, data1));
        data = op(df, data, arg1, arg2);

        StoreU(PromoteTo(dd, LowerHalf(d_half, data)), dd, inout + i);
        StoreU(PromoteTo(dd, UpperHalf(d_half, data)), dd, inout + i + Nd);
    }

    CappedTag<double, 1> d1d;
    CappedTag<float,  1> d1f;
    for (; i < n; i++) {
        auto indices = detail::LoadIndexU(d1f, idx + i);

        auto data = LoadU(d1d, inout + i);

        auto arg1 = PromoteTo(d1d, GatherIndex(d1f, in1, indices));
        auto arg2 = PromoteTo(d1d, GatherIndex(d1f, in2, indices));

        StoreU(op(d1d, data, arg1, arg2), d1d, inout + i);
    }
}


} // namespace BPCells::simd::HWY_NAMESPACE

HWY_AFTER_NAMESPACE();

#endif // BPCELLS_SIMD_TRANSFORM_INL_H_