// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

// Define a second copy of current_target that uses the same compile-time 
// optimization detection. This will allow accurate reporting of the simd set used for
// bp128

#if !defined(NDEBUG) && !defined(__OPTIMIZE__)
#define HWY_COMPILE_ONLY_STATIC 1
#endif

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "simd/bp128/current_target.cpp"
#include <hwy/foreach_target.h>

#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::bp128::HWY_NAMESPACE {

const char* current_target() {
    return hwy::TargetName(HWY_TARGET);
}

} // namespace BPCells::simd::HWY_NAMESPACE
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace BPCells::simd::bp128 {

HWY_EXPORT(current_target);

const char* current_target() {
    return HWY_DYNAMIC_DISPATCH(current_target)();
}

} // namespace BPCells::simd
#endif