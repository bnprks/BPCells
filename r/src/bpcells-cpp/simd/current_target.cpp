// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "simd/current_target.cpp"
#include <hwy/foreach_target.h>

#include <hwy/highway.h>


HWY_BEFORE_NAMESPACE();

namespace BPCells::simd::HWY_NAMESPACE {

const char* current_target() {
    return hwy::TargetName(HWY_TARGET);
}

} // namespace BPCells::simd::HWY_NAMESPACE
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

#include <vector>
#include <string>
#include <stdexcept>
namespace BPCells::simd {

HWY_EXPORT(current_target);

// Return a string with the name of the current SIMD instruction set being used by highway
const char* current_target() {
    return HWY_DYNAMIC_DISPATCH(current_target)();
}

// Return a list of all highway SIMD instruction sets supported on this CPU
std::vector<std::string> supported_targets() {
    std::vector<std::string> res;
    for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
        res.push_back(hwy::TargetName(target));
    }
    return res;
}

// Set the active SIMD instruction set highway should use. 
// This is useful for testing/debugging to force use of an older SIMD instruction set
void set_target(std::string target) {
    for (int64_t t : hwy::SupportedAndGeneratedTargets()) {
        if (target == hwy::TargetName(t)) {
            hwy::SetSupportedTargetsForTest(t);
            return;
        }
    }
    throw std::invalid_argument("set_target(): target '" + target + "' not available");
}


} // namespace BPCells::simd
#endif