// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <vector>
#include <string>

namespace BPCells::simd {

// Return a string with the name of the current SIMD instruction set being used by highway
const char *current_target();

// Return a list of all highway SIMD instruction sets supported on this CPU
std::vector<std::string> supported_targets();

// Set the active SIMD instruction set highway should use. 
// This is useful for testing/debugging to force use of an older SIMD instruction set
void set_target(std::string target);

}