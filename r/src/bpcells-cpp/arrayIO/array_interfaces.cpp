// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "array_interfaces.h"

namespace BPCells {

VecStringReader::VecStringReader(std::vector<std::string> data) : data(data) {}
const char *VecStringReader::get(uint64_t idx) {
    if (idx < data.size()) return data[idx].c_str();
    return NULL;
}
uint64_t VecStringReader::size() { return data.size(); }

} // end namespace BPCells
