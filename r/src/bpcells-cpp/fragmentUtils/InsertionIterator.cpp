// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "InsertionIterator.h"

namespace BPCells {

// I made literally everything else inline, but here we are...
InsertionIterator::InsertionIterator(FragmentLoader &loader) : frags(loader) {}

} // end namespace BPCells
