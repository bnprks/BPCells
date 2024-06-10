// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include "FragmentIterator.h"

namespace BPCells {

// Transform a fragments iterator by subsetting to fragments of a given length range
class LengthSelect : public FragmentLoaderWrapper {
  private:
    uint32_t loaded = 0;
    uint32_t min_len, max_len;

  public:
    // min_len and max_len provide inclusive limits on the size of fragments
    LengthSelect(std::unique_ptr<FragmentLoader> &&loader, uint32_t min_len, uint32_t max_len = UINT32_MAX);

    ~LengthSelect() = default;

    bool load() override;
    uint32_t capacity() const override;
};

} // end namespace BPCells