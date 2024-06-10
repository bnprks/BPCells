// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "LengthSelect.h"

namespace BPCells {

LengthSelect::LengthSelect(std::unique_ptr<FragmentLoader> &&loader, uint32_t min_len, uint32_t max_len)
    : FragmentLoaderWrapper(std::move(loader))
    , min_len(min_len)
    , max_len(max_len) {}

bool LengthSelect::load() {
    loaded = 0;
    // load and filter until we load without filtering out everything
    while (loaded == 0) {
        if (!loader->load()) return false;

        uint32_t *cell = loader->cellData();
        uint32_t *start = loader->startData();
        uint32_t *end = loader->endData();
        uint32_t capacity = loader->capacity();
        for (uint32_t i = 0; i < capacity; i++) {
            cell[loaded] = cell[i];
            start[loaded] = start[i];
            end[loaded] = end[i];
            loaded += (end[i] - start[i]) >= min_len && (end[i] - start[i]) <= max_len;
        }
    }
    return true;
}

uint32_t LengthSelect::capacity() const { return loaded; }

} // end namespace BPCells
