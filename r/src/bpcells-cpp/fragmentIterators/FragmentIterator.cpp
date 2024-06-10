// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include "FragmentIterator.h"

namespace BPCells {

FragmentIterator::FragmentIterator(std::unique_ptr<FragmentLoader> &&loader)
    : FragmentLoaderWrapper(std::move(loader)) {}

bool FragmentIterator::load() {
    if (!loader->load()) {
        current_capacity = 0;
        return false;
    }
    idx = 0;
    current_capacity = loader->capacity();
    current_cell = loader->cellData();
    current_start = loader->startData();
    current_end = loader->endData();
    return true;
}

uint32_t FragmentIterator::capacity() const { return current_capacity; }

FragmentLoaderWrapper::FragmentLoaderWrapper(std::unique_ptr<FragmentLoader> &&loader)
    : loader(std::move(loader)) {}

bool FragmentLoaderWrapper::isSeekable() const { return loader->isSeekable(); }
void FragmentLoaderWrapper::seek(uint32_t chr_id, uint32_t base) { loader->seek(chr_id, base); }

void FragmentLoaderWrapper::restart() { loader->restart(); }

int FragmentLoaderWrapper::chrCount() const { return loader->chrCount(); }
int FragmentLoaderWrapper::cellCount() const { return loader->cellCount(); }

const char *FragmentLoaderWrapper::chrNames(uint32_t chr_id) { return loader->chrNames(chr_id); }
const char *FragmentLoaderWrapper::cellNames(uint32_t cell_id) {
    return loader->cellNames(cell_id);
}

bool FragmentLoaderWrapper::nextChr() { return loader->nextChr(); }
uint32_t FragmentLoaderWrapper::currentChr() const { return loader->currentChr(); }

bool FragmentLoaderWrapper::load() { return loader->load(); }

uint32_t FragmentLoaderWrapper::capacity() const { return loader->capacity(); }

uint32_t *FragmentLoaderWrapper::cellData() { return loader->cellData(); }
uint32_t *FragmentLoaderWrapper::startData() { return loader->startData(); }
uint32_t *FragmentLoaderWrapper::endData() { return loader->endData(); }

} // end namespace BPCells
