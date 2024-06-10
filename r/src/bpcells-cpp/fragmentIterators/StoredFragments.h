// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <atomic>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>

#include "../arrayIO/array_interfaces.h"
#include "../arrayIO/bp128.h"
#include "FragmentIterator.h"

namespace BPCells {

class StoredFragmentsBase : public FragmentLoader {
  protected:
    // cell, start, end = concatenated fragment data across chromosomes
    // end_max[i] = max(ends[chr_start:i*128])
    // chr_ptr[2*i] = start index of chr i in data arrays
    // chr_ptr[2*i + 1] = end index of chr i in data arrays
    UIntReader cell, start, end, end_max;
    ULongReader chr_ptr;
    std::unique_ptr<StringReader> chr_names, cell_names;

    // end_max_buf holds the end_max values for the current chromosome, but to
    // deal with possible overhang from high end values in the previous chromosomes we
    // ensure end_max_buf[0] <= end_max_buf[1]
    std::vector<uint32_t> end_max_buf;
    uint32_t current_chr = UINT32_MAX;
    uint64_t current_idx = UINT64_MAX;
    uint32_t current_capacity = 0;
    uint64_t chr_start_ptr, chr_end_ptr;

    // Read end_max_buf from end_max iterator, making it equal to the values between
    // start_idx and end_idx
    void readEndMaxBuf(uint64_t start_idx, uint64_t end_idx);

  public:
    StoredFragmentsBase(
        UIntReader &&cell,
        UIntReader &&start,
        UIntReader &&end,
        UIntReader &&end_max,
        ULongReader &&chr_ptr,
        std::unique_ptr<StringReader> &&chr_names,
        std::unique_ptr<StringReader> &&cell_names
    );

    StoredFragmentsBase(StoredFragmentsBase&&) = default;

    bool isSeekable() const override;
    void seek(uint32_t chr_id, uint32_t base) override;

    void restart() override;

    int chrCount() const override;
    int cellCount() const override;
    const char *chrNames(uint32_t chr_id) override;
    const char *cellNames(uint32_t cell_id) override;

    bool nextChr() override;
    uint32_t currentChr() const override;

    uint32_t capacity() const override;

    uint32_t *cellData() override;
    uint32_t *startData() override;
    uint32_t *endData() override;
};

class StoredFragments : public StoredFragmentsBase {
  public:
    using StoredFragmentsBase::StoredFragmentsBase;
    static StoredFragments openUnpacked(
        ReaderBuilder &rb,
        std::unique_ptr<StringReader> &&chr_names = nullptr,
        std::unique_ptr<StringReader> &&cell_names = nullptr
    );

    bool load() override;
};

class StoredFragmentsPacked : public StoredFragmentsBase {
  public:
    using StoredFragmentsBase::StoredFragmentsBase;
    static StoredFragmentsPacked openPacked(
        ReaderBuilder &rb,
        uint32_t load_size = 1024,
        std::unique_ptr<StringReader> &&chr_names = nullptr,
        std::unique_ptr<StringReader> &&cell_names = nullptr
    );
    // Just override the methods that load data so we can insert a step to add start+end
    bool load() override;
};

class StoredFragmentsWriter : public FragmentWriter {
  private:
    UIntWriter cell, start, end, end_max;
    ULongWriter chr_ptr;
    std::unique_ptr<StringWriter> chr_names, cell_names;

    bool subtract_start_from_end; // Set to true if writing packed
    void subStartEnd();

  public:
    static StoredFragmentsWriter createUnpacked(WriterBuilder &wb);
    static StoredFragmentsWriter createPacked(WriterBuilder &wb, uint32_t buffer_size = 1024);
    StoredFragmentsWriter(
        UIntWriter &&cell,
        UIntWriter &&start,
        UIntWriter &&end,
        UIntWriter &&end_max,
        ULongWriter &&chr_ptr,
        std::unique_ptr<StringWriter> &&chr_names,
        std::unique_ptr<StringWriter> &&cell_names,
        bool subtract_start_from_end
    );

    void write(FragmentLoader &fragments, std::atomic<bool> *user_interrupt = NULL) override;
};

} // end namespace BPCells