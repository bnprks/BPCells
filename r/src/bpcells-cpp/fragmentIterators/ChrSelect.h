// Copyright 2021 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include "FragmentIterator.h"
#include <algorithm>
#include <unordered_map>

namespace BPCells {

// Transform a fragments loader by renaming or filtering chromosome IDs
class ChrIndexSelect : public FragmentLoaderWrapper {
  private:
    const std::vector<uint32_t> chr_assignments;
    uint32_t current_chr = UINT32_MAX;

  public:
    // chr_assignments -- vector with length <= the number of chromosomes in the input
    //     FragmentLoader. The output chromosome `i` will come from input chromosome
    //     `chr_assignments[i]`. The entries of chr_assignments must be unique
    ChrIndexSelect(std::unique_ptr<FragmentLoader> &&loader, const std::vector<uint32_t> chr_assignments);

    ~ChrIndexSelect() = default;

    void restart() override;

    int chrCount() const override;

    const char *chrNames(uint32_t chr_id) override;

    bool nextChr() override;

    uint32_t currentChr() const override;

    // Move loader to just before fragments which end after "base".
    // It's possible that fragments returned after seek will end before "base",
    // but it's guaranteed that it won't skip over any fragments ending before "base"
    void seek(uint32_t chr_id, uint32_t base) override;
};

// Transform a fragments loader by renaming or filtering chromosome IDs
class ChrNameSelect : public FragmentLoaderWrapper {
  private:
    const std::vector<std::string> chr_names;
    std::unordered_map<std::string, uint32_t> output_index;
    // if iter is seekable, map the output IDs to input IDs. input_index[i] gives
    // the input index for output chromosome i. This is just needed for seeking
    std::vector<uint32_t> input_index;
    uint32_t current_chr = UINT32_MAX;

  public:
    // chr_names -- vector with length <= the number of chromosomes in the input
    //     FragmentsIterator. The output chromosome `i` will come from input chromosome with name
    //     `chr_names[i]`. The entries of chr_names must be unique
    ChrNameSelect(std::unique_ptr<FragmentLoader> &&loader, const std::vector<std::string> chr_names);

    ~ChrNameSelect() = default;

    void restart() override;

    int chrCount() const override;

    const char *chrNames(uint32_t chr_id) override;

    bool nextChr() override;

    uint32_t currentChr() const override;

    void seek(uint32_t chr_id, uint32_t base) override;
};

} // end namespace BPCells