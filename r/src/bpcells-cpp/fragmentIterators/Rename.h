// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include "../arrayIO/array_interfaces.h"
#include "FragmentIterator.h"
#include <algorithm>

namespace BPCells {

// Transform a fragments loader by renaming chromosomes
// Errors may throw at load time if the number of chromosome names provided is not
// equal to the number of chromosome names in the wrapped loader
class RenameChrs : public FragmentLoaderWrapper {
  private:
    const std::unique_ptr<StringReader> chr_names;

  public:
    RenameChrs(std::unique_ptr<FragmentLoader> &&loader, std::unique_ptr<StringReader> &&chr_names);

    const char *chrNames(uint32_t chr_id) override;
};

// Transform a fragments loader by renaming cells
// Errors may throw at load time if the number of cell names provided is not
// equal to the number of cell names in the wrapped loader
class RenameCells : public FragmentLoaderWrapper {
  private:
    const std::unique_ptr<StringReader> cell_names;

  public:
    RenameCells(std::unique_ptr<FragmentLoader> &&loader, std::unique_ptr<StringReader> &&cell_names);

    const char *cellNames(uint32_t cell_id) override;
};

// Tranform a fragments loader by adding a prefix to all the cell names
// This is a limited alternative to RenameCells for situations where
// the cell names are not known ahead-of-time
class PrefixCells : public FragmentLoaderWrapper {
  private:
    std::string prefix;
    std::vector<char> name_buffer;

  public:
    PrefixCells(std::unique_ptr<FragmentLoader> &&loader, std::string prefix);

    const char *cellNames(uint32_t cell_id) override;
};

} // end namespace BPCells