// Copyright 2024 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <atomic>

#include "../arrayIO/array_interfaces.h"
#include "MatrixIterator.h"

namespace BPCells {

// Filter out zero values from a MatrixLoader
// This is useful when reading dense matrices that have many zero values,
// or when performing operations that will cause new zeros to be created (e.g. multiplying a row by zero)
template<typename T>
class FilterZeros : public MatrixLoaderWrapper<T> {
  private:
    size_t capacity_ = 0;
  public:
    FilterZeros(std::unique_ptr<MatrixLoader<T>> &&loader) : MatrixLoaderWrapper<T>(std::move(loader)) {}

    // Return false if there are no more entries to load
    bool load() override {
        capacity_ = 0;
        
        while (capacity_ == 0) {
            if (!this->loader->load()) return false;
            
            uint32_t *row_data = this->loader->rowData();
            T* val_data = this->loader->valData();
            size_t cap = this->loader->capacity();
            
            for (size_t i = 0; i < cap; i++) {
                row_data[capacity_] = row_data[i];
                val_data[capacity_] = val_data[i];
                capacity_ += val_data[capacity_] != 0;
            }
        }
        return true;
    }

    // Number of loaded entries available
    uint32_t capacity() const override {return capacity_;}
};

} // end namespace BPCells