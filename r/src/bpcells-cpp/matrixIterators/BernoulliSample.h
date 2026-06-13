// Copyright 2026 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <cstdint>
#include <memory>
#include <stdexcept>

#include "MatrixIterator.h"

namespace BPCells {

// Keep each non-zero entry independently with probability `prob`, with the
// keep/drop decision a deterministic function of (seed, col, row) via splitmix64.
template <typename T> class BernoulliSample : public MatrixLoaderWrapper<T> {
  private:
    uint64_t seed;
    uint32_t threshold;
    size_t loaded = 0;

    // Hash algorithm
    static inline uint64_t splitmix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

  public:
    BernoulliSample(std::unique_ptr<MatrixLoader<T>> &&loader, double prob, uint64_t seed)
        : MatrixLoaderWrapper<T>(std::move(loader)), seed(seed) {
        if (!(prob > 0.0 && prob <= 1.0))
            throw std::runtime_error("BernoulliSample: prob must be in (0, 1]");
        double scaled = prob * 4294967296.0; // 2^32
        threshold = (prob >= 1.0 || scaled >= 4294967295.5)
            ? UINT32_MAX
            : static_cast<uint32_t>(scaled);
    }

    bool load() override {
        loaded = 0;

        while (loaded == 0) {
            if (!this->loader->load()) return false;
            uint32_t *row_data = this->loader->rowData();
            T *val_data = this->loader->valData();
            size_t cap = this->loader->capacity();

            // Calculate per-column hash
            uint64_t col_mix = splitmix64(seed ^ (uint64_t(this->loader->currentCol()) << 32));

            for (size_t i = 0; i < cap; i++) {
                // Calculate final hash with row number
                uint32_t h = static_cast<uint32_t>(splitmix64(col_mix ^ row_data[i]));
                row_data[loaded] = row_data[i];
                val_data[loaded] = val_data[i];
                loaded += h < threshold;
            }
        }
        return true;
    }

    uint32_t capacity() const override { return loaded; }
};

} // end namespace BPCells
