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
// keep/drop decision a deterministic function of (seed, col, row) via splitmix64 and invariant to storage order.
template <typename T, bool Transpose = false> class BernoulliSample : public MatrixLoaderWrapper<T> {
  private:
    uint64_t seed;
    uint32_t threshold;
    size_t loaded = 0;

    // splitmix64 via Daniel Lemire's testingRNG:
    // https://github.com/lemire/testingRNG/blob/master/source/splitmix64.h
    static constexpr uint64_t GOLDEN_GAMMA = UINT64_C(0x9E3779B97F4A7C15);

    static inline uint64_t splitmix64_r(uint64_t *seed) {
        uint64_t z = (*seed += GOLDEN_GAMMA);
        z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
        return z ^ (z >> 31);
    }

    // returns the value of splitmix64 "offset" steps from seed
    static inline uint64_t splitmix64_stateless(uint64_t seed, uint64_t offset) {
        seed += offset * GOLDEN_GAMMA;
        return splitmix64_r(&seed);
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
        const uint64_t nrow = this->loader->rows();
        const uint64_t ncol = this->loader->cols();

        while (loaded == 0) {
            if (!this->loader->load()) return false;
            uint32_t *row_data = this->loader->rowData();
            T *val_data = this->loader->valData();
            size_t cap = this->loader->capacity();
            uint32_t col = this->loader->currentCol();
            
            for (size_t i = 0; i < cap; i++) {
                uint32_t row = row_data[i];
                uint64_t offset = !Transpose ? (col * nrow + row) : (row * ncol + col);
                uint32_t h = static_cast<uint32_t>(splitmix64_stateless(seed, offset));
                row_data[loaded] = row;
                val_data[loaded] = val_data[i];
                loaded += h < threshold;
            }
        }
        return true;
    }

    uint32_t capacity() const override { return loaded; }
};

} // end namespace BPCells
