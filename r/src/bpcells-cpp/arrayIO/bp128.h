// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once
#include "array_interfaces.h"
#include <cstring>

namespace BPCells {

// Base class for BP128 readers. Derived classes are
// just responsible for defining load128 and seek appropriately
class BP128UIntReader : public UIntBulkReader {
  protected:
    UIntReader data;         // Bitpacked data
    UIntReader idx;          // Index of data for each 128-element chunk
    ULongReader idx_offsets; // Ranges in the idx array that must have UINT32_MAX added to them

    uint64_t pos = 0;  // Logical index in the reader
    uint64_t count;    // Total emlements in the reader
    uint64_t prev_idx; // Physical index in the data
    uint64_t idx_offset =
        0; // Current offset to apply on values read from idx. (Multiple of UINT32_MAX)
    uint64_t prev_offset_boundary,
        next_offset_boundary; // Logical index / 128 of current and next offset chunk

    uint32_t buf[128];

    // Load 128 values from in to out, assuming given bits per element
    virtual void load128(uint32_t *in, uint32_t *out, uint32_t bits);

    // Seek the loader objects to position listed in pos,
    // and update any position-associated variables like prev_idx
    virtual void seekLoaders();

  public:
    BP128UIntReader(UIntReader &&data, UIntReader &&idx, ULongReader &&idx_offsets, uint64_t count);

    // Return total number of integers in the reader
    uint64_t size() const final override;

    // Change the next load to start at index pos
    void seek(uint64_t pos) final override;

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 if count is >0
    // Note: It is the caller's responsibility to ensure there is no data overflow, i.e.
    // that a load does not try to read past size() total elements
    uint64_t load(uint32_t *out, uint64_t count) final override;

    // For testing purposes only -- set the offset multiple so it can be smaller
    // than UINT32_MAX
    static void setOffsetIncrement(uint64_t val);

  private:
    void load128(uint32_t *out);

    static inline uint64_t OFFSET_INCREMENT = UINT32_MAX + 1ULL;
};

// Base class for BP128 writers.
// WARNING: Constructor should propbably write first 0 index pointer, don't forget!
class BP128UIntWriter : public UIntBulkWriter {
  protected:
    UIntWriter data, idx;
    ULongWriter idx_offsets;

    uint64_t pos = 0;
    uint64_t cur_idx = 0;
    uint64_t buf_pos = 0;

    uint32_t buf[128];
    // Pack 128 values from in pointer, using specified number of bits per element
    virtual void pack128(uint32_t *in, uint32_t *out, uint32_t bits);
    // Return the number of bits needed to pack new input data
    virtual uint32_t bits(const uint32_t *in) const;

  public:
    BP128UIntWriter(UIntWriter &&data, UIntWriter &&idx, ULongWriter &&idx_offsets);

    // Write up to `count` integers from `in`, returning the actual number written.
    // Will always write >0, otherwise throwing an exception for an error
    // Note: The writer is allowed to modify the input data, so it might be
    // modified after calling write()
    uint64_t write(uint32_t *in, uint64_t count) final override;

    void finalize() override;

    // For testing purposes only -- set the offset multiple so it can be smaller
    // than UINT32_MAX
    static void setOffsetIncrement(uint64_t val);

  private:
    void pack128(uint32_t *in);

    static inline uint64_t OFFSET_INCREMENT = UINT32_MAX + 1ULL;
};

// ################### BP128 (Vanilla) #########################################

// class BP128UIntReader final : public BP128UIntReaderBase {
//   protected:
//     UIntReader data, idx;
//     uint64_t prev_idx;

//     void load128(uint32_t *out) override;
//     void _seek() override;

//   public:
//     BP128UIntReader(UIntReader &&data, UIntReader &&idx, uint64_t count);
// };

// class BP128UIntWriter final : public BP128UIntWriterBase {
//   protected:
//     UIntWriter data, idx;
//     uint64_t cur_idx = 0;

//     void pack128(uint32_t *in) override;
//     void finalizeWriters() override;

//   public:
//     BP128UIntWriter(UIntWriter &&data, UIntWriter &&idx);
// };

// ######################## BP128 (D1) #########################################

class BP128_D1_UIntReader final : public BP128UIntReader {
  protected:
    UIntReader starts;

    void load128(uint32_t *in, uint32_t *out, uint32_t bits) override;
    void seekLoaders() override;

  public:
    BP128_D1_UIntReader(
        UIntReader &&data,
        UIntReader &&idx,
        ULongReader &&idx_offsets,
        UIntReader &&starts,
        uint64_t count
    );
};

class BP128_D1_UIntWriter final : public BP128UIntWriter {
  protected:
    UIntWriter starts;

    // Pack 128 values from in pointer, using specified number of bits per element
    void pack128(uint32_t *in, uint32_t *out, uint32_t bits) override;
    // Return the number of bits needed to pack new input data
    uint32_t bits(const uint32_t *in) const override;

  public:
    BP128_D1_UIntWriter(
        UIntWriter &&data, UIntWriter &&idx, ULongWriter &&idx_offsets, UIntWriter &&starts
    );

    void finalize() override;
};

// ######################## BP128 (D1Z) #########################################

class BP128_D1Z_UIntReader final : public BP128UIntReader {
  protected:
    UIntReader starts;

    void load128(uint32_t *in, uint32_t *out, uint32_t bits) override;
    void seekLoaders() override;

  public:
    BP128_D1Z_UIntReader(
        UIntReader &&data,
        UIntReader &&idx,
        ULongReader &&idx_offsets,
        UIntReader &&starts,
        uint64_t count
    );
};

class BP128_D1Z_UIntWriter final : public BP128UIntWriter {
  protected:
    UIntWriter starts;

    // Pack 128 values from in pointer, using specified number of bits per element
    void pack128(uint32_t *in, uint32_t *out, uint32_t bits) override;
    // Return the number of bits needed to pack new input data
    uint32_t bits(const uint32_t *in) const override;

  public:
    BP128_D1Z_UIntWriter(
        UIntWriter &&data, UIntWriter &&idx, ULongWriter &&idx_offsets, UIntWriter &&starts
    );

    void finalize() override;
};

// ######################## BP128 (FOR) #########################################
//  This is just FOR encoding using a constant 1 as the frame of reference.
//  (i.e. to encode values of an integer sparse matrix)
class BP128_FOR_UIntReader final : public BP128UIntReader {
  protected:

    void load128(uint32_t *in, uint32_t *out, uint32_t bits) override;
  public:
    BP128_FOR_UIntReader(UIntReader &&data, UIntReader &&idx, ULongReader &&idx_offsets, uint64_t count);
};

class BP128_FOR_UIntWriter final : public BP128UIntWriter {
  protected:
    void pack128(uint32_t *in, uint32_t *out, uint32_t bits) override;
    // Return the number of bits needed to pack new input data
    uint32_t bits(const uint32_t *in) const override;
  public:
    BP128_FOR_UIntWriter(UIntWriter &&data, UIntWriter &&idx, ULongWriter &&idx_offsets);
};

} // end namespace BPCells