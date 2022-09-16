#pragma once
#include "../bitpacking/bp128.h"
#include "array_interfaces.h"
#include <cstring>

namespace BPCells {

// Base class for BP128 readers. Derived classes are
// just responsible for defining load128 and seek appropriately
class BP128UIntReaderBase : public UIntBulkReader {
  protected:
    uint32_t pos = 0;
    uint32_t count;
    std::unique_ptr<uint32_t[]> buf = std::make_unique<uint32_t[]>(128);

    // Load 128 values into out pointer
    virtual void load128(uint32_t *out) = 0;
    // Seek to the position of the pos variable
    virtual void _seek() = 0;

  public:
    BP128UIntReaderBase(uint32_t count);
    virtual ~BP128UIntReaderBase() = default;
    // Return total number of integers in the reader
    uint32_t size() const final override;

    // Change the next load to start at index pos
    void seek(uint32_t pos) final override;

    // Copy up to `count` integers into `out`, returning the actual number copied.
    // Will always load >0 if count is >0
    // Note: It is the caller's responsibility to ensure there is no data overflow, i.e.
    // that a load does not try to read past size() total elements
    uint32_t load(uint32_t *out, uint32_t count) final override;
};

// Base class for BP128 writers.
// WARNING: Constructor should propbably write first 0 index pointer, don't forget!
class BP128UIntWriterBase : public UIntBulkWriter {
  protected:
    std::unique_ptr<uint32_t[]> buf = std::make_unique<uint32_t[]>(128);
    uint32_t buf_pos = 0;

    // Pack 128 values from in pointer
    virtual void pack128(uint32_t *in) = 0;
    // Finalize all the underlying writer objects, after all data has been written
    virtual void finalizeWriters() = 0;

  public:
    virtual ~BP128UIntWriterBase() = default;

    // Write up to `count` integers from `in`, returning the actual number written.
    // Will always write >0, otherwise throwing an exception for an error
    // Note: The writer is allowed to modify the input data, so it might be
    // modified after calling write()
    uint32_t write(uint32_t *in, uint32_t count) override;

    void finalize() final override;
};

//################### BP128 (Vanilla) #########################################

class BP128UIntReader final : public BP128UIntReaderBase {
  protected:
    UIntReader data, idx;
    uint32_t prev_idx;

    void load128(uint32_t *out) override;
    void _seek() override;

  public:
    BP128UIntReader(UIntReader &&data, UIntReader &&idx, uint32_t count);
};

class BP128UIntWriter final : public BP128UIntWriterBase {
  protected:
    UIntWriter data, idx;
    uint32_t cur_idx = 0;

    void pack128(uint32_t *in) override;
    void finalizeWriters() override;

  public:
    BP128UIntWriter(UIntWriter &&data, UIntWriter &&idx);
};

//######################## BP128 (D1) #########################################

class BP128_D1_UIntReader final : public BP128UIntReaderBase {
  protected:
    UIntReader data, idx, starts;
    uint32_t prev_idx;

    void load128(uint32_t *out) override;
    void _seek() override;

  public:
    BP128_D1_UIntReader(UIntReader &&data, UIntReader &&idx, UIntReader &&starts, uint32_t count);
};

class BP128_D1_UIntWriter final : public BP128UIntWriterBase {
  protected:
    UIntWriter data, idx, starts;
    uint32_t cur_idx = 0;

    void pack128(uint32_t *in) override;
    void finalizeWriters() override;

  public:
    BP128_D1_UIntWriter(UIntWriter &&data, UIntWriter &&idx, UIntWriter &&starts);
};

//######################## BP128 (D1Z) #########################################

class BP128_D1Z_UIntReader final : public BP128UIntReaderBase {
  protected:
    UIntReader data, idx, starts;
    uint32_t prev_idx;

    void load128(uint32_t *out) override;
    void _seek() override;

  public:
    BP128_D1Z_UIntReader(UIntReader &&idx, UIntReader &&data, UIntReader &&starts, uint32_t count);
};

class BP128_D1Z_UIntWriter final : public BP128UIntWriterBase {
  protected:
    UIntWriter data, idx, starts;
    uint32_t cur_idx = 0;

    void pack128(uint32_t *in) override;
    void finalizeWriters() override;

  public:
    BP128_D1Z_UIntWriter(UIntWriter &&data, UIntWriter &&idx, UIntWriter &&starts);
};

//######################## BP128 (FOR) #########################################
// This is just FOR encoding using a constant 1 as the frame of reference.
// (i.e. to encode values of an integer sparse matrix)
class BP128_FOR_UIntReader final : public BP128UIntReaderBase {
  protected:
    UIntReader data, idx;
    uint32_t prev_idx;

    void load128(uint32_t *out) override;
    void _seek() override;

  public:
    BP128_FOR_UIntReader(UIntReader &&data, UIntReader &&idx, uint32_t count);
};

class BP128_FOR_UIntWriter final : public BP128UIntWriterBase {
  protected:
    UIntWriter data, idx;
    uint32_t cur_idx = 0;

    void pack128(uint32_t *in) override;
    void finalizeWriters() override;

  public:
    BP128_FOR_UIntWriter(UIntWriter &&data, UIntWriter &&idx);
};

} // end namespace BPCells