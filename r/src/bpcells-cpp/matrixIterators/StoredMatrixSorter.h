// Copyright 2022 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <atomic>
#include <string>

#include "../arrayIO/array_interfaces.h"
#include "../arrayIO/binaryfile.h"
#include "../arrayIO/vector.h"
#include "../../vendor/dary_heap/dary_heap.hpp"
#include "../utils/radix_sort.h"
#include "MatrixIterator.h"
#include "StoredMatrix.h"

namespace BPCells {

// This is a base class used for transposing storage order and for importing matrix market files.
// This handles sorting matrix entries in order of (col, row) using a disk-backed merge sort.
// The subclass is responsible for implementing a function to read the initial matrix entries.
template <typename T> class StoredMatrixSorter {
    const uint64_t load_elements, sort_buffer_elements;
    uint64_t round = 0, total_elements = 0;

    // Load input entries into the provided vectors. The vectors will all be the same length.
    // Don't resize the vectors, but do try to fill the buffer as much as possible.
    // Return 0 at end of input data.
    virtual size_t load_entries(
        std::vector<uint32_t> &row,
        std::vector<uint32_t> &col,
        std::vector<T> &val,
        std::atomic<bool> *user_interrupt
    ) = 0;
    // Return if the data compes pre-sorted by row
    virtual bool row_sorted() const = 0;

  private:
    UIntReader row_reader, col_reader;
    NumReader<T> val_reader;
    UIntWriter row_writer, col_writer;
    NumWriter<T> val_writer;

    WriterBuilder &out_writer_builder;
    FileWriterBuilder file_writer_builder;
    FileReaderBuilder file_reader_builder;

    // Track which round reader/writer were opened so we can close
    // open file handles prior to delete (Windows doesn't like closing files with open handles)
    uint64_t reader_round = SIZE_MAX;
    uint64_t writer_round = SIZE_MAX;

    bool row_major;

    void openWriters(uint64_t round, bool last_round = false);
    void openReaders(uint64_t round, bool last_round = false);
    void deleteWriters(uint64_t round, bool last_round = false);

    template <typename S> NumWriter<S> openWriter(WriterBuilder &wb, std::string name);

    template <typename S> NumWriter<S> openReader(WriterBuilder &wb, std::string name);

    static uint64_t bytesPerElement() { return 2 * sizeof(uint32_t) + sizeof(T); }

    template <typename V>
    static void flushToNumWriter(const V *vals, uint64_t count, NumWriter<V> &writer);

    // Helper class to handle reading from a slice of a NumReader
    // Allow sharing NumReader between multiple SliceReaders, meaning that a
    // seek must be performed before each load
    template <typename V> class SliceReader {
        NumReader<V> &reader;
        std::vector<V> data;
        uint64_t data_idx = 0, num_loaded = 0, reader_idx;
        const uint64_t reader_end_idx;

      public:
        SliceReader(
            NumReader<V> &reader, uint64_t start_offset, uint64_t count, uint64_t chunk_size
        );
        bool advance();
        V value() const;
    };

  public:
    // Here, load_bytes is the smallest read size to use in bytes, and sort_buffer_bytes
    // is the total amount of memory to use for sort buffers
    StoredMatrixSorter(
        WriterBuilder &output,
        const char *tmpdir,
        uint64_t load_bytes,
        uint64_t sort_buffer_bytes,
        bool row_major =
            false // If true, swap metadata in writeValues() to mark outputs as row-major
    )
        : load_elements(load_bytes / sizeof(uint32_t))
        , sort_buffer_elements(sort_buffer_bytes / (bytesPerElement()))
        , out_writer_builder(output)
        , file_writer_builder(tmpdir, 8192)
        , file_reader_builder(tmpdir, 8192)
        , row_major(row_major) {

        if (sort_buffer_elements == 0 || load_elements == 0 ||
            load_elements * 2 > sort_buffer_elements) {
            throw std::invalid_argument(
                "StoredMatrixSorter: invalid/incompatible values for sort_buffer_bytes or "
                "load_bytes"
            );
        }
        if (load_elements < 128) {
            throw std::invalid_argument(
                "StoredMatrixSorter: load_bytes too small (must be big enough to load "
                ">=128 non-zero elements, usually 1536 or 2048)"
            );
        }
    }

    // Write an ouptut stored matrix.
    // Load matrix values using the load_entries() function and sort to be in col-major order
    // row_names, col_names -- row/column names (assuming col-major order)
    // rows, cols -- number of rows/cols (assuming col-major order)
    // row_major -- if true, swap row/col arguments when writing the StoredMatrix format and mark as
    // row-major storage order
    void writeValues(
        std::vector<std::string> &&row_names,
        std::vector<std::string> &&col_names,
        uint32_t rows,
        uint32_t cols,
        std::atomic<bool> *user_interrupt = NULL
    ) {
        if (round != 0)
            throw std::runtime_error(
                "StoredMatrixSorter: can't write more than once to same temporary location"
            );
        std::vector<uint64_t> output_chunk_sizes;
        std::vector<uint64_t> input_chunk_sizes;

        // 1. Read in data + sort in max possible chunks
        {
            openWriters(round);

            // Scope to allocate + free these sorting buffers
            std::vector<uint32_t> row_data(sort_buffer_elements / 2),
                row_buf(sort_buffer_elements / 2), col_data(sort_buffer_elements / 2),
                col_buf(sort_buffer_elements / 2);
            std::vector<T> val_data(sort_buffer_elements / 2), val_buf(sort_buffer_elements / 2);

            // Load input data and write out in sorted chunks
            while (true) {
                size_t loaded = load_entries(row_data, col_data, val_data, user_interrupt);
                if (loaded == 0) break;
                if (!row_sorted()) {
                    lsdRadixSortArrays<uint32_t, uint32_t, T>(
                        loaded, row_data, col_data, val_data, row_buf, col_buf, val_buf
                    );
                }
                lsdRadixSortArrays<uint32_t, uint32_t, T>(
                    loaded, col_data, row_data, val_data, col_buf, row_buf, val_buf
                );
                // Output to tmpStorage
                flushToNumWriter(row_data.data(), loaded, row_writer);
                flushToNumWriter(col_data.data(), loaded, col_writer);
                flushToNumWriter(val_data.data(), loaded, val_writer);
                output_chunk_sizes.push_back(loaded);
            }

            row_writer.finalize();
            col_writer.finalize();
            val_writer.finalize();

            std::swap(output_chunk_sizes, input_chunk_sizes);
        } // row_data, row_buf etc. get freed here
        for (auto &x : input_chunk_sizes)
            total_elements += x;

        // 2. Merge up to (sort_buffer_elements / load_elements) chunks at once into a single sorted
        // chunk
        //    until we have just one sorted chunk encompassing all the data entries
        std::vector<SliceReader<uint32_t>> row_chunks, col_chunks;
        std::vector<SliceReader<T>> val_chunks;

        // heap of which chunk index to draw the next matrix entry from.
        // Format is col, row, chunk_idx such that std::greater will provide the
        // correct lexographic ordering to make a min-heap
        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> heap;

        // Handle cases of empty matrix or a matrix that fits in one chunk
        if (input_chunk_sizes.size() <= 1) {
            // Copy data into the output location
            round += 1;
            openReaders(round - 1);
            openWriters(round, true);

            while (row_reader.requestCapacity()) {
                uint64_t capacity = row_reader.capacity();
                row_writer.ensureCapacity(capacity);
                std::memmove(row_writer.data(), row_reader.data(), capacity * sizeof(uint32_t));
                row_reader.advance(capacity);
                row_writer.advance(capacity);
            }
            row_writer.finalize();

            while (col_reader.requestCapacity()) {
                uint64_t capacity = col_reader.capacity();
                col_writer.ensureCapacity(capacity);
                std::memmove(col_writer.data(), col_reader.data(), capacity * sizeof(uint32_t));
                col_reader.advance(capacity);
                col_writer.advance(capacity);
            }
            col_writer.finalize();

            while (val_reader.requestCapacity()) {
                uint64_t capacity = val_reader.capacity();
                val_writer.ensureCapacity(capacity);
                std::memmove(val_writer.data(), val_reader.data(), capacity * sizeof(T));
                val_reader.advance(capacity);
                val_writer.advance(capacity);
            }
            val_writer.finalize();
        }
        while (input_chunk_sizes.size() > 1) {
            round += 1;
            if (round >= 2) deleteWriters(round - 2);
            openReaders(round - 1);
            if (input_chunk_sizes.size() <= sort_buffer_elements / load_elements) {
                openWriters(round, true); // Open last-round writers with no numeric prefix
            } else {
                openWriters(round);
            }

            output_chunk_sizes.clear();

            uint64_t reader_offset = 0; // Used to calculate offsets for reads

            // Merge up to `row_chunks` sorted chunks at a time, starting at index `chunk`
            for (uint64_t chunk = 0; chunk < input_chunk_sizes.size();
                 chunk += sort_buffer_elements / load_elements) {
                output_chunk_sizes.push_back(0);
                // Set up heap
                heap.clear();
                row_chunks.clear();
                col_chunks.clear();
                val_chunks.clear();
                for (uint64_t i = 0; chunk + i < input_chunk_sizes.size() &&
                                     i < sort_buffer_elements / load_elements;
                     i++) {
                    row_chunks.push_back(SliceReader<uint32_t>(
                        row_reader, reader_offset, input_chunk_sizes[chunk + i], load_elements
                    ));
                    col_chunks.push_back(SliceReader<uint32_t>(
                        col_reader, reader_offset, input_chunk_sizes[chunk + i], load_elements
                    ));
                    val_chunks.push_back(SliceReader<T>(
                        val_reader, reader_offset, input_chunk_sizes[chunk + i], load_elements
                    ));

                    // Initialize the readers with the first batch of data
                    row_chunks[i].advance();
                    col_chunks[i].advance();
                    val_chunks[i].advance();

                    heap.push_back({col_chunks[i].value(), row_chunks[i].value(), i});

                    reader_offset += input_chunk_sizes[chunk + i];
                }
                dary_heap::make_heap<2>(heap.begin(), heap.end(), std::greater<>{});

                // helper to get second to top value in heap
                auto get_second_to_top =
                    [](const std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> &heap) {
                        auto ret = std::make_tuple(UINT32_MAX, UINT32_MAX, UINT32_MAX);
                        uint64_t child_1 = dary_heap::dary_heap_helpers::first_child_index<2>(0);
                        uint64_t child_2 = dary_heap::dary_heap_helpers::last_child_index<2>(0);
                        if (child_1 < heap.size()) ret = std::min(heap[child_1], ret);
                        if (child_2 < heap.size()) ret = std::min(heap[child_2], ret);
                        return ret;
                    };
                auto second_to_top = get_second_to_top(heap);
                while (heap.size() > 0) {
                    if (output_chunk_sizes.back() % (1 << 14) == 0 && user_interrupt != NULL &&
                        *user_interrupt)
                        return;
                    // Output element
                    uint32_t idx = std::get<2>(heap.front());
                    col_writer.write_one(std::get<0>(heap.front()));
                    if (std::get<0>(heap.front()) != col_chunks[idx].value())
                        throw std::runtime_error("Had error with col_chunks matching");
                    row_writer.write_one(std::get<1>(heap.front()));
                    if (std::get<1>(heap.front()) != row_chunks[idx].value())
                        throw std::runtime_error("Had error with row_chunks matching");
                    val_writer.write_one(val_chunks[idx].value());
                    output_chunk_sizes.back() += 1;

                    bool has_more = row_chunks[idx].advance();
                    col_chunks[idx].advance();
                    val_chunks[idx].advance();

                    if (has_more) {
                        std::get<0>(heap.front()) = col_chunks[idx].value();
                        std::get<1>(heap.front()) = row_chunks[idx].value();
                    }

                    // Correctly reposition the head in the heap
                    // Because we will often have repeated elements coming from the same
                    // chunk, do a quick test to see if we even need to adjust the heap.
                    // (Why often repeated? If a chunk spans several adjacent columns of input, then
                    // it will span several adjacent rows of the output)
                    if (!has_more || second_to_top < heap.front()) {
                        dary_heap::pop_heap<2>(heap.begin(), heap.end(), std::greater<>{});
                        if (!has_more) heap.pop_back();
                        else dary_heap::push_heap<2>(heap.begin(), heap.end(), std::greater<>{});
                        second_to_top = get_second_to_top(heap);
                    }
                }
            }
            row_writer.finalize();
            col_writer.finalize();
            val_writer.finalize();
            std::swap(input_chunk_sizes, output_chunk_sizes);
        }

        if (round > 0) deleteWriters(round - 1);

        // 3. Write final metadata to the output writer, such that we can easily create
        // a StoredMatrix object from this later

        openReaders(round, true); // Open just the col_reader as needed

        // Store the col_ptr vector
        ULongWriter col_ptr = out_writer_builder.createULongWriter("idxptr");
        uint32_t cols_seen = 0;
        for (uint64_t i = 0; i < total_elements; i++) {
            uint32_t next_col = col_reader.read_one();
            while (next_col >= cols_seen) {
                col_ptr.write_one(i);
                cols_seen += 1;
            }
        }
        deleteWriters(round, true); // Ready to delete the col_data
        while (cols_seen <= cols) {
            cols_seen += 1;
            col_ptr.write_one(total_elements);
        }
        col_ptr.finalize();

        // Store row and col names. This probably incurs a few extra copies,
        // but it shouldn't matter since writing the actual matrix should dominate cost
        if (row_major) std::swap(row_names, col_names);
        out_writer_builder.createStringWriter("col_names")->write(VecStringReader(col_names));
        out_writer_builder.createStringWriter("row_names")->write(VecStringReader(row_names));

        auto shape = out_writer_builder.createUIntWriter("shape");
        if (row_major) {
            shape.write_one(cols);
            shape.write_one(rows);
        } else {
            shape.write_one(rows);
            shape.write_one(cols);
        }
        shape.finalize();

        out_writer_builder.createStringWriter("storage_order")
            ->write(VecStringReader({row_major ? "row" : "col"}));

        out_writer_builder.writeVersion(StoredMatrix<T>::versionString(true, 2));
    }
};

template <typename T> void StoredMatrixSorter<T>::openWriters(uint64_t round, bool last_round) {
    writer_round = round;

    std::string prefix = last_round ? "" : (std::to_string(round) + "_");
    WriterBuilder &wb = last_round ? out_writer_builder : file_writer_builder;

    // D1Z for rows
    row_writer = UIntWriter(
        std::make_unique<BP128_D1Z_UIntWriter>(
            wb.createUIntWriter(prefix + "index_data"),
            wb.createUIntWriter(prefix + "index_idx"),
            wb.createULongWriter(prefix + "index_idx_offsets"),
            wb.createUIntWriter(prefix + "index_starts")
        ),
        1024
    );
    // D1 for cols. Always write to the temporary output
    col_writer = UIntWriter(
        std::make_unique<BP128_D1_UIntWriter>(
            file_writer_builder.createUIntWriter(std::to_string(round) + "_col_data"),
            file_writer_builder.createUIntWriter(std::to_string(round) + "_col_idx"),
            file_writer_builder.createULongWriter(std::to_string(round) + "_col_idx_offsets"),
            file_writer_builder.createUIntWriter(std::to_string(round) + "_col_starts")
        ),
        1024
    );
    // BP128 if vals is uint32_t, otherwise uncompressed
    if constexpr (std::is_same_v<T, uint32_t>) {
        val_writer = UIntWriter(
            std::make_unique<BP128_FOR_UIntWriter>(
                wb.createUIntWriter(prefix + "val_data"),
                wb.createUIntWriter(prefix + "val_idx"),
                wb.createULongWriter(prefix + "val_idx_offsets")
            ),
            1024
        );
    } else {
        val_writer = wb.create<T>(prefix + "val");
    }
}

template <typename T> void StoredMatrixSorter<T>::openReaders(uint64_t round, bool last_round) {
    reader_round = round;

    std::string prefix = last_round ? "" : (std::to_string(round) + "_");

    col_reader = UIntReader(
        std::make_unique<BP128_D1_UIntReader>(
            file_reader_builder.openUIntReader(std::to_string(round) + "_col_data"),
            file_reader_builder.openUIntReader(std::to_string(round) + "_col_idx"),
            file_reader_builder.openULongReader(std::to_string(round) + "_col_idx_offsets"),
            file_reader_builder.openUIntReader(std::to_string(round) + "_col_starts"),
            total_elements
        ),
        1024,
        1024
    );

    // Only update row_reader and val_reader if we're not on the last round
    if (last_round) return;

    // D1Z for rows
    row_reader = UIntReader(
        std::make_unique<BP128_D1Z_UIntReader>(
            file_reader_builder.openUIntReader(prefix + "index_data"),
            file_reader_builder.openUIntReader(prefix + "index_idx"),
            file_reader_builder.openULongReader(prefix + "index_idx_offsets"),
            file_reader_builder.openUIntReader(prefix + "index_starts"),
            total_elements
        ),
        1024,
        1024
    );
    // D1 for cols
    // BP128 if vals is uint32_t, otherwise uncompressed
    if constexpr (std::is_same_v<T, uint32_t>) {
        val_reader = UIntReader(
            std::make_unique<BP128_FOR_UIntReader>(
                file_reader_builder.openUIntReader(prefix + "val_data"),
                file_reader_builder.openUIntReader(prefix + "val_idx"),
                file_reader_builder.openULongReader(prefix + "val_idx_offsets"),
                total_elements
            ),
            1024,
            1024
        );
    } else {
        val_reader = file_reader_builder.open<T>(prefix + "val");
    }
}

template <typename T> void StoredMatrixSorter<T>::deleteWriters(uint64_t round, bool last_round) {
    // TODO: If I ever change the makeup of BP128 storage to support > 2^34 bytes,
    //  then I'll need to add another deletion here
    std::string prefix = last_round ? "" : (std::to_string(round) + "_");

    // Remove open file handles if needed
    if (reader_round == round) col_reader = UIntReader();
    if (writer_round == round) col_writer = UIntWriter();

    // D1 for cols
    file_writer_builder.deleteWriter(std::to_string(round) + "_col_data");
    file_writer_builder.deleteWriter(std::to_string(round) + "_col_idx");
    file_writer_builder.deleteWriter(std::to_string(round) + "_col_idx_offsets");
    file_writer_builder.deleteWriter(std::to_string(round) + "_col_starts");

    if (last_round) return; // Don't delete the actually needed data

    // Remove open file handles if needed
    if (reader_round == round) {
        row_reader = UIntReader();
        val_reader = NumReader<T>();
    }
    if (writer_round == round) {
        row_writer = UIntWriter();
        val_writer = NumWriter<T>();
    }

    // D1Z for rows
    file_writer_builder.deleteWriter(prefix + "index_data");
    file_writer_builder.deleteWriter(prefix + "index_idx");
    file_writer_builder.deleteWriter(prefix + "index_idx_offsets");
    file_writer_builder.deleteWriter(prefix + "index_starts");

    // BP128 if vals is uint32_t, otherwise uncompressed
    if constexpr (std::is_same_v<T, uint32_t>) {
        file_writer_builder.deleteWriter(prefix + "val_data");
        file_writer_builder.deleteWriter(prefix + "val_idx");
        file_writer_builder.deleteWriter(prefix + "val_idx_offsets");
    } else {
        file_writer_builder.deleteWriter(prefix + "val");
    }
}

template <typename T>
template <typename V>
void StoredMatrixSorter<T>::flushToNumWriter(const V *vals, uint64_t count, NumWriter<V> &writer) {
    uint64_t chunk_size = writer.maxCapacity();
    for (uint64_t i = 0; i < count; i += chunk_size) {
        uint64_t write_size = std::min(chunk_size, count - i);
        writer.ensureCapacity(write_size);
        std::memmove(writer.data(), vals + i, sizeof(V) * write_size);
        writer.advance(write_size);
    }
}

template <typename T>
template <typename V>
StoredMatrixSorter<T>::SliceReader<V>::SliceReader(
    NumReader<V> &reader, uint64_t start_offset, uint64_t count, uint64_t chunk_size
)
    : reader(reader)
    , reader_idx(start_offset)
    , reader_end_idx(start_offset + count) {

    data.resize(chunk_size);
}

template <typename T> template <typename V> bool StoredMatrixSorter<T>::SliceReader<V>::advance() {
    data_idx += 1;
    if (data_idx >= num_loaded) {
        data_idx -= 1;
        if (reader_idx >= reader_end_idx) return false;
        reader.seek(reader_idx);
        num_loaded = 0;
        while (num_loaded < data.size()) {
            if (!reader.requestCapacity()) break;
            uint64_t copy_size = std::min(reader.capacity(), (uint64_t)data.size() - num_loaded);
            std::memmove(data.data() + num_loaded, reader.data(), sizeof(V) * copy_size);
            reader.advance(copy_size);
            num_loaded += copy_size;
        }
        num_loaded = std::min(num_loaded, reader_end_idx - reader_idx);
        reader_idx += num_loaded;
        data_idx = 0;
    }
    return true;
}

template <typename T> template <typename V> V StoredMatrixSorter<T>::SliceReader<V>::value() const {
    return data[data_idx];
}

} // end namespace BPCells