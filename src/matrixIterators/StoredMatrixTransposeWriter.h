#pragma once

#include <string>
#include <filesystem>

#include "../lib/dary_heap.hpp"
#include "../utils/radix_sort.h"
#include "../arrayIO/array_interfaces.h"
#include "../arrayIO/binaryfile.h"
#include "../arrayIO/vector.h"
#include "StoredMatrix.h"
#include "MatrixIterator.h"

namespace BPCells {

template<typename T> 
class StoredMatrixTransposeWriter : public MatrixWriter<T> {
    std::filesystem::path tmpdir;
    const size_t load_elements, sort_buffer_elements;
    size_t round = 0, total_elements=0;
private:
    static size_t bytesPerElement() {return 2*sizeof(uint32_t) + sizeof(T);}

    // I/O routines to read/write data from the temporary folder for a given merge round
    // By default, row will be BP128 D1Z compressed, col will be BP128 D1 compressed, 
    // and val will be BP128 FOR compressed if T is uint32_t, else uncompressed
    class TmpWriters {
    public:
        FileWriterBuilder wb;
        UIntWriter row, col;
        NumWriter<T> val;
        TmpWriters(uint32_t round, const std::filesystem::path &tmpdir, size_t load_elements);
        TmpWriters() = default;
        TmpWriters& operator=(TmpWriters&& other) = default;
    };

    class TmpReaders {
    public:
        FileReaderBuilder rb;
        UIntReader row, col;
        NumReader<T> val;
        TmpReaders(uint32_t round, const std::filesystem::path &tmpdir, size_t load_elements, size_t total_elements);
        TmpReaders() = default;
        TmpReaders& operator=(TmpReaders&& other) = default;
    };

    static void deleteTmpData(uint32_t round, const std::filesystem::path &tmpdir);

    template<typename V>
    static void flushToNumWriter(const V *vals, size_t count,  NumWriter<V> &writer);

    // Helper class to handle reading from a slice of a NumReader
    // Allow sharing NumReader between multiple SliceReaders, meaning that a 
    // seek must be performed before each load
    template<typename V>
    class SliceReader {
        NumReader<V> &reader;
        std::vector<V> data;
        uint32_t data_idx=1, num_loaded=0, reader_idx;
        const uint32_t reader_end_idx;
    public:
        SliceReader(NumReader<V> &reader, uint32_t start_offset, uint32_t count, uint32_t chunk_size);
        bool advance();
        V value() const;
    };

public:
    // Here, load_bytes is the smallest read size to use in bytes, and sort_buffer_bytes
    // is the total amount of memory to use for sort buffers
    StoredMatrixTransposeWriter(const char * tmpdir, size_t load_bytes, size_t sort_buffer_bytes): 
        tmpdir(tmpdir), 
        load_elements(load_bytes/sizeof(uint32_t)), sort_buffer_elements(sort_buffer_bytes/(bytesPerElement())) {
        
        if (sort_buffer_elements == 0 || load_elements == 0 || load_elements * 2 > sort_buffer_elements) {
            throw std::invalid_argument("StoredMatrixTransposeWriter: invalid/incompatible values for sort_buffer_bytes or load_bytes");
        }
        if (load_elements < 128) {
            throw std::invalid_argument("StoredMatrixTransposeWriter: load_bytes too small (must be big enough to load >=128 non-zero elements, usually 1536 or 2048)");
        }
    }

    StoredMatrix<T> read() {
        TmpReaders readers(round, tmpdir, 1024, total_elements);
        return StoredMatrix<T>(
            std::move(readers.row),
            std::move(readers.val),
            readers.rb.openUIntReader("col_ptr"),
            readers.rb.openUIntReader("row_count"),
            readers.rb.openStringReader("row_names"),
            readers.rb.openStringReader("col_names")
        ); 
    }

    void write(MatrixLoader<T> &mat, void (*checkInterrupt)(void) = NULL) override {
        if (round != 0)
            throw std::runtime_error("StoredMatrixTransposeWriter: can't write more than once to same temporary location");
        std::vector<size_t> output_chunk_sizes;
        std::vector<size_t> input_chunk_sizes;

        TmpWriters writers;

        // 1. Read in data + sort in max possible chunks
        { 
            writers = TmpWriters(round, tmpdir, load_elements);
            // Scope to allocate + free these sorting buffers
            std::vector<uint32_t> row_data(sort_buffer_elements / 2), 
                                row_buf(sort_buffer_elements / 2), 
                                col_data(sort_buffer_elements / 2),
                                col_buf(sort_buffer_elements / 2);
            std::vector<T> val_data(sort_buffer_elements / 2), 
                        val_buf(sort_buffer_elements / 2);
            size_t elems = 0;


            while (mat.nextCol()) {
                uint32_t current_col = mat.currentCol();
                while (mat.load()) {
                    uint32_t capacity = mat.capacity();
                    uint32_t read = 0;
                    // Write out sort buffers if they get full
                    while (elems + capacity - read > row_data.size()) {
                        if (checkInterrupt != NULL) checkInterrupt();
                        // Fill out the rest of the buffer size
                        std::memmove(val_data.data() + elems, mat.valData()+read, (row_data.size()-elems) * sizeof(T));
                        std::memmove(col_data.data() + elems, mat.rowData()+read, (row_data.size()-elems) * sizeof(uint32_t));
                        for (uint32_t i = 0; i < (row_data.size()-elems); i++)
                            row_data[elems + i] = current_col;
                        read += (row_data.size()-elems);
                        elems = row_data.size();
                        // Sort by (col, row). Remember though that since the input 
                        // is column-sorted, our output is already row-sorted and 
                        // therefore just needs to be column-sorted
                        // lsdRadixSortArrays<uint32_t, T>(elems, row_data, col_data, val_data, row_buf, col_buf, val_buf);
                        lsdRadixSortArrays<uint32_t, T>(elems, col_data, row_data, val_data, col_buf, row_buf, val_buf);
                        // Output to tmpStorage
                        flushToNumWriter(row_data.data(), elems, writers.row);
                        flushToNumWriter(col_data.data(), elems, writers.col);
                        flushToNumWriter(val_data.data(), elems, writers.val);

                        output_chunk_sizes.push_back(elems);
                        elems = 0;
                    }
                    std::memmove(val_data.data() + elems, mat.valData() + read, (capacity-read) * sizeof(T));
                    std::memmove(col_data.data() + elems, mat.rowData() + read, (capacity-read) * sizeof(uint32_t));
                    for (uint32_t i = 0; i < capacity-read; i++)
                        row_data[elems + i] = current_col;

                    elems += capacity - read;
                }
            }
            // Write out final elements in sort buffers
            if (elems > 0) {
                if (checkInterrupt != NULL) checkInterrupt();
                // Sort by (col, row)
                // lsdRadixSortArrays<uint32_t, T>(elems, row_data, col_data, val_data, row_buf, col_buf, val_buf);
                lsdRadixSortArrays<uint32_t, T>(elems, col_data, row_data, val_data, col_buf, row_buf, val_buf);
                // Output to tmpStorage
                flushToNumWriter(row_data.data(), elems, writers.row);
                flushToNumWriter(col_data.data(), elems, writers.col);
                flushToNumWriter(val_data.data(), elems, writers.val);
                
                output_chunk_sizes.push_back(elems);
            }
            writers.row.finalize();
            writers.col.finalize();
            writers.val.finalize();

            std::swap(output_chunk_sizes, input_chunk_sizes);
        } // row_data, row_buf etc. get freed here
        for (auto x : input_chunk_sizes) total_elements += x;

        // 2. Merge up to (sort_buffer_elements / load_elements) chunks at once into a single sorted chunk
        //    until we have just one sorted chunk encompassing all the data entries
        std::vector<SliceReader<uint32_t>> row_chunks, col_chunks;
        std::vector<SliceReader<T>> val_chunks; 

        // heap of which chunk index to draw the next matrix entry from.
        // Format is col, row, chunk_idx such that std::greater will provide the
        // correct lexographic ordering to make a min-heap
        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> heap; 
        TmpReaders readers = TmpReaders(round, tmpdir, load_elements, total_elements); // dummy initialization because I don't have default constructors
        while (input_chunk_sizes.size() > 1) {
            round += 1;
            if (round >= 2) deleteTmpData(round - 2, tmpdir);
            writers = TmpWriters(round, tmpdir, load_elements);
            readers = TmpReaders(round-1, tmpdir, load_elements, total_elements);

            output_chunk_sizes.clear();

            uint32_t reader_offset = 0; // Used to calculate offsets for reads

            // Merge up to `row_chunks` sorted chunks at a time, starting at index `chunk`
            for (int chunk = 0; chunk < input_chunk_sizes.size(); chunk += sort_buffer_elements / load_elements) {
                output_chunk_sizes.push_back(0);
                // Set up heap
                heap.clear();
                row_chunks.clear(); col_chunks.clear(); val_chunks.clear();
                for (int i = 0; chunk + i < input_chunk_sizes.size() && i < sort_buffer_elements / load_elements; i++) {
                    row_chunks.push_back(SliceReader<uint32_t>(readers.row, reader_offset, input_chunk_sizes[chunk + i], load_elements));
                    col_chunks.push_back(SliceReader<uint32_t>(readers.col, reader_offset, input_chunk_sizes[chunk + i], load_elements));
                    val_chunks.push_back(SliceReader<T>(readers.val, reader_offset, input_chunk_sizes[chunk + i], load_elements));
                    
                    // Initialize the readers with the first batch of data
                    row_chunks[i].advance(); col_chunks[i].advance(); val_chunks[i].advance();
                    
                    heap.push_back({col_chunks[i].value(), row_chunks[i].value(), i});


                    reader_offset += input_chunk_sizes[chunk + i];
                }
                dary_heap::make_heap<2>(heap.begin(), heap.end(), std::greater<>{});
                
                while (heap.size() > 0) {
                    if (output_chunk_sizes.back() %  (1<<14) == 0 && checkInterrupt != NULL) checkInterrupt();
                    // Output element
                    uint32_t idx = std::get<2>(heap.front());
                    writers.col.write_one(std::get<0>(heap.front()));
                    if (std::get<0>(heap.front()) != col_chunks[idx].value()) throw std::runtime_error("Had error with col_chunks matching");
                    writers.row.write_one(std::get<1>(heap.front()));
                    if (std::get<1>(heap.front()) != row_chunks[idx].value()) throw std::runtime_error("Had error with row_chunks matching");
                    writers.val.write_one(val_chunks[idx].value());
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
                    // (Why often repeated? If a chunk spans several adjacent columns of input, then it
                    //  will span several adjacent rows of the output)
                    uint32_t first_child = dary_heap::dary_heap_helpers::first_child_index<2>(0);
                    uint32_t last_child = dary_heap::dary_heap_helpers::last_child_index<2>(0);
                    bool smaller_first = heap.size() > first_child && heap.front() > heap[first_child];
                    bool smaller_last = heap.size() > last_child && heap.front() > heap[last_child];
                    if (smaller_first || smaller_last || !has_more) {
                        dary_heap::pop_heap<2>(heap.begin(), heap.end(), std::greater<>{});
                        if (!has_more) heap.pop_back();
                        else dary_heap::push_heap<2>(heap.begin(), heap.end(), std::greater<>{});
                    }
                }
            }
            writers.row.finalize();
            writers.col.finalize();
            writers.val.finalize();
            std::swap(input_chunk_sizes, output_chunk_sizes);
        }

        if (round > 0) deleteTmpData(round - 1, tmpdir);
        
        // 3. Write final metadata to the output writer, such that we can easily create
        // a StoredMatrix object from this later
        readers = TmpReaders(round, tmpdir, load_elements, total_elements);

        // Store the col_ptr vector
        UIntWriter col_ptr = writers.wb.createUIntWriter("col_ptr");
        uint32_t cols_seen = 0;
        for (size_t i = 0; i < total_elements; i++) {
            uint32_t next_col = readers.col.read_one();
            while (next_col >= cols_seen) {
                col_ptr.write_one(i);
                cols_seen += 1;
            }
        }
        while (cols_seen <= mat.rows()) {
            cols_seen += 1;
            col_ptr.write_one(total_elements);
        }
        col_ptr.finalize();

        // Store row and col names. This probably incurs a few extra copies,
        // but it shouldn't matter since writing the actual matrix should dominate cost
        std::vector<std::string> col_names(0);
        for (int i = 0; ; i++) {
            const char* col_name = mat.colNames(i);
            if (col_name == NULL) break;
            col_names.push_back(std::string(col_name));
        }
        writers.wb.createStringWriter("col_names")->write(VecStringReader(col_names));
                
        std::vector<std::string> row_names(0);
        for (int i = 0; ; i++) {
            const char* row_name = mat.rowNames(i);
            if (row_name == NULL) break;
            row_names.push_back(std::string(row_name));
        }
        writers.wb.createStringWriter("row_names")->write(VecStringReader(row_names));
        
        auto row_count = writers.wb.createUIntWriter("row_count");
        row_count.write_one(mat.cols());
        row_count.finalize();
    }
};


// I sincerely apologize for shoving so much stuff into the initializers here
template<typename T>
StoredMatrixTransposeWriter<T>::TmpWriters::TmpWriters(uint32_t round, const std::filesystem::path &tmpdir, size_t load_elements) 
    : wb(tmpdir / std::filesystem::path(std::to_string(round)), load_elements)
    , row(std::make_unique<BP128_D1Z_UIntWriter>(
        wb.createUIntWriter("row_data"),
        wb.createUIntWriter("row_idx"),
        wb.createUIntWriter("row_starts")
      ), load_elements)
    , col(std::make_unique<BP128_D1_UIntWriter>(
        wb.createUIntWriter("col_data"),
        wb.createUIntWriter("col_idx"),
        wb.createUIntWriter("col_starts")
      ), load_elements)
    , val(std::move(wb.create<T>("val"))) {}

// Specialization for if we can compress val
template<>
StoredMatrixTransposeWriter<uint32_t>::TmpWriters::TmpWriters(uint32_t round, const std::filesystem::path &tmpdir, size_t load_elements) 
    : wb(tmpdir / std::filesystem::path(std::to_string(round)), load_elements)
    , row(std::make_unique<BP128_D1Z_UIntWriter>(
        wb.createUIntWriter("row_data"),
        wb.createUIntWriter("row_idx"),
        wb.createUIntWriter("row_starts")
      ), load_elements)
    , col(std::make_unique<BP128_D1_UIntWriter>(
        wb.createUIntWriter("col_data"),
        wb.createUIntWriter("col_idx"),
        wb.createUIntWriter("col_starts")
      ), load_elements)
    , val(std::make_unique<BP128_FOR_UIntWriter>(
        wb.createUIntWriter("val_data"),
        wb.createUIntWriter("val_idx")
      ), load_elements) {}

template<typename T>
StoredMatrixTransposeWriter<T>::TmpReaders::TmpReaders(uint32_t round, const std::filesystem::path &tmpdir, size_t load_elements, size_t total_elements)
    : rb(tmpdir / std::filesystem::path(std::to_string(round)), load_elements)
    , row(std::make_unique<BP128_D1Z_UIntReader>(
        rb.openUIntReader("row_data"),
        rb.openUIntReader("row_idx"),
        rb.openUIntReader("row_starts"),
        total_elements
      ), load_elements, load_elements)
    , col(std::make_unique<BP128_D1_UIntReader>(
        rb.openUIntReader("col_data"),
        rb.openUIntReader("col_idx"),
        rb.openUIntReader("col_starts"),
        total_elements
      ), load_elements, load_elements)
    , val(std::move(rb.open<T>("val"))) {}

// Specialization for if we can compress val
template<>
StoredMatrixTransposeWriter<uint32_t>::TmpReaders::TmpReaders(uint32_t round, const std::filesystem::path &tmpdir, size_t load_elements, size_t total_elements) 
    : rb(tmpdir / std::filesystem::path(std::to_string(round)), load_elements)
    , row(std::make_unique<BP128_D1Z_UIntReader>(
        rb.openUIntReader("row_data"),
        rb.openUIntReader("row_idx"),
        rb.openUIntReader("row_starts"),
        total_elements
      ), load_elements, load_elements)
    , col(std::make_unique<BP128_D1_UIntReader>(
        rb.openUIntReader("col_data"),
        rb.openUIntReader("col_idx"),
        rb.openUIntReader("col_starts"),
        total_elements
      ), load_elements, load_elements)
    , val(std::make_unique<BP128_FOR_UIntReader>(
        rb.openUIntReader("val_data"),
        rb.openUIntReader("val_idx"),
        total_elements
      ), load_elements, load_elements) {}

template<typename T>
void StoredMatrixTransposeWriter<T>::deleteTmpData(uint32_t round, const std::filesystem::path &tmpdir) {
    std::filesystem::remove_all(tmpdir / std::filesystem::path(std::to_string(round)));
}

template<typename T>
template<typename V>
void StoredMatrixTransposeWriter<T>::flushToNumWriter(const V *vals, size_t count,  NumWriter<V> &writer) {
    size_t chunk_size = writer.maxCapacity();
    for (size_t i = 0; i < count; i += chunk_size) {
        size_t write_size = std::min(chunk_size, count - i);
        writer.ensureCapacity(write_size);
        std::memmove(writer.data(), vals + i, sizeof(V) * write_size);
        writer.advance(write_size);
    }
}

template<typename T>
template<typename V>
StoredMatrixTransposeWriter<T>::SliceReader<V>::SliceReader(NumReader<V> &reader, uint32_t start_offset, uint32_t count, uint32_t chunk_size) :
    reader(reader), reader_idx(start_offset), reader_end_idx(start_offset + count) {

    data.resize(chunk_size);
}

template<typename T>
template<typename V>
bool StoredMatrixTransposeWriter<T>::SliceReader<V>::advance() {
    data_idx += 1;
    if (data_idx >= num_loaded) {
        data_idx -= 1;
        if (reader_idx >= reader_end_idx) return false;
        reader.seek(reader_idx);
        reader.requestCapacity(data.size());
        std::memmove(data.data(), reader.data(), sizeof(V) * reader.capacity());
        num_loaded = std::min(reader.capacity(), reader_end_idx - reader_idx);
        reader_idx += num_loaded;
        data_idx = 0;
    }
    return true;
}
    
template<typename T>
template<typename V>
V StoredMatrixTransposeWriter<T>::SliceReader<V>::value() const {
    return data[data_idx];
}

} // end namespace BPCells