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
    const size_t load_elements, sort_buffer_elements;
    size_t round = 0, total_elements=0;
    
private:
    UIntReader row_reader, col_reader;
    NumReader<T> val_reader;
    UIntWriter row_writer, col_writer;
    NumWriter<T> val_writer;

    FileWriterBuilder file_writer_builder;
    FileReaderBuilder file_reader_builder;

    void openWriters(uint32_t round);
    void openReaders(uint32_t round);
    void deleteWriters(uint32_t round);

    template<typename S>
    NumWriter<S> openWriter(WriterBuilder &wb, std::string name);

    template<typename S>
    NumWriter<S> openReader(WriterBuilder &wb, std::string name);

    static size_t bytesPerElement() {return 2*sizeof(uint32_t) + sizeof(T);}

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
        load_elements(load_bytes/sizeof(uint32_t)), sort_buffer_elements(sort_buffer_bytes/(bytesPerElement())),
        file_writer_builder(tmpdir, 8192), file_reader_builder(tmpdir, 8192) {
        
        if (sort_buffer_elements == 0 || load_elements == 0 || load_elements * 2 > sort_buffer_elements) {
            throw std::invalid_argument("StoredMatrixTransposeWriter: invalid/incompatible values for sort_buffer_bytes or load_bytes");
        }
        if (load_elements < 128) {
            throw std::invalid_argument("StoredMatrixTransposeWriter: load_bytes too small (must be big enough to load >=128 non-zero elements, usually 1536 or 2048)");
        }
    }

    StoredMatrix<T> read() {
        openReaders(round);
        return StoredMatrix<T>(
            std::move(row_reader),
            std::move(val_reader),
            file_reader_builder.openUIntReader("col_ptr"),
            file_reader_builder.openUIntReader("row_count"),
            file_reader_builder.openStringReader("row_names"),
            file_reader_builder.openStringReader("col_names")
        ); 
    }

    void write(MatrixLoader<T> &mat, void (*checkInterrupt)(void) = NULL) override {
        if (round != 0)
            throw std::runtime_error("StoredMatrixTransposeWriter: can't write more than once to same temporary location");
        std::vector<size_t> output_chunk_sizes;
        std::vector<size_t> input_chunk_sizes;

        // 1. Read in data + sort in max possible chunks
        { 
            openWriters(round);
            
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
                        
                        lsdRadixSortArrays<uint32_t, T>(elems, row_data, col_data, val_data, row_buf, col_buf, val_buf);
                        lsdRadixSortArrays<uint32_t, T>(elems, col_data, row_data, val_data, col_buf, row_buf, val_buf);
                        // Output to tmpStorage
                        flushToNumWriter(row_data.data(), elems, row_writer);
                        flushToNumWriter(col_data.data(), elems, col_writer);
                        flushToNumWriter(val_data.data(), elems, val_writer);

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
                lsdRadixSortArrays<uint32_t, T>(elems, row_data, col_data, val_data, row_buf, col_buf, val_buf);
                lsdRadixSortArrays<uint32_t, T>(elems, col_data, row_data, val_data, col_buf, row_buf, val_buf);
                // Output to tmpStorage
                flushToNumWriter(row_data.data(), elems, row_writer);
                flushToNumWriter(col_data.data(), elems, col_writer);
                flushToNumWriter(val_data.data(), elems, val_writer);
                
                output_chunk_sizes.push_back(elems);
            }
            row_writer.finalize();
            col_writer.finalize();
            val_writer.finalize();

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
        
        while (input_chunk_sizes.size() > 1) {
            round += 1;
            if (round >= 2) deleteWriters(round - 2);
            openWriters(round);
            openReaders(round-1);

            output_chunk_sizes.clear();

            uint32_t reader_offset = 0; // Used to calculate offsets for reads

            // Merge up to `row_chunks` sorted chunks at a time, starting at index `chunk`
            for (int chunk = 0; chunk < input_chunk_sizes.size(); chunk += sort_buffer_elements / load_elements) {
                output_chunk_sizes.push_back(0);
                // Set up heap
                heap.clear();
                row_chunks.clear(); col_chunks.clear(); val_chunks.clear();
                for (int i = 0; chunk + i < input_chunk_sizes.size() && i < sort_buffer_elements / load_elements; i++) {
                    row_chunks.push_back(SliceReader<uint32_t>(row_reader, reader_offset, input_chunk_sizes[chunk + i], load_elements));
                    col_chunks.push_back(SliceReader<uint32_t>(col_reader, reader_offset, input_chunk_sizes[chunk + i], load_elements));
                    val_chunks.push_back(SliceReader<T>(val_reader, reader_offset, input_chunk_sizes[chunk + i], load_elements));
                    
                    // Initialize the readers with the first batch of data
                    row_chunks[i].advance(); col_chunks[i].advance(); val_chunks[i].advance();
                    
                    heap.push_back({col_chunks[i].value(), row_chunks[i].value(), i});


                    reader_offset += input_chunk_sizes[chunk + i];
                }
                dary_heap::make_heap<2>(heap.begin(), heap.end(), std::greater<>{});
                
                // helper to get second to top value in heap
                auto get_second_to_top = [](const std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> &heap) {
                    auto ret = std::make_tuple(UINT32_MAX, UINT32_MAX, UINT32_MAX);
                    size_t child_1 = dary_heap::dary_heap_helpers::first_child_index<2>(0);
                    size_t child_2 = dary_heap::dary_heap_helpers::last_child_index<2>(0);
                    if (child_1 < heap.size()) ret = std::min(heap[child_1], ret);
                    if (child_2 < heap.size()) ret = std::min(heap[child_2], ret);
                    return ret;
                };
                auto second_to_top = get_second_to_top(heap);
                while (heap.size() > 0) {
                    if (output_chunk_sizes.back() %  (1<<14) == 0 && checkInterrupt != NULL) checkInterrupt();
                    // Output element
                    uint32_t idx = std::get<2>(heap.front());
                    col_writer.write_one(std::get<0>(heap.front()));
                    if (std::get<0>(heap.front()) != col_chunks[idx].value()) throw std::runtime_error("Had error with col_chunks matching");
                    row_writer.write_one(std::get<1>(heap.front()));
                    if (std::get<1>(heap.front()) != row_chunks[idx].value()) throw std::runtime_error("Had error with row_chunks matching");
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
                    // (Why often repeated? If a chunk spans several adjacent columns of input, then it
                    //  will span several adjacent rows of the output)
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
        openReaders(round);

        // Store the col_ptr vector
        UIntWriter col_ptr = file_writer_builder.createUIntWriter("col_ptr");
        uint32_t cols_seen = 0;
        for (size_t i = 0; i < total_elements; i++) {
            uint32_t next_col = col_reader.read_one();
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
            const char* col_name = mat.rowNames(i);
            if (col_name == NULL) break;
            col_names.push_back(std::string(col_name));
        }
        file_writer_builder.createStringWriter("col_names")->write(VecStringReader(col_names));
                
        std::vector<std::string> row_names(0);
        for (int i = 0; ; i++) {
            const char* row_name = mat.colNames(i);
            if (row_name == NULL) break;
            row_names.push_back(std::string(row_name));
        }
        file_writer_builder.createStringWriter("row_names")->write(VecStringReader(row_names));
        
        auto row_count = file_writer_builder.createUIntWriter("row_count");
        row_count.write_one(mat.cols());
        row_count.finalize();
    }
};


template<typename T>
void StoredMatrixTransposeWriter<T>::openWriters(uint32_t round) {
    // D1Z for rows
    row_writer = UIntWriter(std::make_unique<BP128_D1Z_UIntWriter>(
        file_writer_builder.createUIntWriter(std::to_string(round) + "_row_data"),
        file_writer_builder.createUIntWriter(std::to_string(round) + "_row_idx"),
        file_writer_builder.createUIntWriter(std::to_string(round) + "_row_starts")
    ), 1024);
    // D1 for cols
    col_writer = UIntWriter(std::make_unique<BP128_D1_UIntWriter>(
        file_writer_builder.createUIntWriter(std::to_string(round) + "_col_data"),
        file_writer_builder.createUIntWriter(std::to_string(round) + "_col_idx"),
        file_writer_builder.createUIntWriter(std::to_string(round) + "_col_starts")
    ), 1024);
    // BP128 if vals is uint32_t, otherwise uncompressed
    if constexpr (std::is_same_v<T, uint32_t>) {
        val_writer = UIntWriter(std::make_unique<BP128_FOR_UIntWriter>(
            file_writer_builder.createUIntWriter(std::to_string(round) + "_val_data"),
            file_writer_builder.createUIntWriter(std::to_string(round) + "_val_idx")
        ), 1024);
    } else {
        val_writer = file_writer_builder.create<T>(std::to_string(round) + "_val");
    }
}

template<typename T>
void StoredMatrixTransposeWriter<T>::openReaders(uint32_t round) {
    // D1Z for rows
    row_reader = UIntReader(std::make_unique<BP128_D1Z_UIntReader>(
        file_reader_builder.openUIntReader(std::to_string(round) + "_row_data"),
        file_reader_builder.openUIntReader(std::to_string(round) + "_row_idx"),
        file_reader_builder.openUIntReader(std::to_string(round) + "_row_starts"),
        total_elements
    ), 1024, 1024);
    // D1 for cols
    col_reader = UIntReader(std::make_unique<BP128_D1_UIntReader>(
        file_reader_builder.openUIntReader(std::to_string(round) + "_col_data"),
        file_reader_builder.openUIntReader(std::to_string(round) + "_col_idx"),
        file_reader_builder.openUIntReader(std::to_string(round) + "_col_starts"),
        total_elements
    ), 1024, 1024);
    // BP128 if vals is uint32_t, otherwise uncompressed
    if constexpr (std::is_same_v<T, uint32_t>) {
        val_reader = UIntReader(std::make_unique<BP128_FOR_UIntReader>(
            file_reader_builder.openUIntReader(std::to_string(round) + "_val_data"),
            file_reader_builder.openUIntReader(std::to_string(round) + "_val_idx"),
            total_elements
        ), 1024, 1024);
    } else {
        val_reader = file_reader_builder.open<T>(std::to_string(round) + "_val");
    }
}

template<typename T>
void StoredMatrixTransposeWriter<T>::deleteWriters(uint32_t round) {
    //TODO: If I ever change the makeup of BP128 storage to support > 2^34 bytes,
    // then I'll need to add another deletion here

    // D1Z for rows
    file_writer_builder.deleteWriter(std::to_string(round) + "_row_data");
    file_writer_builder.deleteWriter(std::to_string(round) + "_row_idx");
    file_writer_builder.deleteWriter(std::to_string(round) + "_row_starts");

    // D1 for cols
    file_writer_builder.deleteWriter(std::to_string(round) + "_col_data");
    file_writer_builder.deleteWriter(std::to_string(round) + "_col_idx");
    file_writer_builder.deleteWriter(std::to_string(round) + "_col_starts");

    // BP128 if vals is uint32_t, otherwise uncompressed
    if constexpr (std::is_same_v<T, uint32_t>) {
        file_writer_builder.deleteWriter(std::to_string(round) + "_val_data");
        file_writer_builder.deleteWriter(std::to_string(round) + "_val_idx");
    } else {
        file_writer_builder.deleteWriter(std::to_string(round) + "_val");
    }
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
        num_loaded = 0;
        while (num_loaded < data.size()) {
            if (!reader.requestCapacity()) break;
            uint32_t copy_size = std::min(reader.capacity(), (uint32_t) data.size() - num_loaded);
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
    
template<typename T>
template<typename V>
V StoredMatrixTransposeWriter<T>::SliceReader<V>::value() const {
    return data[data_idx];
}

} // end namespace BPCells