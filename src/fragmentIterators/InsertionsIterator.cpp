#include "InsertionsIterator.h"
namespace BPCells {

// Construct iterator with a given internal buffer size (must be a power of 2 >= 128)
InsertionsIterator::InsertionsIterator(FragmentsLoader &loader, uint32_t buffer_size) : 
    FragmentsLoaderWrapper(loader),
    chunk_capacity(buffer_size), 
    chunk_size(buffer_size), 
    last_start(0), 
    current_chr(UINT32_MAX), 
    start_buf(buffer_size), end_buf(buffer_size), cell_buf(buffer_size) {

    if (buffer_size < 128)
        throw std::invalid_argument("buffer_size must be >= 128");
    else if (buffer_size & (buffer_size - 1)) 
        throw std::invalid_argument("buffer_size must be a power of 2");
    fragments_buf.start = &start_buf[0];
    fragments_buf.end = &end_buf[0];
    fragments_buf.cell = &cell_buf[0];
}

int32_t InsertionsIterator::load(uint32_t count, FragmentArray &buffer) {return loader.load(count, buffer);};

} // end namespace BPCells