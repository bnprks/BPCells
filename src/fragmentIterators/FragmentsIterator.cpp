#include "FragmentsIterator.h"

namespace BPCells {

// Construct iterator with a given internal buffer size (must be a power of 2 >= 128)
FragmentsIterator::FragmentsIterator(FragmentsLoader &loader, uint32_t buffer_size) : 
    FragmentsLoaderWrapper(loader),
    chunk_capacity(buffer_size), 
    chunk_size(buffer_size), 
    idx(buffer_size), 
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

int32_t FragmentsIterator::load(uint32_t count, FragmentArray &buffer) {return loader.load(count, buffer);};

FragmentsLoaderWrapper::FragmentsLoaderWrapper(FragmentsLoader &loader) : loader(loader) {};

bool FragmentsLoaderWrapper::isSeekable() const {return loader.isSeekable();}
void FragmentsLoaderWrapper::seek(uint32_t chr_id, uint32_t base) {loader.seek(chr_id, base);}

void FragmentsLoaderWrapper::restart() {loader.restart();}

int FragmentsLoaderWrapper::chrCount() const {return loader.chrCount();}
int FragmentsLoaderWrapper::cellCount() const {return loader.cellCount();}

const char* FragmentsLoaderWrapper::chrNames(uint32_t chr_id) const {return loader.chrNames(chr_id);}
const char* FragmentsLoaderWrapper::cellNames(uint32_t cell_id) const {return loader.cellNames(cell_id);}

bool FragmentsLoaderWrapper::nextChr() {return loader.nextChr();}
uint32_t FragmentsLoaderWrapper::currentChr() const {return loader.currentChr();}

} // end namespace BPCells