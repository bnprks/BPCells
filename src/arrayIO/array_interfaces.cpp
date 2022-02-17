#include "array_interfaces.h"

namespace BPCells {

UIntReader::UIntReader(std::unique_ptr<UIntBulkReader> &&reader, uint32_t load_size, uint32_t read_size) :
    buffer(load_size), reader(std::move(reader)), total_size(this->reader->size()), read_size(read_size) {}

UIntWriter::UIntWriter(std::unique_ptr<UIntBulkWriter> &&writer, uint32_t write_size) :
    buffer(write_size), writer(std::move(writer)) {}

VecStringReader::VecStringReader(std::vector<std::string> data) : data(data) {}
const char* VecStringReader::get(uint32_t idx) const {
    if (idx < data.size()) return data[idx].c_str();
    return NULL;
}
uint32_t VecStringReader::size() const {return data.size();}


} // end namespace BPCells