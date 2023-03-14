#include "array_interfaces.h"

namespace BPCells {

VecStringReader::VecStringReader(std::vector<std::string> data) : data(data) {}
const char *VecStringReader::get(uint64_t idx) const {
    if (idx < data.size()) return data[idx].c_str();
    return NULL;
}
uint64_t VecStringReader::size() const { return data.size(); }

} // end namespace BPCells
