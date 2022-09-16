#include "LengthSelect.h"

namespace BPCells {

LengthSelect::LengthSelect(FragmentLoader &loader, uint32_t min_len, uint32_t max_len)
    : FragmentLoaderWrapper(loader)
    , min_len(min_len)
    , max_len(max_len) {}

bool LengthSelect::load() {
    loaded = 0;
    // load and filter until we load without filtering out everything
    while (loaded == 0) {
        if (!loader.load()) return false;

        uint32_t *cell = loader.cellData();
        uint32_t *start = loader.startData();
        uint32_t *end = loader.endData();
        uint32_t capacity = loader.capacity();
        for (uint32_t i = 0; i < capacity; i++) {
            cell[loaded] = cell[i];
            start[loaded] = start[i];
            end[loaded] = end[i];
            loaded += (end[i] - start[i]) >= min_len && (end[i] - start[i]) <= max_len;
        }
    }
    return true;
}

uint32_t LengthSelect::capacity() const { return loaded; }

} // end namespace BPCells