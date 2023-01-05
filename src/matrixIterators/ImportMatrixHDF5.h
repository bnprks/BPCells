#pragma once
#include "../arrayIO/array_interfaces.h"
#include "../arrayIO/hdf5.h"
#include "../arrayIO/vector.h"
#include "StoredMatrix.h"
#include "StoredMatrixWriter.h"

namespace BPCells {

// Helper class for when we need to fake reading a single integer from memory
template <class T> class SingletonNumReader : public BulkNumReader<T> {
  private:
    T num;
    bool read = false;

  public:
    SingletonNumReader(T num) : num(num) {}
    uint32_t size() const override { return 1; }
    void seek(uint32_t pos) override { read = pos > 0; }
    uint32_t load(T *out, uint32_t count) override {
        if (read) return 0;
        out[0] = num;
        return 1;
    }
};

// Reader interfaces for 10x and AnnData matrices

StoredMatrix<uint32_t>
open10xFeatureMatrix(std::string file, uint32_t buffer_size, uint32_t read_size = 1024);
StoredMatrixWriter<uint32_t> create10xFeatureMatrix(
    std::string file,
    const StringReader &barcodes,
    const StringReader &feature_ids,
    const StringReader &feature_names,
    const StringReader &feature_types,
    const std::map<std::string, std::unique_ptr<StringReader>> &feature_metadata,
    uint32_t buffer_size,
    uint32_t chunk_size
);

// Read AnnData sparse matrix, with an implicit transpose to CSC format for
// any data stored in CSR format
StoredMatrix<float> openAnnDataMatrix(
    std::string file, std::string group, uint32_t buffer_size, uint32_t read_size = 1024
);

bool isRowOrientedAnnDataMatrix(std::string file, std::string group);

// Read a single member from a struct-type dataset
template <typename T>
void readMember(HighFive::DataSet &&dataset, std::string name, std::vector<T> &out);

} // end namespace BPCells