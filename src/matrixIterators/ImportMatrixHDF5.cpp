#include "ImportMatrixHDF5.h"

namespace BPCells {
    
// Reader interfaces for 10x and AnnData matrices

StoredMatrix<uint32_t> open10xFeatureMatrix(std::string file, std::string group, uint32_t buffer_size, uint32_t read_size) {
    H5ReaderBuilder rb(file, group, buffer_size, read_size);

    return StoredMatrix(
        rb.openULongReader("indices").convert<uint32_t>(),
        rb.openUIntReader("data"),
        rb.openULongReader("indptr").convert<uint32_t>(),
        rb.openUIntReader("shape"),
        rb.openStringReader("features/id"),
        rb.openStringReader("barcodes")
    );
}

// Read AnnData sparse matrix, with an implicit transpose to CSC format for
// any data stored in CSR format
StoredMatrix<double> openAnnDataMatrix(std::string file, std::string group, uint32_t buffer_size, uint32_t read_size) {
    H5ReaderBuilder rb(file, group, buffer_size, read_size);
    
    HighFive::SilenceHDF5 s;
    HighFive::Group &g = rb.getGroup();
    std::string encoding;
    std::vector<uint32_t> dims;
    g.getAttribute("encoding-type").read(encoding);
    g.getAttribute("shape").read(dims);

    std::string col_ids;
    std::string row_ids;
    HighFive::Group root = g.getFile().getGroup("/");
    root.getGroup("var").getAttribute("_index").read(col_ids);
    root.getGroup("obs").getAttribute("_index").read(row_ids);
    col_ids = "var/" + col_ids;
    row_ids = "obs/" + row_ids;

    uint32_t rows;
    if (encoding == "csr_matrix") {
        rows = dims[1];
        std::swap(row_ids, col_ids);
    }
    else if (encoding == "csc_matrix") {
        rows = dims[0];
    }
    else throw std::runtime_error("Unsupported matrix encoding: " + encoding);

    return StoredMatrix<double>(
        rb.openUIntReader("indices"),
        rb.openFloatReader("data").convert<double>(),
        rb.openUIntReader("indptr"),
        UIntReader(std::make_unique<SingletonNumReader<uint32_t>>(rows), 1),
        std::make_unique<H5StringReader>(root, row_ids),
        std::make_unique<H5StringReader>(root, col_ids)
    );
}

bool isRowOrientedAnnDataMatrix(std::string file, std::string group) {
    H5ReaderBuilder rb(file, group, 1024, 1024);
    
    HighFive::SilenceHDF5 s;
    HighFive::Group &g = rb.getGroup();
    std::string encoding;
    
    g.getAttribute("encoding-type").read(encoding);

    if (encoding == "csr_matrix") return true;
    else if (encoding == "csc_matrix") return false;
    else throw std::runtime_error("Unsupported matrix encoding: " + encoding);
}


} // end namespace BPCells