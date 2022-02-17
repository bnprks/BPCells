#include <Rcpp.h>

#include "fragmentIterators/FragmentIterator.h"
#include "fragmentIterators/BedFragments.h"
#include "fragmentIterators/StoredFragments.h"

#include "arrayIO/binaryfile.h"
#include "arrayIO/hdf5.h"
#include "arrayIO/vector.h"

#include "R_array_io.h"

using namespace Rcpp;
using namespace BPCells;


// [[Rcpp::export]]
SEXP load_10x_fragments_cpp(std::string path, std::string comment) {
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new BedFragments(path.c_str(), comment.c_str()))
    );
}

// [[Rcpp::export]]
void write_10x_fragments_cpp(std::string path, SEXP fragments, bool append_5th_column=false) {
    XPtr<FragmentLoader> loader(fragments);
    FragmentIterator iter(*loader);
    BedFragmentsWriter writer(path.c_str(), append_5th_column);
    writer.write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP iterate_packed_fragments_cpp(S4 s4) {
    S4ReaderBuilder rb(s4);
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new StoredFragmentsPacked(StoredFragmentsPacked::openPacked(rb)))
    );
}

// [[Rcpp::export]]
IntegerVector calculate_end_max_cpp(IntegerVector end, IntegerVector chr_ptr) {
    std::vector<uint32_t> end_max;
    uint32_t current_max = 0;
    uint32_t prev_max = 0;
    uint32_t current_chr = 0;
    
    uint32_t i;
    for (i = 0; i < end.size(); i++) {
        if (i % 128 == 127) {
            end_max.push_back(std::max(current_max, prev_max));
            prev_max = 0;
        }
        if (chr_ptr[current_chr*2 + 1] >= i) {
            prev_max = std::max(prev_max, current_max);
            current_max = 0;
        }
        current_max = std::max(current_max, (uint32_t) end[i]);
    }
    if (i % 128 != 0)
        end_max.push_back(std::max(current_max, prev_max));
    IntegerVector ret;
    ret.assign(end_max.begin(), end_max.end());
    return ret;
}

// [[Rcpp::export]]
List write_packed_fragments_cpp(SEXP fragments) {
    XPtr<FragmentLoader> loader(fragments);
    FragmentIterator iter(*loader);
    
    ListWriterBuilder wb;
    StoredFragmentsWriter::createPacked(wb).write(iter, &Rcpp::checkUserInterrupt);

    return wb.getList();
}

// [[Rcpp::export]]
SEXP iterate_unpacked_fragments_cpp(S4 s4) {
    S4ReaderBuilder rb(s4);
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new StoredFragments(StoredFragments::openUnpacked(rb)))
    );
}

// [[Rcpp::export]]
List write_unpacked_fragments_cpp(SEXP fragments) {
    XPtr<FragmentLoader> loader(fragments);
    FragmentIterator iter(*loader);
    
    ListWriterBuilder wb;
    StoredFragmentsWriter::createUnpacked(wb).write(iter, &Rcpp::checkUserInterrupt);

    return wb.getList();
}

// [[Rcpp::export]]
bool is_compressed_fragments_file_cpp(std::string dir, uint32_t buffer_size) {
    FileReaderBuilder rb(dir, buffer_size);
    std::string version = rb.readVersion();
    bool compressed = false;
    if (version == "unpacked-fragments-v1") {
        StoredFragments::openUnpacked(rb);
    } else if (version == "packed-fragments-v1") {
        StoredFragmentsPacked::openPacked(rb);
        compressed = true;
    } else {
        throw std::runtime_error(std::string("Fragments directory has unrecognized version ") + version);
    }
    return compressed;
}

// [[Rcpp::export]]
SEXP iterate_unpacked_fragments_file_cpp(std::string dir, uint32_t buffer_size) {
    FileReaderBuilder rb(dir, buffer_size);
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new StoredFragments(StoredFragments::openUnpacked(rb)))
    );
}

// [[Rcpp::export]]
void write_unpacked_fragments_file_cpp(SEXP fragments, std::string dir, uint32_t buffer_size) {
    XPtr<FragmentLoader> loader(fragments);
    FragmentIterator iter(*loader);
    FileWriterBuilder wb(dir, buffer_size);    
    StoredFragmentsWriter::createUnpacked(wb).write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP iterate_packed_fragments_file_cpp(std::string dir, uint32_t buffer_size) {
    FileReaderBuilder rb(dir, buffer_size);
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new StoredFragmentsPacked(StoredFragmentsPacked::openPacked(rb)))
    );
}

// [[Rcpp::export]]
void write_packed_fragments_file_cpp(SEXP fragments, std::string dir, uint32_t buffer_size) {
    XPtr<FragmentLoader> loader(fragments);
    FragmentIterator iter(*loader);
    FileWriterBuilder wb(dir, buffer_size);    
    StoredFragmentsWriter::createPacked(wb).write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
bool is_compressed_fragments_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    H5ReaderBuilder rb(file, group, buffer_size);
    std::string version = rb.readVersion();
    bool compressed = false;
    if (version == "unpacked-fragments-v1") {
        StoredFragments::openUnpacked(rb);
    } else if (version == "packed-fragments-v1") {
        StoredFragmentsPacked::openPacked(rb);
        compressed = true;
    } else {
        throw std::runtime_error(std::string("Fragments hdf5 has unrecognized version ") + version);
    }
    return compressed;
}

// [[Rcpp::export]]
SEXP iterate_unpacked_fragments_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new StoredFragments(StoredFragments::openUnpacked(rb)))
    );
}

// [[Rcpp::export]]
void write_unpacked_fragments_hdf5_cpp(SEXP fragments, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    XPtr<FragmentLoader> loader(fragments);
    FragmentIterator iter(*loader);

    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    StoredFragmentsWriter::createUnpacked(wb).write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP iterate_packed_fragments_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    H5ReaderBuilder rb(file, group, buffer_size);
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new StoredFragmentsPacked(StoredFragmentsPacked::openPacked(rb)))
    );
}

// [[Rcpp::export]]
void write_packed_fragments_hdf5_cpp(SEXP fragments, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    XPtr<FragmentLoader> loader(fragments);
    FragmentIterator iter(*loader);

    H5WriterBuilder wb(file, group, buffer_size, chunk_size);
    StoredFragmentsWriter::createPacked(wb).write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
int get_bp128_version_cpp() {
    return _SIMDBP128_MODE_;
}

// [[Rcpp::export]]
bool fragments_identical_cpp(SEXP fragments1, SEXP fragments2) {
    XPtr<FragmentLoader> l1(fragments1);
    XPtr<FragmentLoader> l2(fragments2);
    l1->restart(); l2->restart();
    FragmentIterator i1(*l1);
    FragmentIterator i2(*l2);

    while(true) {
        bool res1 = i1.nextChr();
        bool res2 = i2.nextChr();
        if(res1 != res2) return false;
        if (!res1) break;
        if(i1.currentChr() != i2.currentChr()) return false;
        while(true) {
            bool res1 = i1.nextFrag();
            bool res2 = i2.nextFrag();
            if (res1 != res2) return false;
            if (!res1) break;
            if (i1.cell() != i2.cell()) return false;
            if (i1.start() != i2.start()) return false;
            if (i1.end() != i2.end()) return false;
        }
    }
    for (uint32_t i = 0; ;i++) {
        const char* n1 = i1.cellNames(i);
        const char* n2 = i2.cellNames(i);
        if ((n1 == NULL) != (n2 == NULL)) return false;
        if (n1 == NULL) break;
        if (strcmp(n1, n2) != 0) return false;
    }
    for (uint32_t i = 0; ;i++) {
        const char* n1 = i1.chrNames(i);
        const char* n2 = i2.chrNames(i);
        if ((n1 == NULL) != (n2 == NULL)) return false;
        if (n1 == NULL) break;
        if (strcmp(n1, n2) != 0) return false;
    }
    return true;
}