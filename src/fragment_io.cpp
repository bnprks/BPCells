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
SEXP iterate_10x_fragments_cpp(std::string path, std::string comment) {
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

List fragments_get_names(const StoredFragmentsBase &frags) {
    uint32_t chr_count = frags.chrCount();
    uint32_t cell_count = frags.cellCount();
    StringVector chr_names(chr_count);
    StringVector cell_names(cell_count);
    for (uint32_t i = 0; i < chr_count; i++) {
        if (frags.chrNames(i) == NULL) {
            throw std::runtime_error("Error retrieving chrNames from stored fragments");
        }
        chr_names[i] = frags.chrNames(i);
    }
    for (uint32_t i = 0; i < cell_count; i++) {
        if (frags.cellNames(i) == NULL) {
            throw std::runtime_error("Error retrieving cellNames from stored fragments");
        }
        cell_names[i] = frags.cellNames(i);
    }
    return List::create(
        Named("chr_names") = chr_names,
        Named("cell_names") = cell_names
    );
}

List info_fragments_reader_builder(ReaderBuilder &rb) {
    std::string version = rb.readVersion();
    if (version == "unpacked-fragments-v1") {
        List l = fragments_get_names(StoredFragments::openUnpacked(rb));
        l["compressed"] = false;
        return l;
    } else if (version == "packed-fragments-v1") {
        List l = fragments_get_names(StoredFragmentsPacked::openPacked(rb));
        l["compressed"] = true;
        return l;
    } else {
        throw std::runtime_error(std::string("Fragments directory has unrecognized version ") + version);
    }
}

// [[Rcpp::export]]
List info_fragments_file_cpp(std::string dir, uint32_t buffer_size) {
    FileReaderBuilder rb(dir, buffer_size);
    return info_fragments_reader_builder(rb);
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
List info_fragments_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    H5ReaderBuilder rb(file, group, buffer_size);
   return info_fragments_reader_builder(rb);
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
        if(res1 != res2) {
            Rcerr << "Different number of remaining chromosomes." << std::endl;
            return false;
        }
        if (!res1) break;
        if(i1.currentChr() != i2.currentChr()) {
            Rcerr << "Different chromosome ID loaded" << std::endl;
            return false;
        }
        while(true) {
            bool res1 = i1.nextFrag();
            bool res2 = i2.nextFrag();
            if (res1 != res2) {
                Rcerr << "Different number of fragments in chromosome." << std::endl;
                return false;
            }
            if (!res1) break;
            if (i1.cell() != i2.cell() || i1.start() != i2.start() || i1.end() != i2.end()) {
                REprintf("Mismatched fragments: %s(id=%d):%d-%d vs. %s(id=%d):%d-%d\n",
                    i1.cellNames(i1.cell()), i1.cell(), i1.start(), i1.end(),
                    i2.cellNames(i2.cell()), i2.cell(), i2.start(), i2.end()
                );
                return false;
            }
        }
    }
    for (uint32_t i = 0; ;i++) {
        const char* n1 = i1.cellNames(i);
        const char* n2 = i2.cellNames(i);
        if ((n1 == NULL) != (n2 == NULL)) {
            Rcerr << "Mismatched number of cell names" << std::endl;
            return false;
        }
        if (n1 == NULL) break;
        if (strcmp(n1, n2) != 0) {
            Rcerr << "Mismatched cell names" << std::endl;
            return false;
        }
    }
    for (uint32_t i = 0; ;i++) {
        const char* n1 = i1.chrNames(i);
        const char* n2 = i2.chrNames(i);
        if ((n1 == NULL) != (n2 == NULL)) {
            Rcerr << "Mismatched number of chr names" << std::endl;
            return false;
        }
        if (n1 == NULL) break;
        if (strcmp(n1, n2) != 0) {
            Rcerr << "Mismatched  chr names" << std::endl;
            return false;
        }
    }
    return true;
}