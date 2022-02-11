#include <Rcpp.h>

#include "fragmentIterators/FragmentsIterator.h"
#include "fragmentIterators/BedFragments.h"
#include "fragmentIterators/PackedFragments.h"
#include "fragmentIterators/UnpackedFragments.h"

#include "fragmentIterators/UnpackedFragments2.h"
#include "fragmentIterators/PackedFragments2.h"

#include "arrayIO/binaryfile.h"
#include "arrayIO/hdf5.h"
#include "arrayIO/vector.h"

using namespace Rcpp;
using namespace BPCells;

packed_frags R_to_packed_frags(List R_frags);
List packed_frags_to_R(const packed_frags &frags);
vec_uint32 R_to_uint_vector(SEXP s);

// [[Rcpp::export]]
SEXP load_10x_fragments_cpp(std::string path, std::string comment) {
    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new BedFragments(path.c_str(), comment.c_str()))
    );
}


// [[Rcpp::export]]
void write_10x_fragments_cpp(std::string path, SEXP fragments, bool append_5th_column=false) {
    XPtr<FragmentsLoader> loader(fragments);
    FragmentsIterator iter(*loader);
    BedFragmentsWriter writer(path.c_str(), append_5th_column);
    writer.write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP load_packed_fragments_cpp(S4 s4) {
    List R_frags_list = s4.slot("fragments");
    std::vector<std::string> chr_names = s4.slot("chr_names");
    std::vector<std::string> cell_names = s4.slot("cell_names");

    std::vector<packed_frags> frags_vec(R_frags_list.size());
    for(int i = 0; i < R_frags_list.size(); i++) 
        frags_vec[i] = R_to_packed_frags(R_frags_list[i]);

    return Rcpp::wrap( 
        XPtr<FragmentsLoader>(new PackedFragments(frags_vec, cell_names, chr_names))
    );
}

// [[Rcpp::export]]
List write_packed_fragments_cpp(SEXP fragments) {
    XPtr<FragmentsLoader> loader(fragments);
    FragmentsIterator iter(*loader);

    PackedFragmentsWriter writer;
    writer.write(iter, &Rcpp::checkUserInterrupt);
    std::vector<packed_frags> frags = writer.getPackedFrags();

    List fragList(frags.size());
    for (int i = 0; i < frags.size(); i++) {
        fragList[i] = packed_frags_to_R(frags[i]);
    }

    return List::create(
        Named("packed_frags") = fragList,
        Named("chr_names") = writer.getChrNames(),
        Named("cell_names") = writer.getCellNames()
    );
}


// [[Rcpp::export]]
SEXP iterate_unpacked_fragments_cpp(S4 s4) {
    List R_frags_list = s4.slot("fragments");
    std::vector<std::string> chr_names = s4.slot("chr_names");
    std::vector<std::string> cell_names = s4.slot("cell_names");

    std::vector<vec_uint32_t> start, end, cell_id;

    for(int i = 0; i < R_frags_list.size(); i++) {
        List chr_list = R_frags_list[i];
        start.push_back(R_to_uint_vector(chr_list["start"]));
        end.push_back(R_to_uint_vector(chr_list["end"]));
        cell_id.push_back(R_to_uint_vector(chr_list["cell_id"]));
    }
    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new UnpackedFragments(start, end, cell_id, cell_names, chr_names))
    );
}

// [[Rcpp::export]]
List write_unpacked_fragments_cpp(SEXP fragments) {
    XPtr<FragmentsLoader> loader(fragments);
    UnpackedFragmentsWriter writer;

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);

    const std::vector<UnpackedFragmentsWriter::unpacked_frags> *frags = writer.getFrags();
    List fragList(frags->size());
    for (int i = 0; i < frags->size(); i++) {
        fragList[i] = List::create(
            Named("start") = IntegerVector((*frags)[i].start_data.begin(), (*frags)[i].start_data.end()),
            Named("end") = IntegerVector((*frags)[i].end_data.begin(), (*frags)[i].end_data.end()),
            Named("cell_id") = IntegerVector((*frags)[i].cell_data.begin(), (*frags)[i].cell_data.end())
        );
    }
    
    return List::create(
        Named("fragments") = fragList,
        Named("chr_names") = writer.getChrNames(),
        Named("cell_names") = writer.getCellNames()
    );
}

// ######################################## UNPACKED FRAGMENTS ########################################


// [[Rcpp::export]]
SEXP iterate_unpacked_fragments2_cpp(S4 s4) {
    List R_frags_list = s4.slot("fragments");
    
    std::vector<std::string> chr_names = s4.slot("chr_names");
    std::vector<std::string> cell_names = s4.slot("cell_names");

    VecUnpackedFragmentsLoader::Storage storage;
    
    for(int i = 0; i < R_frags_list.size(); i++) {
        List chr_frags = R_frags_list[i];
        IntegerVector start = chr_frags["start"];
        IntegerVector cell_id = chr_frags["cell_id"];
        IntegerVector end = chr_frags["end"];
        IntegerVector end_max = chr_frags["end_max"];
        storage.fragments.push_back(UnpackedFrags<VecUIntReader>{
            VecUIntReader((const uint32_t *) start.cbegin(),start.size()),
            VecUIntReader((const uint32_t *) end.cbegin(),end.size()),
            VecUIntReader((const uint32_t *) cell_id.cbegin(),cell_id.size()),
            VecUIntReader((const uint32_t *) end_max.cbegin(),end_max.size())
        });
    }

    storage.chr_names = chr_names;
    storage.cell_names = cell_names;

    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new UnpackedFragments2<VecUnpackedFragmentsLoader>(VecUnpackedFragmentsLoader(storage)))
    );
}

// [[Rcpp::export]]
List write_unpacked_fragments2_cpp(SEXP fragments) {
    XPtr<FragmentsLoader> loader(fragments);
    VecUnpackedFragmentsSaver::Storage storage;
    
    UnpackedFragmentsWriter2<VecUnpackedFragmentsSaver> writer((VecUnpackedFragmentsSaver(storage)));

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);
    auto &frags = storage.fragments;

    List fragList(frags.size());
    for (int i = 0; i < frags.size(); i++) {
        fragList[i] = List::create(
            Named("start") = IntegerVector(frags[i].start.begin(), frags[i].start.end()),
            Named("cell_id") = IntegerVector(frags[i].cell.begin(), frags[i].cell.end()),
            Named("end") = IntegerVector(frags[i].end.begin(), frags[i].end.end()),
            Named("end_max") = IntegerVector(frags[i].end_max.begin(), frags[i].end_max.end())
        );
    }

    return List::create(
        Named("fragments") = fragList,
        Named("chr_names") = storage.chr_names,
        Named("cell_names") = storage.cell_names
    );
}

// [[Rcpp::export]]
SEXP iterate_unpacked_fragments_file_cpp(std::string dir) {
    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new UnpackedFragments2<FileFragmentsLoader>(FileFragmentsLoader(dir, false)))
    );
}

// [[Rcpp::export]]
void write_unpacked_fragments_file_cpp(SEXP fragments, std::string dir) {
    XPtr<FragmentsLoader> loader(fragments);
    
    UnpackedFragmentsWriter2<FileFragmentsSaver> writer((FileFragmentsSaver(dir)));

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP iterate_unpacked_fragments_hdf5_cpp(std::string file, std::string group) {
    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new UnpackedFragments2<H5FragmentsLoader>(H5FragmentsLoader(file, group, false)))
    );
}

// [[Rcpp::export]]
void write_unpacked_fragments_hdf5_cpp(SEXP fragments, std::string file, std::string group, uint32_t chunk_size) {
    XPtr<FragmentsLoader> loader(fragments);
    
    UnpackedFragmentsWriter2<H5FragmentsSaver> writer(H5FragmentsSaver(file, group, chunk_size));

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);
}

// ######################################## PACKED FRAGMENTS ########################################

// [[Rcpp::export]]
SEXP iterate_packed_fragments2_cpp(S4 s4) {
    List R_frags_list = s4.slot("fragments");
    
    std::vector<std::string> chr_names = s4.slot("chr_names");
    std::vector<std::string> cell_names = s4.slot("cell_names");

    VecPackedFragmentsLoader::Storage storage;

    for(int i = 0; i < R_frags_list.size(); i++) {
        List chr_frags = R_frags_list[i];
        IntegerVector start_data = chr_frags["start_data"];
        IntegerVector start_idx = chr_frags["start_idx"];
        IntegerVector start_starts = chr_frags["start_starts"];
        IntegerVector end_data = chr_frags["end_data"];
        IntegerVector end_idx = chr_frags["end_idx"];
        IntegerVector end_max = chr_frags["end_max"];
        IntegerVector cell_data = chr_frags["cell_data"];
        IntegerVector cell_idx = chr_frags["cell_idx"];
        IntegerVector count = chr_frags["count"];
        storage.fragments.push_back(PackedFrags<VecUIntReader>{
            VecUIntReader((const uint32_t *) start_data.cbegin(), start_data.size()),
            VecUIntReader((const uint32_t *) start_idx.cbegin(), start_idx.size()),
            VecUIntReader((const uint32_t *) start_starts.cbegin(), start_starts.size()),
            VecUIntReader((const uint32_t *) end_data.cbegin(), end_data.size()),
            VecUIntReader((const uint32_t *) end_idx.cbegin(), end_idx.size()),
            VecUIntReader((const uint32_t *) end_max.cbegin(), end_max.size()),
            VecUIntReader((const uint32_t *) cell_data.cbegin(), cell_data.size()),
            VecUIntReader((const uint32_t *) cell_idx.cbegin(), cell_idx.size()),
            VecUIntReader((const uint32_t *) count.cbegin(), count.size())
        });
    }

    storage.chr_names = chr_names;
    storage.cell_names = cell_names;

    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new PackedFragments2<VecPackedFragmentsLoader>(VecPackedFragmentsLoader(storage)))
    );
}

// [[Rcpp::export]]
List write_packed_fragments2_cpp(SEXP fragments) {
    XPtr<FragmentsLoader> loader(fragments);
    VecPackedFragmentsSaver::Storage storage;
    
    PackedFragmentsWriter2<VecPackedFragmentsSaver> writer((VecPackedFragmentsSaver(storage)));

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);
    auto &frags = storage.fragments;

    List fragList(frags.size());
    for (int i = 0; i < frags.size(); i++) {
        fragList[i] = List::create(
            Named("start_data") = IntegerVector(frags[i].start_data.begin(), frags[i].start_data.end()),
            Named("start_idx") = IntegerVector(frags[i].start_idx.begin(), frags[i].start_idx.end()),
            Named("start_starts") = IntegerVector(frags[i].start_starts.begin(), frags[i].start_starts.end()),
            Named("end_data") = IntegerVector(frags[i].end_data.begin(), frags[i].end_data.end()),
            Named("end_idx") = IntegerVector(frags[i].end_idx.begin(), frags[i].end_idx.end()),
            Named("end_max") = IntegerVector(frags[i].end_max.begin(), frags[i].end_max.end()),
            Named("cell_data") = IntegerVector(frags[i].cell_data.begin(), frags[i].cell_data.end()),
            Named("cell_idx") = IntegerVector(frags[i].cell_idx.begin(), frags[i].cell_idx.end()),
            Named("count") = IntegerVector(frags[i].count.begin(), frags[i].count.end())
        );
    }

    return List::create(
        Named("fragments") = fragList,
        Named("chr_names") = storage.chr_names,
        Named("cell_names") = storage.cell_names
    );
}


// [[Rcpp::export]]
SEXP iterate_packed_fragments_file_cpp(std::string dir) {
    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new PackedFragments2<FileFragmentsLoader>(FileFragmentsLoader(dir, true)))
    );
}

// [[Rcpp::export]]
void write_packed_fragments_file_cpp(SEXP fragments, std::string dir) {
    XPtr<FragmentsLoader> loader(fragments);
    
    PackedFragmentsWriter2<FileFragmentsSaver> writer((FileFragmentsSaver(dir)));

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP iterate_packed_fragments_hdf5_cpp(std::string file, std::string group) {
    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new PackedFragments2<H5FragmentsLoader>(H5FragmentsLoader(file, group, true)))
    );
}

// [[Rcpp::export]]
void write_packed_fragments_hdf5_cpp(SEXP fragments, std::string file, std::string group, uint32_t chunk_size) {
    XPtr<FragmentsLoader> loader(fragments);
    
    PackedFragmentsWriter2<H5FragmentsSaver> writer(H5FragmentsSaver(file, group, chunk_size));

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);
}

// ######################################## UNPACKED FRAGMENTS V2 ########################################

class RcppStringReader : public StringReader {
private:
    StringVector data;
public:
    RcppStringReader(StringVector data) : data(data) {}
    const char* get(uint32_t idx) const override {
        if (idx < data.size()) return data[idx];
        return NULL;
    }
    uint32_t size() const override {return data.size();} 
};

// [[Rcpp::export]]
SEXP iterate_unpacked_fragments3_cpp(S4 s4) {    
    IntegerVector cell = s4.slot("cell");
    IntegerVector start = s4.slot("start");
    IntegerVector end = s4.slot("end");
    IntegerVector end_max = s4.slot("end_max");
    IntegerVector chr_ptr = s4.slot("chr_ptr");
    StringVector chr_names = s4.slot("chr_names");
    StringVector cell_names = s4.slot("cell_names");
    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new UnpackedFragments3(
            std::make_unique<ZVecUIntReader>((const uint32_t *) cell.cbegin(), cell.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) start.cbegin(), start.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) end.cbegin(), end.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) end_max.cbegin(), end_max.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) chr_ptr.cbegin(), chr_ptr.size()),
            std::make_unique<RcppStringReader>(chr_names),
            std::make_unique<RcppStringReader>(cell_names)
        ))
    );
}

// [[Rcpp::export]]
List write_unpacked_fragments3_cpp(SEXP fragments) {
    XPtr<FragmentsLoader> loader(fragments);
    
    std::vector<uint32_t> cell, start, end, end_max, chr_ptr;
    std::vector<std::string> chr_names, cell_names;
    
    UnpackedFragmentsWriter3 writer(
        std::make_unique<ZVecUIntWriter>(cell),
        std::make_unique<ZVecUIntWriter>(start),
        std::make_unique<ZVecUIntWriter>(end),
        std::make_unique<ZVecUIntWriter>(end_max),
        std::make_unique<ZVecUIntWriter>(chr_ptr),
        std::make_unique<VecStringWriter>(chr_names),
        std::make_unique<VecStringWriter>(cell_names)
    );

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);
    
    return List::create(
        Named("cell") = IntegerVector(cell.begin(), cell.end()),
        Named("start") = IntegerVector(start.begin(), start.end()),
        Named("end") = IntegerVector(end.begin(), end.end()),
        Named("end_max") = IntegerVector(end_max.begin(), end_max.end()),
        Named("chr_ptr") = IntegerVector(chr_ptr.begin(), chr_ptr.end()),
        Named("chr_names") = StringVector(chr_names.begin(), chr_names.end()),
        Named("cell_names") = StringVector(cell_names.begin(), cell_names.end())
    );
}


// [[Rcpp::export]]
SEXP iterate_packed_fragments3_cpp(S4 s4) {   
    IntegerVector cell_data = s4.slot("cell_data");
    IntegerVector cell_idx = s4.slot("cell_idx");
    IntegerVector start_data = s4.slot("start_data");
    IntegerVector start_idx = s4.slot("start_idx");
    IntegerVector start_starts = s4.slot("start_starts");
    IntegerVector end_data = s4.slot("end_data");
    IntegerVector end_idx = s4.slot("end_idx");
    IntegerVector end_max = s4.slot("end_max");
    IntegerVector chr_ptr = s4.slot("chr_ptr");
    
    StringVector chr_names = s4.slot("chr_names");
    StringVector cell_names = s4.slot("cell_names");
    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new PackedFragments3(
            std::make_unique<ZVecUIntReader>((const uint32_t *) cell_data.cbegin(), cell_data.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) cell_idx.cbegin(), cell_idx.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) start_data.cbegin(), start_data.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) start_idx.cbegin(), start_idx.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) start_starts.cbegin(), start_starts.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) end_data.cbegin(), end_data.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) end_idx.cbegin(), end_idx.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) end_max.cbegin(), end_max.size()),
            std::make_unique<ZVecUIntReader>((const uint32_t *) chr_ptr.cbegin(), chr_ptr.size()),

            std::make_unique<RcppStringReader>(chr_names),
            std::make_unique<RcppStringReader>(cell_names)
        ))
    );
}

// [[Rcpp::export]]
List write_packed_fragments3_cpp(SEXP fragments) {
    XPtr<FragmentsLoader> loader(fragments);
    
    std::vector<uint32_t> cell_data, cell_idx, start_data, start_idx, start_starts, 
        end_data, end_idx, end_max, chr_ptr;
    std::vector<std::string> chr_names, cell_names;
    
    PackedFragmentsWriter3 writer(
        std::make_unique<ZVecUIntWriter>(cell_data),
        std::make_unique<ZVecUIntWriter>(cell_idx),
        std::make_unique<ZVecUIntWriter>(start_data),
        std::make_unique<ZVecUIntWriter>(start_idx),
        std::make_unique<ZVecUIntWriter>(start_starts),
        std::make_unique<ZVecUIntWriter>(end_data),
        std::make_unique<ZVecUIntWriter>(end_idx),
        std::make_unique<ZVecUIntWriter>(end_max),
        std::make_unique<ZVecUIntWriter>(chr_ptr),
        
        std::make_unique<VecStringWriter>(chr_names),
        std::make_unique<VecStringWriter>(cell_names)
    );

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);
    
    return List::create(
        Named("cell_data") = IntegerVector(cell_data.begin(), cell_data.end()),
        Named("cell_idx") = IntegerVector(cell_idx.begin(), cell_idx.end()),
        Named("start_data") = IntegerVector(start_data.begin(), start_data.end()),
        Named("start_idx") = IntegerVector(start_idx.begin(), start_idx.end()),
        Named("start_starts") = IntegerVector(start_starts.begin(), start_starts.end()),
        Named("end_data") = IntegerVector(end_data.begin(), end_data.end()),
        Named("end_idx") = IntegerVector(end_idx.begin(), end_idx.end()),
        Named("end_max") = IntegerVector(end_max.begin(), end_max.end()),
        Named("chr_ptr") = IntegerVector(chr_ptr.begin(), chr_ptr.end()),
        
        Named("chr_names") = StringVector(chr_names.begin(), chr_names.end()),
        Named("cell_names") = StringVector(cell_names.begin(), cell_names.end())
    );
}


// [[Rcpp::export]]
bool is_compressed_fragments_file_cpp(std::string dir, uint32_t buffer_size) {
    std::string version(loadVersionDir(dir));
    bool compressed = false;
    if (version == "v1-unpacked-fragments") {
        UnpackedFragments3 frags(openUnpackedFragmentsDir(dir, buffer_size));
    } else if (version == "v1-packed-fragments") {
        PackedFragments3 frags(openPackedFragmentsDir(dir, buffer_size));
        compressed = true;
    } else {
        throw std::runtime_error(std::string("Fragments directory has unrecognized version ") + version);
    }
    return compressed;
}
// [[Rcpp::export]]
SEXP iterate_unpacked_fragments_file2_cpp(std::string dir, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new UnpackedFragments3(openUnpackedFragmentsDir(dir, buffer_size)))
    );
}

// [[Rcpp::export]]
void write_unpacked_fragments_file2_cpp(SEXP fragments, std::string dir, uint32_t buffer_size) {
    XPtr<FragmentsLoader> loader(fragments);
    
    UnpackedFragmentsWriter3 writer(((createUnpackedFragmentsDir(dir, buffer_size))));

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP iterate_packed_fragments_file2_cpp(std::string dir, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new PackedFragments3(openPackedFragmentsDir(dir, buffer_size)))
    );
}

// [[Rcpp::export]]
void write_packed_fragments_file2_cpp(SEXP fragments, std::string dir, uint32_t buffer_size) {
    XPtr<FragmentsLoader> loader(fragments);
    
    PackedFragmentsWriter3 writer(((createPackedFragmentsDir(dir, buffer_size))));

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);
}


// [[Rcpp::export]]
bool is_compressed_fragments_hdf5_cpp(std::string file, std::string group, uint32_t buffer_size) {
    std::string version(loadVersionH5(file, group));
    bool compressed = false;
    if (version == "v1-unpacked-fragments") {
        UnpackedFragments3 frags(openUnpackedFragmentsH5(file, group, buffer_size));
    } else if (version == "v1-packed-fragments") {
        PackedFragments3 frags(openPackedFragmentsH5(file, group, buffer_size));
        compressed = true;
    } else {
        throw std::runtime_error(std::string("Fragments hdf5 has unrecognized version ") + version);
    }
    return compressed;
}
// [[Rcpp::export]]
SEXP iterate_unpacked_fragments_hdf52_cpp(std::string file, std::string group, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new UnpackedFragments3(openUnpackedFragmentsH5(file, group, buffer_size)))
    );
}

// [[Rcpp::export]]
void write_unpacked_fragments_hdf52_cpp(SEXP fragments, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    XPtr<FragmentsLoader> loader(fragments);
    
    UnpackedFragmentsWriter3 writer(createUnpackedFragmentsH5(file, group, buffer_size, chunk_size));

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
SEXP iterate_packed_fragments_hdf52_cpp(std::string file, std::string group, uint32_t buffer_size) {
    return Rcpp::wrap(
        XPtr<FragmentsLoader>(new PackedFragments3(openPackedFragmentsH5(file, group, buffer_size)))
    );
}

// [[Rcpp::export]]
void write_packed_fragments_hdf52_cpp(SEXP fragments, std::string file, std::string group, uint32_t buffer_size, uint32_t chunk_size) {
    XPtr<FragmentsLoader> loader(fragments);
    
    PackedFragmentsWriter3 writer(createPackedFragmentsH5(file, group, buffer_size, chunk_size));

    FragmentsIterator iter(*loader);

    writer.write(iter, &Rcpp::checkUserInterrupt);
}

// [[Rcpp::export]]
int get_bp128_version_cpp() {
    return _SIMDBP128_MODE_;
}