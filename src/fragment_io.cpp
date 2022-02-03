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

// [[Rcpp::export]]
int get_bp128_version_cpp() {
    return _SIMDBP128_MODE_;
}