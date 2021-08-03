#include <Rcpp.h>

#include "fragmentIterators/FragmentsIterator.h"
#include "fragmentIterators/BedFragments.h"
#include "fragmentIterators/PackedFragments.h"
#include "fragmentIterators/UnpackedFragments.h"

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
