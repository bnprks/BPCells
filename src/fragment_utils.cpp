#include <Rcpp.h>
#include <RcppEigen.h>

#include "R_array_io.h"

#include "fragmentIterators/FragmentIterator.h"
#include "fragmentIterators/ShiftCoords.h"
#include "fragmentIterators/ChrSelect.h"
#include "fragmentIterators/CellSelect.h"
#include "fragmentIterators/MergeFragments.h"
#include "fragmentIterators/RegionSelect.h"
#include "fragmentIterators/Rename.h"
// #include "fragmentIterators/InsertionsIterator2.h"

#include "matrixIterators/PeakMatrix.h"
#include "matrixIterators/TileMatrix.h"
// #include "matrixIterators/PeakMatrix2.h"
// #include "matrixIterators/PeakMatrix3.h"
// #include "matrixIterators/PeakMatrix4.h"
// #include "matrixIterators/PeakMatrix5.h"
// #include "matrixIterators/PeakMatrix6.h"
// #include "matrixIterators/TileMatrix.h"
// #include "matrixIterators/TileMatrix2.h"

using namespace Rcpp;
using namespace BPCells;

// [[Rcpp::export]]
NumericVector scan_fragments_cpp(SEXP fragments) {
    XPtr<FragmentLoader> xptr(fragments);
    
    FragmentIterator iter(*xptr);
    uint64_t start_sum = 0;
    uint64_t end_sum = 0;
    uint64_t cell_sum = 0;
    uint64_t len = 0;

    while (iter.nextChr()) {
        while (iter.nextFrag()) {
            len++;
            start_sum += iter.start();
            end_sum += iter.end();
            cell_sum += iter.cell();
        }
    }
    return {(double) len, (double) start_sum, (double) end_sum, (double) (((start_sum % 104729) + (end_sum % 104729)) % 104729)};
}

// [[Rcpp::export]]
SEXP iterate_peak_matrix_cpp(SEXP fragments,
        std::vector<uint32_t> chr, std::vector<uint32_t> start, std::vector<uint32_t> end,
        StringVector chr_levels) {
    XPtr<FragmentLoader> loader(fragments);
    return Rcpp::wrap(
        XPtr<PeakMatrix>(new PeakMatrix(*loader, chr, start, end, 
            std::make_unique<RcppStringReader>(chr_levels)))
    );
}

// [[Rcpp::export]]
SEXP iterate_tile_matrix_cpp(SEXP fragments,
        std::vector<uint32_t> chr, std::vector<uint32_t> start, std::vector<uint32_t> end, std::vector<uint32_t> width,
        StringVector chr_levels) {
    XPtr<FragmentLoader> loader(fragments);
    return Rcpp::wrap(
        XPtr<TileMatrix>(new TileMatrix(*loader, chr, start, end, width,
            std::make_unique<RcppStringReader>(chr_levels)))
    );
}

//Compute number of fragments like ArchR does for cell stats
// [[Rcpp::export]]
List nucleosome_counts_cpp(SEXP fragments, 
        uint32_t nuc_width=147) {
    FragmentLoader *loader = &(*XPtr<FragmentLoader>(fragments));
    FragmentIterator iter(*loader);

    std::vector<int> subNuc(0);
    std::vector<int> monoNuc(0);
    std::vector<int> multiNuc(0);

    if (iter.cellCount() != -1) {
        subNuc.resize(iter.cellCount(),0);
        monoNuc.resize(iter.cellCount(),0);
        multiNuc.resize(iter.cellCount(),0);
    }

    while (iter.nextChr()) {
        while (iter.nextFrag()) {
            uint32_t cell = iter.cell();
            if (cell >= subNuc.size()) {
                subNuc.resize(cell+1, 0);
                monoNuc.resize(cell+1, 0);
                multiNuc.resize(cell+1, 0);
            }
            uint32_t size = (iter.end() - iter.start())/nuc_width;
            if (size == 0) subNuc[cell]++;
            else if (size == 1) monoNuc[cell]++;
            else multiNuc[cell]++;
        }
    }

    return List::create(
        Named("subNucleosomal") = subNuc,
        Named("monoNucleosomal") = monoNuc,
        Named("multiNucleosomal") = multiNuc
    );
}

// [[Rcpp::export]]
SEXP iterate_shift_cpp(SEXP fragments, int32_t shift_start, int32_t shift_end) {
    XPtr<FragmentLoader> loader(fragments);
    
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new ShiftCoords(*loader, shift_start, shift_end))
    );
} 


// [[Rcpp::export]]
SEXP iterate_chr_index_select_cpp(SEXP fragments, std::vector<uint32_t> chr_selection) {
    XPtr<FragmentLoader> loader(fragments);
    
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new ChrIndexSelect(*loader, chr_selection))
    );
} 

// [[Rcpp::export]]
SEXP iterate_chr_name_select_cpp(SEXP fragments, std::vector<std::string> chr_selection) {
    XPtr<FragmentLoader> loader(fragments);
    
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new ChrNameSelect(*loader, chr_selection))
    );
} 

// [[Rcpp::export]]
SEXP iterate_cell_index_select_cpp(SEXP fragments, std::vector<uint32_t> cell_selection) {
    XPtr<FragmentLoader> loader(fragments);
    
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new CellIndexSelect(*loader, cell_selection))
    );
} 

// [[Rcpp::export]]
SEXP iterate_cell_name_select_cpp(SEXP fragments, std::vector<std::string> cell_selection) {
    XPtr<FragmentLoader> loader(fragments);
    
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new CellNameSelect(*loader, cell_selection))
    );
} 

// [[Rcpp::export]]
SEXP iterate_chr_rename_cpp(SEXP fragments, const StringVector &chr_names) {
    XPtr<FragmentLoader> loader(fragments);
    
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new RenameChrs(*loader, std::make_unique<RcppStringReader>(chr_names)))
    );
} 

// [[Rcpp::export]]
SEXP iterate_cell_rename_cpp(SEXP fragments, const StringVector &cell_names) {
    XPtr<FragmentLoader> loader(fragments);
    
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new RenameCells(*loader, std::make_unique<RcppStringReader>(cell_names)))
    );
} 

// [[Rcpp::export]]
SEXP iterate_region_select_cpp(SEXP fragments,
        std::vector<uint32_t> chr, std::vector<uint32_t> start, std::vector<uint32_t> end,
        StringVector chr_levels, bool invert_selection) {
    XPtr<FragmentLoader> loader(fragments);
    return Rcpp::wrap(
        XPtr<FragmentLoader>(new RegionSelect(*loader, chr, start, end, 
            std::make_unique<RcppStringReader>(chr_levels), invert_selection))
    );
} 

// [[Rcpp::export]]
SEXP iterate_merge_fragments_cpp(SEXP fragments_list) {
    std::vector<FragmentLoader*> fragments_vec;
    List l = fragments_list;
    for (uint32_t i = 0; i < l.size(); i++) {
        XPtr<FragmentLoader> loader(l[i]);
        fragments_vec.push_back(&(*loader));
    }
    
    return Rcpp::wrap(XPtr<FragmentLoader>(new MergeFragments(fragments_vec)));
} 