#include <Rcpp.h>

#include "fragmentIterators/FragmentsIterator.h"
#include "fragmentIterators/InsertionsIterator.h"
#include "fragmentIterators/InsertionsIterator2.h"

using namespace Rcpp;
using namespace BPCells;

// [[Rcpp::export]]
NumericVector scan_insertions_cpp(SEXP fragments) {
    XPtr<FragmentsLoader> xptr(fragments);
    InsertionsIterator iter(*xptr);

    uint64_t len = 0;
    uint64_t coord_sum = 0;
    while (iter.nextChr()) {
        uint32_t last_coord = 0;
        while (iter.nextInsertion()) {
            if (iter.coord() < last_coord) stop("Unsorted coord detected");
            last_coord = iter.coord();
            coord_sum += iter.coord();
            len++;
        }
    }
    return {(double) len, (double) coord_sum, (double) ((coord_sum + len/2) % 104729)};
}

// [[Rcpp::export]]
NumericVector scan_insertions2_cpp(SEXP fragments) {
    XPtr<FragmentsLoader> xptr(fragments);
    InsertionsIterator2 iter(*xptr);

    uint64_t len = 0;
    uint64_t coord_sum = 0;
    while (iter.nextChr()) {
        uint32_t last_coord = 0;
        while (iter.nextInsertion()) {
            if (iter.coord() < last_coord) stop("Unsorted coord detected");
            last_coord = iter.coord();
            coord_sum += iter.coord();
            len++;
        }
    }
    return {(double) len, (double) coord_sum, (double) ((coord_sum + len/2) % 104729)};
}

// [[Rcpp::export]]
NumericVector scan_insertions2_verbose_cpp(SEXP fragments) {
    XPtr<FragmentsLoader> xptr(fragments);
    InsertionsIterator2 iter(*xptr);

    uint64_t len = 0;
    uint64_t coord_sum = 0;
    while (iter.nextChr()) {
        uint32_t last_coord = 0;
        while (iter.nextInsertion()) {
            if (iter.coord() < last_coord) stop("Unsorted coord detected");
            last_coord = iter.coord();
            coord_sum += iter.coord();
            len++;
            printf("%s\t%s\t%d\n", iter.chrNames(iter.currentChr()), iter.cellNames(iter.cell()), iter.coord());
        }
    }
    return {(double) len, (double) coord_sum, (double) ((coord_sum + len/2) % 104729)};
}