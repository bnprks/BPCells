#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>
#include <RcppEigen.h>

#include "utils/radix_sort.h"
using namespace Rcpp;


// Given a list of k nearest neighbors for each cell,
// return a sparse matrix of the number of shared nearest neighbors.
// param neghbor_indices: cells x k matrix, where entry (i,j) lists the
//    1-based index of the jth nearest neighbor of cell i
// param min_neighbors: if a pair (i, j) has less than min_neighbors shared, don't output it.
// Return: List of 3 integer vectors: i, j, and snn_count, marking triplets of (cell 1, cell 2, SNN count).
//    and the values in i are guaranteed to be less or equal to the corresponding
//    values in j. Given in order sorted by (i, j)
// [[Rcpp::export]]
SEXP build_snn_graph_cpp(IntegerMatrix neighbor_indices, int min_neighbors) {
    int cells = neighbor_indices.rows();
    int k = neighbor_indices.cols();

    // For each cell, mark a list of which cells mark it as a direct neighbor
    std::vector<std::vector<int>> reverse_neighbors(cells);

    for (int i = 0; i < cells; i++) {
        for (int j = 0; j < k; j++) {
            reverse_neighbors[neighbor_indices(i, j)-1].push_back(i);
        }
    }

    // For each cell, count up the shared neighbors with lower-index cells (so we get lower-triangle matrix out)
    std::vector<std::vector<int>> shared_neighbors(cells);

    for (int i = 0; i < cells; i++) {
        shared_neighbors[i].insert(shared_neighbors[i].end(), k, i);
        for (size_t j = 0; j < reverse_neighbors[i].size(); j++) {
            for (size_t k = j+1; k < reverse_neighbors[i].size(); k++) {
                int c1 = std::min(reverse_neighbors[i][j], reverse_neighbors[i][k]);
                int c2 = std::max(reverse_neighbors[i][j], reverse_neighbors[i][k]);
                shared_neighbors[c2].push_back(c1);
            }
        }
        reverse_neighbors[i].clear();
    }
    reverse_neighbors.clear();

    std::vector<int> buf;
    std::vector<int> c1, c2, snn_count;

    for (int i = 0; i < cells; i++) {
        if (shared_neighbors[i].size() == 0) continue;
        if (buf.size() < shared_neighbors[i].size()) {
            buf.resize(shared_neighbors[i].size());
        }
        size_t s = shared_neighbors[i].size();
        BPCells::lsdRadixSortArrays(shared_neighbors[i].size(), shared_neighbors[i], buf);
        
        int prev = shared_neighbors[i][0];
        int count = 0;
        for (size_t j = 0; j < s; j++) {
            if (prev == shared_neighbors[i][j]) {
                count++;
            } 
            else {
                if (count >= min_neighbors) {
                    c1.push_back(i);
                    c2.push_back(prev);
                    snn_count.push_back(count);
                }
                prev = shared_neighbors[i][j];
                count = 1;
            }
        }
        if (count >= min_neighbors) {
            c1.push_back(i);
            c2.push_back(prev);
            snn_count.push_back(count);
        }
        
        shared_neighbors[i].clear();
    }

    return List::create(
        Named("i") = c1,
        Named("j") = c2,
        Named("snn_count") = snn_count
    );
}