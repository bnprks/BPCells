
#include <vector>
#include <Rcpp.h>

#include "bitpacking/packing_utils.h"
#include "bitpacking/bp128.h"

using namespace Rcpp;
using namespace BPCells;

IntegerVector uint_vector_to_R(vec_uint32 x) {
    return IntegerVector((int32_t *) &x.data[0], (int32_t *) &x.data[x.size]);
}

vec_uint32 R_to_uint_vector(SEXP s) {
    // This is necessary since otherwise Rcpp converts INTEGER_MIN to 0
    // when casting to uint32_t

    // I suspect this results in a memory copy when compiled down. 
    // This could probably be avoided if I switched my packed_fragments interface
    // to use pointers rather than full-fledged vectors since the type casting
    // would work better. (Or pointers to vectors amenable to reinterpret_cast'ing)
    IntegerVector v1 = as<IntegerVector>(s);
    vec_uint32 v2;
    v2.capacity = v1.size();
    v2.size = v1.size();
    v2.data = (uint32_t *) &v1[0];
    return v2;
}

List packed_frags_to_R(const packed_frags &frags) {    
    return List::create(
        Named("start_data") = uint_vector_to_R(frags.starts.data),
        Named("start_idx") = uint_vector_to_R(frags.starts.idx),
        Named("start_starts") = uint_vector_to_R(frags.starts.starts),

        Named("end_data") = uint_vector_to_R(frags.ends.data),
        Named("end_idx") = uint_vector_to_R(frags.ends.idx),

        Named("cell_id_data") = uint_vector_to_R(frags.cell_ids.data),
        Named("cell_id_idx") = uint_vector_to_R(frags.cell_ids.idx),
        //Named("cell_id_mins") = uint_vector_to_R(frags.cell_ids.mins),

        Named("end_max") = uint_vector_to_R(frags.end_max),
        Named("length") = (int32_t) frags.starts.len
    );
}

packed_frags R_to_packed_frags(List R_frags) {
    packed_frags frags;

    frags.starts.data = R_to_uint_vector(R_frags["start_data"]);
    frags.starts.idx = R_to_uint_vector(R_frags["start_idx"]);
    frags.starts.starts = R_to_uint_vector(R_frags["start_starts"]);
    frags.starts.len = as<uint32_t>(R_frags["length"]);

    frags.ends.data = R_to_uint_vector(R_frags["end_data"]);
    frags.ends.idx = R_to_uint_vector(R_frags["end_idx"]);
    frags.ends.len = frags.starts.len;

    frags.end_max = R_to_uint_vector(R_frags["end_max"]);

    frags.cell_ids.data = R_to_uint_vector(R_frags["cell_id_data"]);
    frags.cell_ids.idx = R_to_uint_vector(R_frags["cell_id_idx"]);
    //frags.cell_ids.mins = R_to_uint_vector(R_frags["cell_id_mins"]);
    frags.cell_ids.len = frags.starts.len;

    return frags;
}
