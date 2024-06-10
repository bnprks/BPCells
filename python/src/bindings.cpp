// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <string>

#include "py_interrupts.hpp"
#include "fragments.hpp"
#include "matrix.hpp"

#include "bpcells-cpp/arrayIO/vector.h"

#include "bpcells-cpp/simd/bp128.h"
#include "bpcells-cpp/simd/current_target.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using namespace BPCells;


PYBIND11_MODULE(cpp, m) {

    m.def("import_10x_fragments", &BPCells::py::import_10x_fragments);
    m.def("cell_names_fragments_dir", &BPCells::py::cell_names_fragments_dir);
    m.def("chr_names_fragments_dir", &BPCells::py::chr_names_fragments_dir);
    m.def("pseudobulk_coverage", &BPCells::py::pseudobulk_coverage);
    m.def("precalculate_pseudobulk_coverage", &BPCells::py::precalculate_pseudobulk_coverage);
    m.def("query_precalculated_pseudobulk_coverage", &BPCells::py::query_precalculated_pseudobulk_coverage);
        
    m.def("write_matrix_dir_from_memory", &BPCells::py::write_matrix_dir_from_memory);
    m.def("write_matrix_dir_from_concat", &BPCells::py::write_matrix_dir_from_concat);
    m.def("write_matrix_dir_from_h5ad", &BPCells::py::write_matrix_dir_from_h5ad);
    
    m.def("load_matrix_dir_subset", &BPCells::py::load_matrix_dir_subset);
    m.def("dims_matrix_dir", &BPCells::py::dims_matrix_dir);

    m.def("load_matrix_dir_to_memory", &BPCells::py::load_matrix_dir_to_memory);
    m.def("load_matrix_memory_subset", &BPCells::py::load_matrix_memory_subset);
    pybind11::class_<VecReaderWriterBuilder, std::shared_ptr<VecReaderWriterBuilder>>(m, "VecReaderWriterBuilder");
    
    m.def("simd_current_target", &BPCells::simd::current_target);
    m.def("simd_current_target_bp128", &BPCells::simd::bp128::current_target);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif

}
