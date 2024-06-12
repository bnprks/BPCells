// Copyright 2024 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.


#pragma once


namespace BPCells::simd {

/// @brief Helper function for denseMultiplyRight inner loop for sparse matrix A * dense matrix B
///   Given sparse matrix entries from the same column of A, update the output matrix (res)
/// @param row_data Row index for each entry (length count)
/// @param val_data Value for each entry (length count)
/// @param count Number of matrix entries
/// @param res Dense matrix in memory for the result (row-major, entries index [c, c+dim) cover row c)
/// @param B_row Values in the relevant row of B (length dim)
/// @param dim Number of columns in the matrix B
void denseMultiplyRightHelper(const uint32_t *row_data, const double *val_data, uint32_t count, double *res, const double *B_row, uint32_t dim);


/// @brief Helper function for denseMultiplyRight inner loop for sparse matrix A * dense matrix B
///   Given sparse matrix entries from the same column of A, update the appropriate column of the output matrix (res_col)
/// @param row_data Row index for each entry (length count)
/// @param val_data Value for each entry (length count)
/// @param count Number of matrix entries
/// @param res_col Values in the relevant col of res (length dim)
/// @param B Dense matrix in memory holding B (column-major; entries index [c, c+dim) cover column c)
/// @param dim Number of rows in the matrix B
void denseMultiplyLeftHelper(const uint32_t *row_data, const double *val_data, uint32_t count, double *res_col, const double *B, uint32_t dim);


}