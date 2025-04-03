# Copyright 2025 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

test_that("assert_is_mat works", {
    mat <- matrix(1:4, nrow = 2)
    mat_dgc <- as(mat, "dgCMatrix")
    mat_iterable <- as(mat, "IterableMatrix")
    expect_no_error(assert_is_mat(mat))
    expect_no_error(assert_is_mat(mat_dgc))
    expect_error(assert_is_mat(c(mat_iterable, mat_iterable)))
    expect_error(assert_is_mat("a"))
    expect_error(assert_is_mat(c("a", "a")))
    expect_error(assert_is_mat(1))
})