# Copyright 2025 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.


test_that("Getting test data works", {
    mat <- get_demo_mat()
    frags <- get_demo_frags()
    expect_true(is(mat, "IterableMatrix"))
    expect_true(is(frags, "IterableFragments"))
    ## Also make sure data can be built
    expect_no_error(prepare_demo_data(file.path(tools::R_user_dir("BPCells", which = "data"), "demo_data")))
    remove_demo_data()
})