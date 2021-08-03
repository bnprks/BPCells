test_that("BP128 pack/unpack, d1, and FOR round trip properly", {
    expect_true(test_bitpacking_cpp())
})