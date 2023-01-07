generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero = 0.5, max_val = 10) {
  m <- matrix(rbinom(nrow * ncol, 1, fraction_nonzero) * sample.int(max_val, nrow * ncol, replace = TRUE), nrow = nrow)
  as(m, "dgCMatrix")
}

test_that("log1p works", {
    m <- generate_sparse_matrix(20, 10, max_val=1e5)
    res <- m %>%
        as("IterableMatrix") %>%
        log1p() %>%
        as("dgCMatrix")
    expect_equal(log1p(m), res, tolerance=10*2^-23)
})

test_that("expm1 works", {
    m <- generate_sparse_matrix(20, 10, max_val=100) * .01
    res <- m %>%
        as("IterableMatrix") %>%
        expm1() %>%
        as("dgCMatrix")
    expect_equal(expm1(m), res, tolerance=10*2^-23)
})

test_that("scale row/col works", {
    m <- generate_sparse_matrix(20, 10, max_val=1e5)
    m2 <- as(m, "IterableMatrix")
    
    scale_row <- runif(20)
    scale_col <- runif(10)

    expect_equal(multiply_rows(m, scale_row), multiply_rows(m2, scale_row) %>% as("dgCMatrix"))
    expect_equal(multiply_cols(m, scale_col), multiply_cols(m2, scale_col) %>% as("dgCMatrix"))
})

test_that("shift row/col works", {
    m <- generate_sparse_matrix(20, 10, max_val=1e5)
    m2 <- as(m, "IterableMatrix")
    
    shift_row <- runif(20)
    shift_col <- runif(10)

    expect_equal(add_rows(m, shift_row) %>% as("CsparseMatrix"), add_rows(m2, shift_row) %>% as("dgCMatrix"))
    expect_equal(add_cols(m, shift_col) %>% as("CsparseMatrix"), add_cols(m2, shift_col) %>% as("dgCMatrix"))
})

test_that("min_scalar works", {
    m <- generate_sparse_matrix(20, 10, max_val=10)
    m2 <- as(m, "IterableMatrix")
    min_val <- 5
    expect_equal(pmin(m, min_val), min_scalar(m2, min_val) %>% as("dgCMatrix"))
})

test_that("pow works", {
    m <- generate_sparse_matrix(20, 10, max_val=10)
    m2 <- as(m, "IterableMatrix")

    expect_equal(m^2, m2^2 %>% as("dgCMatrix"))
    expect_equal(m^3, m2^3 %>% as("dgCMatrix"))
})