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

test_that("min_by_row works", {
    m <- generate_sparse_matrix(20, 10, max_val=10)
    m2 <- as(m, "IterableMatrix")
    min_vals <- sample.int(10, size=nrow(m), replace=TRUE)
    expect_equal(pmin(m, min_vals), min_by_row(m2, min_vals) %>% as("dgCMatrix"))
})

test_that("min_by_col works", {
    m <- generate_sparse_matrix(20, 10, max_val=10)
    m2 <- as(m, "IterableMatrix")
    min_vals <- sample.int(10, size=ncol(m), replace=TRUE)
    expect_equal(t(pmin(t(m), min_vals)), min_by_col(m2, min_vals) %>% as("dgCMatrix"))
})

test_that("pow works", {
    m <- generate_sparse_matrix(20, 10, max_val=10)
    m2 <- as(m, "IterableMatrix")

    expect_equal(m^2, m2^2 %>% as("dgCMatrix"))
    expect_equal(m^3, m2^3 %>% as("dgCMatrix"))
})

test_that("sctransform works", {
    m <- generate_sparse_matrix(200, 100, max_val=10)

    mi <- as(m, "IterableMatrix")
    mi_t <- t(as(t(m), "IterableMatrix"))

    # Make dummy model fit parameters
    cell_counts <- Matrix::colSums(m)
    gene_means <- Matrix::rowSums(m) / sum(m)
    gene_theta <- 1 / (runif(nrow(m)) * 3)
    gene_theta[runif(nrow(m)) < 0.1] <- Inf

    clip_range <- c(-10, 10)
    min_var <- 1/25

    # Calculate desired sctransform result
    mu <- tcrossprod(gene_means, cell_counts)
    var <- mu + mu^2 / gene_theta
    var[var < min_var] <- min_var
    ans <- (as.matrix(m) - mu) / sqrt(var)
    ans[ans < clip_range[1]] <- clip_range[1]
    ans[ans > clip_range[2]] <- clip_range[2]

    vec_left <- runif(nrow(m))
    vec_right <- runif(ncol(m))

    for (slow in c(FALSE, TRUE)) {
        for (mat in c(mi, mi_t)) {
            res1 <- sctransform_pearson(
                mat, gene_theta, gene_means, cell_counts, min_var, clip_range, slow=slow
            )
            res2 <- sctransform_pearson(
                t(mat), gene_theta, gene_means, cell_counts, min_var, clip_range, columns_are_cells=FALSE, slow=slow
            ) %>% t()
            
            # Note: I expect SIMD will have better accuracy than 1e-5 typically, but just leaving lots of margin
            # to avoid test failures
            tol <- ifelse(slow, testthat::testthat_tolerance(), 1e-5)
            # Check matrix conversion
            expect_equal(ans, as.matrix(res1), tolerance=tol)
            expect_equal(ans, as.matrix(res2), tolerance=tol)

            # Check vector products
            expect_equal(ans %*% vec_right, res1 %*% vec_right, tolerance=tol)
            expect_equal(ans %*% vec_right, res2 %*% vec_right, tolerance=tol)

            expect_equal(vec_left %*% ans, vec_left %*% res1, tolerance=tol)
            expect_equal(vec_left %*% ans, vec_left %*% res2, tolerance=tol)
        }
    }
})