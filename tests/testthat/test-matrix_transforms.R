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

test_that("binarize works", {
    v <- c(-1.0, 0, .1, .2, .5, 1.0, 1.5, 2.0)
    m <- matrix(v, nrow=2, byrow=TRUE)
    m2 <- as(as(m, 'dgCMatrix'), 'IterableMatrix')
    m3 <- as(binarize(m2), 'matrix')
    expect_identical(m3, matrix(c(0, 0, 0, 1, 0, 1, 0, 1), nrow=2))
    m4 <- as(binarize(m2, threshold=.5), 'matrix')
    expect_identical(m4, matrix(c(0, 1, 0, 1, 0, 1, 0, 1), nrow=2))
})

test_that("round works", {
    v <- c(0.1, 0.5, 0.7, 1.3, 1.5, 1.9, -0.1, -0.5, -0.7, -1.3, -1.5, -1.9)
    m <- matrix(v, nrow=2, byrow=TRUE)
    vr <- c(0, 0, 1, 1, 2, 2, 0, 0, -1, -1, -2, -2)
    mr <- matrix(vr, nrow=2, byrow=TRUE)
    m2 <- as(as(m, 'dgCMatrix'), 'IterableMatrix')
    m3 <- as(round(m2), 'matrix')
    expect_identical(m3, mr)
})

test_that("sctransform works", {
    withr::local_seed(15123)
    # Generate poisson-distributed data from randomized parameters
    genes <- 201
    cells <- 101
    gene_means <- exp(rnorm(genes))
    gene_means <- gene_means / sum(gene_means)
    cell_counts <- exp(rnorm(cells)) * 1e3

    mat_dense <- matrix(rpois(genes*cells, tcrossprod(gene_means, cell_counts)), nrow=genes)
    m <- as(mat_dense, "dgCMatrix")

    mi <- as(m, "IterableMatrix")
    mi_t <- t(as(t(m), "IterableMatrix"))

    gene_theta <- 1 / (runif(nrow(m)) * 3)
    gene_theta[runif(nrow(m)) < 0.1] <- Inf

    # Set a narrow clip range so we actually hit it on this test
    clip_range <- c(-2, 2)
    min_var <- 1/25

    # Calculate desired sctransform result
    mu <- tcrossprod(gene_means, cell_counts)
    var <- mu + mu^2 / gene_theta
    testthat::expect_true(any(var < min_var))
    var[var < min_var] <- min_var
    ans <- (as.matrix(m) - mu) / sqrt(var)

    # Make sure clipping code will get exercised
    testthat::expect_true(any(ans < min(clip_range)))
    testthat::expect_true(any(ans > max(clip_range)))
    # This condition is the trickiest to cover with clipping -- where
    # we'll need to consider both the top and bottom clip range in a single entry
    testthat::expect_true(any(-mu/sqrt(var) < min(clip_range) & ans > max(clip_range)))
    
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
