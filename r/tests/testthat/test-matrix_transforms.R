# Copyright 2023 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero = 0.5, max_val = 10) {
  m <- matrix(rbinom(nrow * ncol, 1, fraction_nonzero) * sample.int(max_val, nrow * ncol, replace = TRUE), nrow = nrow)
  as(m, "dgCMatrix")
}

test_that("Subsetting transform works", {
    set.seed(125124)
    m1 <- generate_sparse_matrix(20, 10, max_val=1e5)
    m2 <- as(m1, "IterableMatrix")
    m1 <- as.matrix(m1)

    # Apply some transformation layers
    y <- runif(nrow(m1))
    m1 <- multiply_rows(m1, y)
    m2 <- multiply_rows(m2, y)

    y <- runif(ncol(m1))
    m1 <- multiply_cols(m1, y)
    m2 <- multiply_cols(m2, y)

    m1 <- log1p(m1) + 1
    m2 <- log1p_slow(m2) + 1

    expect_equal(as.matrix(m2), m1)
    
    i <- sample.int(nrow(m1))
    j <- sample.int(ncol(m2))
    x <- 8
    expect_equal(as.matrix(m2[,c(7,8)]), m1[,c(7,8),drop=FALSE])


    expect_equal(as.matrix(m2[i,j]), m1[i,j])
    expect_equal(as.matrix(m2[i[1:5],j[1:6]]), m1[i[1:5],j[1:6]])
    expect_equal(as.matrix(m2[sort(i[1:5]),sort(j[1:6])]), m1[sort(i[1:5]),sort(j[1:6])])
})

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

test_that("square works", {
    m <- generate_sparse_matrix(20, 10, max_val=1e5)
    res <- m %>%
        as("IterableMatrix") %>%
        {.^2} %>%
        as("dgCMatrix")
    expect_equal(m^2, res, tolerance=10*2^-23)
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
    expect_error(m2^-2)
})

test_that("binarize works", {
    v <- c(-1.0, 0, .1, .2, .5, 1.0, 1.5, 2.0)
    m <- matrix(v, nrow=2, byrow=TRUE)
    m2 <- as(as(m, 'dgCMatrix'), 'IterableMatrix')

    m3 <- as(binarize(m2), 'matrix')
    ans3 <- matrix(as.integer(c(0, 1, 0, 1, 1, 1, 1, 1)), nrow=2)
    expect_identical(m3, ans3)
    expect_identical(as(m2 > 0L, 'matrix'), ans3)
    expect_identical(as(0L < m2, 'matrix'), ans3)

    m4 <- as(binarize(m2, threshold=.5), 'matrix')
    ans4 <- matrix(as.integer(c(0, 0, 0, 1, 0, 1, 0, 1)), nrow=2)
    expect_identical(m4, ans4)
    expect_identical(as(m2 > 0.5, 'matrix'), ans4)
    expect_identical(as(0.5 < m2, 'matrix'), ans4)

    m5 <- as(binarize(m2, threshold=.5, strict_inequality=FALSE), 'matrix')
    ans5 <- matrix(as.integer(c(0, 1, 0, 1, 0, 1, 0, 1)), nrow=2)
    expect_identical(m5, ans5)
    expect_identical(as(m2 >= 0.5, 'matrix'), ans5)
    expect_identical(as(0.5 <= m2, 'matrix'), ans5)

    short_description(m2 > 0.5)
})

test_that("Issue 43 regression (preserve colnames when cancelling type conversion)",  {
    m <- matrix(1:12, nrow=3) |> as("dgCMatrix")
    rownames(m) <- paste0("row", seq_len(nrow(m)))
    colnames(m) <- paste0("col", seq_len(ncol(m)))

    m2 <- t(m) |> as("IterableMatrix")
    rownames(m2) <- paste0("row", seq_len(nrow(m2)))
    colnames(m2) <- paste0("col", seq_len(ncol(m2)))

    res <- m %*% (m2 >= 1)
    expect_identical(rownames(res), rownames(m))
    expect_identical(colnames(res), colnames(m2))
})

test_that("Multiply cols of transposed TransformScaleShift works", {
    m <- matrix(1:12, nrow=3) |> as("dgCMatrix") |> as("IterableMatrix") |> t()

    res <- (m - seq_len(nrow(m)))/seq_len(nrow(m))
    res <- multiply_cols(res, seq_len(ncol(m)))

    ans <- (t(matrix(1:12, nrow=3)) - seq_len(nrow(m)))/seq_len(nrow(m))
    ans <- multiply_cols(ans, seq_len(ncol(m)))
    expect_equal(as.matrix(res), as.matrix(ans))
})

test_that("Complicated TransformScaleShift works", {
    # Idea: randomly generate a series of scales + shifts, 
    # checking that the whole sequence returns correct results
    # Don't set a seed so if there's a rare bug we have a better chance of
    # hitting it eventually

    # Run this test 5 times to improve chances of catching something weird.
    # We can't just extend the random operations indefinitely without running
    # into precision issues.
    for (j in seq_len(5)) {
        m <- matrix(1:12, nrow=3)
        bp <- as(m, "dgCMatrix") |> as("IterableMatrix")
        bp_t <- as(t(m), "dgCMatrix") |> as("IterableMatrix") |> t()

        # Do 30 random operations and check at the end
        for (i in seq_len(30)) {
            axis <- sample(c("row", "col", "global"), 1)
            op <- sample(c(`+`, `*`), 1)[[1]]
            bp_prev <- bp
            bp_t_prev <- bp_t
            if (axis == "row") {
                y <- sample(c(-2, -1, 1, 2), nrow(m))
                m <- op(m, y)
                bp <- op(bp, y)
                bp_t <- op(bp_t, y)
            } else if (axis == "col") {
                y <- sample(c(-2, -1, 1, 2), ncol(m))
                m <- t(op(t(m), y))
                bp <- t(op(t(bp), y))
                bp_t <- t(op(t(bp_t), y))
            } else {
                y <- sample(c(-2, -1, 1, 2), 1)
                m <- op(m, y)
                bp <- op(bp, y)
                bp_t <- op(bp_t, y)
            }
        }
        expect_identical(m, as.matrix(bp))
        expect_identical(m, as.matrix(bp_t))
    }
    
})

test_that("round works", {
    m <- generate_sparse_matrix(20, 10, max_val=1e5) / 70
    digits <- 0
    res <- m %>%
        as("IterableMatrix") %>%
        round(digits=digits) %>%
        as("dgCMatrix")
    expect_equal(round(m, digits=digits), res, tolerance=1.49e-08)

    digits <- 2
    res <- m %>%
        as("IterableMatrix") %>%
        round(digits=digits) %>%
        as("dgCMatrix")
    expect_equal(round(m, digits=digits), res, tolerance=1.49e-08)
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

test_that("linear regression works", {
    nrow <- 5
    ncol <- 30

    m0 <- generate_sparse_matrix(nrow, ncol, max_val = 10)
    m <- as(m0, "IterableMatrix")
    mt <- t(as(t(m0), "IterableMatrix"))
    latent_data <- data.frame(
        nUMI = colSums(m),
        group = as.factor(sample(1:5, ncol, replace = TRUE)),
        random = runif(ncol),
        gene1 = m0[3, ]
    )
    # Calculate residuals manually
    # This is less efficient than calculating the QR once and using qr.resid, but
    # it's a bit more foolproof for our reference answer
    ans <- matrix(nrow = nrow, ncol = ncol)
    for (i in seq_len(nrow)) {
        regression_data <- cbind(latent_data, m0[i, ])
        colnames(regression_data) <- c(colnames(latent_data), "y")
        fmla <- fmla <- as.formula(paste("y ~", paste(colnames(latent_data), collapse="+")))
        ans[i,] <- lm(fmla, regression_data)$residuals
    }
    m1 <- regress_out(m, latent_data = latent_data)
    m1t <- regress_out(mt, latent_data = latent_data)
    expect_equal(as(m1, "matrix"), ans)
    expect_equal(as(m1t, "matrix"), ans)
    
    # Also test when the predicted variables are in the columns.
    latent_data <- data.frame(
        nUMI = rowSums(m),
        group = as.factor(sample(1:5, nrow, replace = TRUE)),
        random = runif(nrow),
        cell1 = m0[, 4]
    )
    ans <- matrix(nrow = nrow, ncol = ncol)
    for (i in seq_len(ncol)) {
        regression_data <- cbind(latent_data, m0[, i])
        colnames(regression_data) <- c(colnames(latent_data), "y")
        fmla <- fmla <- as.formula(paste("y ~", paste(colnames(latent_data), collapse="+")))
        ans[, i] <- lm(fmla, regression_data)$residuals
    }
    m1 <- regress_out(m, latent_data = latent_data, prediction_axis = "col")
    m1t <- regress_out(mt, latent_data = latent_data, prediction_axis = "col")
    expect_equal(as(m1, "matrix"), ans)
    expect_equal(as(m1t, "matrix"), ans)
})

test_that("tf-idf normalization works", {
    m <- generate_sparse_matrix(5, 5)
    rownames(m) <- paste0("row", seq_len(nrow(m)))
    rev_rownames <- rev(rownames(m))
    # Create tf-idf normalization for dgCMatrix
    res_dgc <- log1p((diag(1/rowMeans(m)) %*% (m %*% diag(1/colSums(m))) %>% as("dgCMatrix")) * 1e4)
    
    rownames(res_dgc) <- rownames(m)
    m2 <- as(m, "IterableMatrix")
    # Check that we can pass in row means as a (named) vector
    row_means <- matrix_stats(m2, row_stats = c("mean"))$row_stats["mean",]
    # Test that row means ordering does not matter as long as names exist
    row_means_shuffled <- row_means[sample(1:length(row_means))]
    # Test that row means can have an extra element as long as all rownames are in the vector
    row_means_plus_one <- c(row_means, row6 = 1)
    # Check when row means has no row names
    row_means_no_names <- unname(row_means)
    

    res <- normalize_tfidf(m2)
    expect_equal(res %>% as("dgCMatrix"), res_dgc, tolerance = 1e-6)
    res_with_row_means <- normalize_tfidf(m2, feature_means = row_means)
    res_with_row_means_partial <- normalize_tfidf(feature_means = row_means)(m2)
    expect_equal(res_with_row_means, res_with_row_means_partial)
    expect_identical(res, res_with_row_means)

    res_with_shuffled_row_means <- normalize_tfidf(m2, feature_means = row_means_shuffled)
    expect_identical(res_with_row_means, res_with_shuffled_row_means, res)

    res_with_row_means_with_extra_element <- normalize_tfidf(m2, feature_means = row_means_plus_one)
    expect_identical(res, res_with_row_means_with_extra_element)
    
    # Check cases where names exists in either row means xor rownames(mat)
    # check where both don't have names
    res_with_unnamed_row_means <- normalize_tfidf(m2, feature_means = row_means_no_names)
    expect_identical(res_with_unnamed_row_means, res)
    rownames(m2) <- NULL
    res_with_unnamed_row_names <- normalize_tfidf(m2, feature_means = row_means)
    res_with_unnamed_row_names_row_means <- normalize_tfidf(m2, feature_means = row_means_no_names)
    rownames(res) <- NULL
    expect_identical(as(res_with_unnamed_row_names, "dgCMatrix"), as(res, "dgCMatrix"))
    expect_identical(as(res_with_unnamed_row_names_row_means, "dgCMatrix"), as(res, "dgCMatrix"))
})

test_that("normalize_log works", {
    m <- generate_sparse_matrix(5, 5)
    res_dgc <- m %*% diag(1/colSums(m)) %>% as("dgCMatrix")
    m2 <- as(m, "IterableMatrix")
    # Test that default params yield the same as log1p on dgCMatrix
    res_1 <- as(normalize_log(m2), "dgCMatrix")
    expect_equal(res_1, log1p(res_dgc*1e4), tolerance = 1e-6)
    
    # Test that changing scale factor works
    res_2 <- as(normalize_log(m2, scale_factor = 1e5), "dgCMatrix")
    res_2_partial <- as(normalize_log(scale_factor = 1e5)(m2), "dgCMatrix")
    expect_equal(res_2, log1p(res_dgc*1e5), tolerance = 1e-6)
})