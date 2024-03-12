generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero = 0.5, max_val = 10) {
  m <- matrix(rbinom(nrow * ncol, 1, fraction_nonzero) * sample.int(max_val, nrow * ncol, replace = TRUE), nrow = nrow)
  as(m, "dgCMatrix")
}
generate_dense_matrix <- function(nrow, ncol) {
  matrix(runif(nrow * ncol), nrow = nrow)
}

to_matrix <- function(x) {
  x <- as.matrix(x)
  attr(x, "dimnames") <- NULL
  x
}
to_vector <- function(x) {
  x <- as.numeric(x)
  attributes(x) <- NULL
  x
}

test_that("Conversion to base types works", {
  for (row_names in c(TRUE, FALSE)) {
    for (col_names in c(TRUE, FALSE)) {
      base_mat <- matrix(as.numeric(1:12), nrow=3)
      if (row_names) rownames(base_mat) <- sprintf("row_%d", 1:3)
      if (col_names) colnames(base_mat) <- sprintf("col_%d", 1:4)
      dgc_mat <- as(base_mat, "dgCMatrix")
      bpcells_mat <- as(dgc_mat, "IterableMatrix")
    
      expect_identical(as(bpcells_mat, "dgCMatrix"), dgc_mat)
      expect_identical(as.matrix(bpcells_mat), base_mat)
      expect_identical(as(bpcells_mat, "matrix"), base_mat)
    }
  }
})

test_that("Conversion to dense matrix keep data types", {
  m <- generate_sparse_matrix(10, 10)
  m_double <- as(m, "IterableMatrix")
  # convert to double matrix
  expect_identical(matrix_type(m_double), "double")
  expect_identical(as.matrix(m_double), as.matrix(m))
  # convert to double matrix
  m_float <- convert_matrix_type(m_double, "float")
  expect_identical(matrix_type(m_float), "float")
  expect_identical(as.matrix(m_float), as.matrix(m))
  # convert to integer matrix
  m_integer <- convert_matrix_type(m_double, "uint32_t")
  expect_identical(matrix_type(m_integer), "uint32_t")
  expect_identical(as.matrix(m_integer), matrix(as.integer(as.matrix(m)), nrow=nrow(m)))
  # value exceed `.Machine$integer.max` return double mode and warn message
  m_large <- generate_dense_matrix(10, 10)
  m_large[1L] <- m_large[1L] + .Machine$integer.max
  m_large_integer <- as(m_large, "dgCMatrix") %>%
    as("IterableMatrix") %>%
    convert_matrix_type("uint32_t")
  expect_identical(matrix_type(m_large_integer), "uint32_t")
  expect_warning(dense_mat <- as.matrix(m_large_integer), "integer.max")
  expect_identical(dense_mat, m_large)
})

test_that("Chained subsetting works", {
  m1 <- generate_dense_matrix(10, 10) %>% as("dgCMatrix")

  m2 <- as(m1, "IterableMatrix")

  # Singleton selections
  expect_identical(m1[c(2, 4, 3, 5), ], m2[c(2, 4, 3, 5), ] %>% as("dgCMatrix"))
  expect_identical(m1[, c(2, 4, 3, 5)], m2[, c(2, 4, 3, 5)] %>% as("dgCMatrix"))

  # Repeated selections on 1 axis
  expect_identical(m1[c(2, 4, 3, 5), ][c(3, 1), ], m2[c(2, 4, 3, 5), ][c(3, 1), ] %>% as("dgCMatrix"))
  expect_identical(m1[, c(2, 4, 3, 5)][, c(3, 1)], m2[, c(2, 4, 3, 5)][, c(3, 1)] %>% as("dgCMatrix"))

  # Repeated selections on 2 axis
  expect_identical(m1[c(2, 4, 3, 5), c(5, 4, 2, 3)][c(3, 1), c(1, 3)], m2[c(2, 4, 3, 5), c(5, 4, 2, 3)][c(3, 1), c(1, 3)] %>% as("dgCMatrix"))

  # Missing one axis in 2nd selection
  expect_identical(m1[c(2, 4, 3, 5), c(5, 4, 2, 3)][c(3, 1), ], m2[c(2, 4, 3, 5), c(5, 4, 2, 3)][c(3, 1), ] %>% as("dgCMatrix"))
  expect_identical(m1[c(2, 4, 3, 5), c(5, 4, 2, 3)][, c(1, 3)], m2[c(2, 4, 3, 5), c(5, 4, 2, 3)][, c(1, 3)] %>% as("dgCMatrix"))

  # Missing one axis in 1st selection
  expect_identical(m1[, c(5, 4, 2, 3)][c(3, 1), c(1, 3)], m2[, c(5, 4, 2, 3)][c(3, 1), c(1, 3)] %>% as("dgCMatrix"))
  expect_identical(m1[c(2, 4, 3, 5), ][c(3, 1), c(1, 3)], m2[c(2, 4, 3, 5), ][c(3, 1), c(1, 3)] %>% as("dgCMatrix"))
})

test_that("Subsetting transpose missing variable works", {
  # This addresses a specific bug and checks to prevent regression
  x <- matrix(1:12, nrow=3) %>%
    as("dgCMatrix") %>%
    as("IterableMatrix")

  # This had a missing variable error
  ans <- x[1:2,] %>%
    t() %>%
    .[,1:2] %>%
    as("dgCMatrix")
  
  expected <- matrix(1:12, nrow=3)[1:2,] %>%
    t() %>%
    as("dgCMatrix")
  expect_identical(ans, expected)
})

test_that("Subsetting matrix multiply works", {
  m1 <- generate_dense_matrix(10, 5) %>% as("dgCMatrix")
  m2 <- generate_dense_matrix(5, 6) %>% as("dgCMatrix")

  m <- as(m1, "IterableMatrix") %*% as(m2, "IterableMatrix")
  res <- m1 %*% m2

  expect_equal(as(m[c(2, 4, 9, 5), ], "dgCMatrix"), res[c(2, 4, 9, 5), ], tolerance=testthat_tolerance())
  expect_s4_class(m[c(2, 4, 9, 5), ], "MatrixSubset")
  expect_s4_class(m[sort(c(2, 4, 9, 5)), ], "MatrixMultiply")

  expect_equal(as(m[, c(4, 3, 6)], "dgCMatrix"), res[, c(4, 3, 6)], tolerance=testthat_tolerance())
  expect_s4_class(m[, sort(c(4, 3, 6))], "MatrixMultiply")
  expect_s4_class(m[, c(4, 3, 6)], "MatrixSubset")

  expect_equal(as(m[c(2, 4, 9, 5), c(4, 3, 6)], "dgCMatrix"), res[c(2, 4, 9, 5), c(4, 3, 6)], tolerance=testthat_tolerance())
  expect_s4_class(m[sort(c(2, 4, 9, 5)), sort(c(4, 3, 6))], "MatrixMultiply")
  expect_s4_class(m[c(2, 4, 9, 5), c(4, 3, 6)], "MatrixSubset")
})

test_that("Subsetting RowBindMatrices/ColBindMatrices works", {
  withr::local_seed(195123) 

  mat_list <- lapply(1:10, function(i) {
    generate_dense_matrix(i, 15)
  })
  
  ref <- do.call(rbind, mat_list)
  bpmat_list1 <- lapply(mat_list, function(x) as(x, "dgCMatrix") |> as("IterableMatrix"))
  bpmat_list2 <- lapply(mat_list, function(x) as(t(x), "dgCMatrix") |> as("IterableMatrix"))
  bpmat1 <- do.call(rbind, bpmat_list1)
  bpmat2 <- t(do.call(cbind, bpmat_list2))

  for (bpmat in list(bpmat1, bpmat2)) {
    j <- sample.int(ncol(ref), ncol(ref) %/% 2, replace=FALSE)
    
    # Try random subsets of 25% of the rows
    for (x in 1:5) {
      i <- sample.int(nrow(ref), nrow(ref) %/% 4, replace=FALSE)
      expect_identical(
        ref[i, j],
        bpmat[i, j] |> as("dgCMatrix") |> as.matrix()
      )

      expect_identical(
        ref[i, j],
        t(bpmat)[j,i] |> as("dgCMatrix") |> as.matrix() |> t()
      )
    }

    # Try random subsets of 3 contiguous rows
    for (x in 1:5) {
      i <- sample.int(nrow(ref) - 3, 1)
      i <- i:(i+3L)
      expect_identical(
        ref[i, j],
        bpmat[i, j] |> as("dgCMatrix") |> as.matrix()
      )

      expect_identical(
        ref[i, j],
        t(bpmat)[j,i] |> as("dgCMatrix") |> as.matrix() |> t()
      )
    }

    # Try just subsetting j 
    expect_identical(
      ref[, j],
      bpmat[, j] |> as("dgCMatrix") |> as.matrix()
    )

    expect_identical(
      ref[, j],
      t(bpmat)[j,] |> as("dgCMatrix") |> as.matrix() |> t()
    )
  }
})

test_that("rbind and cbind check types (#68 regression)", {
  m <- generate_dense_matrix(10, 5) %>% as("dgCMatrix") %>% as("IterableMatrix")

  expect_error(
    rbind(m, convert_matrix_type(m, "uint32_t")), "type"
  )
  expect_error(
    cbind(m, convert_matrix_type(m, "uint32_t")), "type"
  )
})

test_that("Subsetting to 0 dimensions works", {
  m1 <- generate_dense_matrix(10, 5) %>% as("dgCMatrix")
  m2 <- as(m1, "IterableMatrix")
  expect_identical(
    m2[rep_len(FALSE, nrow(m1)),] %>% as("dgCMatrix"),
    m1[rep_len(FALSE, nrow(m1)),]
  )

  expect_identical(
    m2[,rep_len(FALSE, ncol(m1))] %>% as("dgCMatrix"),
    m1[,rep_len(FALSE, ncol(m1))]
  )

  expect_identical(
    m2[1:3,][integer(0),] %>% as("dgCMatrix"),
    m1[1:3,][integer(0),]
  )

  expect_identical(
    m2[1:3,][,integer(0)] %>% as("dgCMatrix"),
    m1[1:3,][,integer(0)]
  )

  expect_identical(
    m2[,1:3][integer(0),] %>% as("dgCMatrix"),
    m1[,1:3][integer(0),]
  )

  expect_identical(
    m2[,1:3][,integer(0)] %>% as("dgCMatrix"),
    m1[,1:3][,integer(0)]
  )
})

test_that("Subset rename dims works (Issue #65)", {
  m1 <- generate_sparse_matrix(20, 20) %>% as("dgCMatrix")
  m2 <- as(m1, "IterableMatrix")

  dimnames(m2) <- dimnames(m1)
  m2 <- m2*2
  m1 <- m1*2
  i <- sample.int(nrow(m1), nrow(m1) - 2)
  j <- sample.int(ncol(m1), nrow(m1) - 2)
  m2 <- m2[i,j]
  m1 <- m1[i,j]
  
  m1 <- cbind(m1, m1)
  m2 <- cbind(m2, m2)

  rownames(m1) <- paste("row", seq_len(nrow(m1)))
  colnames(m1) <- paste("col", seq_len(ncol(m1)))
  dimnames(m2) <- dimnames(m1)

  expect_identical(
    m2[,1:2] %>% as("dgCMatrix"),
    m1[,1:2]
  )

  expect_identical(
    m2[1:2,] %>% as("dgCMatrix"),
    m1[1:2,]
  )
})

test_that("Subset assignment works", {
  m1 <- generate_dense_matrix(10, 20) %>% as("dgCMatrix")
  rownames(m1) <- paste0("m1_row_", seq_len(nrow(m1)))
  colnames(m1) <- paste0("m1_col_", seq_len(ncol(m1)))
  r <- -(1 * generate_dense_matrix(10, 20)) %>% as("dgCMatrix")
  rownames(r) <- paste0("r_row_", seq_len(nrow(r)))
  colnames(r) <- paste0("r_col_", seq_len(ncol(r)))
  r2 <- as(r, "IterableMatrix")

  test_cases <- list(
    list(rows = sample.int(nrow(m1), 5), cols = sample.int(ncol(m1), 10)),
    list(rows = seq_len(nrow(m1)), cols = seq_len(ncol(m1))),
    list(rows = rep_len(FALSE, nrow(m1)), cols=rep_len(FALSE, ncol(m1)))
  )
  for (t in test_cases) {
    # Test subset just rows
    m2 <- as(m1, "IterableMatrix")
    m2[t$rows,] <- r2[t$rows,]
    m1b <- m1
    m1b[t$rows,] <- r[t$rows,]
    expect_identical(as(m2, "dgCMatrix"), m1b)
    expect_error(m2[t$rows,] <- cbind(r2, r2))

    # Test subset just cols
    m2 <- as(m1, "IterableMatrix")
    m2[,t$cols] <- r2[,t$cols]
    m1b <- m1
    m1b[,t$cols] <- r[,t$cols]
    expect_identical(as(m2, "dgCMatrix"), m1b)
    expect_error(m2[,t$cols] <- cbind(r2, r2))

    # Test subset rows + cols
    m2 <- as(m1, "IterableMatrix")
    m2[t$rows, t$cols] <- r2[t$rows, t$cols]
    m1b <- m1
    m1b[t$rows, t$cols] <- r[t$rows, t$cols]
    expect_identical(as(m2, "dgCMatrix"), m1b)
    expect_error(m2[t$rows,t$cols] <- cbind(r2, r2))
  }

  # Test replacing the full matrix with rows + cols missing
  m2 <- as(m1, "IterableMatrix")
  m2[,] <- r2
  m1b <- m1
  m1b[,] <- r
  expect_identical(as(m2, "dgCMatrix"), m1b)

  # Test assigning a dense matrix and a dgCMatrix
  t <- test_cases[[1]]
  m2 <- as(m1, "IterableMatrix")
  m2[t$rows, t$cols] <- as.matrix(r)[t$rows, t$cols]
  m1b <- m1
  m1b[t$rows, t$cols] <- r[t$rows, t$cols]
  expect_identical(as(m2, "dgCMatrix"), m1b)

  m2 <- as(m1, "IterableMatrix")
  m2[t$rows, t$cols] <- r[t$rows, t$cols]
  expect_identical(as(m2, "dgCMatrix"), m1b)

  # Check that rownames + colnames are preserved
  rownames(m1) <- paste0("r", seq_len(nrow(m1)))
  colnames(m1) <- paste0("c", seq_len(ncol(m1)))
  m2 <- as(m1, "IterableMatrix")
  m2[t$rows, t$cols] <- as.matrix(r)[t$rows, t$cols]
  m1b <- m1
  m1b[t$rows, t$cols] <- r[t$rows, t$cols]
  expect_identical(as(m2, "dgCMatrix"), m1b)

  # Additional regression test from #67
  m2 <- as(m1, "IterableMatrix") |> rank_transform("col")
  m1b <- as(m2, "dgCMatrix")
  m2[1:5,1:5] <- as.matrix(r)[1:5,1:5]
  m1b[1:5,1:5] <- r[1:5,1:5]
  expect_identical(as(m2, "dgCMatrix"), m1b)

})

test_that("Transposing a 0 dimension matrix works", {
  m1 <- generate_dense_matrix(10, 5) %>% as("dgCMatrix") 
  rownames(m1) <- paste0("row", seq_len(nrow(m1)))
  colnames(m1) <- paste0("col", seq_len(ncol(m1)))
  m2 <- m1[integer(0),]
  rownames(m2) <- NULL # Don't worry about edge-case of rownames mismatching on NULL vs character(0)
  expect_identical(
    as(m2, "IterableMatrix") %>% transpose_storage_order() %>% as("dgCMatrix"),
    m2
  )

  m3 <- m1[,integer(0)]
  colnames(m3) <- NULL
  expect_identical(
    as(m3, "IterableMatrix") %>% transpose_storage_order() %>% as("dgCMatrix"),
    m3
  )
})

# Test dense multiplication ops given dgCMatrix m1 and equivalent iterable matrix i1
test_dense_multiply_ops <- function(m1, i1) {
  withr::local_seed(195123) 
  m2 <- t(m1)
  i2 <- t(i1)
  b1 <- generate_dense_matrix(nrow(m1), 12)
  b2 <- generate_dense_matrix(ncol(m1), 12)

  expect_identical(to_matrix(t(b1) %*% m1), t(b1) %*% i1)
  expect_identical(to_matrix(m2 %*% b1), i2 %*% b1)

  expect_identical(to_matrix(t(m1) %*% b1), t(i1) %*% b1)
  expect_identical(to_matrix(t(b1) %*% t(m2)), t(b1) %*% t(i2))

  expect_identical(to_matrix(t(b2) %*% m2), t(b2) %*% i2)
  expect_identical(to_matrix(m1 %*% b2), i1 %*% b2)

  expect_identical(to_matrix(t(m2) %*% b2), t(i2) %*% b2)
  expect_identical(to_matrix(t(b2) %*% t(m1)), t(b2) %*% t(i1))

  y1 <- as.numeric(generate_dense_matrix(nrow(m1), 1))
  y2 <- as.numeric(generate_dense_matrix(ncol(m1), 1))
  expect_identical(y1 %*% as.matrix(m1), y1 %*% i1)
  expect_identical(as.matrix(m2) %*% y1, i2 %*% y1)

  expect_identical(y2 %*% as.matrix(m2), y2 %*% i2)
  expect_identical(as.matrix(m1) %*% y2, i1 %*% y2)

  # Check that everything works fine with integers
  i3 <- convert_matrix_type(i1, "uint32_t")
  i4 <- convert_matrix_type(i2, "uint32_t")

  expect_identical(to_matrix(t(b1) %*% m1), t(b1) %*% i3)
  expect_identical(to_matrix(m2 %*% b1), i4 %*% b1)

  expect_identical(to_matrix(t(m1) %*% b1), t(i3) %*% b1)
  expect_identical(to_matrix(t(b1) %*% t(m2)), t(b1) %*% t(i4))

  expect_identical(to_matrix(t(b2) %*% m2), t(b2) %*% i4)
  expect_identical(to_matrix(m1 %*% b2), i3 %*% b2)

  expect_identical(to_matrix(t(m2) %*% b2), t(i4) %*% b2)
  expect_identical(to_matrix(t(b2) %*% t(m1)), t(b2) %*% t(i3))

  expect_identical(y1 %*% as.matrix(m1), y1 %*% i3)
  expect_identical(as.matrix(m2) %*% y1, i4 %*% y1)

  expect_identical(y2 %*% as.matrix(m2), y2 %*% i4)
  expect_identical(as.matrix(m1) %*% y2, i3 %*% y2)
}

test_that("Dense matrix-vector multiply works", {
  withr::local_seed(195123)

  m1 <- generate_sparse_matrix(5, 1000)
  i1 <- as(m1, "IterableMatrix")

  test_dense_multiply_ops(m1, i1)
})

# Test rows/col sum/mean given dgCMatrix m1 and equivalent iterable matrix i1
test_rowsum_colsum_rowmean_colmean <- function(m1, i1) {
  m2 <- t(m1)
  i2 <- t(i1)

  expect_identical(rowSums(i1), rowSums(m1))
  expect_identical(rowSums(i2), rowSums(m2))
  expect_identical(colSums(i1), colSums(m1))
  expect_identical(colSums(i2), colSums(m2))

  expect_identical(rowMeans(i1), rowMeans(m1))
  expect_identical(rowMeans(i2), rowMeans(m2))
  expect_identical(colMeans(i1), colMeans(m1))
  expect_identical(colMeans(i2), colMeans(m2))

  # Check that everything works fine with integers
  i3 <- convert_matrix_type(i1, "uint32_t")
  i4 <- convert_matrix_type(i2, "uint32_t")

  expect_identical(rowSums(i3), rowSums(m1))
  expect_identical(rowSums(i4), rowSums(m2))
  expect_identical(colSums(i3), colSums(m1))
  expect_identical(colSums(i4), colSums(m2))

  expect_identical(rowMeans(i3), rowMeans(m1))
  expect_identical(rowMeans(i4), rowMeans(m2))
  expect_identical(colMeans(i3), colMeans(m1))
  expect_identical(colMeans(i4), colMeans(m2))
}

test_that("Row/Col sum/mean works", { #nolint
  m1 <- generate_sparse_matrix(5, 1000)
  i1 <- as(m1, "IterableMatrix")
  test_rowsum_colsum_rowmean_colmean(m1, i1)
})

test_that("rbind Math works", {
  m1 <- generate_sparse_matrix(1000, 10)

  mat_list <- lapply(1:10, function(i) {
    as(m1, "IterableMatrix")[(i-1)*100 + 1:100,]
  })

  i1 <- do.call(rbind, mat_list)
  test_dense_multiply_ops(m1, i1)
  test_rowsum_colsum_rowmean_colmean(m1, i1)


  # Test for the specialized multiply ops with MatrixSubset of RowBind/ColBind
  row_permute <- sample.int(nrow(m1))
  col_permute <- sample.int(ncol(m1))
  test_dense_multiply_ops(m1[row_permute, col_permute], i1[row_permute, col_permute])

  i1 <- set_threads(i1, 3)
  test_dense_multiply_ops(m1, i1)
  test_rowsum_colsum_rowmean_colmean(m1, i1)
})

test_that("cbind Math works", {
  m1 <- generate_sparse_matrix(10, 1000)

  mat_list <- lapply(1:10, function(i) {
    as(m1, "IterableMatrix")[,(i-1)*100 + 1:100]
  })

  i1 <- do.call(cbind, mat_list)
  test_dense_multiply_ops(m1, i1)
  test_rowsum_colsum_rowmean_colmean(m1, i1)

  # Test for the specialized multiply ops with MatrixSubset of RowBind/ColBind
  row_permute <- sample.int(nrow(m1))
  col_permute <- sample.int(ncol(m1))
  test_dense_multiply_ops(m1[row_permute, col_permute], i1[row_permute, col_permute])

  i1 <- set_threads(i1, 3)
  test_dense_multiply_ops(m1, i1)
  test_rowsum_colsum_rowmean_colmean(m1, i1)
})

test_that("LinearOperator works", {
  m1 <- generate_sparse_matrix(5, 1000)
  op <- linear_operator(as(m1, "IterableMatrix"))

  b <- generate_dense_matrix(5, 12)
  y <- as.numeric(generate_dense_matrix(5, 1))

  expect_identical(to_matrix(t(b) %*% m1), t(b) %*% op)
  expect_identical(as.numeric(y %*% as.matrix(m1)), as.numeric(y %*% op))
})

test_that("Garbage collection between iterate_matrix doesn't mess things up", {
  m <- generate_sparse_matrix(20, 20) 
  # Apply log1p %>% square %>% pow(0.5) %>% expm1
  it <- m %>% 
    as("IterableMatrix") %>% 
    iterate_matrix() %>%
    iterate_matrix_log1p_cpp() %>%
    iterate_matrix_square_cpp() %>%
    iterate_matrix_pow_cpp(0.5) %>%
    iterate_matrix_expm1_cpp()
  # Garbage collect so R will destroy any intermediate data it has on the pointers
  gc()
  # Make sure everything still works
  res <- build_csparse_matrix_double_cpp(it)
  expect_equal(m, res, tolerance=testthat_tolerance())
})

test_that("Adverserial garbage collection of dgCMatrix doesn't mess things up", {
  gc_flags <- list(gc_1=FALSE, gc_2=FALSE)
  test_gc <- function() {
    set_gc_flag <- function(name) {function(e) gc_flags[[name]] <<- TRUE}
    # Run everything in a closure to create new env
    construct_matrix <- function(flag_name) {
      m <- as(matrix(1:12, nrow=3), "dgCMatrix")
      x <- environment()
      reg.finalizer(x, set_gc_flag(flag_name))
      attr(m, "gc_flag") <- x
      m
    }
    m <- construct_matrix("gc_1")
    m2 <- construct_matrix("gc_2") # Positive control for GC independent of BPCells
    
    it <- iterate_matrix(as(m, "IterableMatrix"))
    rm(m)
    rm(m2)
    gc()
    # Don't want GC to happen for the matrix we're iterating on, since the iterator object should maintain a reference
    expect_false(gc_flags[["gc_1"]])
    expect_true(gc_flags[["gc_2"]])
    res <- BPCells:::build_csparse_matrix_double_cpp(it)
    gc()
    # Now we do want GC to happen, since there's no longer a reference to maintain
    expect_true(gc_flags[["gc_1"]])
  }
  test_gc()
})

test_that("Matrix mask works", {
  dims <- list(c(5, 1000), c(5000, 5))
  for (d in dims) {
    m1 <- generate_sparse_matrix(d[1], d[2])
    mask <- generate_sparse_matrix(d[1], d[2], max_val = 1)
    
    res <- mask_matrix(as(m1, "IterableMatrix"), mask)
    res2 <- mask_matrix(t(as(m1, "IterableMatrix")), t(mask))
    res_inv <- mask_matrix(as(m1, "IterableMatrix"), mask, invert=TRUE)

    expect_identical(as(res, "dgCMatrix"), as(m1 * (1 - mask), "dgCMatrix"))
    expect_identical(as(res2, "dgCMatrix"), t(as(m1 * (1- mask), "dgCMatrix")))

    expect_identical(as(res_inv, "dgCMatrix"), as(m1 * (mask), "dgCMatrix"))
  }
})

test_that("Rank transform works", {
  test_mats <- list(generate_sparse_matrix(4000, 10), generate_sparse_matrix(10, 10, max_val=100))
  # Add offsets so we have explicit zeros and negative values
  test_mats[[1]]@x <- test_mats[[1]]@x - 5
  test_mats[[2]]@x <- test_mats[[2]]@x - 50
  
  for (m in test_mats) {
    for (type in c("uint32_t", "float", "double")) {
      if (type == "uint32_t") m@x <- m@x - min(m@x) # Get rid of negative value for uint32_t test
      r_rank <- as.matrix(m)
      for (i in seq_len(ncol(r_rank))) {
        r_rank[,i] <- rank(as.numeric(m[,i]))
      }

      r <- m %>% as("IterableMatrix") %>%
        convert_matrix_type(type) %>%
        rank_transform("col")
        
      bp_rank <- as(r, "dgCMatrix") %>% as.matrix()

      # Rank of 0 entries should map to 0
      expect_true(all(bp_rank[as.matrix(m) == 0] == 0))

      # When we undo the rank offset, we should match R
      for (i in seq_len(ncol(r_rank))) {
        bp_rank[,i] <- bp_rank[,i] + (nrow(r) + 1)/2 - mean(bp_rank[,i])
      }

      expect_identical(r_rank, bp_rank)
    }
  }
})

test_that("Relocating matrix inputs works", {
  # Test case: do a matrix multiply of two concatenated matrices,
  # then swap out the inputs and check it all works

  # X: 4x5
  x1 <- matrix(rnorm(9), nrow=3) 
  x2 <- matrix(rnorm(6), nrow=3)
  x3 <- matrix(rnorm(5), ncol=5)
  x <- rbind(cbind(x1, x2), x3)

  # Y: 5x3
  y1 <- matrix(rnorm(9), ncol=3) 
  y2 <- matrix(rnorm(6), ncol=3)
  y <- rbind(y1, y2)

  z <- x %*% y
  x1_bp <- as(as(x1, "dgCMatrix"), "IterableMatrix")
  x2_bp <- as(as(x2, "dgCMatrix"), "IterableMatrix")
  x3_bp <- as(as(x3, "dgCMatrix"), "IterableMatrix")

  y1_bp <- as(as(y1, "dgCMatrix"), "IterableMatrix")
  y2_bp <- as(as(y2, "dgCMatrix"), "IterableMatrix")

  x_bp <- rbind(cbind(x1_bp, x2_bp), x3_bp)
  y_bp <- rbind(y1_bp, y2_bp)

  z_bp <- x_bp %*% y_bp
  # Add in a length 1 layer which would catch an early error in all_matrix_inputs
  rownames(z_bp) <- rownames(z_bp)

  in_x <- all_matrix_inputs(x_bp)
  in_y <- all_matrix_inputs(y_bp)
  in_z <- all_matrix_inputs(z_bp)

  expect_identical(as.matrix(x_bp), x)
  expect_identical(as.matrix(y_bp), y)
  expect_equal(as.matrix(z_bp), z)

  expect_length(all_matrix_inputs(x_bp), 3)
  expect_length(all_matrix_inputs(y_bp), 2)
  expect_length(all_matrix_inputs(z_bp), 5)

  expect_identical(
      lapply(all_matrix_inputs(z_bp), class) %>% as.character(),
      rep("Iterable_dgCMatrix_wrapper", 5)
  )

  new_inputs <- lapply(all_matrix_inputs(z_bp), function(x) write_matrix_dir(x, tempfile(), compress=FALSE))

  z_bp2 <- z_bp
  all_matrix_inputs(z_bp2) <- new_inputs

  expect_equal(as.matrix(z_bp2), z)
  expect_length(all_matrix_inputs(z_bp2), 5)
  expect_identical(
      lapply(all_matrix_inputs(z_bp2), class) %>% as.character(),
      rep("MatrixDir", 5)
  )

})

test_that("Generic methods work", {
  # Generic methods to test:
  # - dim
  # - dimnames
  # - matrix_type
  # - matrix_is_transform
  # - show
  # - rowSums, colSums,
  m <- generate_sparse_matrix(5, 1000)
  rownames(m) <- paste0("row", seq_len(nrow(m)))
  colnames(m) <- paste0("col", seq_len(ncol(m)))
  mi <- as(m, "IterableMatrix")

  dir <- withr::local_tempdir()

  id_right <- as(Matrix::Diagonal(ncol(m)), "dgCMatrix")
  colnames(id_right) <- colnames(m)
  id_left <- as(Matrix::Diagonal(nrow(m)), "dgCMatrix")
  rownames(id_left) <- rownames(m)
  ident_transforms <- list(
    write_memory_uint32_t = mi %>% convert_matrix_type("uint32_t") %>% write_matrix_memory(compress = TRUE),
    write_memory_unpacked_uint32_t = mi %>% convert_matrix_type("uint32_t") %>% write_matrix_memory(compress = FALSE),
    write_memory_float = mi %>% convert_matrix_type("float") %>% write_matrix_memory(compress = TRUE),
    write_memory_unpacked_float = mi %>% convert_matrix_type("float") %>% write_matrix_memory(compress = FALSE),
    write_memory_double = mi %>% convert_matrix_type("double") %>% write_matrix_memory(compress = TRUE),
    write_memory_unpacked_double = mi %>% convert_matrix_type("double") %>% write_matrix_memory(compress = FALSE),
    shift_scale_1 = {
      t((mi * 1 + 0) * rep_len(1, nrow(m))) * rep_len(1, ncol(m)) + rep_len(0, ncol(m))
    } %>% t(),
    shift_scale_2 = {
      t((mi / 1 - 0) / rep_len(1, nrow(m)) - rep_len(0, nrow(m))) / rep_len(1, ncol(m))
    } %>% t(),
    multiply_right_1 = mi %*% id_right,
    multiply_right_2 = mi %*% as(id_right, "IterableMatrix"),
    multiply_left_1 = id_left %*% mi,
    multiply_types_1 = convert_matrix_type(mi, "uint32_t") %*% convert_matrix_type(as(id_right, "IterableMatrix"), "uint32_t"),
    multiply_types_2 = convert_matrix_type(mi, "float") %*% convert_matrix_type(as(id_right, "IterableMatrix"), "float"),
    mask = mask_matrix(mi, Matrix::sparseMatrix(i=integer(0), j=integer(0), x=integer(0), dims=dim(mi))),
    mask_float = convert_matrix_type(mi, "float") %>% mask_matrix(Matrix::sparseMatrix(i=integer(0), j=integer(0), x=integer(0), dims=dim(mi))),
    mask_int = convert_matrix_type(mi, "uint32_t") %>% mask_matrix(Matrix::sparseMatrix(i=integer(0), j=integer(0), x=integer(0), dims=dim(mi))),
    min_1 = min_scalar(mi, 1e9 + 0.1),
    min_row = min_by_row(mi, rep.int(1e9 + 0.1, nrow(mi))),
    min_col = min_by_col(mi, rep.int(1e9 + 0.1, ncol(mi))),
    round = round(mi),
    subset = mi[seq_len(nrow(m)), ][, seq_len(ncol(m))][seq_len(nrow(m)), seq_len(ncol(m))],
    rbind = rbind2(mi[1:2, ], mi[3:nrow(m)]),
    cbind = cbind2(mi[, 1:2], mi[, 3:ncol(m)]),
    rbind_uint32_t = convert_matrix_type(mi, "uint32_t")  %>% {rbind2(.[1:2, ], .[3:nrow(m)])},
    cbind_uint32_t = convert_matrix_type(mi, "uint32_t")  %>% {cbind2(.[, 1:2], .[, 3:ncol(m)])},
    rbind_float = convert_matrix_type(mi, "float")  %>% {rbind2(.[1:2, ], .[3:nrow(m)])},
    cbind_float = convert_matrix_type(mi, "float")  %>% {cbind2(.[, 1:2], .[, 3:ncol(m)])}
  )

  for (i in names(ident_transforms)) {
    trans <- ident_transforms[[i]]
    expect_no_error(short_description(trans))
    expect_identical(dim(trans), dim(m))
    expect_identical(dimnames(trans), dimnames(m))
    expect_true(matrix_type(trans) %in% c("uint32_t", "float", "double"))
    expect_type(matrix_is_transform(trans), "logical")
    expect_true(storage_order(trans) %in% c("roww", "col"))

    expect_identical(rowSums(trans), rowSums(m))
    expect_identical(colSums(trans), colSums(m))

    if (i %in% c("shift_scale_1", "shift_scale_2")) {
      expect_identical(as.matrix(as(trans, "dgCMatrix")), as.matrix(m))
    } else {
      expect_identical(as(trans, "dgCMatrix"), m)
    }

    # Test that garbage collection after creating the iterator doesn't cause issues
    # Create the C++ iterator
    if (matrix_type(trans) != "double")
      convert_function <- get(sprintf("convert_matrix_%s_%s_cpp", matrix_type(trans), "double"))
    else
      convert_function <- identity
    it <- iterate_matrix(trans) %>%
      convert_function() %>%
      iterate_matrix_square_cpp() %>%
      iterate_matrix_pow_cpp(0.5) 
    expect_type(it, "externalptr")
    # Delete any trailing XPtr references from R
    gc()
    # Check that the C++ iterator still works
    res <- build_csparse_matrix_double_cpp(it)
    res@Dimnames <- m@Dimnames
    if (i %in% c("shift_scale_1", "shift_scale_2")) {
      expect_identical(as.matrix(res), as.matrix(m))
    } else {
      expect_identical(res, m)
    }
  }
})
