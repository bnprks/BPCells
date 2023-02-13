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
  expect_s4_class(m[c(2, 4, 9, 5), ], "MatrixMultiply")

  expect_equal(as(m[, c(4, 3, 6)], "dgCMatrix"), res[, c(4, 3, 6)], tolerance=testthat_tolerance())
  expect_s4_class(m[, c(4, 3, 6)], "MatrixMultiply")

  expect_equal(as(m[c(2, 4, 9, 5), c(4, 3, 6)], "dgCMatrix"), res[c(2, 4, 9, 5), c(4, 3, 6)], tolerance=testthat_tolerance())
  expect_s4_class(m[c(2, 4, 9, 5), c(4, 3, 6)], "MatrixMultiply")
})

test_that("Dense matrix-vector multiply works", {
  withr::local_seed(195123)

  m1 <- generate_sparse_matrix(5, 1000)
  m2 <- t(m1)

  i1 <- as(m1, "IterableMatrix")
  i2 <- t(i1)

  b <- generate_dense_matrix(5, 12)

  expect_identical(to_matrix(t(b) %*% m1), t(b) %*% i1)
  expect_identical(to_matrix(m2 %*% b), i2 %*% b)

  expect_identical(to_matrix(t(m1) %*% b), t(i1) %*% b)
  expect_identical(to_matrix(t(b) %*% t(m2)), t(b) %*% t(i2))

  y <- as.numeric(generate_dense_matrix(5, 1))
  expect_identical(y %*% as.matrix(m1), y %*% i1)
  expect_identical(as.matrix(m2) %*% y, i2 %*% y)

  expect_identical(t(as.matrix(m1)) %*% y, t(i1) %*% y)
  expect_identical(y %*% t(as.matrix(m2)), y %*% t(i2))

  # Check that everything works fine with integers
  i3 <- convert_matrix_type(i1, "uint32_t")
  i4 <- convert_matrix_type(i2, "uint32_t")

  expect_identical(to_matrix(t(b) %*% m1), t(b) %*% i3)
  expect_identical(to_matrix(m2 %*% b), i4 %*% b)

  expect_identical(to_matrix(t(m1) %*% b), t(i3) %*% b)
  expect_identical(to_matrix(t(b) %*% t(m2)), t(b) %*% t(i4))

  expect_identical(y %*% as.matrix(m1), y %*% i3)
  expect_identical(as.matrix(m2) %*% y, i4 %*% y)

  expect_identical(t(as.matrix(m1)) %*% y, t(i3) %*% y)
  expect_identical(y %*% t(as.matrix(m2)), y %*% t(i4))
})

test_that("Row/Col sum/mean works", { #nolint
  m1 <- generate_sparse_matrix(5, 1000)
  m2 <- t(m1)

  i1 <- as(m1, "IterableMatrix")
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
    min_1 = min_scalar(mi, 1e9),
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
    short_description(trans)
    expect_identical(dim(trans), dim(m))
    expect_identical(dimnames(trans), dimnames(m))
    expect_true(matrix_type(trans) %in% c("uint32_t", "float", "double"))
    matrix_is_transform(trans)

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
