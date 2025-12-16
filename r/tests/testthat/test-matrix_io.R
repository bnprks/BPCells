# Copyright 2022 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero = 0.5, max_val = 10) {
  m <- matrix(rbinom(nrow * ncol, 1, fraction_nonzero) * sample.int(max_val, nrow * ncol, replace = TRUE), nrow = nrow)
  rownames(m) <- paste0("row", seq_len(nrow(m)))
  colnames(m) <- paste0("col", seq_len(ncol(m)))
  as(m, "dgCMatrix")
}

test_that("Write 10x matrix to HDF5", {
  dir <- withr::local_tempdir()

  m <- generate_sparse_matrix(10, 9, fraction_nonzero=0.2)
  for (type in c("uint32_t", "float", "double")) {
    m2 <- m %>%
      as("IterableMatrix") %>%
      convert_matrix_type(type)
    if (type != "uint32_t") {
      expect_warning({
        mm1 <- write_matrix_10x_hdf5(m2, file.path(dir, paste0("mm1_", type, 10, ".h5")))
        mm2 <- write_matrix_10x_hdf5(m2, file.path(dir, paste0("mm2_", type, 10, ".h5")))
      })
    } else {
      mm1 <- write_matrix_10x_hdf5(m2, file.path(dir, paste0("mm1_", type, 10, ".h5")))
      mm2 <- write_matrix_10x_hdf5(m2, file.path(dir, paste0("mm2_", type, 10, ".h5")))
    }
    mm3 <- write_matrix_10x_hdf5(m2, file.path(dir, paste0("mm3_", type, 10, ".h5")), type = "auto")
    
    expect_identical(matrix_type(mm1), "uint32_t")
    expect_identical(matrix_type(mm2), "uint32_t")
    expect_identical(matrix_type(mm3), type)
  }
})

test_that("Write 10x matrix uses correct HDF5 data types", {
  dir <- withr::local_tempdir()
  mat_path <- file.path(dir, "10x.h5")
  m <- generate_sparse_matrix(10, 9, fraction_nonzero=0.2) %>%
    as("IterableMatrix") %>%
    convert_matrix_type("uint32_t")
  write_matrix_10x_hdf5(m, mat_path)
  expect_identical(hdf5_storage_type_cpp(mat_path, "matrix/indices"), "int64_t")
  expect_identical(hdf5_storage_type_cpp(mat_path, "matrix/indptr"), "int64_t")
  expect_identical(hdf5_storage_type_cpp(mat_path, "matrix/shape"), "int32_t")
})

test_that("Memory Matrix round-trips", {
  dir <- withr::local_tempdir()
  args <- list(
    list(nrow = 10, ncol = 9, fraction_nonzero = 0.2, max_val = 10),
    list(nrow = 100, ncol = 99, fraction_nonzero = 0.2, max_val = 100),
    list(nrow = 1000, ncol = 999, fraction_nonzero = 0.1, max_val = 1000)
  )
  for (a in args) {
    m <- do.call(generate_sparse_matrix, a)
    for (type in c("uint32_t", "float", "double")) {
      m2 <- m %>%
        as("IterableMatrix") %>%
        convert_matrix_type(type)
      mp <- write_matrix_memory(m2, compress = TRUE)
      mu <- write_matrix_memory(m2, compress = FALSE)

      mp_t <- write_matrix_memory(t(m2), compress = TRUE)
      mu_t <- write_matrix_memory(t(m2), compress = FALSE)

      expect_identical(m, as(mp, "dgCMatrix"))
      expect_identical(m, as(mu, "dgCMatrix"))
      expect_identical(matrix_type(mp), type)
      expect_identical(matrix_type(mu), type)

      expect_identical(t(m), as(mp_t, "dgCMatrix"))
      expect_identical(t(m), as(mu_t, "dgCMatrix"))
      expect_identical(matrix_type(mp_t), type)
      expect_identical(matrix_type(mu_t), type)
      expect_identical(mu, write_matrix_memory(mp, compress = FALSE))
    }
  }
})

test_that("DirPacked Matrix round-trips", {
  dir <- withr::local_tempdir()
  args <- list(
    list(nrow = 10, ncol = 9, fraction_nonzero = 0.2, max_val = 10),
    list(nrow = 100, ncol = 99, fraction_nonzero = 0.2, max_val = 100),
    list(nrow = 1000, ncol = 999, fraction_nonzero = 0.1, max_val = 1000)
  )
  for (a in args) {
    m <- do.call(generate_sparse_matrix, a)
    for (type in c("uint32_t", "float", "double")) {
      m2 <- m %>%
        as("IterableMatrix") %>%
        convert_matrix_type(type)
      mem <- write_matrix_memory(m2, compress = FALSE)
      mp <- write_matrix_dir(m2, file.path(dir, "packed", type, a[["nrow"]]), compress = TRUE)
      mu <- write_matrix_dir(m2, file.path(dir, "unpacked", type, a[["nrow"]]), compress = FALSE)
      expect_identical(m, as(mp, "dgCMatrix"))
      expect_identical(m, as(mu, "dgCMatrix"))
      expect_identical(matrix_type(mp), type)
      expect_identical(matrix_type(mu), type)
      expect_identical(mem, write_matrix_memory(mp, compress = FALSE))
      expect_identical(mem, write_matrix_memory(mu, compress = FALSE))

      mp_t <- write_matrix_dir(t(m2), file.path(dir, "t_packed", type, a[["nrow"]]), compress = TRUE)
      mu_t <- write_matrix_dir(t(m2), file.path(dir, "t_unpacked", type, a[["nrow"]]), compress = FALSE)
      expect_identical(t(m), as(mp_t, "dgCMatrix"))
      expect_identical(t(m), as(mu_t, "dgCMatrix"))
      expect_identical(matrix_type(mp_t), type)
      expect_identical(matrix_type(mu_t), type)
    }
  }
})


test_that("H5Packed Matrix round-trips", {
  dir <- withr::local_tempdir()
  args <- list(
    list(nrow = 10, ncol = 9, fraction_nonzero = 0.2, max_val = 10),
    list(nrow = 100, ncol = 99, fraction_nonzero = 0.2, max_val = 100),
    list(nrow = 1000, ncol = 999, fraction_nonzero = 0.1, max_val = 1000)
  )
  for (a in args) {
    m <- do.call(generate_sparse_matrix, a)
    for (type in c("uint32_t", "float", "double")) {
      m2 <- m %>%
        as("IterableMatrix") %>%
        convert_matrix_type(type)
      mem <- write_matrix_memory(m2, compress = FALSE)
      mu <- write_matrix_hdf5(m2, file.path(dir, "subdir", type, "file.h5"), paste0("packed/", as.character(a[["nrow"]])), compress = FALSE)
      mp <- write_matrix_hdf5(m2, file.path(dir, "subdir", type, "file.h5"), paste0("unpacked/", as.character(a[["nrow"]])), compress = TRUE)

      expect_identical(m, as(mp, "dgCMatrix"))
      expect_identical(m, as(mu, "dgCMatrix"))
      expect_identical(matrix_type(mp), type)
      expect_identical(matrix_type(mu), type)
      expect_identical(mem, write_matrix_memory(mp, compress = FALSE))
      expect_identical(mem, write_matrix_memory(mu, compress = FALSE))

      mu_t <- write_matrix_hdf5(t(m2), file.path(dir, "subdir", type, "file.h5"), paste0("t_packed/", as.character(a[["nrow"]])), compress = FALSE)
      mp_t <- write_matrix_hdf5(t(m2), file.path(dir, "subdir", type, "file.h5"), paste0("t_unpacked/", as.character(a[["nrow"]])), compress = TRUE)
      expect_identical(t(m), as(mp_t, "dgCMatrix"))
      expect_identical(t(m), as(mu_t, "dgCMatrix"))
      expect_identical(matrix_type(mp_t), type)
      expect_identical(matrix_type(mu_t), type)
    }
  }
})

test_that("H5Packed Matrix round-trips with gzip", {
  dir <- withr::local_tempdir()
  args <- list(
    list(nrow = 10, ncol = 9, fraction_nonzero = 0.2, max_val = 10),
    list(nrow = 100, ncol = 99, fraction_nonzero = 0.2, max_val = 100),
    list(nrow = 1000, ncol = 999, fraction_nonzero = 0.1, max_val = 1000)
  )
  for (a in args) {
    m <- do.call(generate_sparse_matrix, a)
    for (type in c("uint32_t", "float", "double")) {
      m2 <- m %>%
        as("IterableMatrix") %>%
        convert_matrix_type(type)
      mem <- write_matrix_memory(m2, compress = FALSE)
      mu <- write_matrix_hdf5(m2, file.path(dir, "subdir", type, "file.h5"), paste0("packed/", as.character(a[["nrow"]])), compress = FALSE, gzip_level = 4L)
      mp <- write_matrix_hdf5(m2, file.path(dir, "subdir", type, "file.h5"), paste0("unpacked/", as.character(a[["nrow"]])), compress = TRUE, gzip_level = 4L)
      
      expect_identical(m, as(mp, "dgCMatrix"))
      expect_identical(m, as(mu, "dgCMatrix"))
      expect_identical(matrix_type(mp), type)
      expect_identical(matrix_type(mu), type)
      expect_identical(mem, write_matrix_memory(mp, compress = FALSE))
      expect_identical(mem, write_matrix_memory(mu, compress = FALSE))
      
      mu_t <- write_matrix_hdf5(t(m2), file.path(dir, "subdir", type, "file.h5"), paste0("t_packed/", as.character(a[["nrow"]])), compress = FALSE, gzip_level = 4L)
      mp_t <- write_matrix_hdf5(t(m2), file.path(dir, "subdir", type, "file.h5"), paste0("t_unpacked/", as.character(a[["nrow"]])), compress = TRUE, gzip_level = 4L)
      expect_identical(t(m), as(mp_t, "dgCMatrix"))
      expect_identical(t(m), as(mu_t, "dgCMatrix"))
      expect_identical(matrix_type(mp_t), type)
      expect_identical(matrix_type(mu_t), type)
    }
  }
})

test_that("Packed Matrix works on all bit widths", {
  # Test matrix design
  # - 512 columns, each column has as many entries in it as the column index
  # - Attempt to use col_idx % 32 bits for storing values and col_idx % 8 bits for row offsets
  #   (I had slowness and crashing issues with dgCMatrices that had too many rows)
  # - Test packing the matrix with columns in-order and with shuffled columns

  withr::local_seed(1258123)

  x <- integer(0)
  row <- integer(0)
  col_ptr <- 0L
  for (col in 0:512) {
    max_bits <- col %% 32
    vals <- sample.int(2^max_bits, col, replace = TRUE)
    rows <- seq.int(0, col * 2^(col %% 8), 2^(col %% 8))[seq_len(col)] + 1

    col_ptr <- c(col_ptr, col)
    row <- c(row, rows)
    x <- c(x, vals)
  }
  col_ptr <- cumsum(col_ptr)
  m <- Matrix::sparseMatrix(i = row, x = x, p = col_ptr)

  mp <- write_matrix_memory(convert_matrix_type(m, "uint32_t"))

  expect_identical(m, as(mp, "dgCMatrix"))

  col_reorder <- sample.int(ncol(m))
  m <- m[, col_reorder]
  mp <- write_matrix_memory(convert_matrix_type(m, "uint32_t"))
  expect_identical(m, as(mp, "dgCMatrix"))
})

test_that("Transpose storage order works", {
  dir <- withr::local_tempdir()
  args <- list(
    list(nrow = 10, ncol = 9, fraction_nonzero = 0.2, max_val = 10),
    list(nrow = 100, ncol = 99, fraction_nonzero = 0.2, max_val = 100),
    list(nrow = 1000, ncol = 999, fraction_nonzero = 0.1, max_val = 1000)
  )
  for (a in args) {
    m <- do.call(generate_sparse_matrix, a)
    for (type in c("uint32_t", "float", "double")) {
      m2 <- m %>%
        as("IterableMatrix") %>%
        convert_matrix_type(type)

      out <- transpose_storage_order(m2)
      out_t <- transpose_storage_order(t(m2))

      expect_identical(out@transpose, TRUE)
      expect_identical(out_t@transpose, FALSE)
      expect_identical(m, as(out, "dgCMatrix"))
      expect_identical(t(m), as(out_t, "dgCMatrix"))
    }
  }
})

test_that("Transpose storage order on dense-transformed matrix (#71 regression test)", {
  mat <- matrix(1337, nrow=1)
  obj <- as(mat, "dgCMatrix") |> as("IterableMatrix")

  res <- BPCells::transpose_storage_order(obj + 1)
  expect_identical(as.matrix(res), mat + 1)
})

test_that("AnnData read backwards compatibility", {
  # Make a copy since apparently reading the test hdf5 file causes modifications that git detects
  dir <- withr::local_tempdir()

  ans <- matrix(c(1., 2., 0., 3., 0., 
                  0., 0., 2., 2., 1., 
                  0., 0., 2., 0., 2.), ncol=3) %>%
         as("dgCMatrix")
  rownames(ans) <- as.character(0:4)
  colnames(ans) <- as.character(0:2)
  obsm_ans <- ans[1:2,]
  rownames(obsm_ans) <- NULL
  varm_ans <- t(ans[,1:2])
  rownames(varm_ans) <- NULL 

  test_files <- c(
    "mini_mat.anndata-v0.6.22.h5ad", 
    "mini_mat.anndata-v0.7.h5ad", 
    "mini_mat.anndata-v0.10.9.h5ad", 
    "mini_mat_uint_16.anndata-v0.10.9.h5ad"
  )
  for (f in test_files) {
    file.copy(file.path("../data", f), file.path(dir, f))
    open_matrix_anndata_hdf5(file.path(dir, f)) %>%
      as("dgCMatrix") %>%
      expect_identical(ans)
    open_matrix_anndata_hdf5(file.path(dir, f), group="layers/transpose") %>%
      as("dgCMatrix") %>%
      expect_identical(ans)
    open_matrix_anndata_hdf5(file.path(dir, f), group="layers/dense") %>%
      as("dgCMatrix") %>%
      expect_identical(ans)
    
    open_matrix_anndata_hdf5(file.path(dir, f), group="obsm/obs_mat") %>%
      as("dgCMatrix") %>%
      expect_identical(obsm_ans)
    open_matrix_anndata_hdf5(file.path(dir, f), group="varm/var_mat") %>%
      as("dgCMatrix") %>%
      expect_identical(varm_ans)
  }
})

test_that("Anndata reads work with 64 bit datasets", {
  test_files <- c(
    "mini_mat_int_64_neg.anndata-v0.12.6.h5ad",
    "mini_mat_int_64.anndata-v0.12.6.h5ad"
  )

  ans_neg <- matrix(
    c(
      -83,  54,  30, -13, -14,  71,
      -83,  39, -60, -82,   5,  95,
       47,  52,  43,  57,   2, -75,
       67, -10,   0, -26, -64,  85,
       56,  28, -20,  64,   9, -12
    ),
    nrow = 6,
    ncol = 5
  ) %>% as("dgCMatrix")
  ans <- matrix(
    c(
      45, 22,  9, 55, 88,  6,
      85, 82, 27, 63, 16, 75,
      70, 35,  6, 97, 44, 89,
      67, 77, 75, 19, 36, 46,
      49,  4, 54, 15, 74, 68
    ),
    nrow = 6,
    ncol = 5
  ) %>% as("dgCMatrix")

  for (f in test_files) {
    mat_x <- open_matrix_anndata_hdf5(file.path("../data", f)) %>% as("dgCMatrix")
    mat_transpose <- open_matrix_anndata_hdf5(file.path("../data", f), group="layers/transpose") %>% as("dgCMatrix")
    mat_dense <- open_matrix_anndata_hdf5(file.path("../data", f), group="layers/dense") %>% as("dgCMatrix")
    mat_obs <- open_matrix_anndata_hdf5(file.path("../data", f), group="obsm/obs_mat") %>% as("dgCMatrix")
    mat_var <- open_matrix_anndata_hdf5(file.path("../data", f), group="varm/var_mat") %>% as("dgCMatrix")
    rownames(mat_obs) <- c("0", "1")
    rownames(mat_var) <- c("0", "1")

    if (grepl("neg", f)) {
      res <- ans_neg
    } else {
      res <- ans
    }
    rownames(res) <- as.character(0:(nrow(res)-1))
    colnames(res) <- as.character(0:(ncol(res)-1))
    expect_identical(mat_x, res)
    expect_identical(mat_transpose, res)
    expect_identical(mat_dense, res)
    expect_identical(mat_obs, res[1:2, ])
    expect_identical(mat_var, t(res[, 1:2]) )
  }
})

test_that("AnnData subset hasn't regressed", {
  dir <- withr::local_tempdir()
  # Make a copy since apparently reading the test hdf5 file causes modifications that git detects
  file.copy("../data/mini_mat.anndata-v0.10.9.h5ad", file.path(dir, "mini_mat.anndata-v0.10.9.h5ad"))
  x <- open_matrix_anndata_hdf5(file.path(dir, "mini_mat.anndata-v0.10.9.h5ad")) %>%
    .[1:nrow(.),1:ncol(.)]

  s <- matrix_stats(x, row_stats="mean", col_stats="mean")
  colnames(s$row_stats) <- NULL
  colnames(s$col_stats) <- NULL
  expect_identical(s$row_stats["mean",], c(1/3,2/3, 4/3, 5/3, 1))
  expect_identical(s$col_stats["mean",], c(1.2, 1, 0.8))
})

test_that("AnnData write works", {
  dir <- withr::local_tempdir()
  mat_1 <- generate_sparse_matrix(10, 15)
  rownames(mat_1) <- paste0("mat1_row", seq_len(nrow(mat_1)))
  colnames(mat_1) <- NULL
  mat_2 <- generate_sparse_matrix(10, 20)
  
  mat_1_res <- write_matrix_anndata_hdf5(as(mat_1, "IterableMatrix"), file.path(dir, "mat.h5ad")) %>%
    as("dgCMatrix")
  mat_2_res <- write_matrix_anndata_hdf5(t(as(t(mat_2), "IterableMatrix")), file.path(dir, "mat.h5ad"), group="varm/mat2") %>%
    as("dgCMatrix")

  expect_identical(rownames(mat_1_res), rownames(mat_1))
  expect_identical(colnames(mat_1_res), as.character(seq_len(ncol(mat_1)) - 1L))

  expect_identical(rownames(mat_2_res), rownames(mat_1))
  expect_identical(colnames(mat_2_res), NULL)

  dimnames(mat_1) <- NULL
  dimnames(mat_2) <- NULL
  dimnames(mat_1_res) <- NULL
  dimnames(mat_2_res) <- NULL
  expect_identical(mat_1, mat_1_res)
  expect_identical(mat_2, mat_2_res)
})

test_that("AnnData write dense matrix works", {
  dir <- withr::local_tempdir()
  mat_1 <- generate_sparse_matrix(10, 15)
  rownames(mat_1) <- paste0("mat1_row", seq_len(nrow(mat_1)))
  colnames(mat_1) <- NULL
  mat_2 <- generate_sparse_matrix(10, 20)

  mat_1_res <- write_matrix_anndata_hdf5_dense(as(mat_1, "IterableMatrix"), file.path(dir, "mat.h5ad")) %>%
    as.matrix()
  mat_2_res <- write_matrix_anndata_hdf5_dense(t(as(t(mat_2), "IterableMatrix")), file.path(dir, "mat.h5ad"), dataset = "varm/mat2") %>%
    as.matrix()
  expect_identical(rownames(mat_1_res), rownames(mat_1))
  expect_identical(colnames(mat_1_res), as.character(seq_len(ncol(mat_1)) - 1L))

  expect_identical(rownames(mat_2_res), rownames(mat_1))
  expect_identical(colnames(mat_2_res), NULL)

  dimnames(mat_1) <- NULL
  dimnames(mat_2) <- NULL
  dimnames(mat_1_res) <- NULL
  dimnames(mat_2_res) <- NULL
  expect_identical(as.matrix(mat_1), mat_1_res)
  expect_identical(as.matrix(mat_2), mat_2_res)

  # Test empty columns
  mat_3 <- generate_sparse_matrix(10, 15)
  mat_3[, 4] <- 0
  mat_3_res <- write_matrix_anndata_hdf5_dense(as(mat_3, "IterableMatrix"), file.path(dir, "mat_3.h5ad")) %>%
    as.matrix()
  expect_identical(as.matrix(mat_3), mat_3_res)

  # Test empty columns
  mat_3 <- generate_sparse_matrix(10, 15)
  mat_3[, 4:6] <- 0
  mat_3_res <- write_matrix_anndata_hdf5_dense(as(mat_3, "IterableMatrix"), file.path(dir, "mat_4.h5ad")) %>%
    as.matrix()
  expect_identical(as.matrix(mat_3), mat_3_res)
  
  m <- matrix(0, nrow = 3, ncol = 4)
  m[2, 2] <- 1
  m[3, 4] <- 1
  rownames(m) <- paste0("row", seq_len(nrow(m)))
  colnames(m) <- paste0("col", seq_len(ncol(m)))
  mat <- m |> as("dgCMatrix") |> as("IterableMatrix")
  ans <- write_matrix_anndata_hdf5_dense(mat, file.path(dir, "zeros.h5"))
  expect_identical(as.matrix(mat), as.matrix(ans))

  # Create a dense IterableMatrix
  mat_3 <- as(mat_1, "IterableMatrix") %>%
    multiply_cols(1 / Matrix::colSums(mat_1)) %>%
    log1p()
  stats <- matrix_stats(mat_3, row_stats = "variance")
  gene_means <- stats$row_stats["mean", ]
  gene_vars <- stats$row_stats["variance", ]
  mat_3 <- (mat_3 - gene_means) / gene_vars
  rownames(mat_3) <- paste0("mat3_row", seq_len(nrow(mat_3)))
  colnames(mat_3) <- paste0("mat3_col", seq_len(ncol(mat_3)))
  mat_3_res <- write_matrix_anndata_hdf5_dense(mat_3, file.path(dir, "mat2.h5ad")) %>%
    as.matrix()
  expect_identical(as.matrix(mat_3), mat_3_res)
  mat_3_res <- write_matrix_anndata_hdf5_dense(t(mat_3), file.path(dir, "mat3.h5ad")) %>%
    as.matrix()
  expect_identical(as.matrix(t(mat_3)), mat_3_res)
})

test_that("AnnData types round-trip", {
  dir <- withr::local_tempdir()
  mat <- generate_sparse_matrix(10, 15)
  for (type in c("uint32_t", "float", "double")) {
    typed_mat <- as(mat, "IterableMatrix") %>% convert_matrix_type(type)
    mat_res <- write_matrix_anndata_hdf5(typed_mat, file.path(dir, "mat.h5ad"), group=paste0("layers/", type))
    expect_identical(matrix_type(mat_res), type)
    expect_identical(as(mat_res, "dgCMatrix"), mat)
  }
})

test_that("AnnData and 10x row/col rename works", {
  dir <- withr::local_tempdir()
  # Make a copy since apparently reading the test hdf5 file causes modifications that git detects
  file.copy("../data/mini_mat.anndata-v0.10.9.h5ad", file.path(dir, "mini_mat.anndata-v0.10.9.h5ad"))

  # Test change row+col names on hdf5
  x <- open_matrix_anndata_hdf5(file.path(dir, "mini_mat.anndata-v0.10.9.h5ad"))
  
  orig_colnames <- colnames(x)
  orig_rownames <- rownames(x)
  rownames(x) <- paste0("row", rownames(x))
  colnames(x) <- paste0("col", colnames(x))
  expect_false(identical(orig_colnames, colnames(x)))
  expect_false(identical(orig_rownames, rownames(x)))

  x2 <- write_matrix_memory(x, compress=FALSE)
  expect_identical(rownames(x2), rownames(x))
  expect_identical(colnames(x2), colnames(x))

  x2_t <- write_matrix_memory(t(x), compress=FALSE)
  expect_identical(colnames(x2_t), rownames(x))
  expect_identical(rownames(x2_t), colnames(x))

  # Test change row+col names on 10x
  x_10x <- open_matrix_anndata_hdf5(file.path(dir, "mini_mat.anndata-v0.10.9.h5ad")) %>%
    convert_matrix_type("uint32_t") %>%
    write_matrix_10x_hdf5(file.path(dir, "mini_mat.h5"))
  rownames(x_10x) <- paste0("2row", rownames(x_10x))
  colnames(x_10x) <- paste0("2col", colnames(x_10x))
  expect_false(identical(orig_colnames, colnames(x_10x)))
  expect_false(identical(orig_rownames, rownames(x_10x)))

  x2 <- write_matrix_memory(x_10x, compress=FALSE)
  expect_identical(rownames(x2), rownames(x_10x))
  expect_identical(colnames(x2), colnames(x_10x))

  x2_t <- write_matrix_memory(t(x_10x), compress=FALSE)
  expect_identical(colnames(x2_t), rownames(x_10x))
  expect_identical(rownames(x2_t), colnames(x_10x))
})

test_that("AnnData write to dir matrix works (#57 regression test)", {
  dir <- withr::local_tempdir()
  # Make a copy since apparently reading the test hdf5 file causes modifications that git detects
  file.copy("../data/mini_mat.anndata-v0.10.9.h5ad", file.path(dir, "mini_mat.anndata-v0.10.9.h5ad"))
  x <- open_matrix_anndata_hdf5(file.path(dir, "mini_mat.anndata-v0.10.9.h5ad"))
  xt <- open_matrix_anndata_hdf5(file.path(dir, "mini_mat.anndata-v0.10.9.h5ad"), group="layers/transpose")
  y <- write_matrix_dir(x, file.path(dir, "mini_mat"))
  yt <- write_matrix_dir(xt, file.path(dir, "mini_mat_t"))
  
  expect_identical(as(x, "dgCMatrix"), as(y, "dgCMatrix"))
  expect_identical(as(xt, "dgCMatrix"), as(yt, "dgCMatrix"))
  expect_identical(as(x, "dgCMatrix"), as(xt, "dgCMatrix"))
})

test_that("Renaming transformed matrix works", {
  x <- matrix(1:12, nrow=3)
  rownames(x) <- paste0("row", seq_len(nrow(x)))
  colnames(x) <- paste0("col", seq_len(ncol(x)))
  
  # Check that dim names changes are properly preserved (or removed)
  x1 <- x %>%
    as("dgCMatrix") %>%
    as("IterableMatrix") %>%
    .[1:3,1:4]
  
  rownames(x1) <- paste0("newrow", seq_len(nrow(x)))
  colnames(x1) <- paste0("newcol", seq_len(ncol(x)))
  expect_identical(
    dimnames(x1), dimnames(as.matrix(x1))
  )

  x2 <- x1
  rownames(x2) <- NULL
  colnames(x2) <- NULL
  expect_identical(
    rownames(x2), rownames(as.matrix(x2))
  )
  expect_identical(
    colnames(x2), colnames(as.matrix(x2))
  )

  # Check that this works across subsets
  x3 <- x1[1:2,3:4]
  expect_s4_class(x3, "RenameDims")
  expect_identical(
    as.matrix(x3), as.matrix(x1)[1:2,3:4]
  )

  # Check that changing dimnames twice coalesces
  x4 <- x3
  rownames(x4) <- paste0("newnewrow", seq_len(nrow(x4)))
  colnames(x4) <- paste0("newnewcol", seq_len(ncol(x4)))
  expect_s4_class(x4, "RenameDims")
  expect_s4_class(x4@matrix, class(x3@matrix))

  # Check that the dimnames are preserved after a transform and a write
  # (Bug re-reported in issue #29)
  expect_identical(
    dimnames(x1), dimnames(write_matrix_memory(log1p(x1), compress=FALSE))
  )
  # Test that this still works with rbind and cbind
  expect_identical(
    dimnames(x1), dimnames(write_matrix_memory(rbind(x1[1:2,],x1[3,]), compress=FALSE))
  )
  expect_identical(
    dimnames(x1), dimnames(write_matrix_memory(cbind(x1[,1:2],x1[,3:4]), compress=FALSE))
  )
})



test_that("Matrix without names works", {
  dir <- withr::local_tempdir()
  m <- generate_sparse_matrix(nrow = 10, ncol = 9, fraction_nonzero = 0.2, max_val = 10)
  rownames(m) <- NULL
  colnames(m) <- NULL
  m1 <- m %>% as("IterableMatrix") %>% write_matrix_hdf5(file.path(dir, "nameless.h5"), "mat")
  m2 <- m %>% as("IterableMatrix") %>% write_matrix_dir(file.path(dir, "nameless-mat"))
  expect_identical(m, as(m1, "dgCMatrix"))
  expect_identical(m, as(m2, "dgCMatrix"))
})

test_that("H5 overwrite works", {
  dir <- withr::local_tempdir()
  m1 <- matrix(1:12, nrow=3)
  m2 <- m1
  m2[1,1] <- 5
  m1 %>% 
    as("dgCMatrix") %>% 
    as("IterableMatrix") %>%
    write_matrix_hdf5(file.path(dir, "overwrite.h5"), "mat")
  # writing without "overwrite" set should result in an error, and no data changed
  expect_error({
    m2 %>% 
      as("dgCMatrix") %>%
      as("IterableMatrix") %>%
      write_matrix_hdf5(file.path(dir, "overwrite.h5"), "mat")
  })
  expect_identical(
    as(m1, "dgCMatrix"),
    open_matrix_hdf5(file.path(dir, "overwrite.h5"), "mat") %>%
      as("dgCMatrix")
  )
  # writing with "overwrite" set should run and result in data change
  rlang::reset_message_verbosity("hdf5_overwrite")
  expect_message({
    m2 %>% 
      as("dgCMatrix") %>%
      as("IterableMatrix") %>%
      write_matrix_hdf5(file.path(dir, "overwrite.h5"), "mat", overwrite=TRUE)
  }, "dataset does not free old storage")
  expect_identical(
    as(m2, "dgCMatrix"),
    open_matrix_hdf5(file.path(dir, "overwrite.h5"), "mat") %>%
      as("dgCMatrix")
  )
  # Should work even when the overwritten h5 doesn't exist
  expect_identical(
    as(m2, "dgCMatrix"),
    m2 %>% 
      as("IterableMatrix") %>%
      write_matrix_hdf5(file.path(dir, "overwrite_new.h5"), "mat", overwrite = TRUE) %>%
      as("dgCMatrix")
  )
  # It should be okay even if we overwrite the same source we're loading from
  expect_identical(
    as(m2, "dgCMatrix")[c(1,3),],
    open_matrix_hdf5(file.path(dir, "overwrite.h5"), "mat") %>% 
      .[c(1,3),] %>%
      write_matrix_hdf5(file.path(dir, "overwrite.h5"), "mat", overwrite=TRUE) %>%
      as("dgCMatrix")
  )

  
})

test_that("Dir overwrite works", {
  dir <- withr::local_tempdir()
  m1 <- matrix(1:12, nrow=3)
  m2 <- m1
  m2[1,1] <- 5
  m1 %>% 
    as("dgCMatrix") %>% 
    as("IterableMatrix") %>%
    write_matrix_dir(file.path(dir, "overwrite-mat"))
  # writing without "overwrite" set should result in an error, and no data changed
  expect_error({
    m2 %>% 
      as("dgCMatrix") %>%
      as("IterableMatrix") %>%
      write_matrix_dir(file.path(dir, "overwrite-mat"))
  }, "exists")
  expect_identical(
    as(m1, "dgCMatrix"),
    open_matrix_dir(file.path(dir, "overwrite-mat")) %>%
      as("dgCMatrix")
  )
  # writing with "overwrite" set should run and result in data change
  m2 %>% 
    as("dgCMatrix") %>%
    as("IterableMatrix") %>%
    write_matrix_dir(file.path(dir, "overwrite-mat"), overwrite=TRUE)
  expect_identical(
    as(m2, "dgCMatrix"),
    open_matrix_dir(file.path(dir, "overwrite-mat")) %>%
      as("dgCMatrix")
  )
  # It should be okay even if we overwrite the same source we're loading from
  expect_identical(
    as(m2, "dgCMatrix")[c(1,3),],
    open_matrix_dir(file.path(dir, "overwrite-mat")) %>% 
      .[c(1,3),] %>%
      write_matrix_dir(file.path(dir, "overwrite-mat"), overwrite=TRUE) %>%
      as("dgCMatrix")
  )
})

test_that("Mtx import works", {
  md <- import_matrix_market("../data/double_mat.mtx")
  md_t <- import_matrix_market("../data/double_mat.mtx", row_major=TRUE)
  mi <- import_matrix_market("../data/int_mat.mtx.gz")
  mi_t <- import_matrix_market("../data/int_mat.mtx.gz", row_major=TRUE)

  expect_identical(storage_order(md), "col")
  expect_identical(storage_order(mi), "col")
  expect_identical(storage_order(md_t), "row")
  expect_identical(storage_order(mi_t), "row")

  expect_identical(md@type, "double")
  expect_identical(md_t@type, "double")
  expect_identical(mi@type, "uint32_t")
  expect_identical(mi_t@type, "uint32_t")

  ans_d <- Matrix::readMM("../data/double_mat.mtx") %>% as("dgCMatrix")
  ans_i <- Matrix::readMM("../data/int_mat.mtx.gz") %>% as("dgCMatrix")
  expect_identical(
    ans_d,
    as(md, "dgCMatrix")
  )
  expect_identical(
    ans_d,
    as(md_t, "dgCMatrix")
  )
  expect_identical(
    ans_i,
    as(mi, "dgCMatrix")
  )
  expect_identical(
    ans_i,
    as(mi_t, "dgCMatrix")
  )
})

test_that("Opening >64 matrices works", {
  skip_on_ci() # TODO: diagnose why github CI on windows fails here. Suggested start by debug printing around `_setmaxstdio()` calls

  # Note: this test is primarily for Windows, where there's a 512 file limit by default
  dir <- withr::local_tempdir()
  
  n_matrices <- 80
  matrices <- lapply(seq_len(n_matrices), function(x) generate_sparse_matrix(10, 20)) %>%
    lapply(function(x) write_matrix_dir(convert_matrix_type(x, "uint32_t"), tempfile("mat", tmpdir=dir)))
  
  matrix_rbind <- do.call(rbind, matrices)
  matrix_cbind <- do.call(cbind, matrices)

  # Test that we can do a basic operation on the matrix
  expect_identical(colSums(matrix_rbind), colSums(as(matrix_rbind, "dgCMatrix")))
  expect_identical(colSums(matrix_cbind), colSums(as(matrix_cbind, "dgCMatrix")))
})

test_that("Regression test for Rcpp out-of-bounds warning", {
  m <- Matrix::sparseMatrix(i=1, j=1, x=1, dims=c(1,1)) %>% as("dgCMatrix")
  expect_no_warning({
    m2 <- m %>%
      as("IterableMatrix") %>%
      write_matrix_memory() %>%
      as("dgCMatrix")
  })
  expect_identical(m, m2)
})

test_that("Matrix reads work from Windows with CRLF/LF", {
  # This is a regression test for issue #253, making sure that
  # matrices with crlf line endings for string arrays open without issues

  # Confirm that we actually have lf and crlf files as expected
  expect_false(any(charToRaw("\r") == readBin("../data/mini_mat_lf/version", "raw", 1000)))
  expect_true(any(charToRaw("\r") == readBin("../data/mini_mat_crlf/version", "raw", 1000)))

  # We should be able to read crlf-containing files without complaint
  expect_no_error(open_matrix_dir("../data/mini_mat_lf"))
  expect_no_error(open_matrix_dir("../data/mini_mat_crlf"))

  # We should output only lf moving forward
  dir <- withr::local_tempdir()
  path <- file.path(dir, "test_mat")
  x <- open_matrix_dir("../data/mini_mat_lf") %>%
    write_matrix_dir(path)
  expect_false(any(charToRaw("\r") == readBin(file.path(path, "version"), "raw", 1000)))
  expect_false(any(charToRaw("\r") == readBin(file.path(path, "row_names"), "raw", 1000)))
  expect_false(any(charToRaw("\r") == readBin(file.path(path, "col_names"), "raw", 1000)))
})