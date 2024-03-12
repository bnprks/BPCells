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

test_that("AnnData subset hasn't regressed", {
  dir <- withr::local_tempdir()
  # Make a copy since apparently reading the test hdf5 file causes modifications that git detects
  file.copy("../data/mini_mat.h5ad", file.path(dir, "mini_mat.h5ad"))
  x <- open_matrix_anndata_hdf5(file.path(dir, "mini_mat.h5ad")) %>%
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
  file.copy("../data/mini_mat.h5ad", file.path(dir, "mini_mat.h5ad"))

  # Test change row+col names on hdf5
  x <- open_matrix_anndata_hdf5(file.path(dir, "mini_mat.h5ad"))
  
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
  x_10x <- open_matrix_anndata_hdf5(file.path(dir, "mini_mat.h5ad")) %>%
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
  file.copy("../data/mini_mat.h5ad", file.path(dir, "mini_mat.h5ad"))
  x <- open_matrix_anndata_hdf5(file.path(dir, "mini_mat.h5ad"))
  xt <- open_matrix_anndata_hdf5(file.path(dir, "mini_mat.h5ad"), group="layers/transpose")
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