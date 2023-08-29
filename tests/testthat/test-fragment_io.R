
context("Fragment IO v2")

test_that("10x Fragments round-trip", {
  dir <- withr::local_tempdir()

  in_path <- "../data/mini_fragments.tsv.gz" #nolint
  out_path1 <- file.path(dir, "fragments_copy1.tsv")
  out_path2 <- file.path(dir, "fragments_copy2.tsv.gz")


  input <- open_fragments_10x(in_path)
  write_fragments_10x(input, out_path1)
  write_fragments_10x(input, out_path2)

  expected <- readr::read_file(in_path)
  res1 <- readr::read_file(out_path1)
  res2 <- readr::read_file(out_path2)
  expect_identical(res1, expected)
  expect_identical(res2, expected)
})

test_that("Packed Fragments example data round-trip", {
  in_path <- "../data/mini_fragments.tsv.gz" #nolint
  raw_fragments <- write_fragments_memory(open_fragments_10x(in_path), compress = FALSE)
  packed_fragments <- write_fragments_memory(raw_fragments)
  raw_fragments2 <- write_fragments_memory(packed_fragments, compress = FALSE)

  expect_identical(raw_fragments, raw_fragments2)
})

test_that("Unpacked Fragments example data round-trip", {
  in_path <- "../data/mini_fragments.tsv.gz" #nolint
  raw_fragments <- write_fragments_memory(open_fragments_10x(in_path), compress = FALSE)
  packed_fragments <- write_fragments_memory(raw_fragments, compress = FALSE)
  raw_fragments2 <- write_fragments_memory(packed_fragments, compress = FALSE)

  expect_identical(raw_fragments, raw_fragments2)
})

test_that("Binary File Fragments example data round-trip", {
  dir <- withr::local_tempdir()

  in_path <- "../data/mini_fragments.tsv.gz" #nolint
  raw_fragments <- write_fragments_memory(open_fragments_10x(in_path), compress = FALSE)

  write_fragments_dir(raw_fragments, file.path(dir, "unpacked"), compress = FALSE)
  write_fragments_dir(raw_fragments, file.path(dir, "packed"), compress = TRUE)

  unpacked <- open_fragments_dir(file.path(dir, "unpacked"))
  packed <- open_fragments_dir(file.path(dir, "packed"))

  expect_identical(raw_fragments, write_fragments_memory(unpacked, compress = FALSE))
  expect_identical(raw_fragments, write_fragments_memory(packed, compress = FALSE))
})


test_that("HDF5 File Fragments example data round-trip", {
  dir <- withr::local_tempdir()

  in_path <- "../data/mini_fragments.tsv.gz" #nolint
  raw_fragments <- write_fragments_memory(open_fragments_10x(in_path), compress = FALSE)

  write_fragments_hdf5(raw_fragments, file.path(dir, "file.h5"), "unpacked", compress = FALSE)
  write_fragments_hdf5(raw_fragments, file.path(dir, "file.h5"), "packed", compress = TRUE)

  unpacked <- open_fragments_hdf5(file.path(dir, "file.h5"), "unpacked")
  packed <- open_fragments_hdf5(file.path(dir, "file.h5"), "packed")

  expect_identical(raw_fragments, write_fragments_memory(unpacked, compress = FALSE))
  expect_identical(raw_fragments, write_fragments_memory(packed, compress = FALSE))
})

test_that("HDF5 File Fragments example data round-trip with gzip", {
  dir <- withr::local_tempdir()
  
  in_path <- "../data/mini_fragments.tsv.gz" #nolint
  raw_fragments <- write_fragments_memory(open_fragments_10x(in_path), compress = FALSE)
  
  write_fragments_hdf5(raw_fragments, file.path(dir, "file.h5"), "unpacked", compress = FALSE, gzip_level = 4L)
  write_fragments_hdf5(raw_fragments, file.path(dir, "file.h5"), "packed", compress = TRUE, gzip_level = 4L)
  
  unpacked <- open_fragments_hdf5(file.path(dir, "file.h5"), "unpacked")
  packed <- open_fragments_hdf5(file.path(dir, "file.h5"), "packed")
  
  expect_identical(raw_fragments, write_fragments_memory(unpacked, compress = FALSE))
  expect_identical(raw_fragments, write_fragments_memory(packed, compress = FALSE))
})

test_that("H5 overwrite works", {
  dir <- withr::local_tempdir()
  frags1 <- tibble::tibble(
    chr = c("chr1", "chr1", "chr2", "chr3", "chr3"),
    start = c(10,20,30,40,50),
    end = start + 5,
    cell_id = 1:5
  )
  frags2 <- tibble::tibble(
    chr = c("chr1", "chr1", "chr2", "chr3", "chr4"),
    start = 10 + c(10,20,30,40,50),
    end = start + 5,
    cell_id = 1:5
  )
  f1 <- as(frags1, "IterableFragments")
  f2 <- as(frags2, "IterableFragments")
  
  write_fragments_hdf5(f1, file.path(dir, "overwrite.h5"), "frags")
  # writing without "overwrite" set should result in an error, and no data changed
  expect_error({
    write_fragments_hdf5(f2, file.path(dir, "overwrite.h5"), "frags")
  })
  expect_identical(
    as.data.frame(f1),
    open_fragments_hdf5(file.path(dir, "overwrite.h5"), "frags") %>%
      as.data.frame()
  )
  # writing with "overwrite" set should run and result in data change
  rlang::reset_message_verbosity("hdf5_overwrite")
  expect_message(
    write_fragments_hdf5(f2, file.path(dir, "overwrite.h5"), "frags", overwrite=TRUE),
    "dataset does not free old storage"
  )

  expect_identical(
    as.data.frame(f2),
    open_fragments_hdf5(file.path(dir, "overwrite.h5"), "frags") %>%
      as.data.frame()
  )
  # Overwriting from the same source data should work
  expect_identical(
    as.data.frame(select_cells(f2, c(1,3,5))),
    open_fragments_hdf5(file.path(dir, "overwrite.h5"), "frags") %>%
      select_cells(c(1,3,5)) %>%
      write_fragments_hdf5(file.path(dir, "overwrite.h5"), "frags", overwrite=TRUE) %>%
      as.data.frame()
  )
})

test_that("Dir overwrite works", {
  dir <- withr::local_tempdir()
  frags1 <- tibble::tibble(
    chr = c("chr1", "chr1", "chr2", "chr3", "chr3"),
    start = c(10,20,30,40,50),
    end = start + 5,
    cell_id = 1:5
  )
  frags2 <- tibble::tibble(
    chr = c("chr1", "chr1", "chr2", "chr3", "chr4"),
    start = 10 + c(10,20,30,40,50),
    end = start + 5,
    cell_id = 1:5
  )
  f1 <- as(frags1, "IterableFragments")
  f2 <- as(frags2, "IterableFragments")
  
  write_fragments_dir(f1, file.path(dir, "overwrite-frags"))
  # writing without "overwrite" set should result in an error, and no data changed
  expect_error({
    write_fragments_dir(f2, file.path(dir, "overwrite-frags"))
  })
  expect_identical(
    as.data.frame(f1),
    open_fragments_dir(file.path(dir, "overwrite-frags")) %>%
      as.data.frame()
  )
  # writing with "overwrite" set should run and result in data change
  write_fragments_dir(f2, file.path(dir, "overwrite-frags"), overwrite=TRUE)
  expect_identical(
    as.data.frame(f2),
    open_fragments_dir(file.path(dir, "overwrite-frags")) %>%
      as.data.frame()
  )
  # Overwriting from the same source data should work
  expect_identical(
    as.data.frame(select_cells(f2, c(1,3,5))),
    open_fragments_dir(file.path(dir, "overwrite-frags")) %>%
      select_cells(c(1,3,5)) %>%
      write_fragments_dir(file.path(dir, "overwrite-frags"), overwrite=TRUE) %>%
      as.data.frame()
  )
})

test_that("Regression test on multiple of 128 fragments", {
  frags_table <- tibble::tibble(
    chr = as.factor("chr1"),
    start = 1:(128*10),
    end = start+5L,
    cell_id = as.factor(seq_along(start) %% 5)
  )
  f <- as(frags_table, "IterableFragments") %>% write_fragments_memory()
  x <- as.data.frame(f) %>% tibble::as_tibble()
  expect_identical(frags_table, x)
})