# Copyright 2021 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

tibble_to_fragments <- function(x, chr_names, cell_names) {
  x %>%
    dplyr::mutate(
      chr = factor(chr_names[chr], levels = chr_names),
      cell_id = factor(cell_names[cell_id], levels = cell_names)
    ) %>%
    convert_to_fragments()
}

test_that("Basic insertion matrix succeeds", {
  # chr1 tests accuracy of start coordinate overlaps alone
  # chr2 tests accuracy of start+end coordinate overlaps combined
  # The ordering tests string matching of chromosome name & correct
  # sorting of output rows

  chr1_list <- list()
  for (i in 1:5) {
    chr1_list[[i]] <- tibble::tibble(
      cell_id = i - 1,
      start = rep(i:5, i),
      end = 1001 + i
    )
  }
  chr1_coords <- dplyr::bind_rows(chr1_list) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(chr = 2, cell_id = cell_id + 1)

  chr2_coords <- tibble::tribble(
    ~cell_id, ~start, ~end,
    0, 9, 21,
    1, 9, 20,
    2, 10, 21,
    3, 10, 20
  ) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(chr = 1, cell_id = cell_id + 1)

  raw_fragments <- tibble_to_fragments(
    dplyr::bind_rows(chr2_coords, chr1_coords),
    chr_names = c("chr2", "chr1"),
    cell_names = sprintf("cell%d", 1:5)
  )

  peaks <- list(
    chr = c("chr2", "chr1", "chr1", "chr1"),
    start = c(10, 3, 1002, 1004),
    end = c(20, 5, 1005, 1006)
  )
  res <- peak_matrix(
    raw_fragments,
    peaks
  ) %>% as("dgCMatrix")

  expect_s4_class(res, "dgCMatrix")

  answer <- matrix(c(
    0, 1, 1, 2, 0,
    2, 4, 6, 4, 0,
    0, 8, 9, 8, 0,
    0, 0, 0, 8, 5
  ), ncol = 4) %>% t()

  my_answer <- as.matrix(res)
  attr(my_answer, "dimnames") <- NULL

  expect_identical(my_answer, answer)

  answer_fragments <- matrix(c(
    0, 1, 1, 1, 0,
    2, 4, 6, 4, 0,
    0, 8, 9, 8, 0,
    0, 0, 0, 8, 5
  ), ncol=4) %>% t()
  res_fragments <- peak_matrix(
    raw_fragments,
    peaks,
    mode="fragments"
  ) %>% as("dgCMatrix") %>% as.matrix()
  attr(res_fragments, "dimnames") <- NULL
  expect_identical(res_fragments, answer_fragments)

  answer_overlaps <- matrix(c(
    1, 1, 1, 1, 0,
    4, 6, 6, 4, 0,
    0, 8, 9, 8, 5,
    0, 0, 0, 8, 5
  ), ncol=4) %>% t()
  res_overlaps <- peak_matrix(
    raw_fragments,
    peaks,
    mode="overlaps"
  ) %>% as("dgCMatrix") %>% as.matrix()
  attr(res_overlaps, "dimnames") <- NULL
  expect_identical(res_overlaps, answer_overlaps)
})


test_that("Out of range peaks work", {
  # chr1 tests having some peaks before the fragments start
  # chr2 tests having some peaks after the end of the fragments
  # chr3 tests having some fragments but no peaks
  # chr4 tests having some fragments, with a duplicated 1bp-wide peak
  chr1_list <- list()
  for (i in 1:5) {
    chr1_list[[i]] <- tibble::tibble(
      cell_id = i - 1,
      start = i:5 + 400,
      end = 1001 + i
    )
  }
  chr1_coords <- dplyr::bind_rows(chr1_list) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(chr = 2, cell_id = cell_id + 1)

  chr2_coords <- tibble::tribble(
    ~cell_id, ~start, ~end,
    0, 9, 21,
    1, 9, 20,
    2, 10, 21,
    3, 10, 20
  ) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(chr = 1, cell_id = cell_id + 1)

  chr3_coords <- dplyr::mutate(chr2_coords, chr = 3)
  chr4_coords <- dplyr::mutate(chr2_coords, chr = 4)

  raw_fragments <- tibble_to_fragments(
    dplyr::bind_rows(chr2_coords, chr1_coords, chr3_coords, chr4_coords),
    chr_names = c("chr2", "chr1", "chr3", "chr4"),
    cell_names = sprintf("cell%d", 1:5)
  )

  res <- peak_matrix(
    raw_fragments,
    list(
      chr = c("chr2", "chr2", "chr1", "chr1", "chr4", "chr4"),
      start = c(10, 1002, 3, 1004, 1, 1),
      end = c(20, 1005, 5, 1006, 2, 2)
    )
  ) %>% as("dgCMatrix")

  expect_s4_class(res, "dgCMatrix")

  answer <- matrix(c(
    0, 1, 1, 2, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 2, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0
  ), ncol = 6) %>% t()

  my_answer <- as.matrix(res)
  attr(my_answer, "dimnames") <- NULL

  expect_identical(my_answer, answer)
})



test_that("Basic tile matrix works", {
  # chr1 tests having some peaks before the fragments start
  # chr2 tests having some peaks after the end of the fragments
  # chr3 tests having some fragments but no peaks
  # chr4 tests having some fragments, with a duplicated 1bp-wide peak
  chr1_list <- list()
  for (i in 1:5) {
    chr1_list[[i]] <- tibble::tibble(
      cell_id = i - 1,
      start = i:5 + 400,
      end = 1001 + i
    )
  }
  chr1_coords <- dplyr::bind_rows(chr1_list) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(chr = 2, cell_id = cell_id + 1)

  chr2_coords <- tibble::tribble(
    ~cell_id, ~start, ~end,
    0, 9, 21,
    1, 9, 20,
    2, 10, 21,
    3, 10, 20
  ) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(chr = 1, cell_id = cell_id + 1)

  chr3_coords <- dplyr::mutate(chr2_coords, chr = 3)
  chr4_coords <- dplyr::mutate(chr2_coords, chr = 4)

  raw_fragments <- tibble_to_fragments(
    dplyr::bind_rows(chr2_coords, chr1_coords, chr3_coords, chr4_coords),
    chr_names = c("chr2", "chr1", "chr3", "chr4"),
    cell_names = sprintf("cell%d", 1:5)
  )

  res <- tile_matrix(
    raw_fragments,
    list(
      chr = c("chr2", "chr1", "chr1", "chr1"),
      start = c(10, 3, 402, 1004),
      end = c(20, 5, 405, 1006),
      tile_width = c(4, 1, 2, 2)
    )
  ) %>% as("dgCMatrix")

  expect_s4_class(res, "dgCMatrix")

  answer <- matrix(c(
    0, 0, 1, 1, 0, # 10-14
    0, 0, 0, 0, 0, # 14-18
    0, 1, 0, 1, 0, # 18-19
    0, 0, 0, 0, 0, # 3
    0, 0, 0, 0, 0, # 4
    2, 2, 1, 0, 0, # 402-403
    1, 1, 1, 1, 0, # 404
    0, 0, 0, 2, 1 # 1004-1006
  ), ncol = 8) %>% t()

  my_answer <- as.matrix(res)
  attr(my_answer, "dimnames") <- NULL

  expect_identical(my_answer, answer)

  res_fragments <- tile_matrix(
    raw_fragments,
    list(
      chr = c("chr2", "chr1", "chr1", "chr1"),
      start = c(10, 3, 402, 1004),
      end = c(20, 5, 405, 1006),
      tile_width = c(10, 1, 2, 2)
    ),
    mode="fragments"
  ) %>% as("dgCMatrix") %>% as.matrix()
  attr(res_fragments, "dimnames") <- NULL

  answer_fragments <- matrix(c(
    0, 1, 1, 1, 0, # 10-19
    0, 0, 0, 0, 0, # 3
    0, 0, 0, 0, 0, # 4
    2, 2, 1, 0, 0, # 402-403
    1, 1, 1, 1, 0, # 404
    0, 0, 0, 2, 1 # 1004-1006
  ), ncol = 6) %>% t()
  expect_identical(res_fragments, answer_fragments)
})


test_that("Subsetting peak matrix works", {
  f <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>% write_fragments_memory()
  n_peaks <- 100
  peak_width <- 1e6
  max_coord <- 180e6
  peaks <- tibble::tibble(
    chr = sample(paste0("chr", 1:5), n_peaks, replace=TRUE),
    start = as.integer(runif(n_peaks) * max_coord),
    end = start + max_coord
  )
  peaks <- peaks[order_ranges(peaks, chrNames(f)),]

  m <- peak_matrix(f, peaks)

  cols <- sample.int(ncol(m), 1000)
  rows <- sample.int(nrow(m), 50)

  expect_identical(
    m[rows, cols] %>% as("dgCMatrix"),
    as(m, "dgCMatrix")[rows, cols]
  )

  expect_identical(
    t(m)[cols, rows] %>% as("dgCMatrix"),
    t(as(m, "dgCMatrix")[rows, cols])
  )
})

# From https://stackoverflow.com/a/3791284
prime_sieve <- function(n)
{
   n <- as.integer(n)
   if(n > 1e8) stop("n too large")
   primes <- rep(TRUE, n)
   primes[1] <- FALSE
   last.prime <- 2L
   fsqr <- floor(sqrt(n))
   while (last.prime <= fsqr)
   {
      primes[seq.int(2L*last.prime, n, last.prime)] <- FALSE
      sel <- which(primes[(last.prime+1):(fsqr+1)])
      if(any(sel)){
        last.prime <- last.prime + min(sel)
      }else last.prime <- fsqr+1
   }
   return(primes)
}

test_that("Subsetting tile matrix works", {
  f <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>% write_fragments_memory()
  
  tile_width <- 1e6
  max_coord <- 180e6
  tiles <- tibble::tibble(
    chr = paste0("chr", 1:5),
    start = 0,
    end = max_coord,
    tile_width = tile_width
  )

  m <- tile_matrix(f, tiles, explicit_tile_names=TRUE)

  cols <- sample.int(ncol(m), 1000)
  rows <- sample.int(nrow(m), 50)

  # Test non-contiguous subset
  expect_identical(
    m[rows, cols] %>% as("dgCMatrix"),
    as(m, "dgCMatrix")[rows, cols]
  )
  expect_s4_class(m[rows, cols]@matrix, "PeakMatrix")

  expect_identical(
    t(m)[cols, rows] %>% as("dgCMatrix"),
    t(as(m, "dgCMatrix")[rows, cols])
  )
  expect_s4_class(t(m)[cols, rows]@matrix, "PeakMatrix")

  # Test contiguous subset (remove some prime-numbered tiles)
  primes <- matrix(which(prime_sieve(nrow(m))), nrow=2)[1,]
  rows_cont <- seq_len(nrow(m))[-primes]
  expect_identical(
    m[rows_cont, cols] %>% as("dgCMatrix"),
    as(m, "dgCMatrix")[rows_cont, cols]
  )
  expect_s4_class(m[rows_cont, cols], "TileMatrix")
  expect_gt(length(m[rows_cont, cols]@start), length(m@start))
  expect_lt(length(m[rows_cont, cols]@start), 0.2*length(rows_cont))

  expect_identical(
    t(m)[cols, rows_cont] %>% as("dgCMatrix"),
    t(as(m, "dgCMatrix")[rows_cont, cols])
  )
  expect_s4_class(t(m)[cols, rows_cont], "TileMatrix")
})