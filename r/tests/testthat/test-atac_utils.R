# Copyright 2022 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.


test_that("range_overlaps works", {
  a <- tibble::tibble(
    chr = c("chr1", "chr2", "chr1"),
    start = c(1, 5, 10),
    end = c(20, 20, 20)
  )
  b <- tibble::tibble(
    chr = c("chr2", "chr2", "chr1", "chr1", "chr1", "chr1"),
    start = c(1, 25, 2, 11, 0, 2),
    end = c(10, 30, 20, 25, 25, 4)
  )

  expected <- tibble::tibble(
    from = as.integer(c(1, 1, 1, 1, 2, 3, 3, 3)),
    to = as.integer(c(3, 4, 5, 6, 1, 3, 4, 5))
  )

  expect_identical(range_overlaps(a, b), expected)
})

test_that("tile_ranges works", {
  frags <- convert_to_fragments(tibble::tibble(
    chr = paste0("chr", 1:10),
    start = 0,
    end = 5,
    cell_id = "cell1"
  ))
  chr_sizes <- tibble::tibble(
    chr = paste0("chr", 1:10),
    start = 0,
    end = 1000 * 1:10,
    tile_width = 300
  )
  expect_warning(
    tile_mat <- tile_matrix(frags, chr_sizes, explicit_tile_names = TRUE),
    "Tiles given out of order"
  )

  selection1 <- sample.int(nrow(tile_mat))
  expect_identical(
    normalize_ranges(rownames(tile_mat)[selection1]) %>% dplyr::mutate(chr=as.character(chr)),
    tile_ranges(tile_mat, selection1) %>% dplyr::mutate(chr=as.character(chr))
  )

  selection2 <- as.logical(rbinom(nrow(tile_mat), 1, 0.5))
  expect_identical(
    normalize_ranges(rownames(tile_mat)[selection2]) %>% dplyr::mutate(chr=as.character(chr)),
    tile_ranges(tile_mat, selection2) %>% dplyr::mutate(chr=as.character(chr))
  )
})

test_that("merge_peaks works", {
  peaks1 <- tibble::tibble(
    chr = "chr1",
    start = as.integer(1:10),
    end = start + 2L
  )
  res1 <- tibble::tibble(
    chr = "chr1",
    start = as.integer(c(1, 3, 5, 7, 9)),
    end = start + 2L
  )
  expect_identical(merge_peaks_iterative(peaks1), res1)

  peaks2 <- tibble::tibble(
    chr = "chr1",
    start = as.integer(c(1, 22, 11, 22, 2, 3, 13, 1, 4, 21, 12)),
    end = start + 2L
  )
  res2 <- tibble::tibble(
    chr = "chr1",
    start = as.integer(c(1, 22, 11, 3, 13)),
    end = start + 2L
  )
  expect_identical(merge_peaks_iterative(peaks2), res2)
})

test_that("write_insertion_bedgraph works", {
  dir <- withr::local_tempdir()

  # These parameters should make it virtually certain we have
  # counts of 0, 1, and >1 represented
  chr1 <- tibble::tibble(
    chr = "chr1",
    start = sort(sample.int(1000, 1000, replace=TRUE)),
    end = start + sample.int(150, 1000, replace=TRUE),
    cell_id = sample(LETTERS, 1000, replace=TRUE)
  )
  chr2 <- tibble::tibble(
    chr = "chr2",
    start = sort(sample.int(1000, 1000, replace=TRUE)),
    end = start + sample.int(150, 1000, replace=TRUE),
    cell_id = sample(LETTERS, 1000, replace=TRUE)
  )
  frags <- convert_to_fragments(dplyr::bind_rows(chr1, chr2))

  cell_groups <- dplyr::if_else(cellNames(frags) %in% c("A", "E", "I", "O", "U"), "vowel", "consonant")

  frag_table <- dplyr::bind_rows(chr1, chr2) %>%
    dplyr::mutate(cell_group = dplyr::if_else(cell_id %in% c("A", "E", "I", "O", "U"), "vowel", "consonant"))

  coverage_start <- frag_table %>%
    dplyr::mutate(end=start+1L) %>%
    dplyr::group_by(chr, start, end, cell_group) %>%
    dplyr::summarize(value = dplyr::n(), .groups="drop")
  
  coverage_end <- frag_table %>%
    dplyr::mutate(start=end-1L) %>%
    dplyr::group_by(chr, start, end, cell_group) %>%
    dplyr::summarize(value = dplyr::n(), .groups="drop")
  
  coverage <- dplyr::bind_rows(coverage_start, coverage_end) %>%
    dplyr::group_by(chr, start, end, cell_group) %>%
    dplyr::summarize(value = sum(value), .groups="drop")

  answers <- list(
    "start_only" = coverage_start,
    "end_only" = coverage_end,
    "both" = coverage
  )


  # Test start + end insertions
  for (mode in c("start_only", "end_only", "both")) {
    write_insertion_bedgraph(
      frags,
      c("vowel"=file.path(dir, "vowel.bg"), "consonant"=file.path(dir, "consonant.bg.gz")),
      cell_groups,
      mode
    )

    # Check for gzip encoding
    expect_false(readBin(file.path(dir, "vowel.bg"), "integer", size=2, endian="big") == 0x1f8b)
    expect_true(readBin(file.path(dir, "consonant.bg.gz"), "integer", size=2, endian="big") == 0x1f8b)

    result_vowel <- readr::read_tsv(file.path(dir, "vowel.bg"), col_names=c("chr", "start", "end", "value"), col_types="ciii")
    result_consonant <- readr::read_tsv(file.path(dir, "consonant.bg.gz"), col_names=c("chr", "start", "end", "value"), col_types="ciii")

    expect_identical(
      dplyr::filter(answers[[mode]], cell_group == "vowel") %>% dplyr::select(-c(cell_group)) %>% data.frame(), 
      data.frame(result_vowel)
    )
    expect_identical(
      dplyr::filter(answers[[mode]], cell_group == "consonant") %>% dplyr::select(-c(cell_group)) %>% data.frame(), 
      data.frame(result_consonant)
    )
  }
})