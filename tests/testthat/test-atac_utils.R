
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
