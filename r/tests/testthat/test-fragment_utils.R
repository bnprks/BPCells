# Copyright 2021 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

raw_fragments_to_tibble <- function(raw_fragments) {
  dplyr::bind_rows(raw_fragments$fragments, .id = "chr")
}

tibble_to_fragments <- function(x, chr_names, cell_names) {
  x %>%
    dplyr::mutate(
      chr = factor(chr_names[chr], levels = chr_names),
      cell_id = factor(cell_names[cell_id], levels = cell_names)
    ) %>%
    convert_to_fragments()
}

test_that("Chromosome name/index select works", { #nolint
  # Test lots of short chromosomes with index and name selection
  withr::local_seed(1258123)

  nfrags <- 10e3
  short_chrs_tibble <- tibble::tibble(
    chr = sample.int(500, nfrags, replace = TRUE),
    start = sample.int(10000, nfrags, replace = TRUE),
    end = start + sample.int(500, nfrags, replace = TRUE),
    cell_id = sample.int(500, nfrags, replace = TRUE)
  ) %>% dplyr::arrange(chr, start)

  short_chrs <- tibble_to_fragments(
    short_chrs_tibble,
    paste0("chr", 1:500),
    paste0("cell", 1:500)
  )

  chr_selection <- sort(sample.int(500, 250, replace = FALSE))

  short_chrs_ans <- as(short_chrs, "data.frame")
  short_chrs_ans <- short_chrs_ans[as.character(short_chrs_ans$chr) %in% paste0("chr", chr_selection),] %>%
    dplyr::mutate(chr = factor(as.character(chr), levels=paste0("chr", chr_selection)))

  short_chrs_2 <- select_chromosomes(short_chrs, chr_selection) %>%
    as("data.frame")

  short_chrs_3 <- select_chromosomes(short_chrs, paste0("chr", chr_selection)) %>%
    as("data.frame")

  expect_equal(short_chrs_ans, short_chrs_2, check.attributes = FALSE)
  expect_equal(short_chrs_ans, short_chrs_3, check.attributes = FALSE)
})

test_that("Chromosome select reordering works", {
  x1 <- open_fragments_10x("../data/mini_fragments.tsv.gz")
  x2 <- write_fragments_memory(x1)

  chr_order <- c("chr2", "chr1", "chr3")
  y1 <- select_chromosomes(x1, chr_order)
  y2 <- select_chromosomes(x2, chr_order)

  expect_identical(
    as.data.frame(y1) %>% dplyr::arrange(chr),
    as.data.frame(y2) %>% dplyr::arrange(chr)
  )
  expect_identical(
    as.data.frame(y2)$chr %>% as.character() %>% rle() %>% {.$values},
    chr_order
  )

  z1 <- select_chromosomes(x1, c(2,1,3))
  z2 <- select_chromosomes(x2, c(2,1,3))

  expect_identical(
    as.data.frame(z1) %>% dplyr::arrange(chr),
    as.data.frame(z2) %>% dplyr::arrange(chr)
  )
  expect_identical(
    as.data.frame(z2)$chr %>% as.character() %>% rle() %>% {.$values},
    chr_order
  )
})

test_that("Cell name/index select works", { #nolint
  # Test cells with index and name selection, and try explicitly making
  # chr 2 only contain cells that get filtered out, and
  # chr 3 only contain cells that get kept
  withr::local_seed(1258123)

  nfrags <- 1e4
  nchrs <- 20
  ncells <- 100
  cell_selection <- sample.int(ncells, ncells / 2, replace = FALSE)

  frags_tibble <- tibble::tibble(
    chr = sample.int(nchrs, nfrags, replace = TRUE),
    start = sample.int(10000, nfrags, replace = TRUE),
    end = start + sample.int(500, nfrags, replace = TRUE),
    cell_id = sample.int(ncells, nfrags, replace = TRUE)
  ) %>%
    dplyr::arrange(chr, start) %>%
    dplyr::filter(
      (chr == 2 & !(cell_id %in% cell_selection)) |
        (chr == 3 & cell_id %in% cell_selection) |
        !(chr %in% c(2, 3))
    )

  frags <- tibble_to_fragments(
    frags_tibble,
    paste0("chr", seq_len(nchrs)),
    paste0("cell", seq_len(ncells))
  )

  frags_ans <- tibble_to_fragments(
    frags_tibble %>%
      dplyr::filter(cell_id %in% cell_selection) %>%
      dplyr::mutate(cell_id = match(cell_id, cell_selection)),
    frags@chr_names,
    frags@cell_names[cell_selection]
  ) %>%
    write_fragments_memory() %>%
    as("data.frame")

  frags_2 <- select_cells(frags, cell_selection) %>%
    as("data.frame")
  frags_3 <- select_cells(frags, paste0("cell", cell_selection)) %>%
    as("data.frame")

  expect_identical(frags_ans, frags_2)
  expect_identical(frags_ans, frags_3)
})

test_that("Merge cells works", {
  nfrags <- 1e4
  nchrs <- 20
  ncells <- 100
  ngroups <- 10
  cell_groups <- sample.int(ngroups, ncells, replace = TRUE)

  frags_in <- tibble::tibble(
    chr = sample.int(nchrs, nfrags, replace = TRUE),
    start = sample.int(10000, nfrags, replace = TRUE),
    end = start + sample.int(500, nfrags, replace = TRUE),
    cell_id = sample.int(ncells, nfrags, replace = TRUE)
  ) %>%
    dplyr::arrange(chr, start, end, cell_id)

  frags_out <- tibble_to_fragments(
    frags_in,
    paste0("chr", seq_len(nchrs)),
    paste0("cell", seq_len(ncells))
  ) %>%
    merge_cells(as.factor(cell_groups)) %>%
    as("data.frame")

  frags_ans <- dplyr::mutate(frags_in, cell_id = cell_groups[cell_id]) %>%
    tibble_to_fragments(paste0("chr", seq_len(nchrs)), as.character(seq_len(ngroups))) %>%
    as("data.frame")
  expect_identical(frags_out, frags_ans)
})

test_that("GRanges conversion round-trips", {
  skip_if_not_installed("GenomicRanges")
  frags <- open_fragments_10x("../data/mini_fragments.tsv.gz") #nolint
  raw <- write_fragments_memory(frags)
  ranges <- as(raw, "GRanges")
  ranges2 <- as(frags, "GRanges")
  expect_identical(ranges, ranges2)
  raw2 <- write_fragments_memory(convert_to_fragments(ranges))
  expect_identical(raw, raw2)
})

test_that("Region select works", {
  skip_if_not_installed("GenomicRanges")
  frags <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>% #nolint
    write_fragments_memory()
  gfrags <- as(frags, "GRanges")

  # Include some overlapping regions for good measure
  regions <- tibble::tibble(
    chr = c("chr1", "chr5", "chr5", "chr2", "chr2"),
    start = c(1000000, 6373800, 6373699, 224956193, 225551362),
    end = c(1100000, 6375000, 15144614, 227996287, 230681382)
  )

  inclusive <- frags %>% select_regions(regions, invert_selection = FALSE)
  exclusive <- frags %>% select_regions(regions, invert_selection = TRUE)

  gregions <- as(regions, "GRanges")
  GenomicRanges::end(gregions) <- regions$end - 1

  ans_inclusive <- IRanges::subsetByOverlaps(gfrags, gregions)
  ans_exclusive <- IRanges::subsetByOverlaps(gfrags, gregions, invert = TRUE)

  expect_identical(as(inclusive, "GRanges"), ans_inclusive)
  expect_identical(as(exclusive, "GRanges"), ans_exclusive)
})

test_that("Region select regressions", {
  # To hit the bug, we need to do region select where:
  # - 1st we have a load with no hits (causes a seek)
  # - 2nd we have a load with a hit
  # - 3rd we have a load with no hit, but it's still in our region (causes repeat seek)

  # Separately, there was a bug in conversion from data.frame to IterableFragments that
  # was hit with this same input type. (The end_max calculation was completely wrong)
  load_1 <- tibble::tibble(
    chr = "chr1",
    start = rep.int(1L, 1050),
    end = 3L,
    cell_id = paste0("cell1.", seq_len(1050))
  )

  load_2 <- tibble::tibble(
    chr = "chr1",
    start = rep.int(10L, 100),
    end = 20L,
    cell_id = paste0("cell2.", seq_len(100))
  )

  load_3 <- tibble::tibble(
    chr = "chr1",
    start = rep.int(11L, 2000),
    end = 13L,
    cell_id = paste0("cell3.", seq_len(2000))
  )

  frags <- dplyr::bind_rows(load_1, load_2, load_3) %>% as("IterableFragments")

  expect_identical(
    frags,
    write_fragments_memory(frags, compress=FALSE)
  )

  expect_identical(
    select_regions(frags, "chr1:15-25") %>% tibble::as_tibble() %>% dplyr::mutate(
      chr = as.character(chr),
      cell_id = as.character(cell_id)
    ),
    load_2
  )
})

test_that("subset_lengths works", {
  frags <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>% #nolint
    write_fragments_memory()
  gfrags <- as(frags, "data.frame")

  min_size <- 20
  max_size <- 100

  lengths <- gfrags$end - gfrags$start
  expect_lt(min(lengths), min_size)
  expect_gt(max(lengths), max_size)

  ans <- gfrags[lengths >= min_size & lengths <= max_size,]
  row.names(ans) <- NULL

  subset <- subset_lengths(frags, min_size, max_size)

  expect_identical(as(subset, "data.frame"), ans)
})

test_that("Name prefix works", {
  frags1 <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>% #nolint
    write_fragments_memory()
  prefix <- "MyFancyPrefix##"
  frags2 <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>% #nolint
    prefix_cell_names(prefix) %>%
    write_fragments_memory()
  expect_identical(paste0(prefix, cellNames(frags1)), cellNames(frags2))
})

test_that("footprint works", {
  # Footprint in bases 21-27, making a pattern of 1,2,3,4,5,6,7
  # Note: Doesn't thoroughly test overlapping regions
  motif_positions <- list(chr = c("chr1", "chr2"), start = c(24, 2), end = c(60, 25), strand = c("+", "-"))
  frags1 <- tibble::tibble(
    chr = 1,
    #       fix start;  fix end;   weighted; finish
    start = c(rep(10, 7), 22:27, 24, 26, 23, 27),
    end = c(22:28, rep(40, 6), 26, 28, 26, 40),
    cell_id = c(rep(1, 7), rep(3, 6), 5, 7, 9, 11)
  )
  frags2 <- frags1 %>%
    dplyr::mutate(chr = 2, cell_id = cell_id + 1)

  frags <- tibble_to_fragments(
    dplyr::bind_rows(frags1, frags2) %>% dplyr::arrange(chr, start),
    c("chr1", "chr2"),
    paste0("c", 1:12)
  )

  expected <- tibble::tibble(
    group = c(rep("grp1", 9), rep("grp2", 9)),
    position = c(-4:4, 4:-4),
    count = rep(c(0, 1:7, 0), 2),
    enrichment = NA
  ) %>% dplyr::arrange(group, position)

  for (i in 0:4) {
    res <- footprint(
      frags,
      motif_positions,
      cell_groups = rep(c("grp1", "grp2"), 6),
      cell_weights = c(1, 1, 1, 1, 2, 2, 4, 4, 1, 1, 1, 1),
      flank = i,
      normalization_width = 0L
    )
    expect_identical(
      res %>% dplyr::arrange(group, position),
      expected %>% dplyr::filter(abs(position) <= i)
    )
  }
})

test_that("Concatenate seek works", {
  # Regression test against incorrect tile/peak matrix subset with concatenated fragments
  frags_a <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>%
    prefix_cell_names("A") %>%
    write_fragments_memory()
  frags_b <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>%
    prefix_cell_names("B") %>%
    write_fragments_memory()

  frags_merge <- c(frags_a, frags_b)

  fragment_counts <- nucleosome_counts(frags_merge)$nFrags
  keeper_cells <- order(-fragment_counts)[1:2]
  frags_filt <- select_cells(frags_merge, keeper_cells)

  tiles <- tile_matrix(frags_filt, list(chr="chr1", start=0L, end=248956422L, tile_width=100000), explicit_tile_names=TRUE)
  tile_counts <- matrix_stats(tiles, row_stats = "nonzero")$row_stats
  keeper_tiles <- order(-tile_counts)[1:3]
  x <- as(tiles, "dgCMatrix")[keeper_tiles,]
  y <- as(tiles[keeper_tiles,], "dgCMatrix")
  expect_identical(x, y)
})

test_that("Concatenate/merge different chromosome sets works", {
  frags_a <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>%
    prefix_cell_names("A") %>%
    select_chromosomes(c("chr1", "chr3", "chr2", "chr5")) %>%
    write_fragments_memory()
  frags_b <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>%
    prefix_cell_names("B") %>%
    select_chromosomes(c("chr1", "chr2", "chr5", "chr4", "chr9")) %>%
    write_fragments_memory()

  frags_merge <- c(frags_a, frags_b)

  expect_identical(
    as.data.frame(frags_merge) ,
    dplyr::bind_rows(as.data.frame(frags_a), as.data.frame(frags_b)) %>%
      dplyr::arrange(chr, start, cell_id)
  )
})

test_that("Generic methods work", {
  frags <- open_fragments_10x("../data/mini_fragments.tsv.gz") %>% #nolint
    write_fragments_memory()


  first_half_cells <- sort(sample.int(length(cellNames(frags)), length(cellNames(frags)) / 2))
  second_half_cells <- seq_along(cellNames(frags))[-first_half_cells]

  dir <- withr::local_tempdir()

  ident_transforms <- list(
    write_memory = write_fragments_memory(frags, compress = TRUE),
    write_memory_unpacked = write_fragments_memory(frags, compress = FALSE),
    write_dir = write_fragments_dir(frags, file.path(dir, "fragments_identical")),
    write_h5 = write_fragments_hdf5(frags, file.path(dir, "fragments_identical.h5")),
    # write_bed = write_fragments_10x(frags, file.path(dir, "fragments_identical.tsv")), # Ignored because this doesn't have cell/chr Names
    shift = shift_fragments(frags),
    lengthSelect = subset_lengths(frags),
    chrSelectName = select_chromosomes(frags, chrNames(frags)),
    chrSelectIdx = select_chromosomes(frags, seq_along(chrNames(frags))),
    cellSelectName = select_cells(frags, cellNames(frags)),
    cellSelectIdx = select_cells(frags, seq_along(cellNames(frags))),
    cellMerge = merge_cells(frags, factor(cellNames(frags), levels = cellNames(frags))),
    chrRename = new("ChrRename", fragments = frags, chr_names = chrNames(frags)),
    cellRename = new("CellRename", fragments = frags, cell_names = cellNames(frags)),
    cellPrefix = prefix_cell_names(frags, ""),
    regionSelect = select_regions(frags, list(chr = chrNames(frags), start = rep_len(0, length(chrNames(frags))), end = rep_len(1e9, length(chrNames(frags))))),
    merge = select_cells(c(select_cells(frags, first_half_cells), select_cells(frags, second_half_cells)), cellNames(frags))
  )

  for (i in seq_along(ident_transforms)) {
    trans <- ident_transforms[[i]]
    
    # Test that basic operations work
    short_description(trans)
    expect_identical(chrNames(trans), chrNames(frags))
    expect_identical(cellNames(trans), cellNames(frags))
    expect_identical(write_fragments_memory(trans), frags)

    # Test that selecting a non-existent chromosome is fine
    expect_identical(
      as.data.frame(frags) %>% dplyr::mutate(chr=as.character(chr)), 
      trans %>% 
        select_chromosomes(c(chrNames(frags)[1:5], "NOT_A_CHROMOSOME", chrNames(frags)[-(1:5)])) %>%
        as.data.frame() %>%
        dplyr::mutate(chr=as.character(chr))
    )

    # Test that renaming chromosomes/cells works
    expected_rename <- as.data.frame(trans)
    levels(expected_rename$chr) <- paste0("new-", levels(expected_rename$chr))
    levels(expected_rename$cell_id) <- paste0("new-", levels(expected_rename$cell_id))
    chrNames(trans) <- paste0("new-", chrNames(frags))
    cellNames(trans) <- paste0("new-", cellNames(frags))
    expect_identical(chrNames(trans), paste0("new-", chrNames(frags)))
    expect_identical(cellNames(trans), paste0("new-", cellNames(frags)))
    n <- write_fragments_memory(trans)
    expect_identical(chrNames(n), paste0("new-", chrNames(frags)))
    expect_identical(cellNames(n), paste0("new-", cellNames(frags)))
    expect_identical(as.data.frame(n), expected_rename)
    # Reset cell + chr names
    chrNames(trans) <- chrNames(frags)
    cellNames(trans) <- cellNames(frags)

    # Test that garbage collection after creating the iterator doesn't cause issues
    # Create the C++ iterator
    it <- iterate_fragments(trans) %>%
      iterate_chr_index_select_cpp(seq_along(chrNames(frags)) - 1L) %>%
      iterate_cell_index_select_cpp(seq_along(cellNames(frags)) - 1L) %>%
      iterate_shift_cpp(0L, 0L)
    expect_type(it, "externalptr")
    # Delete any trailing XPtr references from R
    gc()
    # Check that the C++ iterator still works
    res <- do.call(new, c("PackedMemFragments", write_packed_fragments_cpp(it)))
    expect_identical(res, frags)
  }
})
