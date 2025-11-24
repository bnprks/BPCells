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

test_that("range_distance_to_nearest works", {
  frags <- tibble::tibble(
    chr = "chr1",
    start = seq(10, 410, 100),
    end = start + 50,
    strand = "+"
  )
  res <- tibble::tibble(
    upstream = c(Inf, rep(51, 4)),
    downstream = c(rep(51, 4), Inf)
  )
  expect_identical(
    range_distance_to_nearest(frags),
    res
  )
  frags_with_nested <- frags %>% 
    tibble::add_row(chr = "chr1", start = 11, end = 20, strand = "+")
  res_with_nested <- res %>%
    tibble::add_row(upstream = 0, downstream = 0)
  expect_identical(
    range_distance_to_nearest(frags_with_nested),
    res_with_nested
  )
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

  coverage_start_single_group <- dplyr::bind_rows(chr1, chr2) %>%
    dplyr::mutate(end=start+1L) %>%
    dplyr::group_by(chr, start, end) %>%
    dplyr::summarize(value = dplyr::n(), .groups="drop")

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
    result_vowel <- readr::read_tsv(file.path(dir, "vowel.bg"), col_names=c("chr", "start", "end", "value"), col_types="ciid")
    result_consonant <- readr::read_tsv(file.path(dir, "consonant.bg.gz"), col_names=c("chr", "start", "end", "value"), col_types="ciid")

    expect_equal(
      dplyr::filter(answers[[mode]], cell_group == "vowel") %>% dplyr::select(-c(cell_group)) %>% data.frame(), 
      data.frame(result_vowel)
    )
    expect_equal(
      dplyr::filter(answers[[mode]], cell_group == "consonant") %>% dplyr::select(-c(cell_group)) %>% data.frame(), 
      data.frame(result_consonant)
    )
  }
  # Test single group
  write_insertion_bedgraph(frags, file.path(dir, "all.bg"), insertion_mode = "start")
  expect_equal(
    readr::read_tsv(file.path(dir, "all.bg"), col_names = c("chr", "start", "end", "value"), col_types = "ciid") %>% 
      data.frame(),
    coverage_start_single_group %>% data.frame()
  )
  # Define chromosome sizes for tiling
  chrom_sizes <- tibble::tibble(
    chr = c("chr1", "chr2"),
    start = c(0, 0),
    end = c(400, 400)
  )
  
  # Test basic functionality with tile_width = 100
  tile_width <- 100
  write_insertion_bedgraph(
    frags,
    c("vowel" = file.path(dir, "vowel_tiled.bg"), "consonant" = file.path(dir, "consonant_tiled.bg.gz")),
    cell_groups,
    insertion_mode = "both",
    tile_width = tile_width,
    normalization_method = "none",
    chrom_sizes = chrom_sizes
  )
  
  # Check files were created
  expect_true(file.exists(file.path(dir, "vowel_tiled.bg")))
  expect_true(file.exists(file.path(dir, "consonant_tiled.bg.gz")))
  
  # Read and validate basic output format - appears to use bedgraph format (chr, start, end, value)
  result_vowel <- readr::read_tsv(file.path(dir, "vowel_tiled.bg"), col_names=c("chr", "start", "end", "value"), col_types="ciid", show_col_types=FALSE)
  result_consonant <- readr::read_tsv(file.path(dir, "consonant_tiled.bg.gz"), col_names=c("chr", "start", "end", "value"), col_types="ciid", show_col_types=FALSE)
  
  # Check basic properties (files should have content and valid numeric values)
  expect_true(nrow(result_vowel) > 0)
  expect_true(nrow(result_consonant) > 0)
  expect_true(all(is.finite(result_vowel$value)))
  expect_true(all(is.finite(result_consonant$value)))
  
  # Test different insertion modes
  
  for (mode in c("start_only", "end_only", "both")) {
    vowel_cells <- which(cell_groups == "vowel")
    frags_vowel <- select_cells(frags, vowel_cells)
    write_insertion_bedgraph(
      frags_vowel,
      file.path(dir, paste0("vowel_", mode, ".bg")),
      insertion_mode = mode,
      tile_width = tile_width,
      normalization_method = "none",
      chrom_sizes = chrom_sizes
    )
    
    result <- readr::read_tsv(
      file.path(dir, paste0("vowel_", mode, ".bg")), 
      col_names=c("chr", "start", "end", "value"), 
      col_types="ciid", show_col_types=FALSE
    )
    expect_true(nrow(result) > 0)
    expect_true(all(is.finite(result$value)))
  }
  
  # Test different normalization methods
  for (norm_method in c("none", "n_cells", "cpm")) {
    write_insertion_bedgraph(
      frags,
      file.path(dir, paste0("norm_", norm_method, ".bg")),
      normalization_method = norm_method,
      tile_width = tile_width,
      chrom_sizes = chrom_sizes
    )
    
    result <- readr::read_tsv(
      file.path(dir, paste0("norm_", norm_method, ".bg")), 
      col_names=c("chr", "start", "end", "value"), 
      col_types="ciid", 
      show_col_types=FALSE
    )
    expect_true(nrow(result) > 0) 
    expect_true(all(is.finite(result$value)))
  }
  
  # Test single group (no cell_groups specified)
  write_insertion_bedgraph(
    frags, 
    file.path(dir, "single_group.bg"),
    tile_width = tile_width,
    chrom_sizes = chrom_sizes
  )
  
  result_single <- readr::read_tsv(
    file.path(dir, "single_group.bg"), 
    col_names=c("chr", "start", "end", "value"), 
    col_types="ciid", show_col_types=FALSE
    )
  expect_true(nrow(result_single) > 0)
  expect_true(all(is.finite(result_single$value)))
  
  # Test error when chrom_sizes is NULL and tile_width != 1
  expect_error(
    write_insertion_bedgraph(frags, file.path(dir, "error_test.bg"), tile_width = 100, chrom_sizes = NULL),
    "If chrom_sizes is NULL, then tile_width must be 1"
  )
  
  # Test with tile_width = 1 and NULL chrom_sizes  
  write_insertion_bedgraph(
    frags,
    file.path(dir, "tile_width_1.bg"),
    tile_width = 1,
    chrom_sizes = NULL
  )
})

test_that("write_insertion_bed works", {
  dir <- withr::local_tempdir()
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
  frags_vowel <- convert_to_fragments(dplyr::filter(dplyr::bind_rows(chr1, chr2), cell_id %in% c("A", "E", "I", "O", "U")))
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
  write_insertion_bed(
    fragments = frags_vowel,
    path = c("vowel" = file.path(dir, "vowel_only.bed")),
    cell_groups = rep("vowel", length(cellNames(frags_vowel))),
    insertion_mode = "start_only",
    verbose = FALSE,
    threads = 1
  )
  write_insertion_bed(
    fragments = frags_vowel,
    path = c("vowel" = file.path(dir, "no_labels.bed")),
    insertion_mode = "start_only",
    verbose = FALSE,
    threads = 1
  )
  ## check with all insertions in same group
  expect_identical(
    dplyr::filter(coverage_start, cell_group == "vowel") %>% tidyr::uncount(value) %>% dplyr::select(-c(cell_group))  %>% data.frame(), 
    data.frame(readr::read_tsv(file.path(dir, "vowel_only.bed"), col_names = c("chr", "start", "end"), col_types = "cii")),
  )
  ## check with all insertions in same group without labels
  expect_identical(
    dplyr::filter(coverage_start, cell_group == "vowel") %>% tidyr::uncount(value) %>% dplyr::select(-c(cell_group)) %>% data.frame(), 
    data.frame(readr::read_tsv(file.path(dir, "no_labels.bed"), col_names = c("chr", "start", "end"), col_types = "cii"))
  )
  # test across all modes and single/multiple thread counts
  for (mode in c("start_only", "end_only", "both")) {
    for (threads in c(1, 2)) {
      write_insertion_bed(
        fragments = frags,
        path = c("vowel" = file.path(dir, "vowel.bed"), "consonant" = file.path(dir, "consonant.bed.gz")),
        cell_groups = cell_groups,
        insertion_mode = mode,
        verbose = FALSE,
        threads = threads
      )
      # read in the bed.gz files
      result <- list(
        "vowel" = readr::read_tsv(file.path(dir, "vowel.bed"), col_names = c("chr", "start", "end"), col_types="cii"),
        "consonant" = readr::read_tsv(file.path(dir, "consonant.bed.gz"), col_names = c("chr", "start", "end"), col_types="cii")
      )
      expect_identical(
        dplyr::filter(answers[[mode]], cell_group == "vowel") %>% tidyr::uncount(value) %>% 
          dplyr::select(-c(cell_group)) %>% data.frame(), 
        data.frame(result[["vowel"]])
      )
      expect_identical(
        dplyr::filter(answers[[mode]], cell_group == "consonant") %>% tidyr::uncount(value) %>% 
          dplyr::select(-c(cell_group)) %>% data.frame(), 
        data.frame(result[["consonant"]])
      )
    }
  }
})


test_that("CPM normalization scales each value by 1e6 / total_fragments_in_group", {
  dir <- withr::local_tempdir()
  frag_tbl <- tibble::tibble(
    chr     = "chr1",
    start   = c(0L, 1L, 0L, 2L, 2L),   # three frags in grpA, two in grpB
    end     = start + 1L,
    cell_id = c("A",  "A",  "A", "B", "B")
  )
  frags <- convert_to_fragments(frag_tbl)

  cell_groups <- dplyr::if_else(cellNames(frags) == "A", "grpA", "grpB")

  # Write CPM-normalised bedGraphs
  write_insertion_bedgraph(
    frags,
    c("grpA" = file.path(dir, "grpA.bg"),
      "grpB" = file.path(dir, "grpB.bg")),
    cell_groups       = cell_groups,
    insertion_mode    = "start_only",
    normalization_method = "cpm",
    tile_width        = 1,
    chrom_sizes       = NULL
  )

  # Expected values
  counts <- frag_tbl %>%
    dplyr::mutate(group = dplyr::if_else(cell_id == "A", "grpA", "grpB"),
                  end   = start + 1L) %>%
    dplyr::group_by(chr, start, end, group) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  totals <- counts %>% dplyr::count(group, wt = n, name = "tot")

  exp_tbl <- counts %>%
    dplyr::left_join(totals, by = "group") %>%
    dplyr::mutate(value = n * 1e6 / tot) %>%
    dplyr::select(chr, start, end, value, group) %>%
    dplyr::arrange(group, start)

  res_tbl <- dplyr::bind_rows(
    readr::read_tsv(file.path(dir, "grpA.bg"),
      col_names = c("chr", "start", "end", "value"), col_types = "ciid") %>%
      dplyr::mutate(group = "grpA"),
    readr::read_tsv(file.path(dir, "grpB.bg"),
      col_names = c("chr", "start", "end", "value"), col_types = "ciid") %>%
      dplyr::mutate(group = "grpB")
  ) %>% dplyr::arrange(group, start)

  expect_equal(res_tbl, exp_tbl, tolerance = 1e-6)
})

test_that("NCells normalization divides by number of cells in each group", {
  dir <- withr::local_tempdir()
  frag_tbl <- tibble::tibble(
    chr     = "chr1",
    start   = c(0L, 0L, 1L, 2L, 2L),
    end     = start + 1L,
    cell_id = c("A1", "A2", "A1", "B1", "B2")
  )
  frags <- convert_to_fragments(frag_tbl)

  cell_groups <- dplyr::case_when(
    grepl("^A", cellNames(frags)) ~ "grpA",
    TRUE                          ~ "grpB"
  )

  write_insertion_bedgraph(
    frags,
    c("grpA" = file.path(dir, "grpA.bg"),
      "grpB" = file.path(dir, "grpB.bg")),
    cell_groups = cell_groups,
    insertion_mode = "start_only",
    normalization_method = "n_cells",
    tile_width = 1,
    chrom_sizes = NULL
  )

  # Expected values
  cell_map <- tibble::tibble(cell_id = cellNames(frags),
                             n_cells = dplyr::n_distinct(cellNames(frags)),
                             group   = cell_groups) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

  counts <- frag_tbl %>%
    dplyr::mutate(group = dplyr::if_else(grepl("^A", cell_id), "grpA", "grpB"),
                  end   = start + 1L) %>%
    dplyr::group_by(chr, start, end, group) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::left_join(cell_map, by = "group") %>%
    dplyr::mutate(value = n / n_cells) %>%
    dplyr::select(chr, start, end, value, group) %>%
    dplyr::arrange(group, start)

  res_tbl <- dplyr::bind_rows(
    readr::read_tsv(file.path(dir, "grpA.bg"),
      col_names = c("chr", "start", "end", "value"), col_types = "ciid") %>%
      dplyr::mutate(group = "grpA"),
    readr::read_tsv(file.path(dir, "grpB.bg"),
      col_names = c("chr", "start", "end", "value"), col_types = "ciid") %>%
      dplyr::mutate(group = "grpB")
  ) %>% dplyr::arrange(group, start)

  expect_equal(res_tbl, counts, tolerance = 1e-6)
})

test_that("Tiles are clipped to chromosome end when chrom_sizes < tile_width * n_tiles", {
  dir <- withr::local_tempdir()

  frag_tbl <- tibble::tibble(
    chr     = "chr1",
    start   = c(0L, 4L),
    end     = start + 1L,
    cell_id = "A"
  )
  frags <- convert_to_fragments(frag_tbl)

  chrom_sizes <- tibble::tibble(chr = "chr1", start = 0L, end = 5L)
  tile_width  <- 3 

  write_insertion_bedgraph(
    frags,
    file.path(dir, "out.bg"),
    insertion_mode    = "start_only",
    tile_width        = tile_width,
    normalization_method = "none",
    chrom_sizes       = chrom_sizes
  )

  res <- readr::read_tsv(
    file.path(dir, "out.bg"),
    col_names = c("chr", "start", "end", "value"),
    col_types = "ciid",
    show_col_types = FALSE
  ) %>% dplyr::arrange(start)

  exp <- tibble::tibble(
    chr   = "chr1",
    start = c(0L, 3L),
    end   = c(3L, 5L),
    value = c(1L, 1L)
  )

  expect_equal(as.data.frame(res), as.data.frame(exp))
})

test_that("macs_e2e_works", {
  # only run if macs2 or macs3 is installed
  if ((suppressWarnings((system2("macs3", args = "--version", stdout = FALSE, stderr = FALSE) != 0)) ||
      (!grepl("macs", system2("macs3", args = "--version", stdout = TRUE))))) {
    if ((suppressWarnings((system2("macs2", args = "--version", stdout = FALSE, stderr = FALSE) != 0)) ||
        (!grepl("macs", system2("macs2", args = "--version", stdout = TRUE))))) {
      skip("macs2 not installed")
    } else macs_executable <- "macs2"
  } else macs_executable <- "macs3"
  dir <- withr::local_tempdir()
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
  cell_groups <- dplyr::if_else(cellNames(frags) %in% c("A", "E", "I", "O", "U"), "v o w e l", "consonant;")
  # call each step seperately and hold in memory
  macs_prep <- call_peaks_macs(
    fragments = frags,
    cell_groups = cell_groups,
    effective_genome_size = 2.9e9,
    path = dir,
    insertion_mode = "both",
    step = "prep-inputs",
    macs_executable = macs_executable,
    verbose = FALSE,
    threads = 2
  )
  # Check to see if bed/shell files are created
  for (cluster in c("v o w e l", "consonant;")) {
    expect_true(file.exists(file.path(dir, "input", paste0(cluster, ".bed.gz"))))
    expect_true(file.exists(file.path(dir, "input", paste0(cluster, ".sh"))))
  }
  
  # call macs using the prepared inputs
  macs_call <- call_peaks_macs(
    fragments = frags,
    cell_groups = cell_groups,
    effective_genome_size = 2.9e9,
    path = dir,
    insertion_mode = "both",
    step = "run-macs",
    macs_executable = macs_executable,
    verbose = FALSE,
    threads = 2
  )
  # Check to see if the output files are created
  for (cluster in c("v o w e l", "consonant;")) {
    expect_true(file.exists(file.path(dir, "output", cluster, paste0(cluster, "_peaks.narrowPeak"))))
    expect_true(file.exists(file.path(dir, "output", cluster, paste0(cluster, "_peaks.xls"))))
    expect_true(file.exists(file.path(dir, "output", cluster, paste0(cluster, "_summits.bed"))))
  }
  # Read in the outputs
  macs_read <- call_peaks_macs(
    fragments = frags,
    cell_groups = cell_groups,
    effective_genome_size = 2.9e9,
    path = dir,
    insertion_mode = "both",
    step = "read-outputs",
    macs_executable = macs_executable,
    verbose = FALSE,
    threads = 2
  )
  # Check length to see if the same number of clusters are returned
  expect_equal(length(unique(macs_read$group)), length(unique(cell_groups)))
  macs_read_full_pipeline <- call_peaks_macs(
    fragments = frags,
    cell_groups = cell_groups,
    effective_genome_size = 2.9e9,
    path = dir,
    insertion_mode = "both",
    step = "all",
    macs_executable = macs_executable,
    verbose = FALSE,
    threads = 2
  )
  # Make sure the outputs are the same
  expect_equal(macs_read, macs_read_full_pipeline)
})

test_that("macs errors print when running in parallel", {
  # The dummy macs script only works on a unix setup, and windows currently doesn't support
  # multi-threading anyhow
  skip_on_os("windows")
  dir <- withr::local_tempdir()

  # Make a dummy macs script that will register as valid but won't actually run
  bad_macs_path <- file.path(dir, "bad_macs.sh")
  writeLines(c(
    "#!/bin/bash",
    "if [[ $1 == '--version' ]]; then",
    "    echo 'Bad macs demo'",
    "else",
    "    exit 1",
    "fi"
  ), bad_macs_path)
  Sys.chmod(bad_macs_path, "0700")

  frags <- tibble::tibble(chr="chr1", start=1:10, end=start+5, cell_id=rep(c("a","b"), 5)) %>% convert_to_fragments() 
  call_peaks_macs(
    frags, 
    path=file.path(dir, "macs-test"), 
    cell_groups=c(a="a_group", b="b_group"), 
    step="prep-inputs",
    macs_executable=bad_macs_path, 
    threads=2
  )
  expect_error({
    call_peaks_macs(
      frags,
      path=file.path(dir, "macs-test"), 
      cell_groups=c(a="a_group", b="b_group"), 
      step="run-macs",
      macs_executable=bad_macs_path, 
      threads=2
    )
  }, "MACS calls encountered .* failures")
})



test_that("Regression test for gene_score_archr() Issues 185 + 188", {
  # Test setup:
  #  - Cell1 overlaps gene1 and gene3
  #  - cell2 overlaps just gene1
  #  - cell3 overlaps just gene3
  #  - All fragments overlap gene body and gene lengths are the same, so 
  #    we don't need to worry about distance parameters affecting things
  fragments <- tibble::tibble(
    chr = c("chr1", "chr1", "chr2", "chr2"),
    start = 10,
    end = start + 20,
    cell_id = c("cell1", "cell2", "cell1", "cell3")
  )
  genes <- tibble::tibble(
    chr = c("chr1", "chr1", "chr2"),
    start = c(0, 10000, 15),
    end = start + 500,
    strand = "+",
    gene_id = c("gene1", "gene2", "gene3")
  )
  chromosome_sizes <- tibble::tibble(
    chr = c("chr1", "chr2"),
    start = 0,
    end = c(20000, 500)
  )

  bp_frags <- convert_to_fragments(fragments)

  # Check that tile_width is getting handled properly in `gene_score_weights_archr`
  expected_tiles <- as.integer(sum((chromosome_sizes$end + 99) %/% 100))
  expect_identical(ncol(gene_score_weights_archr(genes, chromosome_sizes, tile_width=100)), expected_tiles)

  ans <- matrix(c(
    5000, 10000, 0,
    0, 0, 0,
    5000, 0, 10000
  ), nrow=3, byrow=TRUE, dimnames=list(genes$gene_id, unique(fragments$cell_id))) %>%
    as("dgCMatrix")
  
  bp_ans <- gene_score_archr(bp_frags, genes, chromosome_sizes)
  expect_identical(as(bp_ans, "dgCMatrix"), ans)

  # Try with out-of-order chromosome_sizes
  bp_ans <- gene_score_archr(bp_frags, genes, chromosome_sizes[c(2,1),])
  expect_identical(as(bp_ans, "dgCMatrix"), ans)

  # Check that tile_width changing doesn't cause crashes in `gene_score_archr`
  bp_ans <- gene_score_archr(bp_frags, genes, chromosome_sizes, tile_width=1000)
  expect_identical(as(bp_ans, "dgCMatrix"), ans)
})

test_that("qc_scATAC handles genes near the start of a chromosome", {
  frag_tbl <- tibble::tibble(
    chr = "chr1",
    start = c(0L, 25L, 50L, 75L),
    end = start + 20L,
    cell_id = rep(c("cell1", "cell2"), each = 2)
  )
  fragments <- convert_to_fragments(frag_tbl)

  genes <- tibble::tibble(
    chr = "chr1",
    start = c(0L, 25L),
    end = c(75L, 80L),
    strand = c("+", "-")
  )
  blacklist <- tibble::tibble(
    chr = "chr1",
    start = 500L,
    end = 550L
  )

  qc <- qc_scATAC(fragments, genes, blacklist)
  expect_identical(sort(qc$cellName), sort(cellNames(fragments)))
})