# Copyright 2022 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

#' Count fragments by nucleosomal size
#' @param fragments Fragments object
#' @param nucleosome_width Integer cutoff to use as nucleosome width
#' @return List with names `subNucleosomal`, `monoNucleosomal`, `multiNucleosomal`, and `nFrags`, containing the
#'         count vectors of fragments in each class per cell. 
#' @details
#' Shorter than `nucleosome_width` is `subNucleosomal`,
#' `nucleosome_width` to `2*nucleosome_width-1` is `monoNucleosomal`, and anything longer is `multiNucleosomal`.
#' The sum of all fragments is given as `nFrags`
#'         
#' @export
nucleosome_counts <- function(fragments, nucleosome_width = 147) {
  assert_is(fragments, "IterableFragments")
  assert_is_wholenumber(nucleosome_width)
  assert_len(nucleosome_width, 1)

  iter <- iterate_fragments(fragments)
  res <- nucleosome_counts_cpp(iter, nucleosome_width)
  res[["nFrags"]] <- res[[1]] + res[[2]] + res[[3]]
  return(res)
}

#' Get footprints around a set of genomic coordinates
#'
#' @param fragments IterableFragments object
#' @param ranges `r document_granges("Footprint centers", strand="default")`
#' 
#' "+" strand ranges will footprint around the start coordinate, and "-" strand ranges
#' around the end coordinate.
#'
#' @inheritParams normalize_ranges
#' @param cell_groups Character or factor assigning a group to each cell, in order of
#'   `cellNames(fragments)`
#' @param cell_weights Numeric vector assigning weight factors
#'   (e.g. inverse of total reads) to each cell, in order of `cellNames(fragments)`
#' @param flank Number of flanking basepairs to include on either side of the motif
#' @param normalization_width Number of basepairs at the upstream + downstream
#'   extremes to use for calculating enrichment
#'
#' @return `tibble::tibble()` with columns `group`, `position`, and `count`, `enrichment`
#' @export
footprint <- function(fragments, ranges, zero_based_coords = !is(ranges, "GRanges"),
                      cell_groups = rlang::rep_along(cellNames(fragments), "all"),
                      cell_weights = rlang::rep_along(cell_groups, 1),
                      flank = 125L, normalization_width = flank %/% 10L) {
  assert_is(fragments, "IterableFragments")
  ranges <- normalize_ranges(ranges, metadata_cols = "strand", zero_based_coords = zero_based_coords)
  assert_is(cell_groups, c("character", "factor"))
  assert_len(cell_groups, length(cellNames(fragments)))

  assert_is(cell_weights, c("numeric"))
  assert_len(cell_weights, length(cellNames(fragments)))
  assert_is_wholenumber(flank)

  chr <- as.integer(factor(ranges$chr, chrNames(fragments))) - 1
  cell_groups <- as.factor(cell_groups)

  iter <- iterate_fragments(fragments)
  mat <- footprint_matrix_cpp(
    iter,
    chr,
    ifelse(ranges$strand, ranges$start, ranges$end - 1),
    -1 + 2 * ranges$strand,
    as.integer(flank),
    chrNames(fragments),
    as.integer(cell_groups) - 1,
    cell_weights
  )

  if (normalization_width > 0) {
    flank_indices <- c(seq_len(normalization_width), ncol(mat) + 1 - seq_len(normalization_width))
    mat_norm <- mat / rowMeans(mat[, flank_indices, drop = FALSE])
  } else {
    mat_norm <- NA
  }

  rownames(mat) <- levels(cell_groups)
  data <- tibble::tibble(
    group = rep(rownames(mat), ncol(mat)),
    position = rep(-flank:flank, each = nrow(mat)),
    count = as.vector(mat),
    enrichment = as.vector(mat_norm)
  )
  return(data)
}



#' Calculate ArchR-compatible per-cell QC statistics
#' @param fragments IterableFragments object
#' @param genes `r document_granges("Gene coordinates")`
#' @param blacklist `r document_granges("Blacklisted regions")`
#' @return data.frame with QC data
#' @details
#' This implementation mimics ArchR's default parameters. For uses requiring more flexibility to tweak default parameters,
#' the best option is to re-implement this function with required changes.
#' Output columns of data.frame:
#'  - `cellName`: cell name for each cell
#'  - `nFrags`: number of fragments per cell
#'  - `subNucleosomal`, `monoNucleosomal`, `multiNucleosomal`: number of fragments of size 1-146bp, 147-254bp, and 255bp + respectively.
#'    equivalent to ArchR's nMonoFrags, nDiFrags, nMultiFrags respectively
#'  - `TSSEnrichment`: `AvgInsertInTSS / max(AvgInsertFlankingTSS, 0.1)`, where `AvgInsertInTSS` is `ReadsInTSS / 101` (window size),
#'    and `AvgInsertFlankingTSS` is `ReadsFlankingTSS / (100*2)` (window size). The `max(0.1)` ensures that very low-read cells
#'    do not get assigned spuriously high TSSEnrichment.
#'  - `ReadsInPromoter`: Number of reads from 2000bp upstream of TSS to 101bp downstream of TSS
#'  - `ReadsInBlacklist`: Number of reads in the provided blacklist region
#'  - `ReadsInTSS`: Number of reads overlapping the 101bp centered around each TSS
#'  - `ReadsFlankingTSS`: Number of reads overlapping 1901-2000bp +/- each TSS
#'
#' Differences from ArchR:
#' Note that ArchR by default uses a different set of annotations to derive TSS sites and promoter sites.
#' This function uses just one annotation for gene start+end sites, so must be called twice to exactly
#' re-calculate the ArchR QC stats.
#'
#' ArchR's `PromoterRatio` and `BlacklistRatio` are not included in the output, as they can be easily calculated
#' from `ReadsInPromoter / nFrags` and  `ReadsInBlacklist / nFrags`. Similarly, ArchR's `NucleosomeRatio` can be calculated
#' as `(monoNucleosomal + multiNucleosomal) / subNucleosomal`.
#'
#' @export
qc_scATAC <- function(fragments, genes, blacklist) {
  assert_is(fragments, "IterableFragments")
  genes <- normalize_ranges(genes, metadata_cols = "strand")
  blacklist <- normalize_ranges(blacklist)

  # Things to check: standard chromosomes?
  # standard_chr <- grep("^chr[0-9XY]+$", chrNames(fragments))
  nucleosome_qc <- nucleosome_counts(fragments)

  tss <- genes %>%
    dplyr::mutate(start = dplyr::if_else(strand, start, end - 1L)) %>%
    dplyr::distinct(chr, start)

  # Compute signal & background regions for TSSEnrichment calculation
  tss_window_width <- 101L
  tss_window <- tss %>% dplyr::mutate(start = start - 50L, end = start + tss_window_width)

  tss_flank_width <- 100L
  tss_flank <- dplyr::bind_rows(
    dplyr::mutate(tss, start = start + 1901L, end = start + tss_flank_width),
    dplyr::mutate(tss, start = start - 2000L, end = start + tss_flank_width)
  )

  promoters <- genes %>%
    dplyr::mutate(
      start = dplyr::if_else(strand, start - 2000L, end - 101L),
      end = start + 2000L + 101L
    )

  # Calculate overlaps with all regions in one pass
  regions <- dplyr::bind_rows(tss_window, tss_flank, promoters, blacklist, .id = "origin") %>%
    dplyr::mutate(origin = as.integer(origin)) %>%
    dplyr::filter(chr %in% chrNames(fragments))

  regions <- regions[order_ranges(regions, chrNames(fragments)), ]

  membership_mat <- matrix(0, ncol = nrow(regions), nrow = 4)
  membership_mat[matrix(c(regions$origin, seq_len(nrow(regions))), ncol = 2)] <- 1

  overlap_sums <- membership_mat %*% peak_matrix(fragments, regions)
  rownames(overlap_sums) <- c("tss_window", "tss_flank", "promoters", "blacklist")

  tibble::tibble(
    cellName = cellNames(fragments),
    TSSEnrichment = overlap_sums["tss_window", ] / tss_window_width /
      pmax(overlap_sums["tss_flank", ] / (2 * tss_flank_width), 0.1),
    nFrags = nucleosome_qc[[1]] + nucleosome_qc[[2]] + nucleosome_qc[[3]],
    subNucleosomal = nucleosome_qc$subNucleosomal,
    monoNucleosomal = nucleosome_qc$monoNucleosomal,
    multiNucleosomal = nucleosome_qc$multiNucleosomal,
    ReadsInTSS = overlap_sums["tss_window", ],
    ReadsFlankingTSS = overlap_sums["tss_flank", ],
    ReadsInPromoter = overlap_sums["promoters", ],
    ReadsInBlacklist = overlap_sums["blacklist", ]
  )
}

#' Merge peaks
#'
#' Merge peaks according to ArchR's iterative merging algorithm. More
#' details on the [ArchR website](https://www.archrproject.com/bookdown/the-iterative-overlap-peak-merging-procedure.html)
#'
#' Properties of merged peaks:
#'   - No peaks in the merged set overlap
#'   - Peaks are prioritized according to their order in the original input
#'   - The output peaks are a subset of the input peaks, with no peak boundaries
#'     changed
#'
#' @param peaks `r document_granges("Peaks")`  
#'
#'  Must be ordered by priority and have columns chr, start, end.
#' @return `tibble::tibble()` with a nonoverlapping subset of the rows in peaks. All metadata
#'  columns are preserved
merge_peaks_iterative <- function(peaks) {
  assert_is(peaks, c("data.frame", "list"))
  assert_has_names(peaks, c("chr", "start", "end"))

  # Filter out identical peaks as a quick first pass
  peaks <- tibble::as_tibble(peaks) %>%
    dplyr::distinct(chr, start, end, .keep_all = TRUE)

  overlaps <- range_overlaps(peaks, peaks)
  # Initialize keeper set with any non-overlapping peaks
  keeper_sets <- list(setdiff(seq_len(nrow(peaks)), overlaps$from))

  # Maintain invariant: overlaps contains only overlap pairs where we have not
  #    permanently kept or discarded either of the elements in the pair
  overlaps <- dplyr::filter(overlaps, from < to)
  while (nrow(overlaps) > 0) {
    # Add peaks with no higher-ranked overlap to the keeper set
    keeper_set <- setdiff(overlaps$from, overlaps$to)
    keeper_sets <- c(keeper_sets, list(keeper_set))
    # Mark peaks that overlap directly with the keeper set
    discard_set <- overlaps$to[overlaps$from %in% keeper_set] %>% unique()
    # Discard all overlap information on the keeper and discard sets
    overlaps <- overlaps %>%
      dplyr::filter(
        !(to %in% discard_set), !(from %in% discard_set),
        !(from %in% keeper_set)
      )
  }
  return(peaks[sort(unlist(keeper_sets)), ])
}

#' Call peaks from tiles
#'
#' Calling peaks from a pre-set list of tiles can be much faster than using
#' dedicated peak-calling software like `macs3`. The resulting peaks are less
#' precise in terms of exact coordinates, but should be sufficient for most
#' analyses.
#'
#' @param fragments IterableFragments object
#' @param chromosome_sizes `r document_granges("Chromosome start and end coordinates")`  
#' 
#'   See `read_ucsc_chrom_sizes()`.
#' @param cell_groups Grouping vector with one entry per cell in fragments, e.g.
#'    cluster IDs
#' @param effective_genome_size (Optional) effective genome size for poisson
#'    background rate estimation. See [deeptools](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html)
#'    for values for common genomes. Defaults to sum of chromosome sizes, which
#'    overestimates peak significance
#' @param peak_width Width of candidate peaks
#' @param peak_tiling Number of candidate peaks overlapping each base of genome.
#'    E.g. peak_width = 300 and peak_tiling = 3 results in candidate peaks of
#'    300bp spaced 100bp apart
#' @param fdr_cutoff Adjusted p-value significance cutoff
#' @param merge_peaks How to merge significant peaks with `merge_peaks_iterative()`
#'
#' - `"all"` Merge the full set of peaks
#' - `"group"` Merge peaks within each group
#' - `"none"` Don't perform any merging
#'
#' @return tibble with peak calls and the following columns:
#'
#' - `chr`, `start`, `end`: genome coordinates
#' - `group`: group ID that this peak was identified in
#' - `p_val`, `q_val`: Poission p-value and BH-corrected p-value
#' - `enrichment`: Enrichment of counts in this peak compared to a genome-wide
#'    background
#'
#' @details Peak calling steps:
#'
#' 1. Estimate the genome-wide expected insertions per tile based on
#'  `peak_width`, `effective_genome_size`, and per-group read counts
#' 2. Tile the genome with nonoverlapping tiles of size peak_width
#' 3. For each tile and group, calculate p_value based on a Poisson model
#' 4. Compute adjusted p-values using BH method and using the total number of
#' tiles as the number of hypotheses tested.
#' 5. Repeat steps 2-4 `peak_tiling` times, with evenly spaced offsets
#' 6. If `merge_peaks` is "all" or "group": use `merge_peaks_iterative()` within each group to keep only the most
#' significant of the overlapping candidate peaks
#' 7. If `merge_peaks` is "all", perform a final round of `merge_peaks_iterative()`,
#' prioritizing each peak by its within-group significance rank
#' @export
call_peaks_tile <- function(fragments, chromosome_sizes, cell_groups = rep.int("all", length(cellNames(fragments))),
                            effective_genome_size = NULL,
                            peak_width = 200, peak_tiling = 3, fdr_cutoff = 0.01,
                            merge_peaks = c("all", "group", "none")) {
  assert_is(fragments, "IterableFragments")
  assert_not_null(chrNames(fragments))
  assert_not_null(cellNames(fragments))
  assert_true(length(cell_groups) == length(cellNames(fragments)))
  assert_is_wholenumber(peak_width)
  assert_is_wholenumber(peak_tiling)
  assert_is_numeric(fdr_cutoff)

  merge_peaks <- match.arg(merge_peaks)

  fragments <- merge_cells(fragments, cell_groups)

  ranges <- normalize_ranges(chromosome_sizes)
  ranges$tile_width <- peak_width
  ranges <- ranges[order_ranges(ranges, chrNames(fragments)), ]
  if (is.null(effective_genome_size)) {
    effective_genome_size <- sum(ranges$end - ranges$start)
  } else {
    assert_is_numeric(effective_genome_size)
  }

  group_counts <- peak_matrix(fragments, ranges) %>% colSums()
  background_rate <- group_counts / effective_genome_size * peak_width
  min_cutoffs <- qpois(1 - fdr_cutoff, background_rate)

  offsets <- (peak_width * seq_len(peak_tiling)) %/% peak_tiling
  peak_list <- list()
  for (offset in offsets) {
    ranges$start <- offset

    tile_mat <- tile_matrix(fragments, ranges)
    tile_counts <- as(tile_mat, "dgCMatrix")

    total_tiles <- nrow(tile_counts)

    peaks <- tile_counts %>%
      {
        tibble::tibble(
          tile = .@i + 1,
          group = rep.int(seq_len(ncol(.)), diff(.@p)),
          counts = .@x
        )
      } %>%
      dplyr::filter(counts > min_cutoffs[group]) %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(
        p_val = ppois(background_rate[group], counts),
        q_val = p.adjust(p_val, method = "BH", n = total_tiles)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(q_val < fdr_cutoff) %>%
      dplyr::arrange(tile)
    peak_coords <- tile_ranges(tile_mat, peaks$tile)
    peaks <- peaks %>%
      dplyr::transmute(
        chr = peak_coords$chr,
        start = peak_coords$start,
        end = peak_coords$end,
        group = cellNames(fragments)[group],
        p_val = p_val,
        q_val = q_val,
        enrichment = counts / background_rate[group]
      )
    peak_list <- c(peak_list, list(peaks))
  }

  peaks <- dplyr::bind_rows(peak_list) %>%
    dplyr::arrange(dplyr::desc(enrichment))
  if (merge_peaks == "none") {
    return(peaks)
  } else if (merge_peaks == "group") {
    peaks <- peaks %>%
      dplyr::group_by(group) %>%
      dplyr::summarize(
        data = merge_peaks_iterative(dplyr::cur_data_all())
      ) %>%
      dplyr::ungroup()
    return(dplyr::bind_rows(peaks$data))
  } else {
    peaks <- peaks %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(group_rank = rank(dplyr::desc(enrichment))) %>%
      dplyr::arrange(group_rank) %>%
      dplyr::select(!group_rank) %>%
      merge_peaks_iterative()

    return(peaks)
  }
}

#' Write insertion counts to bedgraph file
#' 
#' Write insertion counts data for one or more pseudobulks to bedgraph format. This reports the total
#' number insertions at each basepair for each group listed in `cell_groups`.
#'
#' @param fragments IterableFragments object
#' @param path Path(s) to save bedgraph to, optionally ending in ".gz" to add gzip compression. If `cell_groups` is provided,
#'   `path` must be a named character vector, with one name for each level in `cell_groups`
#' @param insertion_mode Which fragment ends to use for insertion counts calculation. One of "both", "start_only", or "end_only"
#' @inheritParams footprint
#' @export
write_insertion_bedgraph <- function(fragments, path, cell_groups = NULL, insertion_mode=c("both", "start_only", "end_only")) {
  assert_is(fragments, "IterableFragments")
  assert_is_character(path)
  insertion_mode <- match.arg(insertion_mode)
  
  if (is.null(cell_groups)) {
    cell_groups <- rep.int(0L, length(cellNames(fragments)))
    assert_len(path, 1)
  } else {
    assert_is(cell_groups, c("character", "factor"))
    assert_len(cell_groups, length(cellNames(fragments)))
    cell_groups <- as.factor(cell_groups)
    assert_has_names(path, levels(cell_groups))
    path <- path[levels(cell_groups)]
    cell_groups <- as.integer(cell_groups) - 1L
  }
  write_insertion_bedgraph_cpp(iterate_fragments(fragments), cell_groups, path, insertion_mode)
}

#' Create bed files for MACS2/3 input
#' @param fragments IterableFragments object
#' @param cell_names Character vector of cluster assignments for each cell. If is null, all cells are treated as one group.
#' @param path Path to save MACS2/3 output files. If `cell_names` is provided, this must be a character vector with one name for each level in `cell_names` 
#' Else, this must be a character vector of length 1.
#' @param insertion_mode Which fragment ends to use for insertion counts calculation. One of "both", "start_only", or "end_only"
#' @param keep_dups Whether to keep duplicate reads with the same chrom, start, end and cell. Default is TRUE.
#' @param threads Number of threads to use
#' @return NULL
#' @export
prep_macs_inputs <- function(fragments, cell_names,
                             path, insertion_mode = c("both", "start_only", "end_only"),
                             keep_dups = TRUE,
                             threads=1) {
  assert_is(fragments, "IterableFragments")
  assert_is(keep_dups, "logical")
  assert_is(cell_names, c("character", "factor", "NULL"))
  assert_is_character(path)
  assert_is_wholenumber(threads)
  insertion_mode <- match.arg(insertion_mode)

  # Prep inputs
  if (is.null(cell_names)) {
    cell_names <- as.factor(rep.int("all", length(cellNames(fragments))))
    names(path) <- "all"
  } else {
    assert_len(cell_names, length(cellNames(fragments)))
    cell_names <- as.factor(cell_names)
  }
  cell_names_int <- as.integer(cell_names)
  cluster_name_mapping <- levels(cell_names)
  # Parallelize writing bed inputs into MACS
  parallel::mclapply(seq_along(cluster_name_mapping), function(i) {
    write_insertion_bed_by_pseudobulk_cpp(
      iterate_fragments(fragments),
      cell_names_int,
      i,
      path[cluster_name_mapping[[i]]],
      insertion_mode,
      keep_dups
    )
  }, mc.cores = threads, mc.preschedule = FALSE)
}


#' Call MACS2/3 peaks using fragments and cluster assignments.
#' First creates the input bedfiles for each cluster for input to MACS2/3, then runs MACS2/3, and finally reads the outputs into tibbles.
#' Also can be used to run only one of these steps.
#' @param fragments IterableFragments object. Only used if step is "prep-inputs" or "all"
#' @param cell_names Character vector of cluster assignments for each cell.
#' @param genome_size Numeric of length 1 representing Effective genome size.
#' @param path Character vector of length 1 representing parent directory to store MACS inputs and outputs.
#' Inputs are stored in `<path>/input/` and outputs in `<path>/output/<cluster>/`.
#' @param insertion_mode Which fragment ends to use for insertion counts calculation. One of "both", "start_only", or "end_only"
#' @param keep_dups Character describing Whether to keep duplicate reads with the same chrom, start, end and cell in bedfile creation. Default is TRUE.
#' @param step Character describing which step to run. One of  "all", "prep-inputs", "run-macs", "read-outputs".  If "prep-inputs", create the input bedfiles for macs,
#' and provides a shell script per cluster with the command to run macs.  If run-macs, also run bash scripts to execute macs.
#' If read-outputs, read the outputs into tibbles.
#' @param macs_version Character describing which version of MACS to use. One of "macs2" or "macs3"
#' @param additional_params Character describing additional parameters to pass to MACS2/3.  Default is "--nomodel --nolambda"
#' @param quiet  boolean of whether to suppress output from MACS. Only used if step is "run-macs" or "all".  Default is FALSE.
#' @param use_gz boolean of whether to use gzip compression for bed files. Only used if step is "prep-inputs" or "all".  User should
#' expect that not using gzip would result in sizes of bed files that are 2-3x larger.  Default is FALSE.
#' @param threads Number of threads to use. Default is 1.
#' @return If step is "prep-inputs", script paths for each cluster given as a character vector.
#' If step is "run-macs", return NULL
#' If step is "read-outputs" or "all", returns a list of tibbles with peaks for each cluster.
#' @export
call_macs_peaks <- function(fragments = NULL, cell_names, genome_size = 2.7e9,
                            path, insertion_mode = c("both", "start_only", "end_only"),
                            keep_dups = TRUE,
                            step = c("prep-inputs", "run-macs", "read-outputs", "all"), 
                            macs_version = c("macs2", "macs3"),
                            additionalParams = "--keep-dup all --nomodel --nolambda",
                            quiet = FALSE,
                            use_gz = FALSE,
                            threads=1) {
  assert_is(fragments, c("IterableFragments", "NULL"))
  assert_is(cell_names, c("character", "factor"))
  assert_is_numeric(genome_size)
  assert_is_character(path)
  assert_is(keep_dups, "logical")
  assert_is(use_gz, "logical")
  assert_is_wholenumber(threads)
  insertion_mode <- match.arg(insertion_mode)
  step <- match.arg(step)
  macs_version <- match.arg(macs_version)
  cell_names <- as.factor(cell_names)
  if (use_gz) {
    path_bed_input <- paste0(path, "/input/", levels(cell_names), ".bed.gz")
  } else {
    path_bed_input <- paste0(path, "/input/", levels(cell_names), ".bed")
  }
  names(path_bed_input) <- levels(cell_names)
  path_macs_output <- paste0(path, "/output/", levels(cell_names))
  # prep macs call
  macs_call_template <- c("%s callpeak -g %s --name %s --treatment %s",
                          "--outdir %s --format BED --call-summits",
                          "%s")
  macs_call_template <- paste(macs_call_template, collapse = " ")
  macs_call <- sprintf(macs_call_template,
                       macs_version, genome_size, levels(cell_names),
                       path_bed_input, path_macs_output, additionalParams)
  if (step %in% c("prep-inputs", "all")) {
    dir.create(file.path(path, "input"), showWarnings = FALSE, recursive = TRUE)
    prep_macs_inputs(fragments, cell_names, path_bed_input, insertion_mode, keep_dups, threads)
    shell_paths <- paste0(path, "/input/", levels(cell_names), ".sh")
    for (cluster_idx in seq_along(levels(cell_names))) {
      writeLines(macs_call[cluster_idx], con = shell_paths[cluster_idx])
    }
    if (step == "prep-inputs") {
      # Create a shell file for each cluster
      return(shell_paths)
    }
  }
  # Run macs
  if (step %in% c("run-macs", "all")) {
    # create output dirs
    dir.create(file.path(path, "output"), showWarnings = FALSE, recursive = TRUE)
    for (cluster in levels(cell_names)) {
      dir.create(file.path(path, "output", cluster), showWarnings = FALSE, recursive = TRUE)
    }
    # Run macs through the shell files created in the previous step
    file_names <- list.files(file.path(path, "input"), pattern = "\\.sh$", full.names = TRUE)
    if (length(file_names) != length(levels(cell_names))) {
      warning("Number of shell files does not match number of clusters")
    }
    for (shell_file in file_names) {
      system2("bash", shell_file, stdout = !quiet, stderr = !quiet)
    }
  }
  # Read outputs
  if (step %in%  c("read-outputs", "all")) {
    peaks <- list()
    output_dirs <- list.dirs(file.path(path, "output"), full.names = FALSE, recursive = FALSE)
    for (cluster in output_dirs) {
      peak_path <- paste0(path, "/output/", cluster, "/", cluster, "_peaks.narrowPeak")
      peaks[[cluster]] <- readr::read_tsv(peak_path, 
                                          col_names=c("chr", "start", "end", "name", 
                                                      "score", "strand", "signalValue", 
                                                      "pValue", "qValue", "pointSource"),
                                          show_col_types = FALSE)
      # We want to treat users files as the ground truth, so we give a warning if this gives different information than what 
    }
    if (length(peaks) != length(levels(cell_names))) warning("Number of output files does not match number of clusters")
    return(peaks)
  }
}

range_overlaps <- function(a, b) {
  a <- normalize_ranges(a)
  b <- normalize_ranges(b)

  a$cell_id <- seq_len(nrow(a))
  order_a <- order_ranges(a, levels(a$chr), sort_by_end = FALSE)
  order_b <- order_ranges(b, levels(a$chr))

  peak_matrix(
    convert_to_fragments(a[order_a, ]),
    b[order_b, ],
    explicit_peak_names = FALSE,
    mode = "overlaps"
  ) %>%
    as("dgCMatrix") %>%
    {
      tibble::tibble(
          from = rep.int(seq_len(ncol(.)), diff(.@p)),
          to = order_b[.@i +1]
      )
    } %>%
    dplyr::arrange(from, to)
}
