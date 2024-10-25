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
#' @param path (character vector) Path(s) to save bedgraph to, optionally ending in ".gz" to add gzip compression. If `cell_groups` is provided,
#'   `path` must be a named character vector, with one name for each level in `cell_groups`
#' @param insertion_mode (string) Which fragment ends to use for coverage calculation. One of "both", "start_only", or "end_only"
#' @inheritParams footprint
#' @export
write_insertion_bedgraph <- function(fragments, path, cell_groups = rlang::rep_along(cellNames(fragments), "all"), insertion_mode=c("both", "start_only", "end_only")) {
  assert_is(fragments, "IterableFragments")
  assert_is_character(path)
  assert_is(cell_groups, c("character", "factor"))
  insertion_mode <- match.arg(insertion_mode)
  path_names <- names(path)
  path <- suppressWarnings(normalizePath(path))
  names(path) <- path_names
  cell_groups <- as.factor(cell_groups)
  assert_len(path, length(levels(cell_groups)))
  
  if (length(levels(cell_groups)) == 1) {
    names(path) <- "all"
  } else {
    assert_len(cell_groups, length(cellNames(fragments)))
    assert_has_names(path, levels(cell_groups))
    path <- path[levels(cell_groups)]
    cell_groups <- as.integer(cell_groups) - 1L
  }
  write_insertion_bedgraph_cpp(iterate_fragments(fragments), cell_groups, path, insertion_mode)
}

#' Create bed files from fragments split by cell group.
#' @param path (character vector) Path to save bed files. If `cell_groups` is provided, this must be a character vector with one name for each level in `cell_groups` 
#' Else, this must be a character vector of length 1.
#' @param cell_groups (character vector or factor) Cluster assignments for each cell.
#' @param threads (int) Number of threads to use.
#' @param verbose (bool) Whether to provide verbose progress output to console.
#' @return `NULL`
#' @inheritParams write_insertion_bedgraph
#' @keywords internal
write_insertion_bed <- function(fragments, path,
                                cell_groups = rlang::rep_along(cellNames(fragments), "all"),
                                insertion_mode = c("start_only", "both", "end_only"),
                                verbose = FALSE,
                                threads = 1) {
  assert_is(fragments, "IterableFragments")
  assert_is(cell_groups, c("character", "factor"))
  assert_is_character(path)
  assert_is_wholenumber(threads)
  assert_is(verbose, "logical")
  insertion_mode <- match.arg(insertion_mode)
  path_names <- names(path)
  path <- suppressWarnings(normalizePath(path))
  names(path) <- path_names
  cell_groups <- as.factor(cell_groups)
  assert_len(path, length(levels(cell_groups)))

  # Prep inputs
  if (length(levels(cell_groups)) == 1) {
    names(path) <- levels(cell_groups)
  } else {
    assert_len(cell_groups, length(cellNames(fragments)))
    assert_has_names(path, levels(cell_groups))
  }
  cell_groups_int <- as.integer(cell_groups)
  cluster_name_mapping <- levels(cell_groups)
  # Parallelize writing bed inputs into MACS
  parallel::mclapply(seq_along(cluster_name_mapping), function(i) {
    if (verbose) log_progress(paste0("Writing bed file for cluster: ", cluster_name_mapping[[i]]))
    fragments_by_cluster <- select_cells(fragments, cell_groups_int == i)
    write_insertion_bed_cpp(
      iterate_fragments(fragments_by_cluster),
      path[cluster_name_mapping[[i]]],
      insertion_mode
    )
    if (verbose) log_progress(paste0("Bed file for cluster: ", cluster_name_mapping[[i]],
                                     " written to: ", path[cluster_name_mapping[[i]]]))
  }, mc.cores = threads, mc.preschedule = FALSE)
  if (verbose) {
    message(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " Finished writing bed files"))
  }
}


#' Call peaks using MACS2/3
#' 
#' Export pseudobulk bed files as input for MACS, then run MACS and read the output peaks as a tibble.
#' Each step can can be run independently, allowing for quickly re-loading the results of an already completed call,
#' or running MACS externally (e.g. via cluster job submisison) for increased parallelization. See details for more information.
#'
#' @param path (string) Parent directory to store MACS inputs and outputs.
#' Inputs are stored in `<path>/input/` and outputs in `<path>/output/<group>/`. See "File format" in details
#' @param effective_genome_size (numeric) Effective genome size for MACS. Default is `2.9e9` following MACS default for GRCh38. See [deeptools](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html)
#'    for values for other common genomes.
#' @param insertion_mode (string) Which fragment ends to use for coverage calculation. One of `both`, `start_only`, or `end_only`.
#' @param step (string) Which step to run. One of  `all`, `prep-inputs`, `run-macs`, `read-outputs`.  If `prep-inputs`, create the input bed files for macs,
#' and provides a shell script per cell group with the command to run macs.  If `run-macs`, also run bash scripts to execute macs.
#' If `read-outputs`, read the outputs into tibbles.
#' @param macs_executable (string) Path to either MACS2/3 executable. Default (`NULL`) will autodetect from PATH.
#' @param additional_params (string) Additional parameters to pass to MACS2/3. 
#' @param verbose (bool) Whether to provide verbose output from MACS. Only used if step is `run-macs` or `all`.
#' @param threads (int) Number of threads to use.
#' @return
#'  - If step is `prep-inputs`, return script paths for each cell group given as a character vector.
#'  - If step is `run-macs`, return `NULL`.
#'  - If step is `read-outputs` or `all`, returns a tibble with all the peaks from each cell group concatenated.
#' Columnns are `chr`, `start`, `end`, `group`, `name`, `score`, `strand`, `fold_enrichment`, `log10_pvalue`, `log10_qvalue`, `summit_offset`
#' @details 
#' **File format**:
#'  - Inputs are written such that a bed file used as input into MACS, 
#' as well as a shell file containing a call to MACS are written for each cell group.
#'  - Bed files containing `chr`, `start`, and `end` coordinates of insertions are written at `<path>/input/<group>.bed.gz`.
#'  - Shell commands to run MACS manually are written at `<path>/input/<group>.sh`.
#'
#' Outputs are written to an output directory with a subdirectory for each cell group.
#' Each cell group's output directory contains a file for narrowPeaks, peaks, and summits.
#'  - NarrowPeaks are written at `<path>/output/<group>/<group>_peaks.narrowPeak`.
#'  - Peaks are written at `<path>/output/<group>/<group>_peaks.xls`.
#'  - Summits are written at `<path>/output/<group>/<group>_summits.bed`.
#'
#' Only the narrowPeaks file is read into a tibble and returned.
#' For more information on outputs from MACS, visit the [MACS docs](https://macs3-project.github.io/MACS/docs/callpeak.html)
#'
#' **Performance**:
#' 
#' Running on a 2600 cell dataset and taking both start and end insertions into account, written input bedfiles and MACS outputs 
#' used 364 MB and 158 MB of space respectively.  With 4 threads, running this function end to end took 74 seconds, with 61 of those seconds spent on running MACS.
#'
#' **Running MACS manually**:
#'
#' To run MACS manually, you will first run `call_peaks_macs()` with `step="prep-inputs`. Then, manually run all of the
#' shell scripts generated at `<path>/input/<group>.sh`. Finally, run `call_peaks_macs()` again with the same original arguments, but
#' setting `step="read-outputs"`.
#' @inheritParams call_peaks_tile
#' @export
call_peaks_macs <- function(fragments, path,
                            cell_groups = rlang::rep_along(cellNames(fragments), "all"), effective_genome_size = 2.9e9,
                            insertion_mode = c("both", "start_only", "end_only"),
                            step = c("all", "prep-inputs", "run-macs", "read-outputs"),
                            macs_executable = NULL,
                            additional_params = "--call-summits --keep-dup all --shift -75 --extsize 150 --nomodel --nolambda",
                            verbose = FALSE,
                            threads = 1) {
  assert_is(fragments, c("IterableFragments", "NULL"))
  assert_is(cell_groups, c("character", "factor"))
  assert_is_numeric(effective_genome_size)
  assert_is_character(path)
  assert_is_wholenumber(threads)
  insertion_mode <- match.arg(insertion_mode)
  step <- match.arg(step)
  cell_groups <- as.factor(cell_groups)
  levels(cell_groups) <- normalize_unique_file_names(levels(cell_groups))
  # Create paths
  dir.create(file.path(path, "input"), showWarnings = FALSE, recursive = TRUE)
  path <- normalizePath(path)
  path_bed_input <- paste0(path, "/input/", levels(cell_groups), ".bed.gz")
  names(path_bed_input) <- levels(cell_groups)
  path_macs_output <- paste0(path, "/output/", levels(cell_groups))
  # Check if MACS can be run
  if (!(step %in% c("read-outputs"))) {
    macs_executable <- macs_path_is_valid(macs_executable)
  }
  # Write bed files as input into MACS
  if (step %in% c("prep-inputs", "all")) {
    write_insertion_bed(fragments, path_bed_input, cell_groups, insertion_mode, threads, verbose = verbose)
  }

  # prep macs call
  macs_call_template <- c('"%s" callpeak -g %s --name "%s" --treatment "%s"',
                          '--outdir "%s" --format BED',
                          '%s')
  macs_call_template <- paste(macs_call_template, collapse = " ")
  macs_call <- sprintf(macs_call_template,
                       macs_executable, effective_genome_size, levels(cell_groups),
                       path_bed_input, path_macs_output, additional_params)

  if (step %in% c("prep-inputs", "run-macs", "all")) {
    shell_paths <- paste0(path, "/input/", levels(cell_groups), ".sh")
    for (cluster_idx in seq_along(levels(cell_groups))) { 
      writeLines(macs_call[cluster_idx], con = shell_paths[cluster_idx])
    }
  }
  if (step == "prep-inputs") return(shell_paths)
  # Run macs
  if (step %in% c("run-macs", "all")) {
    # Run macs through the shell files created in the previous step
    dir.create(file.path(path, "output"), showWarnings = FALSE, recursive = TRUE)
    file_names <- list.files(file.path(path, "input"), pattern = "\\.sh$", full.names = FALSE)
    if (length(file_names) != length(levels(cell_groups))) {
      warning("Number of shell files does not match number of clusters")
    }
    parallel::mclapply(file_names, function(shell_file) {
      cluster <- gsub(".sh$", "", shell_file)
      if (verbose) log_progress(paste0("Running MACS for cluster: ", cluster))
      dir.create(file.path(path, "output", cluster), showWarnings = FALSE, recursive = TRUE)
      log_file <- paste0(path, "/output/", cluster, "/log.txt")
      macs_message <- system2("bash", sprintf("'%s'", file.path(path, "input", shell_file)),stdout = log_file, stderr = log_file, env = c("OMP_NUM_THREADS=1"))
      # Try detecting if macs failed before writing that cluster is finished
      if (macs_message != 0) {
        stop(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " Error running MACS for cluster: ", cluster, "\n",
                    "MACS log file written to: ", log_file))
      } else if (verbose) {
        log_progress(paste0(" Finished running MACS for cluster: ", cluster))
        log_progress(paste0(" MACS log file written to: ", log_file))
      }
    }, mc.cores = threads, mc.preschedule = FALSE)
  }
  # Read outputs
  if (step %in%  c("read-outputs", "all")) {
    peaks <- list()
    output_dirs <- list.dirs(file.path(path, "output"), full.names = FALSE, recursive = FALSE)
    for (cluster in output_dirs) {
      peak_path <- paste0(path, "/output/", cluster, "/", cluster, "_peaks.narrowPeak")
      peaks[[cluster]] <- readr::read_tsv(peak_path,
                                          col_names=c("chr", "start", "end", "name", 
                                                      "score", "strand", "fold_enrichment", 
                                                      "log10_pvalue", "log10_qvalue", "summit_offset"),
                                          show_col_types = FALSE)
      peaks[[cluster]]$group <- cluster
      # set cluster column as the fourth column
      peaks[[cluster]] <- peaks[[cluster]][, c(1:3, 11, 4:10)]
    }
    # We want to treat users files as the ground truth, so we give a warning if this gives different information than what we 
    # expect given cell_groups
    if (length(peaks) != length(levels(cell_groups))) warning("Number of output files does not match number of clusters")
    #  combine all peaks together into a single dataframe
    peaks <- dplyr::bind_rows(peaks)
    return(peaks)
  }
}

#' Call peaks using MACS2/3
#'
#' @description
#' `r lifecycle::badge("deprecated")` 
#'
#' This function has been renamed to `call_peaks_macs()`
#' @export
#' @keywords internal
call_macs_peaks <- function(...) {
  lifecycle::deprecate_warn("0.2.0", "call_macs_peaks()", "call_peaks_macs()")
  return(call_peaks_macs(...))
}


#' Test if MACS executable is valid.
#' If macs_executable is NULL, this function will try to auto-detect MACS from PATH, with preference for MACS3 over MACS2.
#' If macs_executable is provided, this function will check if MACS can be called.
#' @return MACS executable path.
#' @inheritParams call_peaks_macs
#' @keywords internal
macs_path_is_valid <- function(macs_executable) {
  if (is.null(macs_executable)) {
    # Check if can call version on macs and that macs is indeed being called
    if ((suppressWarnings((system2("macs3", args = "--version", stdout = FALSE, stderr = FALSE) == 0)) &&
        (grepl("macs", system2("macs3", args = "--version", stdout = TRUE))))) {
      macs_executable <- "macs3"
    } else if ((suppressWarnings((system2("macs2", args = "--version", stdout = FALSE, stderr = FALSE) == 0)) && 
               (grepl("macs", system2("macs2", args = "--version", stdout = TRUE))))) {
      macs_executable <- "macs2"
    } else {
      stop(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
                  paste0(" MACS not found. Please install MACS3 or MACS2")))
    }
  # Only run if macs executable is provided and the executable is indeed macs
  } else if (!(suppressWarnings((system2(macs_executable, args = "--version", stdout = FALSE, stderr = FALSE) == 0)) &&
              (grepl("macs", system2(macs_executable, args = "--version", stdout = TRUE))))) {
    stop(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
                sprintf(" MACS not found for MACS executable: %s \nPlease install MACS3 or MACS2", macs_executable)))
  }
  return(macs_executable)
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
