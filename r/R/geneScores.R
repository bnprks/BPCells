# Copyright 2022 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

#' Find signed distance to nearest genomic ranges
#'
#' Given a set of genomic ranges, find the distance to the nearest neighbors both
#' upstream and downstream.
#' @param ranges `r document_granges(strand="default")`
#' @inheritParams normalize_ranges
#' @param addArchRBug boolean to reproduce ArchR's bug that incorrectly handles nested genes
#' @return A 2-column data.frame with columns upstream and downstream, containing
#' the distances to the nearest neighbor in the respective directions.
#' For ranges on `+` or `*` strand, distance is calculated as:
#'   - upstream = `max(start(range) - end(upstreamNeighbor), 0)`
#'   - downstream = `max(start(downstreamNeighbor) - end(range), 0)`
#'
#' For ranges on `-` strand, the definition of upstream and downstream is flipped.
#' Note that this definition of distance is one off from
#' `GenomicRanges::distance()`, as ranges that neighbor but don't overlap are given
#' a distance of 1 rather than 0.
#' @export
range_distance_to_nearest <- function(ranges, addArchRBug = FALSE, zero_based_coords = !is(ranges, "GRanges")) {
  ranges <- normalize_ranges(ranges, metadata_cols = "strand", zero_based_coords = zero_based_coords)
  # Label indices of the ranges incase they aren't pre-sorted
  ranges$idx <- seq_len(nrow(ranges))

  # Default to distance 0
  upstream <- rep_len(0, nrow(ranges))
  downstream <- rep_len(0, nrow(ranges))

  # Filter down to ranges that aren't fully nested within others (including identical start and/or end)
  # Fully nested genes have distance 0 to neighbors.
  # Once we filter them out, we know that starts + ends are unique, and
  # sorting by starts implies also getting sorted ends.
  if (!addArchRBug) {
    overlaps <- range_overlaps(ranges, ranges) %>%
      dplyr::filter(from != to) %>%
      dplyr::filter(
        ranges$start[from] >= ranges$start[to],
        ranges$end[from] <= ranges$end[to]
      ) %>%
      dplyr::pull(from) %>%
      unique()
    ranges <- ranges[-overlaps, ]
  }

  dists <- ranges %>%
    dplyr::group_by(chr) %>%
    dplyr::arrange(start, end) %>%
    dplyr::mutate(
      neighbor_dist = c(Inf, start[-1] + 1 - end[-length(end)]),
      start_dist = pmax(0, neighbor_dist),
      end_dist = pmax(0, c(neighbor_dist[-1], Inf)),
      upstream = dplyr::if_else(strand, start_dist, end_dist),
      downstream = dplyr::if_else(strand, end_dist, start_dist)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(idx, upstream, downstream)
  upstream[dists$idx] <- dists$upstream
  downstream[dists$idx] <- dists$downstream
  tibble::tibble(
    upstream = upstream,
    downstream = downstream
  )
}

#' Extend genome ranges in a strand-aware fashion.
#' @inheritParams normalize_ranges
#' @param upstream Number of bases to extend each range upstream (negative to shrink width)
#' @param downstream Number of bases to extend each range downstream (negative to shrink width)
#' @param chromosome_sizes (optional) Size of chromosomes as a [genomic-ranges] object
#' @details Note that ranges will be blocked from extending past the beginning of the chromosome (base 0),
#' and if `chromosome_sizes` is given then they will also be blocked from extending past the end of the chromosome
#'
#' @export
extend_ranges <- function(ranges, upstream = 0, downstream = 0, metadata_cols = c("strand"),
                          chromosome_sizes = NULL, zero_based_coords = !is(ranges, "GRanges")) {
  metadata_cols <- union(metadata_cols, "strand")
  ranges <- normalize_ranges(ranges, metadata_cols = metadata_cols, zero_based_coords = zero_based_coords)
  ranges$start <- pmax(
    0,
    ranges$start - dplyr::if_else(ranges$strand, upstream, downstream)
  )
  ranges$end <- pmax(
    0,
    ranges$end + dplyr::if_else(ranges$strand, downstream, upstream)
  )
  if (!is.null(chromosome_sizes)) {
    chr_size_lookup <- dplyr::pull(chromosome_sizes, end, name = chr)
    if (length(setdiff(ranges$chr, names(chr_size_lookup))) > 0) {
      rlang::abort("chromosome_sizes does not contain all chromosomes present in ranges")
    }
    ranges$start <- pmin(chr_size_lookup[as.character(ranges$chr)], ranges$start)
    ranges$end <- pmin(chr_size_lookup[as.character(ranges$chr)], ranges$end)
  }
  ranges
}

#' Calculate gene-tile distances for ArchR gene activities
#'
#' ArchR-style gene activity scores are based on a weighted sum of each tile
#' according to the signed distance from the tile to a gene body. This function
#' calculates the signed distances according to ArchR's default parameters.
#'
#' ArchR's tile distance algorithm works as follows
#' 1. Genes are extended 5kb upstream
#' 2. Genes are linked to any tiles 1kb-100kb upstream + downstream, but tiles
#'    beyond a neighboring gene are not considered
#' @param genes `r document_granges("Gene coordinates", strand="default")`
#' @param tile_width Size of tiles to consider
#' @inheritParams extend_ranges
#' @param addArchRBug Replicate ArchR bug in handling nested genes
#'
#' @return Tibble with one range per tile, with additional metadata
#' columns gene_idx (row index of the gene this tile corresponds to) and
#' distance.
#'
#' Distance is a signed distance calculated such that if the tile has a smaller
#' start coordinate than the gene and the gene is on the + strand, distance will
#' be negative. The distance of adjacent but non-overlapping regions is 1bp, counting
#' up from there.
#' @export
gene_score_tiles_archr <- function(genes, chromosome_sizes = NULL, tile_width = 500, addArchRBug = FALSE) {
  assert_is_wholenumber(tile_width)
  # Extend upstream by 5kb
  genes <- extend_ranges(genes, upstream = 5000, chromosome_sizes = chromosome_sizes)

  # Extend up 1kb - 100kb upstream + downstream, but stop if we hit a neighboring gene.
  extension_lengths <- range_distance_to_nearest(genes, addArchRBug = addArchRBug)
  extension_lengths$upstream <- pmin(pmax(extension_lengths$upstream - tile_width, 1000), 100000)
  extension_lengths$downstream <- pmin(pmax(extension_lengths$downstream - tile_width, 1000), 100000)
  extended <- extend_ranges(genes, upstream = extension_lengths$upstream, downstream = extension_lengths$downstream, chromosome_sizes = chromosome_sizes)
  # 1. Get a list of tile start coords
  first_tile <- extended$start %/% tile_width
  last_tile <- (extended$end - 1) %/% tile_width
  tile_count <- last_tile - first_tile + 1
  # 2. Set up for vectorized operation: tile_idx counts up from 1 to tile_count[1], then restarts and counts
  #    from 1 to tile_count[2], etc. There is 1 element for each tile which has a non-zero weight for a gene score
  region_idx <- rep.int(seq_len(nrow(genes)), times = tile_count)
  tile_idx <- seq_len(sum(tile_count)) - rep.int(cumsum(tile_count) - tile_count, times = tile_count)
  # 3. Calculate distances to the gene "core" region
  tile_start <- (tile_idx - 1 + rep(first_tile, times = tile_count)) * tile_width

  distance <- dplyr::if_else(
    tile_start >= genes$end[region_idx],
    tile_start + 1 - genes$end[region_idx],
    pmin(0, tile_start + tile_width - 1 - genes$start[region_idx])
  )
  distance <- distance * dplyr::if_else(genes$strand[region_idx], 1, -1)
  tiles <- tibble::tibble(
    chr = genes$chr[region_idx],
    start = tile_start,
    end = tile_start + tile_width,
    gene_idx = region_idx,
    distance = distance
  )
  return(tiles)
}

#' Calculate GeneActivityScores
#'
#' Gene activity scores can be calculated as a distance-weighted sum of per-tile accessibility.
#' The tile weights for each gene can be represented as a sparse matrix of dimension genes x tiles.
#' If we multiply this weight matrix by a corresponding tile matrix (tiles x cells), then we can 
#' get a gene activity score matrix of genes x cells. `gene_score_weights_archr()` calculates the 
#' weight matrix (best if you have a pre-computed tile matrix), while `gene_score_archr()` provides
#' a easy-to-use wrapper.
#' 
#' **gene_score_weights_archr**:
#' 
#' Given a set of tile coordinates and distances returned by `gene_score_tiles_archr()`,
#' calculate a weight matrix of dimensions genes x tiles. This matrix can be
#' multiplied with a tile matrix to obtain ArchR-compatible gene activity scores.
#'
#' @inheritParams call_peaks_tile
#' @inheritParams gene_score_tiles_archr
#' @param blacklist `r document_granges("Regions to exclude from calculations,")`
#' @param gene_name_column If not NULL, a column name of `genes` to use as row names
#' @return **gene_score_weights_archr**
#' 
#' Weight matrix of dimension genes x tiles
#' @export
#' @rdname gene_scores
gene_score_weights_archr <- function(genes, chromosome_sizes, blacklist = NULL, tile_width = 500, gene_name_column="gene_id", addArchRBug = FALSE) {
  if (!is.null(gene_name_column)) {
    assert_is_character(gene_name_column)
    gene_names <- normalize_ranges(genes, metadata_cols = c("strand", gene_name_column))[[gene_name_column]]
    if(anyNA(gene_names)) {
      rlang::abort("Gene score weights found NA gene names. Check gene_name_column argument")
    }
  }
  genes <- normalize_ranges(genes, metadata_cols = "strand")

  chromosome_sizes <- normalize_ranges(chromosome_sizes)
  if (!is.null(blacklist)) {
    blacklist <- normalize_ranges(blacklist)
  }
  assert_is_wholenumber(tile_width)
  assert_is(addArchRBug, "logical")

  tiles <- gene_score_tiles_archr(genes, chromosome_sizes = chromosome_sizes, addArchRBug = addArchRBug)
  # Correct distance to match the GenomicRanges convention
  tiles$distance <- tiles$distance - sign(tiles$distance)
  # Filter out blacklist if requested
  if (!is.null(blacklist)) {
    blacklist_tiles <- unique(range_overlaps(tiles, blacklist)$from)
    tiles <- tiles[-blacklist_tiles, ]
  }

  tile_width <- tiles$end[1] - tiles$start[1]
  tile_count <- chromosome_sizes %>%
    dplyr::mutate(tile_count = (end - start + tile_width - 1) %/% tile_width) %>%
    dplyr::pull(tile_count, name = chr)
  start_tile <- cumsum(tile_count) - tile_count

  res <- Matrix::sparseMatrix(
    i = tiles$gene_idx,
    j = start_tile[as.character(tiles$chr)] + tiles$start %/% tile_width + 1,
    x = exp(-abs(tiles$distance) / 5000) + exp(-1),
    dims = c(max(tiles$gene_idx), sum(tile_count))
  )
  if (!is.null(gene_name_column)) {
    rownames(res) <- gene_names
  }

  gene_widths <- genes$end - genes$start + 5000
  gene_scale_factor <- 5
  gene_weights <- 1 + (1 / gene_widths) * (gene_scale_factor - 1) / (max(1 / gene_widths) - min(1 / gene_widths))
  res <- res * gene_weights
  return(res)
}

#' @param fragments Input fragments object
#' @param tile_matrix_path Path of a directory where the intermediate tile matrix will be saved
#' @param tile_max_count Maximum value in the tile counts matrix. If not null, tile counts higher than this will be clipped to `tile_max_count`.
#'   Equivalent to `ceiling` argument of `ArchR::addGeneScoreMatrix()`
#' @param scale_factor If not null, counts for each cell will be scaled to sum to `scale_factor`. Equivalent to `scaleTo` argument of `ArchR::addGeneScoreMatrix()`
#' @return **gene_score_archr**
#' 
#' Gene score matrix of dimension genes x cells.
#' @export
#' @rdname gene_scores
gene_score_archr <- function(fragments, genes, chromosome_sizes, blacklist=NULL, tile_width = 500, gene_name_column="gene_id", addArchRBug = FALSE,
                             tile_max_count=4, scale_factor=10000, tile_matrix_path=tempfile(pattern = "gene_score_tile_mat")) {
  assert_is(fragments, "IterableFragments")
  if (!is.null(tile_max_count)) assert_is_wholenumber(tile_max_count)
  if (!is.null(scale_factor)) assert_is_numeric(scale_factor)
  gene_weights <- gene_score_weights_archr(genes=genes, chromosome_sizes=chromosome_sizes, blacklist=blacklist, tile_width=tile_width, gene_name_column=gene_name_column, addArchRBug=addArchRBug)
  
  tile_coords <- normalize_ranges(chromosome_sizes) %>%
    dplyr::mutate(tile_width = tile_width)
  
  tile_mat <- tile_matrix(fragments, tile_coords) %>%
    write_matrix_dir(tile_matrix_path)

  if (!is.null(tile_max_count)) {
    tile_mat <- min_scalar(tile_mat, tile_max_count)
  }

  res <- gene_weights %*% tile_mat

  if (!is.null(scale_factor)) {
    res <- multiply_cols(res, scale_factor / colSums(res))
  }

  return(res)
}