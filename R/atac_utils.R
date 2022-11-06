

#' Count fragments by nucleosomal size
#' @param fragments Fragments object
#' @param nucleosome_width Integer cutoff to use as nucleosome width
#' @return List with names subNucleosomal, monoNucleosomal, multiNucleosomal containing the
#'         count vectors of fragments in each class per cell.
#' @export
nucleosome_counts <- function(fragments, nucleosome_width = 147) {
    assert_is(fragments, "IterableFragments")
    assert_wholenumber(nucleosome_width)
    assert_len(nucleosome_width, 1)

    iter <- iterate_fragments(fragments)
    res <- nucleosome_counts_cpp(ptr(iter), nucleosome_width)
    res[["nFrags"]] <- res[[1]] + res[[2]] + res[[3]]
    return(res)
}

#' Count fragments by size
#' @param fragments Fragments object
#' @return Numeric vector where index i contans the number of length-i fragments
#' @export
length_distribution <- function(fragments) {
    assert_is(fragments, "IterableFragments")

    iter <- iterate_fragments(fragments)
    res <- fragment_lengths_cpp(ptr(iter))
    return(res)
}

#' Get footprints around a set of genomic coordinates
#'
#' @param fragments IterableFragments object
#' @param ranges GRanges object with the positions to footprint, or
#'     list/data frame with columns chr, start, & end. For list/data frame, must
#'     include strand information as character vector of "+"/"-", or TRUE/FALSE for
#'     positive/negative strand. "+" strand motifs will footprint around the start
#'     coordinate, and "-" strand motifs will footprint around the end coordinate
#' @inheritParams normalize_ranges
#' @param cell_groups Character or factor assigning a group to each cell, in order of
#'   `cellNames(fragments)`
#' @param cell_weights Numeric vector assigning weight factors
#'   (e.g. inverse of total reads) to each cell, in order of `cellNames(fragments)`
#' @param flank Number of flanking basepairs to include on either side of the motif
#' @param normalization_width Number of basepairs at the upstream + downstream
#'   extremes to use for normalizing, or 0 to skip normalization
#'
#' @return tibble with columns "group", "position", and "value"
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
    assert_wholenumber(flank)

    chr <- as.integer(factor(ranges$chr, chrNames(fragments))) - 1
    cell_groups <- as.factor(cell_groups)

    iter <- iterate_fragments(fragments)
    mat <- footprint_matrix_cpp(
        ptr(iter),
        chr,
        ifelse(ranges$strand, ranges$start, ranges$end - 1),
        -1 + 2 * ranges$strand,
        as.integer(flank),
        chrNames(fragments),
        as.integer(cell_groups) - 1,
        cell_weights
    )
    if (normalization_width != 0) {
        flank_indices <- c(1:normalization_width, ncol(mat) + 1 - (1:normalization_width))
        mat <- mat / rowMeans(mat[, flank_indices])
    }
    rownames(mat) <- levels(cell_groups)
    res <- tibble::as_tibble(t(mat)) %>%
        dplyr::mutate(pos = -flank:flank) %>%
        tidyr::pivot_longer(!pos, names_to = "group", values_to = "value")
    return(res)
}



#' Calculate ArchR-compatible per-cell QC statistics
#' @param fragments IterableFragments object
#' @param genes GenomicRanges object or list/data.frame with columns chr, start, end
#' @param blacklist GenomicRanges object or list/data.frame with columns chr, start, end
#' @return data.frame with QC data
#' @details
#' This implementation mimics ArchR's default parameters. For uses requiring more flexibility to tweak default parameters,
#' the best option is to re-implement this function with required changes.
#' Output columns of data.frame:
#'  - cellName: cell name for each cell
#'  - nFrags: number of fragments per cell
#'  - subNucleosomal, monoNucleosomal, multiNucleosomal: number of fragments of size 1-146bp, 147-254bp, and 255bp + respectively.
#'    equivalent to ArchR's nMonoFrags, nDiFrags, nMultiFrags respectively
#'  - TSSEnrichment: AvgInsertInTSS / max(AvgInsertFlankingTSS, 0.1), where AvgInsertInTSS ReadsInTSS / 101 (window size),
#'    and AvgInsertFlankingTSS is ReadsFlankingTSS / (100*2) (window size). The max(0.1) ensures that very low-read cells
#'    do not get assigned spuriously high TSSEnrichment.
#'  - ReadsInPromoter: Number of reads from 2000bp upstream of TSS to 101bp downstream of TSS
#'  - ReadsInBlacklist: Number of reads in the provided blacklist region
#'  - ReadsInTSS: Number of reads overlapping the 101bp centered around each TSS
#'  - ReadsFlankingTSS: Number of reads overlapping 1901-2000bp +/- each TSS
#'
#' Differences from ArchR:
#' Note that ArchR by default uses a different set of annotations to derive TSS sites and promoter sites.
#' This function uses just one annotation for gene start+end sites, so must be called twice to exactly
#' re-calculate the ArchR QC stats.
#'
#' ArchR's PromoterRatio and BlacklistRatio are not included in the output, as they can be easily calculated
#' from ReadsInPromoter / nFrags and  ReadsInBlacklist / nFrags. Similarly, ArchR's NucleosomeRatio can be calculated
#' as (monoNucleosomal + multiNucleosomal) / subNucleosomal.
#'
#' @export
qc_scATAC <- function(fragments, genes, blacklist) {
    assert_is(fragments, "IterableFragments")
    genes <- normalize_ranges(genes, metadata_cols = "strand") %>%
        tibble::as_tibble()
    blacklist <- normalize_ranges(blacklist) %>%
        tibble::as_tibble()

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
range_overlaps <- function(a, b) {
    a <- normalize_ranges(a)
    b <- normalize_ranges(b)

    a$cell_id <- seq_len(nrow(a))
    order_a <- order_ranges(a, levels(a$chr), sort_by_end=FALSE)
    order_b <- order_ranges(b, levels(a$chr))

    peak_matrix(
        convert_to_fragments(a[order_a,]),
        b[order_b,],
        explicit_peak_names = FALSE,
        mode = "overlaps"
    ) %>%
        as("dgCMatrix") %>%
        as("dgTMatrix") %>% 
        {tibble::tibble(
            from = .@j + 1,
            to = order_b[.@i + 1]
        )} %>%
        dplyr::arrange(from, to)
}