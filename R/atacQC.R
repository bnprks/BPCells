

#' Count fragments by nucleosomal size
#' @param fragments Fragments object
#' @param nucleosome_width Integer cutoff to use as nucleosome width
#' @return List with names subNucleosomal, monoNucleosomal, multiNucleosomal containing the 
#'         count vectors of fragments in each class per cell.
#' @export
nucleosome_counts <- function(fragments, nucleosome_width=147) {
    assert_is(fragments, "IterableFragments")
    assert_wholenumber(nucleosome_width)
    assert_len(nucleosome_width, 1)

    iter <- iterate_fragments(fragments)
    nucleosome_counts_cpp(ptr(iter), nucleosome_width)
}


#' Calculate ArchR-compatible per-cell QC statistics according to ArchR's default parameters
#' @param fragments IterableFragments object
#' @param genes GenomicRanges object with genes coordinates (e.g. ArchR::geneAnnoHg38$TSS)
#' @param blacklist GenomicRanges object with blacklist regions (e.g. ArchR::genomeAnnoHg38$blacklist)
#' @inheritParams peakMatrix
#' @return data.frame with QC data
#' @details 
#' For uses requiring more flexibility to tweak default parameters, the best option is to 
#' re-implement this function with required changes.
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
qc_scATAC <- function(fragments, genes, blacklist, zero_based_coords=TRUE) {
    assert_is(fragments, "IterableFragments")
    assert_is(genes, "GenomicRanges")
    assert_is(blacklist, "GenomicRanges")
    assert_is(zero_based_coords, "logical")

    # Things to check: standard chromosomes?
    # standard_chr <- grep("^chr[0-9XY]+$", chrNames(fragments))
    nucleosome_qc <- nucleosome_counts(fragments)

    tss <- GenomicRanges::resize(genes, 1)
    GenomicRanges::strand(tss) <- "*"
    tss <- unique(tss)

    # Compute signal & background regions for TSSEnrichment calculation
    tss_window <- GenomicRanges::resize(tss, 101, fix="center")
    tss_flank <- c(
        #Positive Flank
        GenomicRanges::GRanges(
            GenomicRanges::seqnames(tss), 
            IRanges::IRanges(GenomicRanges::start(tss) + 1901, GenomicRanges::start(tss) + 2000)),
        #Negative Flank
        GenomicRanges::GRanges(
            GenomicRanges::seqnames(tss), 
            IRanges::IRanges(GenomicRanges::start(tss) - 2000, GenomicRanges::start(tss) - 1901))
    )

    promoters <- GenomicRanges::promoters(genes, upstream=2000, downstream=101)

    # Calculate overlaps with all regions in one pass
    regions <- c(tss_window, tss_flank, promoters, blacklist)
    GenomicRanges::mcols(regions) <- NULL
    regions$origin <- c(
        rep(1, length(tss_window)), 
        rep(2, length(tss_flank)), 
        rep(3, length(promoters)), 
        rep(4, length(blacklist))
    )
    regions <- regions[as.character(GenomicRanges::seqnames(regions)) %in% chrNames(fragments)]

    regions <- data.frame(
        chr = as.character(GenomicRanges::seqnames(regions)),
        start =  GenomicRanges::start(regions),
        end =  GenomicRanges::end(regions),
        origin = regions$origin
    )
    
    regions <- regions[order(match(regions$chr, chrNames(fragments),), regions$end, regions$start),]
    rownames(regions) <- NULL

    membership_mat <- matrix(0, nrow=nrow(regions), ncol=4)
    membership_mat[matrix(c(seq_len(nrow(regions)), regions$origin), ncol=2)] <- 1

    overlap_sums <- peakMatrix(fragments, regions, zero_based_coords) %*% membership_mat
    colnames(overlap_sums) <- c("tss_window", "tss_flank", "promoters", "blacklist")

    data.frame(
        cellName = cellNames(fragments),
        
        TSSEnrichment = overlap_sums[,"tss_window"] / GenomicRanges::width(tss_window[1]) /
            pmax(overlap_sums[,"tss_flank"] / (2*GenomicRanges::width(tss_flank)[1]), 0.1),
        
        nFrags = nucleosome_qc[[1]] + nucleosome_qc[[2]] + nucleosome_qc[[3]],
        subNucleosomal = nucleosome_qc$subNucleosomal,
        monoNucleosomal = nucleosome_qc$monoNucleosomal,
        multiNucleosomal = nucleosome_qc$multiNucleosomal,

        ReadsInTSS = overlap_sums[,"tss_window"],
        ReadsFlankingTSS = overlap_sums[,"tss_flank"],
        ReadsInPromoter = overlap_sums[,"promoters"],
        ReadsInBlacklist = overlap_sums[,"blacklist"]
    )
}