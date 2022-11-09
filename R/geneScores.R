#' Find signed distance to nearest genomic ranges
#' 
#' Given a set of genomic ranges, find the distance to the nearest neighbors both
#' upstream and downstream.
#' @param ranges GenomicRanges object
#' @param addArchRBug boolean to reproduce ArchR's bug that incorrectly handles nested genes
#' @return A 2-column data.frame with columns upstream and downstream, containing 
#' the distances to the nearest neighbor in the respective directions.
#' For ranges on + or * strand, distance is calculated as:
#'   - upstream = max(start(range) - end(upstreamNeighbor), 0)
#'   - downstream = max(start(downstreamNeighbor) - end(range), 0)
#' 
#' For ranges on - strand, the definition of upstream and downstream is flipped.
#' Note that this definition of distance is slightly different from 
#' `GRanges::distance`, as it effectively assumes end-coordinates are not included
#' in a range.
#' @keywords internal
distanceToNearestDirectional <- function(ranges, addArchRBug = FALSE) {
    # Label indices of the ranges incase they aren't pre-sorted
    mcols(ranges) <- NULL
    ranges$idx <- seq_along(ranges)

    # Default to distance 0
    upstream <- rep_along(ranges, 0)
    downstream <- rep_along(ranges, 0)
    
    # Filter down to ranges that aren't fully nested within others (including identical start and/or end)
    # Fully nested genes have distance 0 to neighbors.
    # Once we filter them out, we know that starts + ends are unique, and
    # sorting by starts implies also getting sorted ends.
    if (!addArchRBug)
        ranges <- ranges[GenomicRanges::countOverlaps(ranges, type="within", ignore.strand=TRUE) == 1]

    ranges_list <- split(ranges, GenomicRanges::seqnames(ranges))

    for (chr in names(ranges_list)) {
        ranges <- GenomicRanges::sort(ranges_list[[chr]], ignore.strand=TRUE)
        
        neighbor_dist <- GenomicRanges::start(ranges)[-1] - GenomicRanges::end(ranges)[-length(ranges)] 
        start_dist <- pmax(0, c(Inf, neighbor_dist))
        end_dist <- pmax(0, c(neighbor_dist, Inf))

        upstream[ranges$idx] <- ifelse(strand(ranges) == "-", end_dist, start_dist)
        downstream[ranges$idx] <- ifelse(strand(ranges) == "-", start_dist, end_dist)
    }
    data.frame(upstream=upstream, downstream=downstream)
}

#' Extend GenomicRanges in a strand-aware fashion.
#' @param ranges GenomicRanges object
#' @param upstream Number of bases to extend each range upstream (negative to shrink width)
#' @param downstream Number of bases to extend each range downstream (negative to shrink width)
#' @details Note that ranges will be blocked from extending past base 1, or the maximum
#' @keywords internal
extendGenomicRanges <- function(ranges, upstream=0, downstream=0) {
    GenomicRanges::start(ranges) <- pmax(
        1, 
        GenomicRanges::start(ranges) - ifelse(GenomicRanges::strand(ranges) == "-", downstream, upstream)
    )
    GenomicRanges::end(ranges) <- pmin(
        seqlengths(ranges)[as.character(seqnames(ranges))],
        GenomicRanges::end(ranges) + ifelse(GenomicRanges::strand(ranges) == "-", upstream, downstream),
        na.rm=TRUE
    )
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
#' @param genes GRanges object with gene start, end, and strand
#' @param tile_size Size of tiles to consider
#' @param addArchRBug Replicate ArchR bug in handling nested genes
#' @details Note: assumes the 1-based, end inclusive coordinate convention used
#' by GRanges, so 500bp tiles run from bases 1-500, 501-1000, etc.  
#' @return GRanges object with one range per tile, with additional metadata
#' columns gene_idx (row index of the gene this tile corresponds to) and
#' distance.  
#' 
#' Distance is a signed distance calculated such that if the tile has a smaller
#' start coordinate than the gene and the gene is on the + strand, distance will
#' be negative and calculated as min(0, end(tile) - start(gene))
#' @keywords internal
geneScoreTilesArchR <- function(genes, tile_size = 500, addArchRBug = FALSE) {
    # Extend upstream by 5kb
    genes <- extendGenomicRanges(genes, upstream=5000)
    GenomicRanges::mcols(genes) <- NULL
    # Extend up 1kb - 100kb upstream + downstream, but stop if we hit a neighboring gene.
    extension_lengths <- distanceToNearestDirectional(genes, addArchRBug)
    extension_lengths$upstream <- pmin(pmax(extension_lengths$upstream - tile_size, 1000), 100000)
    extension_lengths$downstream <- pmin(pmax(extension_lengths$downstream - tile_size, 1000), 100000)
    extended <- extendGenomicRanges(genes, upstream=extension_lengths$upstream, downstream=extension_lengths$downstream)
    # 1. Get a list of tile start coords 
    first_tile <- GenomicRanges::start(extended) %/% tile_size
    last_tile <- GenomicRanges::end(extended) %/% tile_size
    tile_count <- last_tile - first_tile + 1
    # 2. Set up for vectorized operation: tile_idx counts up from 1 to tile_count[1], then restarts and counts
    #    from 1 to tile_count[2], etc. There is 1 element for each tile which has a non-zero weight for a gene score
    region_idx <- rep.int(seq_along(genes), times=tile_count)
    tile_idx <- seq_len(sum(tile_count)) - rep.int(cumsum(tile_count) - tile_count, times=tile_count)
    # 3. Calculate distances to the gene "core" region
    tile_start <- (tile_idx + rep(first_tile, times=tile_count) - 1) * tile_size + 1

    distance <- ifelse(
        tile_start > GenomicRanges::end(genes)[region_idx], 
        tile_start - GenomicRanges::end(genes)[region_idx],
        pmin(0, tile_start + tile_size - 1 - GenomicRanges::start(genes)[region_idx])
    ) 
    distance <- distance * ifelse(GenomicRanges::strand(genes)[region_idx] == "-", -1, 1)
    tiles <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(genes)[region_idx],
        ranges = IRanges::IRanges(
            start = tile_start,
            width = 1
        ),
        seqinfo = seqinfo(genes),
        gene_idx = region_idx,
        distance = distance
    )
    return(extendGenomicRanges(tiles, downstream=tile_size - 1))
}

#' Calculate GeneActivityScore weights matrix
#' 
#' Given a set of tile coordinates and distances returned by `geneScoreTilesArchR`,
#' calculate a weight matrix of dimensions genes x tiles. This matrix can be 
#' multiplied with a tile matrix to obtain ArchR-compatible gene activity scores.
#' 
#' @param tiles Tile coordinates as returned by `geneScoreTilesArchR`
#' @param chr_sizes Named vector of chromosome sizes with chromosmes in the order they are listed
#' in the tile matrix (used to calculate the index of each tile)
#' @return Weight matrix of dimension genes x tiles
#' @keywords internal
geneScoreWeightsArchR <- function(tiles, chr_sizes = GenomeInfoDb::seqlengths(tiles)) {
    tile_size <- width(tiles[1])
    tile_count <- chr_sizes %/% tile_size + 1
    start_tile <- cumsum(tile_count) - tile_count
    sparseMatrix(
        i = tiles$gene_idx,
        j = start_tile[as.character(GenomicRanges::seqnames(tiles))] + GenomicRanges::start(tiles) %/% tile_size + 1,
        x = exp(-abs(tiles$distance)/5000) + exp(-1),
        dims = c(max(tiles$gene_idx), sum(tile_count))
    )
}

