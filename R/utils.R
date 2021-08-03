
#' Calculate footprinting matrix 
#' 
#' @param fragments IterableFragments object
#' @param motif_positions GRanges object with motif positions. Footprinting will
#'   be centered around the middle of these ranges
#' @param cell_groups Character or factor assigning a group to each cell, in order of
#'   `cellNames(fragments)`
#' @param cell_normalization_factors Numeric vector assigning normalization factors
#'   (e.g. total reads) to each cell, in order of `cellNames(fragments)`
#' @param flank Number of flanking basepairs to include on either side of the motif
#' 
#' @return Matrix of dimensions cell_groups x basepairs (2*flank + 1). Each entry
#'   in the matrix contains the average normalized insertions for a given cell_group
#'   at a given position
#' @export
getFootprints <- function(fragments, motif_positions, cell_groups, 
                         cell_normalization_factors = rep_along(cell_groups, 1), flank=125) {
    assert_is(fragments, "IterableFragments")
    assert_is(motif_positions, "GRanges")
    assert_is(cell_groups, c("character", "factor"))
    assert_len(cell_groups, length(cellNames(fragments)))

    assert_is(cell_normalization_factors, c("numeric"))
    assert_len(cell_normalization_factors, length(cellNames(fragments)))
    assert_wholenumber(flank)

    cell_groups <- as.factor(cell_groups)
    motif_positions <- GenomicRanges::resize(motif_positions, 2*flank + 1, fix="center")

    groups_matrix <- matrix(0, nrow=length(levels(cell_groups)), ncol=length(cell_groups))
    groups_indices <- matrix(ncol=2, c(
        as.integer(cell_groups),
        seq_along(cell_groups)
    ))
    groups_matrix[groups_indices] <- cell_normalization_factors[groups_indices[,2]]

    # Get a matrix of cell_groups x all positions
    motif_bp <- unlist(GenomicRanges::slidingWindows(motif_positions, 1))
    mat_all <- groups_matrix %*% overlapMatrix(fragments, motif_bp)
    mat_all <- mat_all / table(cell_groups)

    # Aggregate positions together: first reshape the matrix so each
    # column is all cell_groups and bp from within a motif position,
    # average across the columns, then re-aggregate
    mat_res <- matrix(
        colSums(matrix(mat_all, ncol=length(motif_positions))),
        ncol=2*flank+1
    )
    rownames(mat_res) <- levels(cell_groups)
    colnames(mat_res) <- as.character(-flank:flank)
    mat_res
}
