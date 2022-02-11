
#' IterableMatrix methods
#'
#' Methods for IterableMatrix objects
#'
#' @name IterableMatrix-methods
#' @rdname IterableMatrix-methods
NULL

setClass("IterableMatrix", 
    slots = c(
        dim = "integer",
        transpose = "logical"
    ),
    prototype = list(
        dim = integer(2),
        transpose = FALSE
    )
)

setGeneric("iterate_matrix", function(x) standardGeneric("iterate_matrix"))
setMethod("iterate_matrix", "externalptr", function(x) { return(x)} )

setMethod("short_description", "IterableMatrix", function(x) {
    character(0)
})


#' @describeIn IterableMatrix-methods Display an IterableMatrix
#' @param object IterableMatrix object
setMethod("show", "IterableMatrix", function(object) {
    cat(sprintf("%d x %d IterableMatrix object with class %s\n", nrow(object), ncol(object), class(object)))

    cat("\n")
    description <- short_description(object)
    if (length(description) > 0) cat("Queued Operations:\n")
    for (i in seq_along(description)) {
        cat(sprintf("%d. %s\n", i, description[i]))
    }
})


#' @describeIn IterableMatrix-methods Transpose an IterableMatrix
#' @param x IterableMatrix object
#' @return * `t()` Transposed object
setMethod("t", signature(x = "IterableMatrix"), function(x) {
    x@transpose <- !x@transpose
    x@dim <- c(x@dim[2L], x@dim[1L])
    return(x)
})


setClass("mat_uint32_t",
    contains = "IterableMatrix"
)
setClass("mat_double",
    contains = "IterableMatrix"
)

# Dense multiply operators (sparse*dense_mat and sparse*dense_vec)

#' @param x IterableMatrix object
#' @param y matrix
#' @describeIn IterableMatrix-methods Multiply by a dense matrix
#' @return * `x %*% y`: dense matrix result
setMethod("%*%", signature(x="mat_double", y="matrix"), function(x, y) {
    iter <- iterate_matrix(x)
    if(x@transpose) {
        return(t(dense_multiply_left_cpp(iter, t(y))))
    } else {
        return(dense_multiply_right_cpp(iter, y))
    }
})

setMethod("%*%", signature(x="matrix", y="mat_double"), function(x, y) {
    iter <- iterate_matrix(y)
    if(y@transpose) {
        return(t(dense_multiply_right_cpp(iter, t(x))))
    } else {
        return(dense_multiply_left_cpp(iter, x))
    }
})

setMethod("%*%", signature(x="mat_double", y="numeric"), function(x, y) {
    iter <- iterate_matrix(x)
    if(x@transpose) {
        return(vec_multiply_left_cpp(iter, y))
    } else {
        return(vec_multiply_right_cpp(iter, y))
    }
})

setMethod("%*%", signature(x="numeric", y="mat_double"), function(x, y) {
    iter <- iterate_matrix(y)
    if(y@transpose) {
        return(vec_multiply_right_cpp(iter, x))
    } else {
        return(vec_multiply_left_cpp(iter, x))
    }
})

setMethod("%*%", signature(x="mat_uint32_t", y="matrix"), function(x, y) as(x, "mat_double") %*% y)
setMethod("%*%", signature(x="matrix", y="mat_uint32_t"), function(x, y) x %*% as(y, "mat_double"))
setMethod("%*%", signature(x="mat_uint32_t", y="numeric"), function(x, y) as(x, "mat_double") %*% y)
setMethod("%*%", signature(x="numeric", y="mat_uint32_t"), function(x, y) x %*% as(y, "mat_double"))

# Row sums and row means

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate rowSums
#' @return * `rowSums()`: vector of row sums
setMethod("rowSums", signature(x="mat_double"), function(x) {
    iter <- iterate_matrix(x)
    if (x@transpose) 
        col_sums_cpp(iter)
    else
        row_sums_cpp(iter)
})

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate colSums
#' @return * `colSums()`: vector of col sums
setMethod("colSums", signature(x="mat_double"), function(x) {
    iter <- iterate_matrix(x)
    if (x@transpose) 
        row_sums_cpp(iter)
    else
        col_sums_cpp(iter)
})

setMethod("rowSums", signature(x="mat_uint32_t"), function(x) rowSums(as(x, "mat_double")))
setMethod("colSums", signature(x="mat_uint32_t"), function(x) colSums(as(x, "mat_double")))

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate rowMeans
#' @return * `rowMeans()`: vector of row means
setMethod("rowMeans", signature(x="IterableMatrix"), function(x) rowSums(x)/ncol(x))

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate colMeans
#' @return * `colMeans()`: vector of col means
setMethod("colMeans", signature(x="IterableMatrix"), function(x) colSums(x)/nrow(x))

# Packed integer matrix
setClass("PackedMatrixMem",
    contains = "mat_uint32_t",
    slots = c(
        val_data = "integer",   
        val_idx = "integer",    
        row_data = "integer",   
        row_starts = "integer", 
        row_idx = "integer",
        col_ptr = "integer",   
        row_count = "integer"
    ),
    prototype = list(
        val_data = integer(0),   
        val_idx = integer(0),    
        row_data = integer(0),   
        row_starts = integer(0), 
        row_idx = integer(0),
        col_ptr = integer(0),   
        row_count = integer(0)
    )
)
setMethod("iterate_matrix", "PackedMatrixMem", function(x) {
    iterate_packed_matrix_cpp(x)
})
setMethod("short_description", "PackedMatrixMem", function(x) {
    "Compressed integer matrix in memory"
})

setClass("UnpackedMatrixMem",
    contains = "mat_uint32_t",
    slots = c(
        val = "integer",   
        row = "integer",    
        col_ptr = "integer",   
        row_count = "integer"
    ),
    prototype = list(
        val = integer(0),   
        row = integer(0),   
        col_ptr = integer(0),   
        row_count = integer(0)
    )
)
setMethod("iterate_matrix", "UnpackedMatrixMem", function(x) {
    iterate_unpacked_matrix_cpp(x)
})
setMethod("short_description", "UnpackedMatrixMem", function(x) {
    "Uncompressed integer matrix in memory"
})

write_matrix_memory <- function(mat, compress=TRUE) {
    assert_is(mat, c("IterableMatrix", "dgCMatrix"))
    if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")
    if (is(mat, "mat_double")) mat <- as(mat, "mat_uint32_t")
    assert_true(mat@transpose == FALSE)

    if (compress) {
        m <- write_packed_matrix_cpp(iterate_matrix(mat))
        res <- do.call(new, c("PackedMatrixMem", m))
        res@dim <- mat@dim
    } else {
        m <- write_unpacked_matrix_cpp(iterate_matrix(mat))
        res <- do.call(new, c("UnpackedMatrixMem", m))
        res@dim <- mat@dim
    }
    res
}

setClass("MatrixDir",
    contains = "mat_uint32_t",
    slots = c(
        dir = "character",
        compressed = "logical",
        buffer_size = "integer"
    ),
    prototype = list(
        dir = character(0),
        compressed = logical(0),
        buffer_size = integer(0)
    )
)
setMethod("iterate_matrix", "MatrixDir", function(x) {
    if (x@compressed)
        iterate_packed_matrix_file_cpp(x@dir, x@buffer_size)
    else
        iterate_unpacked_matrix_file_cpp(x@dir, x@buffer_size)
})
setMethod("short_description", "MatrixDir", function(x) {
    sprintf("%s integer matrix in directory %s", 
        if(x@compressed) "compressed" else "uncompressed",
        x@dir    
    )
})
write_matrix_dir <- function(mat, dir, compress=TRUE, buffer_size=8192L) {
    assert_is(mat, c("IterableMatrix", "dgCMatrix"))
    if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")
    if (is(mat, "mat_double")) mat <- as(mat, "mat_uint32_t")
    assert_true(mat@transpose == FALSE)
    assert_is(dir, "character")
    assert_is(compress, "logical")
    assert_is(buffer_size, "integer")
    
    dir <- path.expand(dir)

    if (compress)
        write_packed_matrix_file_cpp(iterate_matrix(mat), dir, buffer_size)
    else 
        write_unpacked_matrix_file_cpp(iterate_matrix(mat), dir, buffer_size)
    
    open_matrix_dir(dir, buffer_size)
}
open_matrix_dir <- function(dir, buffer_size=8192L) {
    assert_is_file(dir)
    assert_is(buffer_size, "integer")
    
    dir <- path.expand(dir)
    info <- dims_matrix_file_cpp(dir, buffer_size)
    new("MatrixDir", dir=dir, dim=info$dims, compressed=info$compressed, buffer_size=buffer_size)
}

setClass("MatrixH5",
    contains = "mat_uint32_t",
    slots = c(
        path = "character",
        group = "character",
        compressed = "logical",
        buffer_size = "integer"
    ),
    prototype = list(
        path = character(0),
        group = "",
        compressed = logical(0),
        buffer_size = integer(0)
    )
)
setMethod("iterate_matrix", "MatrixH5", function(x) {
    if (x@compressed)
        iterate_packed_matrix_hdf5_cpp(x@path, x@group, x@buffer_size)
    else
        iterate_unpacked_matrix_hdf5_cpp(x@path, x@group, x@buffer_size)
})
setMethod("short_description", "MatrixH5", function(x) {
    sprintf("%s integer matrix in hdf5 file %s, group %s", 
        if(x@compressed) "Compressed" else "Uncompressed",
        x@path, 
        x@group
    )
})
write_matrix_hdf5 <- function(mat, path, group, compress=TRUE, buffer_size=8192L, chunk_size=1024L) {
    assert_is(mat, c("IterableMatrix", "dgCMatrix"))
    if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")
    if (is(mat, "mat_double")) mat <- as(mat, "mat_uint32_t")
    assert_true(mat@transpose == FALSE)
    assert_is(path, "character")
    assert_is(group, "character")
    assert_is(compress, "logical")
    assert_is(buffer_size, "integer")
    assert_is(chunk_size, "integer")

    path <- path.expand(path)

    if (compress)
        write_packed_matrix_hdf5_cpp(iterate_matrix(mat), path, group, buffer_size, chunk_size)
    else 
        write_unpacked_matrix_hdf5_cpp(iterate_matrix(mat), path, group, buffer_size, chunk_size)
    
    open_matrix_hdf5(path, group, buffer_size)
}

open_matrix_hdf5 <- function(path, group, buffer_size=8192L) {
    assert_is_file(path)
    assert_is(group, "character")
    assert_is(buffer_size, "integer")

    info <- dims_matrix_hdf5_cpp(path.expand(path), group, buffer_size)
    new("MatrixH5", path=path.expand(path), group=group, dim=info$dims, compressed=info$compressed, buffer_size=buffer_size)
}

setClass("10xMatrixH5",
    contains = "mat_uint32_t",
    slots = c(
        path = "character",
        buffer_size = "integer"
    ),
    prototype = list(
        path = character(0),
        buffer_size = integer(0)
    )
)
setMethod("iterate_matrix", "10xMatrixH5", function(x) {
    iterate_matrix_10x_hdf5_cpp(x@path, x@buffer_size)
})
setMethod("short_description", "10xMatrixH5", function(x) {
    sprintf("10x HDF5 feature matrix in file %s", 
        x@path
    )
})

open_matrix_10x_hdf5 <- function(path, buffer_size=8192L) {
    assert_is_file(path)
    assert_is(buffer_size, "integer")

    info <- dims_matrix_10x_hdf5_cpp(path.expand(path), buffer_size)
    new("10xMatrixH5", path=path.expand(path), dim=info$dims, buffer_size=buffer_size)
}


# Overlap matrix from fragments
setClass("overlapMatrix",
    contains = "mat_uint32_t",
    slots = c(
        fragments = "IterableFragments",
        chr = "character",
        start = "integer",
        end = "integer",
        version = "integer"
    ),
    prototype = list(
        fragments = NULL,
        chr = character(0),
        start = integer(0),
        end = integer(0),
        version = NA_integer_
    )
)
#' Calculate cell x ranges overlap matrix
#' @param fragments Input fragments object
#' @param ranges GRanges object with the ranges to overlap
#' @param convert_to_0_based_coords Whether to convert the ranges from a 1-based
#'        coordinate system to a 0-based coordinate system. 
#'        (see http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/)
#' @return Iterable matrix object with dimension cell x ranges
#' @export
overlapMatrix <- function(fragments, ranges, tile_width, convert_to_0_based_coords=TRUE, version = 1) {
    assert_is(fragments, "IterableFragments")
    assert_is(ranges, c("GRanges", "list", "data.frame"))
    
    assert_is(convert_to_0_based_coords, "logical")
    
    assert_not_na(fragments@cell_names)

    if(is.list(ranges)) {
        assert_has_names(ranges, c("chr", "start", "end"))
        res <- new("overlapMatrix", fragments=fragments, 
            chr=as.character(ranges$chr),
            start=as.integer(ranges$start) - convert_to_0_based_coords,
            end=as.integer(ranges$end))
    } else {
        res <- new("overlapMatrix", fragments=fragments, 
            chr=as.character(GenomicRanges::seqnames(ranges)),
            start=GenomicRanges::start(ranges) - convert_to_0_based_coords,
            end=GenomicRanges::end(ranges))
    }
    res@dim <- c(length(fragments@cell_names), length(ranges))
    res@version <- as.integer(version)
    return(res)
}

setMethod("iterate_matrix", "overlapMatrix", function(x) {
    chrs <- as.factor(x@chr)
    if(x@version == 1) {
        iterate_overlap_matrix_cpp(
            iterate_fragments(x@fragments),
            as.integer(chrs)-1,
            x@start, x@end,
            levels(chrs)
        )
    } else if (x@version == 2) {
        iterate_overlap_matrix2_cpp(
            iterate_fragments(x@fragments),
            as.integer(chrs)-1,
            x@start, x@end,
            levels(chrs)
        )
    } else if (x@version == 3) {
        iterate_overlap_matrix3_cpp(
            iterate_fragments(x@fragments),
            as.integer(chrs)-1,
            x@start, x@end,
            levels(chrs)
        )
    } else if (x@version == 4) {
        iterate_overlap_matrix4_cpp(
            iterate_fragments(x@fragments),
            as.integer(chrs)-1,
            x@start, x@end,
            levels(chrs)
        )
    } else if (x@version == 5) {
       iterate_overlap_matrix5_cpp(
            iterate_fragments(x@fragments),
            as.integer(chrs)-1,
            x@start, x@end,
            levels(chrs)
        ) 
    } else {
        iterate_overlap_matrix6_cpp(
            iterate_fragments(x@fragments),
            as.integer(chrs)-1,
            x@start, x@end,
            levels(chrs)
        ) 
    }
})
setMethod("short_description", "overlapMatrix", function(x) {
    # Subset strings first to avoid a very slow string concatenation process
    indices <- c(head(seq_along(x@chr), 3), tail(seq_along(x@chr), 1))
    labels <- paste0(x@chr[indices], ":", x@start[indices]+1, "-", x@end[indices])
    c(
        short_description(x@fragments),
        sprintf("Calculate overlaps with %d ranges%s", length(x@chr), 
            pretty_print_vector(labels, prefix=": ", max_len=2)
        )
    )
})

# Tile matrix from fragments
setClass("tileMatrix",
    contains = "mat_uint32_t",
    slots = c(
        fragments = "IterableFragments",
        chr = "character",
        start = "integer",
        end = "integer",
        tile_width = "integer",
        version = "integer"
    ),
    prototype = list(
        fragments = NULL,
        chr = character(0),
        start = integer(0),
        end = integer(0),
        tile_width = integer(0)
    )
)
#' Calculate cell x tiles overlap matrix
#' @param fragments Input fragments object
#' @param ranges GRanges object with the ranges to overlap
#' @param convert_to_0_based_coords Whether to convert the ranges from a 1-based
#'        coordinate system to a 0-based coordinate system. 
#'        (see http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/)
#' @return Iterable matrix object with dimension cells x tiles
#' @export
tileMatrix <- function(fragments, ranges, tile_width, convert_to_0_based_coords=TRUE, version=1L) {
    assert_is(fragments, "IterableFragments")
    assert_not_na(fragments@cell_names)

    assert_is(convert_to_0_based_coords, "logical")

    ranges <- normalize_ranges(ranges)

    assert_wholenumber(tile_width)
    assert_greater_than_zero(tile_width)
    tile_width <- normalize_length(tile_width, length(ranges$start))

    
    res <- new("tileMatrix", fragments=fragments,
        chr = as.character(ranges$chr),
        start = as.integer(ranges$start) - convert_to_0_based_coords,
        end = as.integer(ranges$end),
        tile_width = as.integer(tile_width),
        version = as.integer(version)
    )
    
    tiles <- ceiling((res@end - res@start) / res@tile_width)
    res@dim <- c(length(fragments@cell_names), as.integer(sum(tiles)))
    return(res)
}

setMethod("iterate_matrix", "tileMatrix", function(x) {
    chrs <- as.factor(x@chr)

    if (x@version == 1) {
        iterate_tile_matrix_cpp(
            iterate_fragments(x@fragments),
            as.integer(chrs)-1,
            x@start, x@end,
            x@tile_width,
            levels(chrs)
        )
    } else {
        iterate_tile_matrix2_cpp(
            iterate_fragments(x@fragments),
            as.integer(chrs)-1,
            x@start, x@end,
            x@tile_width,
            levels(chrs)
        )
    }
})

setMethod("short_description", "tileMatrix", function(x) {
    # Subset strings first to avoid a very slow string concatenation process
    indices <- c(head(seq_along(x@chr), 3), tail(seq_along(x@chr), 1))
    labels <- paste0(x@chr[indices], ":", x@start[indices]+1, "-", x@end[indices], " (", x@tile_width[indices], "bp)")
    c(
        short_description(x@fragments),
        sprintf("Calculate %d tiles over %d ranges%s", ncol(x), length(x@chr), 
            pretty_print_vector(labels, prefix=": ", max_len=2)
        )
    )
})

# Conversions between integer and double IterableMatrix
setClass("convertMatrixIntDouble",
    contains = "mat_double",
    slots = c(
        mat = "mat_uint32_t"
    ),
    prototype = list(
        mat = NULL
    )
)
setMethod("iterate_matrix", "convertMatrixIntDouble", function(x) {
    convert_matrix_uint32_t_double_cpp(iterate_matrix(x@mat))
})
setMethod("short_description", "convertMatrixIntDouble", function(x) {
    c(
        short_description(x@mat),
        "Convert integers to doubles"
    )
})
setAs("mat_uint32_t", "mat_double", function(from) {
    new("convertMatrixIntDouble", dim=from@dim, transpose=from@transpose, mat=from)
})

setClass("convertMatrixDoubleInt",
    contains = "mat_uint32_t",
    slots = c(
        mat = "mat_double"
    ),
    prototype = list(
        mat = NULL
    )
)
setMethod("iterate_matrix", "convertMatrixDoubleInt", function(x) {
    convert_matrix_double_uint32_t_cpp(iterate_matrix(x@mat))
})
setMethod("short_description", "convertMatrixDoubleInt", function(x) {
    c(
        short_description(x@mat),
        "Convert doubles to integers"
    )
})
setAs("mat_double", "mat_uint32_t", function(from) {
    new("convertMatrixDoubleInt", dim=from@dim, transpose=from@transpose, mat=from)
})

# Conversions with dgCMatrix
setClass("Iterable_dgCMatrix_wrapper",
    contains = "mat_double",
    slots = c(
        mat = "dgCMatrix"
    ),
    prototype = list(
        mat=NULL
    )
)
setAs("dgCMatrix", "IterableMatrix", function(from) {
    new("Iterable_dgCMatrix_wrapper", dim=dim(from), transpose=FALSE, mat=from)
})
setMethod("iterate_matrix", "Iterable_dgCMatrix_wrapper", function(x) {
    iterate_csparse_matrix_cpp(x@mat)
})
setMethod("short_description", "convertMatrixDoubleInt", function(x) {
    "Load dgCMatrix from memory"
})

setMethod("iterate_matrix", "dgCMatrix", function(x) {
    iterate_csparse_matrix_cpp(x)
})

setAs("mat_uint32_t", "dgCMatrix", function(from) {
    iter <- iterate_matrix(from)
    iter <- convert_matrix_uint32_t_double_cpp(iter)
    res <- build_csparse_matrix_double_cpp(iter)
    if (from@transpose) {
        res <- t(res)
    }
    return(res)
})

setAs("mat_double", "dgCMatrix", function(from) {
    iter <- iterate_matrix(from)
    res <- build_csparse_matrix_double_cpp(iter)
    if (from@transpose) {
        res <- t(res)
    }
    return(res)
})

#' Calculate matrix stats
#' @param matrix Input matrix object
#' @param row_stats Which row statistics to compute 
#' @param col_stats Which col statistics to compute 
#' @return List of row_stats: matrix of n_stats x n_rows, 
#'          col_stats: matrix of n_stats x n_cols
#' @details The statistics will be calculated in a single pass over the matrix,
#' so this method is desirable to use for efficiency purposes compared to 
#' the more standard rowMeans or colMeans if multiple statistics are needed.
#' If variance is calculated, then mean and nonzero count will be included in the
#' output, and if mean is calculated then nonzero count will be included in the output.
#' @export
matrixStats <- function(matrix, 
                        row_stats = c("none", "nonzero", "mean", "variance"), 
                        col_stats = c("none", "nonzero", "mean", "variance")) {

    stat_options <- c("none", "nonzero", "mean", "variance")
    row_stats <- match.arg(row_stats)
    col_stats <- match.arg(col_stats)

    if (matrix@transpose) {
        tmp <- row_stats 
        row_stats <- col_stats
        col_stats <- tmp
    }

    row_stats_number <- match(row_stats, stat_options) - 1
    col_stats_number <- match(col_stats, stat_options) - 1
    
    res <- matrix_stats_cpp(iterate_matrix(matrix), row_stats_number, col_stats_number)
    rownames(res$row_stats) <- stat_options[seq_len(row_stats_number) + 1]
    rownames(res$col_stats) <- stat_options[seq_len(col_stats_number) + 1]
    
    if (matrix@transpose) {
        res <- list(
            row_stats = res$col_stats,
            col_stats = res$row_stats
        )
    }
    return(res)
}