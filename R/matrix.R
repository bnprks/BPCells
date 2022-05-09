
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
        transpose = "logical",
        dimnames = "list"
    ),
    prototype = list(
        dim = integer(2),
        transpose = FALSE,

        dimnames = list(NULL, NULL)
    )
)

#' Construct an S4 matrix object wrapping another matrix object, avoiding
#' duplicate storage of dimnames
wrapMatrix <- function(class, m, ...) {
    dimnames <- dimnames(m)
    m@dimnames <- list(NULL, NULL)
    new(class, matrix=m, transpose=m@transpose, dim=m@dim, dimnames=dimnames, ...)
}

#' Helper function to set dimnames to NULL instead of 0-length character vectors
normalized_dimnames <- function(row_names, col_names) {
    list(
        if(length(row_names) == 0) NULL else row_names,
        if(length(col_names) == 0) NULL else col_names
    )
}
denormalize_dimnames <- function(dimnames) {
    if (is.null(dimnames[[1]])) dimnames[[1]] <- character(0)
    if (is.null(dimnames[[2]])) dimnames[[2]] <- character(0)
    dimnames
}

setGeneric("iterate_matrix", function(x) standardGeneric("iterate_matrix"))
setMethod("iterate_matrix", "externalptr", function(x) { return(x)} )

setGeneric("matrix_type", function(x) standardGeneric("matrix_type"))

setMethod("short_description", "IterableMatrix", function(x) {
    character(0)
})


#' @describeIn IterableMatrix-methods Display an IterableMatrix
#' @param object IterableMatrix object
setMethod("show", "IterableMatrix", function(object) {
    cat(sprintf("%d x %d IterableMatrix object with class %s\n", nrow(object), ncol(object), class(object)))

    cat("\n")

    if(!is.null(rownames(object)))
        cat(sprintf("Row names: %s\n",
            pretty_print_vector(rownames(object), empty="unknown names")
        ))

    if(!is.null(rownames(object)))
        cat(sprintf("Col names: %s\n",
            pretty_print_vector(colnames(object), empty="unknown names")
        ))

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
    x@dimnames <- list(x@dimnames[[2]], x@dimnames[[1]])
    return(x)
})


setClass("mat_uint32_t",
    contains = "IterableMatrix"
)
setMethod("matrix_type", signature(x="mat_uint32_t"), function(x) "mat_uint32_t")
setClass("mat_double",
    contains = "IterableMatrix"
)
setMethod("matrix_type", signature(x="mat_double"), function(x) "mat_double")

# Dense multiply operators (sparse*dense_mat and sparse*dense_vec)

#' @param x IterableMatrix object
#' @param y matrix
#' @describeIn IterableMatrix-methods Multiply by a dense matrix
#' @return * `x %*% y`: dense matrix result
setMethod("%*%", signature(x="IterableMatrix", y="matrix"), function(x, y) {
    iter <- iterate_matrix(cast_matrix_double(x))
    if(x@transpose) {
        return(t(dense_multiply_left_cpp(iter, t(y))))
    } else {
        return(dense_multiply_right_cpp(iter, y))
    }
})

setMethod("%*%", signature(x="matrix", y="IterableMatrix"), function(x, y) {
    iter <- iterate_matrix(cast_matrix_double(y))
    if(y@transpose) {
        return(t(dense_multiply_right_cpp(iter, t(x))))
    } else {
        return(dense_multiply_left_cpp(iter, x))
    }
})

setMethod("%*%", signature(x="IterableMatrix", y="numeric"), function(x, y) {
    iter <- iterate_matrix(cast_matrix_double(x))
    if(x@transpose) {
        return(vec_multiply_left_cpp(iter, y))
    } else {
        return(vec_multiply_right_cpp(iter, y))
    }
})

setMethod("%*%", signature(x="numeric", y="IterableMatrix"), function(x, y) {
    iter <- iterate_matrix(cast_matrix_double(y))
    if(y@transpose) {
        return(vec_multiply_right_cpp(iter, x))
    } else {
        return(vec_multiply_left_cpp(iter, x))
    }
})

# Row sums and row means

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate rowSums
#' @return * `rowSums()`: vector of row sums
setMethod("rowSums", signature(x="IterableMatrix"), function(x) {
    iter <- iterate_matrix(cast_matrix_double(x))
    if (x@transpose) 
        res <- col_sums_double_cpp(iter)
    else
        res <- row_sums_double_cpp(iter)
    names(res) <- rownames(x)
    res
})

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate colSums
#' @return * `colSums()`: vector of col sums
setMethod("colSums", signature(x="IterableMatrix"), function(x) {
    iter <- iterate_matrix(cast_matrix_double(x))
    if (x@transpose) 
        res <- row_sums_double_cpp(iter)
    else
        res <- col_sums_double_cpp(iter)
    names(res) <- colnames(x)
    res
})

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate rowMeans
#' @return * `rowMeans()`: vector of row means
setMethod("rowMeans", signature(x="IterableMatrix"), function(x) rowSums(x)/ncol(x))

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate colMeans
#' @return * `colMeans()`: vector of col means
setMethod("colMeans", signature(x="IterableMatrix"), function(x) colSums(x)/nrow(x))


# Index subsetting
setClass("MatrixSubset",
    contains = "IterableMatrix",
    slots = c(
        matrix = "IterableMatrix",
        row_selection = "integer",
        col_selection = "integer"
    ),
    prototype = list(
        row_selection = integer(0),
        col_selection = integer(0)
    )
)
setMethod("matrix_type", signature(x="MatrixSubset"), function(x) matrix_type(x@matrix))

setMethod("[", "IterableMatrix", function(x, i, j, ...) {
    if (missing(x)) stop("x is missing in matrix selection")
    ret <- wrapMatrix("MatrixSubset", x)

    if(!missing(i)) {
        indices <- seq_len(nrow(x))
        names(indices) <- rownames(x)
        ret@row_selection <- vctrs::vec_slice(indices, i)
        if(!is.null(ret@dimnames[[1]])) ret@dimnames[[1]] <- ret@dimnames[[1]][ret@row_selection]
        if (!is.logical(i)) assert_distinct(i)
        ret@dim[1] <- length(ret@row_selection)
    }
    if (!missing(j)) {
        indices <- seq_len(ncol(x))
        names(indices) <- colnames(x)
        ret@col_selection <- vctrs::vec_slice(indices, j)
        if(!is.null(ret@dimnames[[2]])) ret@dimnames[[2]] <- ret@dimnames[[2]][ret@col_selection]
        assert_distinct(ret@col_selection)
        if (!is.logical(j)) assert_distinct(j)
        ret@dim[2] <- length(ret@col_selection)
    }
    ret
})

# Handle chained subsets by just modifying the selections
setMethod("[", "MatrixSubset", function(x, i, j, ...) {
    if (missing(x)) stop("x is missing in matrix selection")
    ret <- wrapMatrix("MatrixSubset", x@matrix)
    ret@dim <- dim(x)

    if(!missing(i)) {
        indices <- seq_len(nrow(x))
        names(indices) <- rownames(x)
        selection <- vctrs::vec_slice(indices, i)
        ret@row_selection <- x@row_selection[selection]
        if(!is.null(ret@dimnames[[1]])) ret@dimnames[[1]] <- ret@dimnames[[1]][selection]
        if (!is.logical(i)) assert_distinct(i)
        ret@dim[1] <- length(ret@row_selection)
    }
    if (!missing(j)) {
        indices <- seq_len(ncol(x))
        names(indices) <- colnames(x)
        selection <- vctrs::vec_slice(indices, j)
        ret@col_selection <- x@col_selection[selection]
        if(!is.null(ret@dimnames[[2]])) ret@dimnames[[2]] <- ret@dimnames[[2]][selection]
        assert_distinct(ret@col_selection)
        if (!is.logical(j)) assert_distinct(j)
        ret@dim[2] <- length(ret@col_selection)
    }
    ret
})

setMethod("iterate_matrix", "MatrixSubset", function(x) {
    ret <- iterate_matrix(x@matrix)

    if (matrix_type(x) == "mat_double") {
        if (length(x@row_selection) != 0) ret <- iterate_matrix_row_select_double_cpp(ret, x@row_selection-1L)
        if (length(x@col_selection) != 0) ret <- iterate_matrix_col_select_double_cpp(ret, x@col_selection-1L)
    } else if (matrix_type(x) == "mat_uint32_t") {
        if (length(x@row_selection) != 0) ret <- iterate_matrix_row_select_uint32_t_cpp(ret, x@row_selection-1L)
        if (length(x@col_selection) != 0) ret <- iterate_matrix_col_select_uint32_t_cpp(ret, x@col_selection-1L)
    } else {
        stop("Unrecognized matrix_type in MatrixSubset")
    }
    ret
})

setMethod("short_description", "MatrixSubset", function(x) {
    c(
        short_description(x@matrix),
        sprintf("Select rows: %s and cols: %s",
            pretty_print_vector(x@row_selection, empty="all"),
            pretty_print_vector(x@col_selection, empty="all"))
    )
})

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
        row_count = "integer",

        version = "character"
    ),
    prototype = list(
        val_data = integer(0),   
        val_idx = integer(0),    
        row_data = integer(0),   
        row_starts = integer(0), 
        row_idx = integer(0),
        col_ptr = integer(0),   
        row_count = integer(0),

        version = character(0)
    )
)
setMethod("iterate_matrix", "PackedMatrixMem", function(x) {
    x@dimnames <- denormalize_dimnames(x@dimnames)
    iterate_packed_matrix_cpp(x, x@dimnames[[1]], x@dimnames[[2]])
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
        row_count = "integer",

        version="character"
    ),
    prototype = list(
        val = integer(0),   
        row = integer(0),   
        col_ptr = integer(0),   
        row_count = integer(0),

        version=character(0)
    )
)
setMethod("iterate_matrix", "UnpackedMatrixMem", function(x) {
    x@dimnames <- denormalize_dimnames(x@dimnames)
    iterate_unpacked_matrix_cpp(x, x@dimnames[[1]], x@dimnames[[2]])
})
setMethod("short_description", "UnpackedMatrixMem", function(x) {
    "Uncompressed integer matrix in memory"
})

#' Make an integer sparse matrix object in memory. 
#' @param matrix Input matrix, either IterableMatrix or dgCMatrix
#' @param compress Whether or not to compress the data. 
#' @return MatrixMem object
#' @details This function will convert non-integer numbers to integers, so should only
#' be used on integer-valued matrices. With compression, memory usage 
#' for a typical RNA counts matrix should be about 6-8x smaller than an R dgCMatrix,
#' and 4-6x smaller than a scipy csc_matrix. The compression will be more effective when
#' the count values in the matrix are small, and when the rows of the matrix are 
#' sorted by rowMeans. In tests on RNA-seq data optimal ordering could save up to
#' 40% of storage space.
#' @export
write_matrix_memory <- function(mat, compress=TRUE) {
    assert_is(mat, c("IterableMatrix", "dgCMatrix"))
    if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")
    mat <- cast_matrix_int(mat)
    assert_true(mat@transpose == FALSE)


    if (compress) {
        m <- write_packed_matrix_cpp(iterate_matrix(mat))
        class <- "PackedMatrixMem" 
    } else {
        m <- write_unpacked_matrix_cpp(iterate_matrix(mat))
        class <- "UnpackedMatrixMem"
    }
    m[["dimnames"]] <- normalized_dimnames(m$row_names, m$col_names)
    m$row_names <- NULL
    m$col_names <- NULL
    res <- do.call(new, c(class, m))
    res@dim <- mat@dim
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
    x@dimnames <- denormalize_dimnames(x@dimnames)
    if (x@compressed)
        iterate_packed_matrix_file_cpp(x@dir, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]])
    else
        iterate_unpacked_matrix_file_cpp(x@dir, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]])
})
setMethod("short_description", "MatrixDir", function(x) {
    sprintf("%s integer matrix in directory %s", 
        if(x@compressed) "compressed" else "uncompressed",
        x@dir    
    )
})


#' Make an integer sparse matrix object in binary files within a directory
#' @inheritParams write_matrix_memory
#' @param dir Directory to save the data into
#' @param buffer_size For performance tuning only. The number of items to be buffered
#' in memory before calling writes to disk.
#' @inherit write_matrix_memory details
#' @return MatrixDir object
#' @export
write_matrix_dir <- function(mat, dir, compress=TRUE, buffer_size=8192L) {
    assert_is(mat, c("IterableMatrix", "dgCMatrix"))
    if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")
    mat <- cast_matrix_int(mat)
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

#' Open an existing integer sparse matrix object written with [write_matrix_dir]
#' @param dir Directory to read the data from
#' @param buffer_size For performance tuning only. The number of items to be buffered
#' in memory before calling writes to disk.
#' @export
open_matrix_dir <- function(dir, buffer_size=8192L) {
    assert_is_file(dir)
    assert_is(buffer_size, "integer")
    
    dir <- path.expand(dir)
    info <- dims_matrix_file_cpp(dir, buffer_size)
    new("MatrixDir", dir=dir, dim=info$dims, compressed=info$compressed, buffer_size=buffer_size,
        dimnames=normalized_dimnames(info$row_names, info$col_names))
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
    x@dimnames <- denormalize_dimnames(x@dimnames)
    if (x@compressed)
        iterate_packed_matrix_hdf5_cpp(x@path, x@group, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]])
    else
        iterate_unpacked_matrix_hdf5_cpp(x@path, x@group, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]])
})
setMethod("short_description", "MatrixH5", function(x) {
    sprintf("%s integer matrix in hdf5 file %s, group %s", 
        if(x@compressed) "Compressed" else "Uncompressed",
        x@path, 
        x@group
    )
})

#' Write a sparse integer matrix to an hdf5 file
#' @inheritParams write_fragments_hdf5
#' @inheritParams write_matrix_dir
#' @inherit write_matrix_memory details
#' @return MatrixH5 object
#' @export
write_matrix_hdf5 <- function(mat, path, group, compress=TRUE, buffer_size=8192L, chunk_size=1024L) {
    assert_is(mat, c("IterableMatrix", "dgCMatrix"))
    if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")
    mat <- cast_matrix_int(mat)
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

#' Read a sparse integer matrix from an hdf5 file
#' @inheritParams open_matrix_dir
#' @inheritParams open_fragments_hdf5
#' @return MatrixH5 object
#' @export
open_matrix_hdf5 <- function(path, group, buffer_size=16384L) {
    assert_is_file(path)
    assert_is(group, "character")
    assert_is(buffer_size, "integer")

    info <- dims_matrix_hdf5_cpp(path.expand(path), group, buffer_size)
    new("MatrixH5", path=path.expand(path), group=group, dim=info$dims, compressed=info$compressed, buffer_size=buffer_size,
        dimnames=normalized_dimnames(info$row_names, info$col_names))
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
    sprintf("10x HDF5 feature matrix in file %s", x@path)
})

#' Read a sparse integer matrix from a 10x formatted hdf5 feature matrix file
#' @inheritParams open_matrix_hdf5
#' @return 10xMatrixH5 object
#' @details The 10x format makes use of gzip compression for the matrix data,
#' which can slow down read performance. Consider writing into another format
#' if the read performance is important to you.
#' @export
open_matrix_10x_hdf5 <- function(path, buffer_size=16384L) {
    assert_is_file(path)
    assert_is(buffer_size, "integer")

    info <- dims_matrix_10x_hdf5_cpp(path.expand(path), buffer_size)
    new("10xMatrixH5", path=path.expand(path), dim=info$dims, buffer_size=buffer_size,
        dimnames=normalized_dimnames(info$row_names, info$col_names))
}


setClass("AnnDataMatrixH5",
    contains = "mat_uint32_t",
    slots = c(
        path = "character",
        group = "character",
        buffer_size = "integer"
    ),
    prototype = list(
        path = character(0),
        group = "matrix",
        buffer_size = integer(0)
    )
)
setMethod("iterate_matrix", "AnnDataMatrixH5", function(x) {
    iterate_matrix_anndata_hdf5_cpp(x@path, x@group, x@buffer_size)
})
setMethod("short_description", "AnnDataMatrixH5", function(x) {
    sprintf("AnnData HDF5 matrix in file %s, group %s", 
        x@path, x@group
    )
})

#' Read a sparse integer matrix from an anndata matrix in an hdf5 file
#' @inheritParams open_matrix_hdf5
#' @return 10xMatrixH5 object
#' @export
open_matrix_anndata_hdf5 <- function(path, group="X", buffer_size=16384L) {
    assert_is_file(path)
    assert_is(buffer_size, "integer")

    info <- dims_matrix_anndata_hdf5_cpp(path.expand(path), group, buffer_size)
    res <- new("AnnDataMatrixH5", path=path.expand(path), dim=info$dims, buffer_size=buffer_size,
        group=group, dimnames=normalized_dimnames(info$row_names, info$col_names))
    
    if (info[["transpose"]]) {
        return(t(res))
    } else {
        return(res)
    }
}



# Overlap matrix from fragments
setClass("PeakMatrix",
    contains = "mat_uint32_t",
    slots = c(
        fragments = "IterableFragments",
        chr_id = "integer",
        start = "integer",
        end = "integer",
        chr_levels = "character"
    ),
    prototype = list(
        fragments = NULL,
        chr_id = integer(0),
        start = integer(0),
        end = integer(0),
        chr_levels = character(0)
    )
)
#' Calculate cell x ranges overlap matrix
#' @param fragments Input fragments object. Must have cell names and chromosome names defined
#' @param ranges GRanges object with the ranges to overlap, or list/data frame with columns chr, start, & end.
#'     Must be sorted in order (chr, end, start) where chromosomes are ordered according to the chromosome names of `fragments`.
#' @param zero_based_coords Boolean for whether the input ranges are in a 0-based 
#'        or a 1-based coordinate system. (1-based will be converted to 0-based)
#'        (see http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/)
#' @note When calculating the matrix directly from a fragments tsv, it's necessary to first call `select_chromosomes` in order to 
#'     provide the ordering of chromosomes to expect while reading the tsv.
#' @return Iterable matrix object with dimension cell x ranges
#' @export
peakMatrix <- function(fragments, ranges, zero_based_coords=TRUE) {
    assert_is(fragments, "IterableFragments")
    assert_is(ranges, c("GRanges", "list", "data.frame"))
    
    assert_is(zero_based_coords, "logical")
    
    assert_not_na(cellNames(fragments))
    assert_not_na(chrNames(fragments))
    
    if(is.list(ranges)) {
        assert_has_names(ranges, c("chr", "start", "end"))
        chr_id <- as.integer(factor(as.character(ranges$chr), chrNames(fragments))) - 1L
        start <- as.integer(ranges$start) - !zero_based_coords
        end <- as.integer(ranges$end)
    } else {
        chr_id <- as.integer(factor(as.character(GenomicRanges::seqnames(ranges)), chrNames(fragments))) - 1L
        start <- GenomicRanges::start(ranges) - !zero_based_coords
        end <- GenomicRanges::end(ranges)
        
    }
    if(!all(order(chr_id, end, start) == seq_along(chr_id))) {
        stop("Peaks must be given in order sorted by (chr, end, start)")
    }
    res <- new("PeakMatrix", fragments=fragments, chr_id=chr_id, start=start, end=end)
    res@chr_levels <-chrNames(fragments)
    res@dim <- c(length(fragments@cell_names), length(res@chr_id))
    res@dimnames[[1]] <- fragments@cell_names
    return(res)
}

setMethod("iterate_matrix", "PeakMatrix", function(x) {iterate_peak_matrix_cpp(
    iterate_fragments(x@fragments), x@chr_id, x@start, x@end, x@chr_levels)
})


# Overlap matrix from fragments
setClass("TileMatrix",
    contains = "mat_uint32_t",
    slots = c(
        fragments = "IterableFragments",
        chr_id = "integer",
        start = "integer",
        end = "integer",
        tile_width = "integer", 
        chr_levels = "character"
    ),
    prototype = list(
        fragments = NULL,
        chr_id = integer(0),
        start = integer(0),
        end = integer(0),
        tile_width = integer(0),
        chr_levels = character(0)
    )
)
#' Calculate cell x ranges tile overlap matrix
#' @param fragments Input fragments object
#' @param ranges GRanges object with the ranges to overlap including a metadata column tile_width, 
#'  or a list/data frame with columns chr, start, end, and tile_width. Must be non-overlapping and sorted by
#'  (chr, start), with chromosomes ordered according to the chromosome names of `fragments`
#' @param zero_based_coords Boolean for whether the input ranges are in a 0-based 
#'        or a 1-based coordinate system. (1-based will be converted to 0-based)
#'        (see http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/)
#' @note When calculating the matrix directly from a fragments tsv, it's necessary to first call `select_chromosomes` in order to 
#'     provide the ordering of chromosomes to expect while reading the tsv.
#' @return Iterable matrix object with dimension cell x ranges
#' @export
tileMatrix <- function(fragments, ranges, zero_based_coords=TRUE) {
    assert_is(fragments, "IterableFragments")
    assert_is(ranges, c("GRanges", "list", "data.frame"))
    
    assert_is(zero_based_coords, "logical")
    
    assert_not_na(cellNames(fragments))
    assert_not_na(chrNames(fragments))
    
    if(is.list(ranges)) {
        assert_has_names(ranges, c("chr", "start", "end", "tile_width"))
        chr_id <- as.integer(factor(as.character(ranges$chr), chrNames(fragments))) - 1L
        start <- as.integer(ranges$start) - !zero_based_coords
        end <- as.integer(ranges$end)
        tile_width <- as.integer(ranges$tile_width)
    } else {
        assert_has_names(mcols(ranges), c("tile_width"))
        chr_id <- as.integer(factor(as.character(GenomicRanges::seqnames(ranges)), chrNames(fragments))) - 1L
        start <- GenomicRanges::start(ranges) - !zero_based_coords
        end <- GenomicRanges::end(ranges)
        tile_width <- mcols(ranges)$tile_width
    }   
    if(!all(order(chr_id, start, end) == seq_along(chr_id))) {
        stop("Tile regions must be given in order sorted by (chr, start, end)")
    }
    # Check to make sure tiles are non-overlapping
    pair_1 <- seq_len(length(chr_id)-1)
    pair_2 <- pair_1 + 1
    if (any(chr_id[pair_1] == chr_id[pair_2] & end[pair_1] >= start[pair_2])) {
        stop("Tile regions must be non-overlapping")
    }

    # Construct tile matrix
    res <- new("TileMatrix", fragments=fragments, chr_id=chr_id, start=start, end=end, tile_width=tile_width)
    res@chr_levels <- chrNames(fragments)
    
    tiles <- as.integer(ceiling((end - start) / tile_width))
    res@dim <- c(length(fragments@cell_names), sum(tiles))
    res@dimnames[[1]] <- fragments@cell_names
    return(res)
}

setMethod("iterate_matrix", "TileMatrix", function(x) {iterate_tile_matrix_cpp(
    iterate_fragments(x@fragments), x@chr_id, x@start, x@end, x@tile_width, x@chr_levels)
})

setMethod("short_description", "TileMatrix", function(x) {
    # Subset strings first to avoid a very slow string concatenation process
    indices <- c(head(seq_along(x@chr_id), 3), tail(seq_along(x@chr_id), 1))
    labels <- paste0(x@chr_levels[1+x@chr_id[indices]], ":", x@start[indices]+1, "-", x@end[indices], " (", x@tile_width[indices], "bp)")
    c(
        short_description(x@fragments),
        sprintf("Calculate %d tiles over %d ranges%s", ncol(x), length(x@chr_id), 
            pretty_print_vector(labels, prefix=": ", max_len=2)
        )
    )
})

# Conversions between integer and double IterableMatrix
setClass("convertMatrixIntDouble",
    contains = "mat_double",
    slots = c(
        matrix = "IterableMatrix"
    ),
    prototype = list(
        matrix = NULL
    )
)
setMethod("iterate_matrix", "convertMatrixIntDouble", function(x) {
    stopifnot(matrix_type(x@matrix) == "mat_uint32_t")
    convert_matrix_uint32_t_double_cpp(iterate_matrix(x@matrix))
})
setMethod("short_description", "convertMatrixIntDouble", function(x) {
    c(
        short_description(x@matrix),
        "Convert integers to doubles"
    )
})
cast_matrix_double <- function(from) {
    if (matrix_type(from) == "mat_double") 
        from
    else if (matrix_type(from) == "mat_uint32_t")
        wrapMatrix("convertMatrixIntDouble", from)
    else
        stop(sprintf("Unrecognized matrix type: %s", matrix_type(from)))
}

setClass("convertMatrixDoubleInt",
    contains = "mat_uint32_t",
    slots = c(
        matrix = "IterableMatrix"
    ),
    prototype = list(
        matrix = NULL
    )
)
setMethod("iterate_matrix", "convertMatrixDoubleInt", function(x) {
    stopifnot(matrix_type(x@matrix) == "mat_double")
    convert_matrix_double_uint32_t_cpp(iterate_matrix(x@matrix))
})
setMethod("short_description", "convertMatrixDoubleInt", function(x) {
    c(
        short_description(x@matrix),
        "Convert doubles to integers"
    )
})
cast_matrix_int <- function(from) {
    if (matrix_type(from) == "mat_double") 
        wrapMatrix("convertMatrixDoubleInt", from)
    else if (matrix_type(from) == "mat_uint32_t")
        from
    else
        stop(sprintf("Unrecognized matrix type: %s", matrix_type(from)))
}

# Conversions with dgCMatrix
setClass("Iterable_dgCMatrix_wrapper",
    contains = "mat_double",
    slots = c(
        mat = "dgCMatrix"
    ),
    prototype = list(
        mat = NULL
    )
)
setAs("dgCMatrix", "IterableMatrix", function(from) {
    new("Iterable_dgCMatrix_wrapper", dim=dim(from), dimnames=dimnames(from), transpose=FALSE, mat=from)
})
setMethod("iterate_matrix", "Iterable_dgCMatrix_wrapper", function(x) {
    iterate_csparse_matrix_cpp(x@mat)
})
setMethod("short_description", "Iterable_dgCMatrix_wrapper", function(x) {
    "Load dgCMatrix from memory"
})

setMethod("iterate_matrix", "dgCMatrix", function(x) {
    iterate_csparse_matrix_cpp(x)
})

setAs("IterableMatrix", "dgCMatrix", function(from) {
    iter <- iterate_matrix(from)
    if (matrix_type(from) == "mat_uint32_t") iter <- convert_matrix_uint32_t_double_cpp(iter)
    res <- build_csparse_matrix_double_cpp(iter)
    if (from@transpose) {
        res <- t(res)
    }
    res@Dimnames <- from@dimnames
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

    colnames(res$row_stats) <- rownames(matrix)
    colnames(res$col_stats) <- colnames(matrix)
    
    if (matrix@transpose) {
        res <- list(
            row_stats = res$col_stats,
            col_stats = res$row_stats
        )
    }
    return(res)
}