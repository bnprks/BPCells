
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
setMethod("dimnames<-", signature(x="IterableMatrix", value="list"), function(x, value) {
    d <- dim(x)
    has_error <- FALSE
    if (!is.list(value) || length(value != 2)) has_error <- TRUE
    if (!is.null(value[[1]]) && length(value[[1]]) != d[1]) has_error <- TRUE
    if (!is.null(value[[2]]) && length(value[[2]]) != d[2]) has_error <- TRUE
    if (has_error) stop("Invalid dimnames supplied")

    if (!is.null(value[[1]])) x@dimnames[[1]] <- as.character(value[[1]])
    else x@dimnames[[1]] <- NULL
    if (!is.null(value[[2]])) x@dimnames[[2]] <- as.character(value[[2]])
    else x@dimnames[[2]] <- NULL
    x
})
setMethod("dimnames<-", signature(x="IterableMatrix", value="NULL"), function(x, value) {
    x@dimnames <- list(NULL, NULL)
    x
})

#' Construct an S4 matrix object wrapping another matrix object, avoiding
#' duplicate storage of dimnames
wrapMatrix <- function(class, m, ...) {
    dimnames <- dimnames(m)
    if (matrix_is_transform(m)) m@dimnames <- list(NULL, NULL)
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

#' Get a wrapped pointer to the iterable matrix, in the form of an XPtrList (see fragments.R)
setGeneric("iterate_matrix", function(x) standardGeneric("iterate_matrix"))
setMethod("iterate_matrix", "XPtrList", function(x) {
    stopifnot(x@type == "mat_uint32_t" || x@type == "mat_double")
    x
})

#' Get the matrix data type (mat_uint32_t or mat_double for now)
setGeneric("matrix_type", function(x) standardGeneric("matrix_type"))
setMethod("matrix_type", "XPtrList", function(x) x@type)

#' Return boolean whether the matrix is a transform (i.e. is loading data from another
#' matrix R object rather than a data source it owns). This is used primarily to know when it is safe
#' to clear dimnames from intermediate transformed matrices. C++ relies on the base matrices (non-transform) to have
#' dimnames, while R relies on the outermost matrix (transform) to have dimnames.
setGeneric("matrix_is_transform", function(x) standardGeneric("matrix_is_transform"))
setMethod("matrix_is_transform", "XPtrList", function(x) FALSE)
# As a reasonable default convention, a matrix is a transform if it has a slot named matrix
setMethod("matrix_is_transform", "IterableMatrix", function(x) {.hasSlot(x, "matrix")})

setMethod("short_description", "IterableMatrix", function(x) {
    character(0)
})


#' @describeIn IterableMatrix-methods Display an IterableMatrix
#' @param object IterableMatrix object
setMethod("show", "IterableMatrix", function(object) {
    cat(sprintf("%d x %d IterableMatrix object with class %s\n", nrow(object), ncol(object), class(object)))

    cat("\n")

    
    cat(sprintf("Row names: %s\n",
        pretty_print_vector(rownames(object), empty="unknown names")
    ))

    
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
        return(Matrix::t(dense_multiply_left_cpp(ptr(iter), Matrix::t(y))))
    } else {
        return(dense_multiply_right_cpp(ptr(iter), y))
    }
})

setMethod("%*%", signature(x="matrix", y="IterableMatrix"), function(x, y) {
    iter <- iterate_matrix(cast_matrix_double(y))
    if(y@transpose) {
        return(Matrix::t(dense_multiply_right_cpp(ptr(iter), Matrix::t(x))))
    } else {
        return(dense_multiply_left_cpp(ptr(iter), x))
    }
})

setMethod("%*%", signature(x="IterableMatrix", y="numeric"), function(x, y) {
    iter <- iterate_matrix(cast_matrix_double(x))
    if(x@transpose) {
        return(vec_multiply_left_cpp(ptr(iter), y))
    } else {
        return(vec_multiply_right_cpp(ptr(iter), y))
    }
})

setMethod("%*%", signature(x="numeric", y="IterableMatrix"), function(x, y) {
    iter <- iterate_matrix(cast_matrix_double(y))
    if(y@transpose) {
        return(vec_multiply_right_cpp(ptr(iter), x))
    } else {
        return(vec_multiply_left_cpp(ptr(iter), x))
    }
})

# Row sums and row means

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate rowSums
#' @return * `rowSums()`: vector of row sums
setMethod("rowSums", signature(x="IterableMatrix"), function(x) {
    iter <- iterate_matrix(cast_matrix_double(x))
    if (x@transpose) 
        res <- col_sums_double_cpp(ptr(iter))
    else
        res <- row_sums_double_cpp(ptr(iter))
    names(res) <- rownames(x)
    res
})

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate colSums
#' @return * `colSums()`: vector of col sums
setMethod("colSums", signature(x="IterableMatrix"), function(x) {
    iter <- iterate_matrix(cast_matrix_double(x))
    if (x@transpose) 
        res <- row_sums_double_cpp(ptr(iter))
    else
        res <- col_sums_double_cpp(ptr(iter))
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

    # Handle chained subsets by collapsing repeated selections
    if (is(x, "MatrixSubset")) {
        if (length(x@row_selection) == 0) x@row_selection <- ret@row_selection
        else if (length(ret@row_selection) != 0) x@row_selection <- x@row_selection[ret@row_selection]

        if (length(x@col_selection) == 0) x@col_selection <- ret@col_selection
        else if (length(ret@col_selection) != 0) x@col_selection <- x@col_selection[ret@col_selection]
        
        x@dim <- ret@dim
        x@dimnames <- ret@dimnames
        ret <- x
    }
    ret
})

setMethod("iterate_matrix", "MatrixSubset", function(x) {
    ret <- iterate_matrix(x@matrix)

    if (matrix_type(x) == "mat_double") {
        if (length(x@row_selection) != 0) ret <- wrapMatDouble(iterate_matrix_row_select_double_cpp(ptr(ret), x@row_selection-1L), ret)
        if (length(x@col_selection) != 0) ret <- wrapMatDouble(iterate_matrix_col_select_double_cpp(ptr(ret), x@col_selection-1L), ret)
    } else if (matrix_type(x) == "mat_uint32_t") {
        if (length(x@row_selection) != 0) ret <- wrapMatUInt(iterate_matrix_row_select_uint32_t_cpp(ptr(ret), x@row_selection-1L), ret)
        if (length(x@col_selection) != 0) ret <- wrapMatUInt(iterate_matrix_col_select_uint32_t_cpp(ptr(ret), x@col_selection-1L), ret)
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


# Concatenating matrices by row or by column

#' Helper function for rbind/cbind concatenating dimnames
concat_dimnames <- function(x, y, len_x, len_y, warning_prefix, dim_type) {
    if (is.null(x) && is.null(y)) return(NULL)
    if (!is.null(x) && !is.null(y)) return(c(x, y))
    warning(sprintf("%s: %s names presenent on some but not all matrices. Setting missing names to \"\"",
        warning_prefix, dim_type
    ), call.=FALSE)
    if (is.null(x)) x <- rep_len("", len_x)
    if (is.null(y)) y <- rep_len("", len_y)
    c(x, y)
}

#' Helper function for rbind/cbind merging dimnames
merge_dimnames <- function(x, y, warning_prefix, dim_type) {
    if (!is.null(x) && !is.null(y) && !all(x == y)) {
        warning(sprintf("%s: %s names are mismatched. Setting names to match first matrix", warning_prefix, dim_type), call.=FALSE)
    }
    if (!is.null(x)) return(x)
    if (!is.null(y)) return(y)
    return(NULL)
}

setClass("RowBindMatrices",
    contains = "IterableMatrix",
    slots = c(
        matrix_list = "list"
    ),
    prototype = list(
        matrix_list = list()
    )
)
setMethod("matrix_type", signature(x="RowBindMatrices"), function(x) matrix_type(x@matrix_list[[1]]))

setMethod("iterate_matrix", "RowBindMatrices", function(x) {
    iterators <- lapply(x@matrix_list, iterate_matrix)
    if (matrix_type(x) == "mat_uint32_t")
        wrapMatUInt(iterate_matrix_row_bind_uint32_t_cpp(lapply(iterators, ptr)), new("XPtrList", pointers=iterators))
    else if (matrix_type(x) == "mat_double")
        wrapMatDouble(iterate_matrix_row_bind_double_cpp(lapply(iterators, ptr)), new("XPtrList", pointers=iterators))
})

setMethod("matrix_is_transform", "RowBindMatrices", function(x) TRUE)

setMethod("short_description", "RowBindMatrices", function(x) {
    sprintf("Concatenate rows of %d matrix objects with classes%s",
        length(x@matrix_list),
        pretty_print_vector(vapply(x@matrix_list, class, character(1)), prefix=": ", max_len=3)
    )
})

setMethod("rbind2", signature(x="IterableMatrix", y="IterableMatrix"), function(x, y, ...) {
    if (x@transpose != y@transpose) stop("Cannot merge matrices with different interal transpose states.\nPlease convert one matrix to dgCMatrix, transpose, then re-convert to IterableMatrix.")
    if (x@transpose) return(Matrix::t(cbind2(Matrix::t(x), Matrix::t(y))))

    if (ncol(x) != ncol(y)) stop("Error in rbind: matrices must have equal number of columns")
    # Handle dimnames
    col_names <- merge_dimnames(colnames(x), colnames(y), "rbind", "column")
    row_names <- concat_dimnames(rownames(x), rownames(y), nrow(x), nrow(y), "rbind", "row")
    if (matrix_is_transform(x)) x@dimnames <- list(NULL, NULL)
    if (matrix_is_transform(y)) y@dimnames <- list(NULL, NULL)

    matrix_list <- list()
    if (is(x, "RowBindMatrices")) matrix_list <- c(matrix_list, x@matrix_list)
    else matrix_list <- c(matrix_list, x)
    if (is(y, "RowBindMatrices")) matrix_list <- c(matrix_list, y@matrix_list)
    else matrix_list <- c(matrix_list, y)

    new("RowBindMatrices", matrix_list=matrix_list, dim=c(nrow(x)+nrow(y), ncol(x)), dimnames=list(row_names, col_names), transpose=FALSE)
})

setClass("ColBindMatrices",
    contains = "IterableMatrix",
    slots = c(
        matrix_list = "list"
    ),
    prototype = list(
        matrix_list = list()
    )
)
setMethod("matrix_type", signature(x="ColBindMatrices"), function(x) matrix_type(x@matrix_list[[1]]))

setMethod("iterate_matrix", "ColBindMatrices", function(x) {
    iterators <- lapply(x@matrix_list, iterate_matrix)
    if (matrix_type(x) == "mat_uint32_t")
        wrapMatUInt(iterate_matrix_col_bind_uint32_t_cpp(lapply(iterators, ptr)), new("XPtrList", pointers=iterators))
    else if (matrix_type(x) == "mat_double")
        wrapMatDouble(iterate_matrix_col_bind_double_cpp(lapply(iterators, ptr)), new("XPtrList", pointers=iterators))
})

setMethod("matrix_is_transform", "ColBindMatrices", function(x) TRUE)

setMethod("short_description", "ColBindMatrices", function(x) {
    sprintf("Concatenate columns of %d matrix objects with classes%s",
        length(x@matrix_list),
        pretty_print_vector(vapply(x@matrix_list, class, character(1)), prefix=": ", max_len=3)
    )
})

setMethod("cbind2", signature(x="IterableMatrix", y="IterableMatrix"), function(x, y, ...) {
    if (x@transpose != y@transpose) stop("Cannot merge matrices with different interal transpose states.\nPlease convert one matrix to dgCMatrix, transpose, then re-convert to IterableMatrix.")
    if (x@transpose) return(Matrix::t(rbind2(Matrix::t(x), Matrix::t(y))))

    if (nrow(x) != nrow(y)) stop("Error in cbind: matrices must have equal number of columns")
    # Handle dimnames
    row_names <- merge_dimnames(rownames(x), rownames(y), "cbind", "row")
    col_names <- concat_dimnames(colnames(x), colnames(y), ncol(x), ncol(y), "cbind", "column")
    if (matrix_is_transform(x)) x@dimnames <- list(NULL, NULL)
    if (matrix_is_transform(y)) y@dimnames <- list(NULL, NULL)

    matrix_list <- list()
    if (is(x, "ColBindMatrices")) matrix_list <- c(matrix_list, x@matrix_list)
    else matrix_list <- c(matrix_list, x)
    if (is(y, "ColBindMatrices")) matrix_list <- c(matrix_list, y@matrix_list)
    else matrix_list <- c(matrix_list, y)

    new("ColBindMatrices", matrix_list=matrix_list, dim=c(nrow(x), ncol(x)+ncol(y)), dimnames=list(row_names, col_names), transpose=FALSE)
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
    wrapMatUInt(iterate_packed_matrix_cpp(x, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
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
    wrapMatUInt(iterate_unpacked_matrix_cpp(x, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
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

    it <- iterate_matrix(mat)
    if (compress) {
        m <- write_packed_matrix_cpp(ptr(it))
        class <- "PackedMatrixMem" 
    } else {
        m <- write_unpacked_matrix_cpp(ptr(it))
        class <- "UnpackedMatrixMem"
    }
    m[["dimnames"]] <- normalized_dimnames(m$row_names, m$col_names)
    m$row_names <- NULL
    m$col_names <- NULL
    res <- do.call(new, c(class, m))
    res@dim <- mat@dim
    res@dimnames <- m$dimnames
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
        wrapMatUInt(iterate_packed_matrix_file_cpp(x@dir, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
    else
        wrapMatUInt(iterate_unpacked_matrix_file_cpp(x@dir, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
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

    it <- iterate_matrix(mat)
    if (compress)
        write_packed_matrix_file_cpp(ptr(it), dir, buffer_size)
    else 
        write_unpacked_matrix_file_cpp(ptr(it), dir, buffer_size)
    
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
        wrapMatUInt(iterate_packed_matrix_hdf5_cpp(x@path, x@group, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
    else
        wrapMatUInt(iterate_unpacked_matrix_hdf5_cpp(x@path, x@group, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
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

    it <- iterate_matrix(mat)
    if (compress)
        write_packed_matrix_hdf5_cpp(ptr(it), path, group, buffer_size, chunk_size)
    else 
        write_unpacked_matrix_hdf5_cpp(ptr(it), path, group, buffer_size, chunk_size)
    
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
    wrapMatUInt(iterate_matrix_10x_hdf5_cpp(x@path, x@buffer_size), new("XPtrList"))
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
    contains = "mat_double",
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
    wrapMatDouble(iterate_matrix_anndata_hdf5_cpp(x@path, x@group, x@buffer_size), new("XPtrList"))
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
        return(Matrix::t(res))
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
    res@dim <- c(length(cellNames(fragments)), length(res@chr_id))
    res@dimnames[[1]] <- cellNames(fragments)
    return(res)
}

setMethod("iterate_matrix", "PeakMatrix", function(x) {
    it <- iterate_fragments(x@fragments)
    wrapMatUInt(iterate_peak_matrix_cpp(ptr(it), x@chr_id, x@start, x@end, x@chr_levels), it)
})
setMethod("short_description", "PeakMatrix", function(x) {
    # Subset strings first to avoid a very slow string concatenation process
    indices <- c(head(seq_along(x@chr_id), 3), tail(seq_along(x@chr_id), 1))
    labels <- paste0(x@chr_levels[1+x@chr_id[indices]], ":", x@start[indices]+1, "-", x@end[indices])
    c(
        short_description(x@fragments),
        sprintf("Calculate %d peaks over %d ranges%s", ncol(x), length(x@chr_id), 
            pretty_print_vector(labels, prefix=": ", max_len=2)
        )
    )
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
    res@dim <- c(length(cellNames(fragments)), sum(tiles))
    res@dimnames[[1]] <- cellNames(fragments)
    return(res)
}

setMethod("iterate_matrix", "TileMatrix", function(x) {
    it <- iterate_fragments(x@fragments)
    wrapMatUInt(iterate_tile_matrix_cpp(ptr(it), x@chr_id, x@start, x@end, x@tile_width, x@chr_levels), it)
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
    it <- iterate_matrix(x@matrix)
    wrapMatDouble(convert_matrix_uint32_t_double_cpp(ptr(it)), it)
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
    it <- iterate_matrix(x@matrix)
    wrapMatUInt(convert_matrix_double_uint32_t_cpp(ptr(it)), it)
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
    wrapMatDouble(iterate_csparse_matrix_cpp(x@mat), new("XPtrList"))
})
setMethod("short_description", "Iterable_dgCMatrix_wrapper", function(x) {
    "Load dgCMatrix from memory"
})

setMethod("iterate_matrix", "dgCMatrix", function(x) {
    wrapMatDouble(iterate_csparse_matrix_cpp(x), new("XPtrList"))
})

setAs("IterableMatrix", "dgCMatrix", function(from) {
    iter <- iterate_matrix(from)
    if (matrix_type(from) == "mat_uint32_t") iter <- wrapMatDouble(convert_matrix_uint32_t_double_cpp(ptr(iter)), iter)
    res <- build_csparse_matrix_double_cpp(ptr(iter))
    if (from@transpose) {
        res <- Matrix::t(res)
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
    
    it <- iterate_matrix(matrix)
    res <- matrix_stats_cpp(ptr(it), row_stats_number, col_stats_number)
    rownames(res$row_stats) <- stat_options[seq_len(row_stats_number) + 1]
    rownames(res$col_stats) <- stat_options[seq_len(col_stats_number) + 1]
    if (matrix@transpose) {
        tmp <- res$row_stats
        res$row_stats <- res$col_stats
        res$col_stats <- tmp
    }

    colnames(res$row_stats) <- rownames(matrix)
    colnames(res$col_stats) <- colnames(matrix)
  
    return(res)
}