
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
    stopifnot(x@type == "mat_uint32_t" || x@type == "mat_double" || x@type == "mat_float")
    x
})

#' Get the matrix data type (mat_uint32_t, mat_float, or mat_double for now)
setGeneric("matrix_type", function(x) standardGeneric("matrix_type"))
setMethod("matrix_type", "XPtrList", function(x) {
    # Strip off the "mat_" prefix
    assert_true(substr(x@type, 1, 4) == "mat_")
    substr(x@type, 5L, 10000L)
})


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
    cat(sprintf("Data type: %s\n", matrix_type(object)))
    cat(sprintf("Storage order: %s major\n", ifelse(object@transpose, "row", "column")))
    
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
#' @export
setMethod("t", signature(x = "IterableMatrix"), function(x) {
    x@transpose <- !x@transpose
    x@dim <- c(x@dim[2L], x@dim[1L])
    x@dimnames <- list(x@dimnames[[2]], x@dimnames[[1]])
    return(x)
})

# Dense multiply operators (sparse*dense_mat and sparse*dense_vec)

#' @param x IterableMatrix object
#' @param y matrix
#' @describeIn IterableMatrix-methods Multiply by a dense matrix
#' @return * `x %*% y`: dense matrix result
setMethod("%*%", signature(x="IterableMatrix", y="matrix"), function(x, y) {
    iter <- iterate_matrix(convert_matrix_type(x, "double"))
    if(x@transpose) {
        return(Matrix::t(dense_multiply_left_cpp(ptr(iter), Matrix::t(y))))
    } else {
        return(dense_multiply_right_cpp(ptr(iter), y))
    }
})

setMethod("%*%", signature(x="matrix", y="IterableMatrix"), function(x, y) {
    iter <- iterate_matrix(convert_matrix_type(y, "double"))
    if(y@transpose) {
        return(Matrix::t(dense_multiply_right_cpp(ptr(iter), Matrix::t(x))))
    } else {
        return(dense_multiply_left_cpp(ptr(iter), x))
    }
})

setMethod("%*%", signature(x="IterableMatrix", y="numeric"), function(x, y) {
    iter <- iterate_matrix(convert_matrix_type(x, "double"))
    if(x@transpose) {
        return(vec_multiply_left_cpp(ptr(iter), y))
    } else {
        return(vec_multiply_right_cpp(ptr(iter), y))
    }
})

setMethod("%*%", signature(x="numeric", y="IterableMatrix"), function(x, y) {
    iter <- iterate_matrix(convert_matrix_type(y, "double"))
    if(y@transpose) {
        return(vec_multiply_right_cpp(ptr(iter), x))
    } else {
        return(vec_multiply_left_cpp(ptr(iter), x))
    }
})

#' LinearOperator class for performing sparse matrix - dense vector products, e.g.
#' for downstream matrix solvers. This is to avoid having to call iterate_matrix
#' repeatedly for a possible efficiency gain
setClass("LinearOperator",
    slots = c(
        dim = "integer",
        xptrlist = "XPtrList",
        transpose = "logical"
    ),
    prototype = list(
        dim = integer(2),
        transpose = FALSE
    )
)

#' Construct a C++ matrix object and save the pointer to use for repeated matrix-vector products
#' A bit experimental still so for internal use
linear_operator <- function(mat) {
    assert_is(mat, "IterableMatrix")
    new("LinearOperator", dim=dim(mat), xptrlist=iterate_matrix(convert_matrix_type(mat, "double")), transpose=mat@transpose)
}

setMethod("%*%", signature(x="LinearOperator", y="matrix"), function(x, y) {
    if(x@transpose) {
        return(Matrix::t(dense_multiply_left_cpp(ptr(x@xptrlist), Matrix::t(y))))
    } else {
        return(dense_multiply_right_cpp(ptr(x@xptrlist), y))
    }
})

setMethod("%*%", signature(x="matrix", y="LinearOperator"), function(x, y) {
    if(y@transpose) {
        return(Matrix::t(dense_multiply_right_cpp(ptr(y@xptrlist), Matrix::t(x))))
    } else {
        return(dense_multiply_left_cpp(ptr(y@xptrlist), x))
    }
})

setMethod("%*%", signature(x="LinearOperator", y="numeric"), function(x, y) {
    if(x@transpose) {
        return(vec_multiply_left_cpp(ptr(x@xptrlist), y))
    } else {
        return(vec_multiply_right_cpp(ptr(x@xptrlist), y))
    }
})

setMethod("%*%", signature(x="numeric", y="LinearOperator"), function(x, y) {
    if(y@transpose) {
        return(vec_multiply_right_cpp(ptr(y@xptrlist), x))
    } else {
        return(vec_multiply_left_cpp(ptr(y@xptrlist), x))
    }
})


# Sparse matrix multiply
setClass("MatrixMultiply",
    contains = "IterableMatrix",
    slots = c(
        left = "IterableMatrix",
        right = "IterableMatrix"
    ),
    prototype = list(
        left = NULL,
        right = NULL
    )
)
setMethod("matrix_type", signature(x="MatrixMultiply"), function(x) matrix_type(x@left))
setMethod("matrix_is_transform", "MatrixMultiply", function(x) TRUE)

setMethod("iterate_matrix", "MatrixMultiply", function(x) {
    left <- iterate_matrix(x@left)
    right <- iterate_matrix(x@right)
    inner_iterator <- new("XPtrList", pointers=c(left@pointers, right@pointers))
    if (matrix_type(x) == "uint32_t")
        wrapMat_uint32_t(iterate_matrix_multiply_uint32_t_cpp(ptr(left), ptr(right)), inner_iterator)
    else if (matrix_type(x) == "double")
        wrapMat_double(iterate_matrix_multiply_double_cpp(ptr(left), ptr(right)), inner_iterator)
})

setMethod("short_description", "MatrixMultiply", function(x) {
    sprintf("Multiply sparse matrices: %s (%dx%d) * %s (%dx%d)",
        class(x@left), nrow(x@left), ncol(x@left),
        class(x@right), nrow(x@right), ncol(x@right)
    )
})

setMethod("%*%", signature(x="IterableMatrix", y="IterableMatrix"), function(x, y) {
    if (x@transpose != y@transpose) stop("Cannot multiply matrices with different interal transpose states.\nPlease convert one matrix to dgCMatrix, transpose, then re-convert to IterableMatrix.")
    if (x@transpose) return(Matrix::t(Matrix::t(y) %*% Matrix::t(x)))

    assert_true(ncol(x) == nrow(y))

    # If types are mismatched, default to double precision for both
    type_x <- matrix_type(x)
    type_y <- matrix_type(y)
    if (type_x != type_y && type_x != "double") x <- convert_matrix_type(x, "double")
    if (type_x != type_y && type_y != "double") y <- convert_matrix_type(y, "double")

    dim <- c(nrow(x), ncol(y))
    dimnames <- list(rownames(x), colnames(y))
    new("MatrixMultiply", left=x, right=y, transpose=FALSE, dim=dim, dimnames=dimnames)
})

setMethod("%*%", signature(x="IterableMatrix", y="dgCMatrix"), function(x, y) {
    if (x@transpose)
        t(as(t(y), "IterableMatrix") %*% t(x))
    else
        x %*% as(y, "IterableMatrix")
})

setMethod("%*%", signature(x="dgCMatrix", y="IterableMatrix"), function(x, y) {
    if (y@transpose)
        t(t(y) %*% as(t(x), "IterableMatrix"))
    else
        as(x, "IterableMatrix") %*% y
})

# Row sums and row means

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate rowSums
#' @return * `rowSums()`: vector of row sums
setMethod("rowSums", signature(x="IterableMatrix"), function(x) {
    iter <- iterate_matrix(convert_matrix_type(x, "double"))
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
    iter <- iterate_matrix(convert_matrix_type(x, "double"))
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

    if (matrix_type(x) == "double") {
        if (length(x@row_selection) != 0) ret <- wrapMat_double(iterate_matrix_row_select_double_cpp(ptr(ret), x@row_selection-1L), ret)
        if (length(x@col_selection) != 0) ret <- wrapMat_double(iterate_matrix_col_select_double_cpp(ptr(ret), x@col_selection-1L), ret)
    } else if (matrix_type(x) == "uint32_t") {
        if (length(x@row_selection) != 0) ret <- wrapMat_uint32_t(iterate_matrix_row_select_uint32_t_cpp(ptr(ret), x@row_selection-1L), ret)
        if (length(x@col_selection) != 0) ret <- wrapMat_uint32_t(iterate_matrix_col_select_uint32_t_cpp(ptr(ret), x@col_selection-1L), ret)
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
    inner_iterator <- new("XPtrList", pointers=do.call(c, lapply(iterators, function(i) i@pointers)))
    if (matrix_type(x) == "uint32_t")
        wrapMat_uint32_t(iterate_matrix_row_bind_uint32_t_cpp(lapply(iterators, ptr)), inner_iterator)
    else if (matrix_type(x) == "double")
        wrapMat_double(iterate_matrix_row_bind_double_cpp(lapply(iterators, ptr)), inner_iterator)
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
    inner_iterator <- new("XPtrList", pointers=do.call(c, lapply(iterators, function(i) i@pointers)))
    if (matrix_type(x) == "uint32_t")
        wrapMat_uint32_t(iterate_matrix_col_bind_uint32_t_cpp(lapply(iterators, ptr)), inner_iterator)
    else if (matrix_type(x) == "double")
        wrapMat_double(iterate_matrix_col_bind_double_cpp(lapply(iterators, ptr)), inner_iterator)
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
setClass("PackedMatrixMemBase",
    contains = "IterableMatrix",
    slots = c(
        # Leave out val storage since it's datatype-dependent
        row_data = "integer",   
        row_starts = "integer", 
        row_idx = "integer",
        col_ptr = "integer",   
        row_count = "integer",

        version = "character"
    ),
    prototype = list(
        # Leave out val storage since it's datatype-dependent
        row_data = integer(0),   
        row_starts = integer(0), 
        row_idx = integer(0),
        col_ptr = integer(0),   
        row_count = integer(0),

        version = character(0)
    )
)
setMethod("short_description", "PackedMatrixMemBase", function(x) {
    "Load compressed matrix from memory"
})

setClass("PackedMatrixMem_uint32_t",
    contains = "PackedMatrixMemBase",
    slots = c(
        val_data = "integer",   
        val_idx = "integer"  
    ),
    prototype = list(
        val_data = integer(0),   
        val_idx = integer(0)   
    )
)
setMethod("matrix_type", "PackedMatrixMem_uint32_t", function(x) "uint32_t")
setMethod("iterate_matrix", "PackedMatrixMem_uint32_t", function(x) {
    x@dimnames <- denormalize_dimnames(x@dimnames)
    wrapMat_uint32_t(iterate_packed_matrix_mem_uint32_t_cpp(x, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
})

setClass("PackedMatrixMem_float",
    contains = "PackedMatrixMemBase",
    slots = c(val = "integer"),
    prototype = list(val = integer(0))
)
setMethod("matrix_type", "PackedMatrixMem_float", function(x) "float")
setMethod("iterate_matrix", "PackedMatrixMem_float", function(x) {
    x@dimnames <- denormalize_dimnames(x@dimnames)
    wrapMat_float(iterate_packed_matrix_mem_float_cpp(x, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
})

setClass("PackedMatrixMem_double",
    contains = "PackedMatrixMemBase",
    slots = c(val = "numeric"),
    prototype = list(val = numeric(0))
)
setMethod("matrix_type", "PackedMatrixMem_double", function(x) "double")
setMethod("iterate_matrix", "PackedMatrixMem_double", function(x) {
    x@dimnames <- denormalize_dimnames(x@dimnames)
    wrapMat_double(iterate_packed_matrix_mem_double_cpp(x, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
})

setClass("UnpackedMatrixMemBase",
    contains = "IterableMatrix",
    slots = c(
        # Leave out val storage since it's data-type dependent
        row = "integer",    
        col_ptr = "integer",   
        row_count = "integer",

        version="character"
    ),
    prototype = list(
        row = integer(0),   
        col_ptr = integer(0),   
        row_count = integer(0),

        version=character(0)
    )
)
setMethod("short_description", "UnpackedMatrixMemBase", function(x) {
    "Load uncompressed matrix from memory"
})

setClass("UnpackedMatrixMem_uint32_t",
    contains = "UnpackedMatrixMemBase",
    slots = c(val = "integer"),
    prototype = list(val = integer())
)
setMethod("matrix_type", "UnpackedMatrixMem_uint32_t", function(x) "uint32_t")
setMethod("iterate_matrix", "UnpackedMatrixMem_uint32_t", function(x) {
    x@dimnames <- denormalize_dimnames(x@dimnames)
    wrapMat_uint32_t(iterate_unpacked_matrix_mem_uint32_t_cpp(x, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
})

setClass("UnpackedMatrixMem_float",
    contains = "UnpackedMatrixMemBase",
    slots = c(val = "integer"),
    prototype = list(val = integer(0))
)
setMethod("matrix_type", "UnpackedMatrixMem_float", function(x) "float")
setMethod("iterate_matrix", "UnpackedMatrixMem_float", function(x) {
    x@dimnames <- denormalize_dimnames(x@dimnames)
    wrapMat_float(iterate_unpacked_matrix_mem_float_cpp(x, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
})

setClass("UnpackedMatrixMem_double",
    contains = "UnpackedMatrixMemBase",
    slots = c(val = "numeric"),
    prototype = list(val = numeric(0))
)
setMethod("matrix_type", "UnpackedMatrixMem_double", function(x) "double")
setMethod("iterate_matrix", "UnpackedMatrixMem_double", function(x) {
    x@dimnames <- denormalize_dimnames(x@dimnames)
    wrapMat_double(iterate_unpacked_matrix_mem_double_cpp(x, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
})


#' load_bytes = 4MiB, sort_bytes = 1GiB, for a uint32_t matrix implies ~85GB of data
#' can be sorted with two passes through the data, and ~7.3TB of data can be 
#' sorted in three passes through the data.
write_matrix_transpose <- function(matrix, tmpdir, load_bytes, sort_bytes) {

    assert_true(matrix_type(matrix) %in% c("uint32_t", "float", "double"))

    write_function <- get(sprintf("write_matrix_transpose_%s_cpp", matrix_type(matrix)))
    wrap_function <- get(sprintf("wrapMat_%s", matrix_type(matrix)))

    it <- iterate_matrix(matrix)
    res <- write_function(ptr(it), tmpdir, load_bytes, sort_bytes)

    wrap_function(res, new("XPtrList"))
}

#' Write a sparse matrix object to memory.
#' @param matrix Input matrix, either IterableMatrix or dgCMatrix
#' @param compress Whether or not to compress the data. 
#' @return IterableMatrix object
#' @details This function will convert row major storage orders to column major,
#' which will usually require 2 extra passes over the data. 
#' 
#' For typical RNA counts matrices, which are unsigned integer matrices 
#' (type uint32_t), this function will result in 6-8x less space than an R dgCMatrix,
#' and 4-6x smaller than a scipy csc_matrix. The compression will be more effective when
#' the count values in the matrix are small, and when the rows of the matrix are 
#' sorted by rowMeans. In tests on RNA-seq data optimal ordering could save up to
#' 40% of storage space. On non-integer data only the row indices are compressed,
#' not the values themselves so space savings will be smaller.
#' 
#' Note: this function will not perform automatic conversions to non-negative
#' integer matrices, so use the function `convert_matrix_type_uint32_t` first if high-efficiency
#' compression is desired. 
#' 
#' @export
write_matrix_memory <- function(mat, compress=TRUE, transpose_tmpdir=tempdir()) {
    assert_is(mat, c("IterableMatrix", "dgCMatrix"))
    if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")

    assert_true(matrix_type(mat) %in% c("uint32_t", "float", "double"))
    if (compress && matrix_type(mat) != "uint32_t") {
        rlang::inform(c(
            "Warning: Matrix compression performs poorly with non-integers.",
            "Consider calling convert_matrix_type if a compressed integer matrix is intended."
        ), .frequency="regularly", .frequency_id="matrix_compress_non_integer")
    }

    write_function <- get(sprintf("write_%s_matrix_mem_%s_cpp", ifelse(compress, "packed", "unpacked"), matrix_type(mat)))
    class <- sprintf("%sMatrixMem_%s", ifelse(compress, "Packed", "Unpacked"), matrix_type(mat))

    if (mat@transpose) {
        tmp_path <- tempfile(pattern="transpose")
        it <- write_matrix_transpose(t(mat), tmp_path, load_bytes=4194304L, sort_bytes=1073741824L)
    }
    else it <- iterate_matrix(mat)
    
    m <- write_function(ptr(it))
    
    if (mat@transpose) unlink(tmp_path, recursive=TRUE, expand=FALSE)

    m[["dimnames"]] <- normalized_dimnames(m$row_names, m$col_names)
    m$row_names <- NULL
    m$col_names <- NULL
    res <- do.call(new, c(class, m))
    res@dim <- mat@dim
    res@dimnames <- m$dimnames
    res
}

setClass("MatrixDir",
    contains = "IterableMatrix",
    slots = c(
        dir = "character",
        compressed = "logical",
        buffer_size = "integer",
        type = "character"
    ),
    prototype = list(
        dir = character(0),
        compressed = logical(0),
        buffer_size = integer(0),
        type = character(0)
    )
)
setMethod("matrix_type", "MatrixDir", function(x) x@type)

setMethod("iterate_matrix", "MatrixDir", function(x) {
    x@dimnames <- denormalize_dimnames(x@dimnames)

    iter_function <- get(sprintf("iterate_%s_matrix_file_%s_cpp", ifelse(x@compressed, "packed", "unpacked"), x@type))
    wrap_function <- get(sprintf("wrapMat_%s", x@type))

    wrap_function(iter_function(x@dir, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
})

setMethod("short_description", "MatrixDir", function(x) {
    sprintf("Load %s matrix from directory %s", 
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

    assert_is(dir, "character")
    assert_is(compress, "logical")
    assert_is(buffer_size, "integer")

    assert_true(matrix_type(mat) %in% c("uint32_t", "float", "double"))
    if (compress && matrix_type(mat) != "uint32_t") {
        rlang::inform(c(
            "Warning: Matrix compression performs poorly with non-integers.",
            "Consider calling convert_matrix_type if a compressed integer matrix is intended."
        ), .frequency="regularly", .frequency_id="matrix_compress_non_integer")
    }
    
    dir <- path.expand(dir)
    
    if (mat@transpose) {
        tmp_path <- tempfile(pattern="transpose")
        it <- write_matrix_transpose(t(mat), tmp_path, load_bytes=4194304L, sort_bytes=1073741824L)
    }
    else it <- iterate_matrix(mat)
    
    write_function <- get(sprintf("write_%s_matrix_file_%s_cpp", ifelse(compress, "packed", "unpacked"), matrix_type(mat)))
    write_function(ptr(it), dir, buffer_size)
    
    if (mat@transpose) unlink(tmp_path, recursive=TRUE, expand=FALSE)    

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
        dimnames=normalized_dimnames(info$row_names, info$col_names), type=info$type)
}

setClass("MatrixH5",
    contains = "IterableMatrix",
    slots = c(
        path = "character",
        group = "character",
        compressed = "logical",
        buffer_size = "integer",
        type = "character"
    ),
    prototype = list(
        path = character(0),
        group = "",
        compressed = logical(0),
        buffer_size = integer(0),
        type = character(0)
    )
)
setMethod("matrix_type", "MatrixH5", function(x) x@type)

setMethod("iterate_matrix", "MatrixH5", function(x) {
    x@dimnames <- denormalize_dimnames(x@dimnames)

    iter_function <- get(sprintf("iterate_%s_matrix_hdf5_%s_cpp", ifelse(x@compressed, "packed", "unpacked"), x@type))
    wrap_function <- get(sprintf("wrapMat_%s", x@type))

    wrap_function(iter_function(x@path, x@group, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
})

setMethod("short_description", "MatrixH5", function(x) {
    sprintf("Load %s matrix in hdf5 file %s, group %s", 
        if(x@compressed) "compressed" else "uncompressed",
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

    assert_is(path, "character")
    assert_is(group, "character")
    assert_is(compress, "logical")
    assert_is(buffer_size, "integer")
    assert_is(chunk_size, "integer")

    assert_true(matrix_type(mat) %in% c("uint32_t", "float", "double"))
    if (compress && matrix_type(mat) != "uint32_t") {
        rlang::inform(c(
            "Warning: Matrix compression performs poorly with non-integers.",
            "Consider calling convert_matrix_type if a compressed integer matrix is intended."
        ), .frequency="regularly", .frequency_id="matrix_compress_non_integer")
    }

    path <- path.expand(path)

    if (mat@transpose) {
        tmp_path <- tempfile(pattern="transpose")
        it <- write_matrix_transpose(t(mat), tmp_path, load_bytes=4194304L, sort_bytes=1073741824L)
    }
    else it <- iterate_matrix(mat)
    
    write_function <- get(sprintf("write_%s_matrix_hdf5_%s_cpp", ifelse(compress, "packed", "unpacked"), matrix_type(mat)))
    write_function(ptr(it), path, group, buffer_size, chunk_size)
    
    if (mat@transpose) unlink(tmp_path, recursive=TRUE, expand=FALSE)    
    
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
        dimnames=normalized_dimnames(info$row_names, info$col_names), type=info$type)
}

setClass("10xMatrixH5",
    contains = "IterableMatrix",
    slots = c(
        path = "character",
        buffer_size = "integer"
    ),
    prototype = list(
        path = character(0),
        buffer_size = integer(0)
    )
)
setMethod("matrix_type", "10xMatrixH5", function(x) "uint32_t")
setMethod("iterate_matrix", "10xMatrixH5", function(x) {
    wrapMat_uint32_t(iterate_matrix_10x_hdf5_cpp(x@path, x@buffer_size), new("XPtrList"))
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
    contains = "IterableMatrix",
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
setMethod("matrix_type", "AnnDataMatrixH5", function(x) "float")
setMethod("iterate_matrix", "AnnDataMatrixH5", function(x) {
    wrapMat_double(iterate_matrix_anndata_hdf5_cpp(x@path, x@group, x@buffer_size), new("XPtrList"))
})
setMethod("short_description", "AnnDataMatrixH5", function(x) {
    sprintf("AnnData HDF5 matrix in file %s, group %s", 
        x@path, x@group
    )
})

#' Read a sparse integer matrix from an anndata matrix in an hdf5 file
#' @inheritParams open_matrix_hdf5
#' @return 10xMatrixH5 object, with cells as the columns.
#' @export
open_matrix_anndata_hdf5 <- function(path, group="X", buffer_size=16384L) {
    assert_is_file(path)
    assert_is(buffer_size, "integer")

    info <- dims_matrix_anndata_hdf5_cpp(path.expand(path), group, buffer_size)
    res <- new("AnnDataMatrixH5", path=path.expand(path), dim=info$dims, buffer_size=buffer_size,
        group=group, dimnames=normalized_dimnames(info$row_names, info$col_names))
    
    # We do the reverse of what transpose says, because anndata files are usually
    # stored row-major with cells as rows, whereas BPCells will work best with
    # col-major with cells as cols.
    if (info[["transpose"]]) {
        return(res)
    } else {
        return(t(res))
    }
}



# Overlap matrix from fragments
setClass("PeakMatrix",
    contains = "IterableMatrix",
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
setMethod("matrix_type", "PeakMatrix", function(x) "uint32_t")
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
    
    assert_not_null(cellNames(fragments))
    assert_not_na(cellNames(fragments))
    assert_not_null(chrNames(fragments))
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
    wrapMat_uint32_t(iterate_peak_matrix_cpp(ptr(it), x@chr_id, x@start, x@end, x@chr_levels), it)
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
    contains = "IterableMatrix",
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
setMethod("matrix_type", "TileMatrix", function(x) "uint32_t")

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
    
    assert_not_null(cellNames(fragments))
    assert_not_na(cellNames(fragments))
    assert_not_null(chrNames(fragments))
    assert_not_na(chrNames(fragments))
    
    if(is.list(ranges)) {
        assert_has_names(ranges, c("chr", "start", "end", "tile_width"))
        assert_true(all(as.character(ranges$chr) %in% chrNames(fragments)))

        chr_id <- as.integer(factor(as.character(ranges$chr), chrNames(fragments))) - 1L
        start <- as.integer(ranges$start) - !zero_based_coords
        end <- as.integer(ranges$end)
        tile_width <- as.integer(ranges$tile_width)
    } else {
        assert_has_names(mcols(ranges), c("tile_width"))
        assert_true(all(as.character(GenomicRanges::seqnames(ranges)) %in% chrNames(fragments)))
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
    wrapMat_uint32_t(iterate_tile_matrix_cpp(ptr(it), x@chr_id, x@start, x@end, x@tile_width, x@chr_levels), it)
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

# Convert matrix types
setClass("ConvertMatrixType",
    contains = "IterableMatrix",
    slots = c(
        matrix = "IterableMatrix",
        type = "character"
    ),
    prototype = list(
        matrix = NULL,
        type = character(0)
    )
)
setMethod("matrix_type", signature(x="ConvertMatrixType"), function(x) x@type)
setMethod("iterate_matrix", "ConvertMatrixType", function(x) {
    iter_function <- get(sprintf("convert_matrix_%s_%s_cpp", matrix_type(x@matrix), x@type))
    wrap_function <- get(sprintf("wrapMat_%s", x@type))

    it <- iterate_matrix(x@matrix)
    wrap_function(iter_function(ptr(it)), it)
})
setMethod("short_description", "ConvertMatrixType", function(x) {
    c(
        short_description(x@matrix),
        sprintf("Convert type from %s to %s", matrix_type(x@matrix), x@type)
    )
})

convert_matrix_type <- function(matrix, type=c("uint32_t", "double", "float")) {
    assert_is(matrix, c("dgCMatrix", "IterableMatrix"))
    type <- match.arg(type)
    if (is(matrix, "dgCMatrix")) matrix <- as(matrix, "IterableMatrix")
    if (matrix_type(matrix) == type)
        return(matrix)
    wrapMatrix("ConvertMatrixType", matrix, type=type)
}

# Conversions with dgCMatrix
setClass("Iterable_dgCMatrix_wrapper",
    contains = "IterableMatrix",
    slots = c(
        mat = "dgCMatrix"
    ),
    prototype = list(
        mat = NULL
    )
)
setMethod("matrix_type", signature(x="Iterable_dgCMatrix_wrapper"), function(x) "double")

setAs("dgCMatrix", "IterableMatrix", function(from) {
    new("Iterable_dgCMatrix_wrapper", dim=dim(from), dimnames=dimnames(from), transpose=FALSE, mat=from)
})
setMethod("iterate_matrix", "Iterable_dgCMatrix_wrapper", function(x) {
    x@dimnames <- denormalize_dimnames(x@dimnames)
    wrapMat_double(iterate_csparse_matrix_cpp(x@mat, x@dimnames[[1]], x@dimnames[[2]]), new("XPtrList"))
})
setMethod("short_description", "Iterable_dgCMatrix_wrapper", function(x) {
    "Load dgCMatrix from memory"
})

setMethod("iterate_matrix", "dgCMatrix", function(x) {
    wrapMat_double(iterate_csparse_matrix_cpp(x), new("XPtrList"))
})

setAs("IterableMatrix", "dgCMatrix", function(from) {
    iter <- iterate_matrix(convert_matrix_type(from, "double"))
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