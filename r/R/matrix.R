# Copyright 2021 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

#' IterableMatrix methods
#'
#' Generic methods and built-in functions for IterableMatrix objects
#'
#' @name IterableMatrix-methods
#' @rdname IterableMatrix-methods
NULL

setClass("IterableMatrix",
  slots = c(
    dim = "numeric",
    transpose = "logical",
    dimnames = "list"
  ),
  prototype = list(
    dim = numeric(2),
    transpose = FALSE,
    dimnames = list(NULL, NULL)
  )
)

#' Construct an S4 matrix object wrapping another matrix object
#'
#' Helps to avoid duplicate storage of dimnames
#' @keywords internal
wrapMatrix <- function(class, m, ...) {
  dimnames <- dimnames(m)
  if (matrix_is_transform(m) && !is(m, "RenameDims")) m@dimnames <- list(NULL, NULL)
  new(class, matrix = m, transpose = m@transpose, dim = m@dim, dimnames = dimnames, ...)
}

#' Helper function to set dimnames to NULL instead of 0-length character vectors
#' @keywords internal
normalized_dimnames <- function(row_names, col_names) {
  list(
    if (length(row_names) == 0) NULL else row_names,
    if (length(col_names) == 0) NULL else col_names
  )
}

denormalize_dimnames <- function(dimnames) {
  if (is.null(dimnames[[1]])) dimnames[[1]] <- character(0)
  if (is.null(dimnames[[2]])) dimnames[[2]] <- character(0)
  dimnames
}

#' Get a wrapped pointer to the iterable matrix
#' @keywords internal
setGeneric("iterate_matrix", function(x) standardGeneric("iterate_matrix"))

#' @describeIn IterableMatrix-methods Get the matrix data type (mat_uint32_t, mat_float, or mat_double for now)
#' @examples
#' ## Prep data
#' mat <- matrix(1:25, nrow = 5) %>% as("dgCMatrix")
#' mat
#' mat <- as(mat, "IterableMatrix")
#' mat
#' 
#' 
#' #######################################################################
#' ## matrix_type() example
#' #######################################################################
#' matrix_type(mat)
#' 
#' 
#' @export
setGeneric("matrix_type", function(x) standardGeneric("matrix_type"))

#' @describeIn IterableMatrix-methods Get the matrix storage order ("row" or "col")
#' @examples
#' #######################################################################
#' ## storage_order() example
#' #######################################################################
#' storage_order(mat)
#' 
#' 
#' @export
setGeneric("storage_order", function(x) standardGeneric("storage_order"))

setMethod("storage_order", "IterableMatrix", function(x) if(x@transpose) "row" else "col")

#' Return a list of input matrices to the current matrix (experimental)
#'
#' File objects have 0 inputs. Most transforms have 1 input. Some transforms
#' (e.g. matrix multiplication or matrix concatenation) can have multiple
#' This is used primarily to know when it is safe to clear dimnames from intermediate transformed matrices.
#' C++ relies on the base matrices (non-transform) to have dimnames, while R relies on the outermost matrix (transform) to have dimnames.
#' @keywords internal
setGeneric("matrix_inputs", function(x) standardGeneric("matrix_inputs"))
matrix_is_transform <- function(x) length(matrix_inputs(x)) != 0

# As a reasonable default convention, a matrix is a transform if it has a slot named matrix
setMethod("matrix_inputs", "IterableMatrix", function(x) {
  if (.hasSlot(x, "matrix")) {
    return(list(x@matrix))
  }
  rlang::abort(sprintf("matrix_inputs not defined for inputs of type %s", class(x)))
})

setGeneric("matrix_inputs<-", function(x, ..., value) standardGeneric("matrix_inputs<-"))
setMethod("matrix_inputs<-", "IterableMatrix", function(x, ..., value) {
  if (.hasSlot(x, "matrix")) {
    assert_is(value, "list")
    assert_len(value, 1)
    assert_is(value[[1]], "IterableMatrix")
    x@matrix <- value[[1]]
    return(x)
  }
  rlang::abort(sprintf("matrix_inputs not defined for inputs of type %s", class(x)))
})

#' Get/set inputs to a matrix transform
#' 
#' A matrix object can either be an input (i.e. a file on disk or a raw matrix in memory),
#' or it can represent a delayed operation on one or more matrices. The `all_matrix_inputs()`
#' getter and setter functions allow accessing the base-level input matrices as a list, and
#' changing them. This is useful if you want to re-locate data on disk without losing your
#' transformed BPCells matrix. (Note: experimental API; potentially subject to revisions).
#'
#' 
#'
#' @param x IterableMatrix
#' @param value List of IterableMatrix objects
#' @return List of IterableMatrix objects. If a matrix `m` is itself an input object, then 
#'   `all_matrix_inputs(m)` will return `list(m)`.
#' @export
all_matrix_inputs <- function(x) {
  assert_is(x, "IterableMatrix")
  inputs <- matrix_inputs(x)
  if (length(inputs) == 0) return(list(x))
  do.call(c, lapply(inputs, all_matrix_inputs))
}

#' @rdname all_matrix_inputs
#' @export
`all_matrix_inputs<-` <- function(x, value) {
  # Error checking:
  # - value should be a list of iterable matrices
  # - value should be the same length as the existing matrix inputs
  # - Each iterable matrix should have the same dimensions as the one it's replacing
  assert_is(value, "list")
  prev_inputs <- all_matrix_inputs(x)
  if (length(prev_inputs) != length(value)) {
    rlang::abort("Error assigning matrix inputs: length mismatch of input list")
  }
  for (i in seq_along(value)) {
    if (!is(value[[i]], "IterableMatrix")) {
      rlang::abort(sprintf("Error assigning matrix inputs: Entry %d is not an IterableMatrix", i))
    }
    if (!all(dim(prev_inputs[[i]]) == dim(value[[i]]))) {
      rlang::warn(sprintf("Warning assigning matrix inputs: entry %d has mismatched dimensions with the matrix it replaces", i))
    }
  }

  # Actual work -- recursively reset the inputs
  direct_inputs <- matrix_inputs(x)
  if (length(direct_inputs) == 0) {
    stopifnot(length(value) == 1)
    return(value[[1]])
  }
  subtree_counts <- lapply(direct_inputs, function(x) length(all_matrix_inputs(x)))
  stopifnot(sum(as.integer(subtree_counts)) == length(value))
  start <- 1L
  for (i in seq_along(direct_inputs)) {
    end <- start + subtree_counts[[i]] - 1L
    all_matrix_inputs(direct_inputs[[i]]) <- value[start:end]
    start <- end + 1L
  }
  matrix_inputs(x) <- direct_inputs
  x
}

setMethod("short_description", "IterableMatrix", function(x) {
  character(0)
})


#' @describeIn IterableMatrix-methods Display an IterableMatrix
#' @param object IterableMatrix object
#' @examples
#' #######################################################################
#' ## show() example
#' #######################################################################
#' show(mat)
#' 
#' 
setMethod("show", "IterableMatrix", function(object) {
  cat(sprintf("%d x %d IterableMatrix object with class %s\n", nrow(object), ncol(object), class(object)))

  cat("\n")
  cat(sprintf(
    "Row names: %s\n",
    pretty_print_vector(rownames(object), empty = "unknown names")
  ))
  cat(sprintf(
    "Col names: %s\n",
    pretty_print_vector(colnames(object), empty = "unknown names")
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
#' @examples
#' #######################################################################
#' ## t() example
#' #######################################################################
#' t(mat)
#' 
#' 
#' @export
setMethod("t", signature(x = "IterableMatrix"), function(x) {
  x@transpose <- !x@transpose
  x@dim <- c(x@dim[2L], x@dim[1L])
  x@dimnames <- list(x@dimnames[[2]], x@dimnames[[1]])
  if (matrix_is_transform(x)) {
    matrix_inputs(x) <- lapply(matrix_inputs(x), t)
  }
  return(x)
})

# Dense multiply operators (sparse*dense_mat and sparse*dense_vec)

#' @param x IterableMatrix object
#' @param y matrix
#' @describeIn IterableMatrix-methods Multiply by a dense matrix
#' @return * `x %*% y`: dense matrix result
#' @examples
#' #######################################################################
#' ## `x %*% y` example
#' #######################################################################
#' mat %*% as(matrix(1:50, nrow = 5), "dgCMatrix")
#' 
#' 
setMethod("%*%", signature(x = "IterableMatrix", y = "matrix"), function(x, y) {
  iter <- iterate_matrix(convert_matrix_type(x, "double"))
  if (x@transpose) {
    res <- (t(dense_multiply_left_cpp(iter, t(y))))
  } else {
    res <- (dense_multiply_right_cpp(iter, y))
  }
  rownames(res) <- rownames(x)
  colnames(res) <- colnames(y)
  res
})

setMethod("%*%", signature(x = "matrix", y = "IterableMatrix"), function(x, y) {
  iter <- iterate_matrix(convert_matrix_type(y, "double"))
  if (y@transpose) {
    res <- (t(dense_multiply_right_cpp(iter, t(x))))
  } else {
    res <- (dense_multiply_left_cpp(iter, x))
  }
  rownames(res) <- rownames(x)
  colnames(res) <- colnames(y)
  res
})

setMethod("%*%", signature(x = "IterableMatrix", y = "numeric"), function(x, y) {
  iter <- iterate_matrix(convert_matrix_type(x, "double"))
  if (x@transpose) {
    res <- (vec_multiply_left_cpp(iter, y))
  } else {
    res <- (vec_multiply_right_cpp(iter, y))
  }
  res <- matrix(res, ncol=1)
  rownames(res) <- rownames(x)
  res
})

setMethod("%*%", signature(x = "numeric", y = "IterableMatrix"), function(x, y) {
  iter <- iterate_matrix(convert_matrix_type(y, "double"))
  if (y@transpose) {
    res <- (vec_multiply_right_cpp(iter, x))
  } else {
    res <- (vec_multiply_left_cpp(iter, x))
  }
  res <- matrix(res, nrow=1)
  colnames(res) <- colnames(y)
  res
})

#' Represent a sparse matrix-vector product operation
#'
#' LinearOperators perform sparse matrix-vector product operations for
#' for downstream matrix solvers. They avoid repeatedly calling iterate_matrix
#' from an SVD solver for a possible efficiency gain
#' @keywords internal
setClass("LinearOperator",
  slots = c(
    dim = "integer",
    xptr = "externalptr",
    transpose = "logical"
  ),
  prototype = list(
    dim = integer(2),
    transpose = FALSE
  )
)

#' Construct a LinearOperator object
#'
#' Constructs a C++ matrix object and save the pointer to use for repeated matrix-vector products
#' A bit experimental still so for internal use
#' @keywords internal
linear_operator <- function(mat) {
  assert_is(mat, "IterableMatrix")
  new("LinearOperator", dim = dim(mat), xptr = iterate_matrix(convert_matrix_type(mat, "double")), transpose = mat@transpose)
}

setMethod("%*%", signature(x = "LinearOperator", y = "matrix"), function(x, y) {
  if (x@transpose) {
    return(t(dense_multiply_left_preserve_loader_cpp(x@xptr, t(y))))
  } else {
    return(dense_multiply_right_preserve_loader_cpp(x@xptr, y))
  }
})

setMethod("%*%", signature(x = "matrix", y = "LinearOperator"), function(x, y) {
  if (y@transpose) {
    return(t(dense_multiply_right_preserve_loader_cpp(y@xptr, t(x))))
  } else {
    return(dense_multiply_left_preserve_loader_cpp(y@xptr, x))
  }
})

setMethod("%*%", signature(x = "LinearOperator", y = "numeric"), function(x, y) {
  if (x@transpose) {
    return(vec_multiply_left_preserve_loader_cpp(x@xptr, y))
  } else {
    return(vec_multiply_right_preserve_loader_cpp(x@xptr, y))
  }
})

setMethod("%*%", signature(x = "numeric", y = "LinearOperator"), function(x, y) {
  if (y@transpose) {
    return(vec_multiply_right_preserve_loader_cpp(y@xptr, x))
  } else {
    return(vec_multiply_left_preserve_loader_cpp(y@xptr, x))
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
setMethod("matrix_type", signature(x = "MatrixMultiply"), function(x) matrix_type(x@left))
setMethod("matrix_inputs", "MatrixMultiply", function(x) list(x@left, x@right))
setMethod("matrix_inputs<-", "MatrixMultiply", function(x, ..., value) {
  assert_is(value, "list")
  assert_len(value, 2)
  for (v in value) assert_is(v, "IterableMatrix")
  x@left <- value[[1]]
  x@right <- value[[2]]
  x
})

setMethod("iterate_matrix", "MatrixMultiply", function(x) {
  iter_function <- get(sprintf("iterate_matrix_multiply_%s_cpp", matrix_type(x)))
  iter_function(iterate_matrix(x@left), iterate_matrix(x@right))
})

setMethod("short_description", "MatrixMultiply", function(x) {
  if (x@transpose) {
    # Flip the display order for transposed case
    sprintf(
      "Multiply sparse matrices: %s (%dx%d) * %s (%dx%d)",
      class(x@right), ncol(x@right), nrow(x@right),
      class(x@left), ncol(x@left), nrow(x@left)
    )
  } else {
    sprintf(
      "Multiply sparse matrices: %s (%dx%d) * %s (%dx%d)",
      class(x@left), nrow(x@left), ncol(x@left),
      class(x@right), nrow(x@right), ncol(x@right)
    )
  }
})

setMethod("%*%", signature(x = "IterableMatrix", y = "IterableMatrix"), function(x, y) {
  if (x@transpose != y@transpose) stop("Cannot multiply matrices with different internal transpose states.\nPlease use transpose_storage_order().")
  if (x@transpose) {
    return(t(t(y) %*% t(x)))
  }

  assert_true(ncol(x) == nrow(y))

  # If types are mismatched, default to double precision for both
  type_x <- matrix_type(x)
  type_y <- matrix_type(y)
  if (type_x != type_y && type_x != "double") x <- convert_matrix_type(x, "double")
  if (type_x != type_y && type_y != "double") y <- convert_matrix_type(y, "double")

  dim <- c(nrow(x), ncol(y))
  dimnames <- list(rownames(x), colnames(y))
  new("MatrixMultiply", left = x, right = y, transpose = FALSE, dim = dim, dimnames = dimnames)
})

setMethod("%*%", signature(x = "IterableMatrix", y = "dgCMatrix"), function(x, y) {
  if (x@transpose) {
    t(as(t(y), "IterableMatrix") %*% t(x))
  } else {
    x %*% as(y, "IterableMatrix")
  }
})

setMethod("%*%", signature(x = "dgCMatrix", y = "IterableMatrix"), function(x, y) {
  if (y@transpose) {
    t(t(y) %*% as(t(x), "IterableMatrix"))
  } else {
    as(x, "IterableMatrix") %*% y
  }
})


# Subsetting on MatrixMultiply
setMethod("[", "MatrixMultiply", function(x, i, j, ...) {
  if (missing(x)) stop("x is missing in matrix selection")
  # Handle transpose via recursive call
  if (x@transpose) {
    return(t(t(x)[rlang::maybe_missing(j), rlang::maybe_missing(i)]))
  }

  i <- split_selection_index(i, nrow(x), rownames(x))
  j <- split_selection_index(j, ncol(x), colnames(x))
  # If we're just reordering rows/cols, do a standard matrix selection
  if (rlang::is_missing(i$subset) && rlang::is_missing(j$subset)) {
    return(callNextMethod(x, unsplit_selection(i), unsplit_selection(j)))
  }
  x <- selection_fix_dims(x, rlang::maybe_missing(i$subset), rlang::maybe_missing(j$subset))

  # Selection will be a no-op if i or j is missing
  x@left <- x@left[rlang::maybe_missing(i$subset), ]
  x@right <- x@right[, rlang::maybe_missing(j$subset)]

  x[rlang::maybe_missing(i$reorder),rlang::maybe_missing(j$reorder)]
})

setClass("MatrixMask",
  contains = "IterableMatrix",
  slots = c(
    matrix = "IterableMatrix",
    mask = "IterableMatrix",
    invert = "logical"
  ),
  prototype = list(
    matrix = NULL,
    mask = NULL,
    invert = FALSE
  )
)
setMethod("matrix_type", signature(x = "MatrixMask"), function(x) matrix_type(x@matrix))
setMethod("matrix_inputs", "MatrixMask", function(x) list(x@matrix, x@mask))
setMethod("matrix_inputs<-", "MatrixMask", function(x, ..., value) {
  assert_is(value, "list")
  assert_len(value, 2)
  for (v in value) assert_is(v, "IterableMatrix")
  x@matrix <- value[[1]]
  x@mask <- value[[2]]
  x
})

setMethod("iterate_matrix", "MatrixMask", function(x) {
  iter_function <- get(sprintf("iterate_matrix_mask_%s_cpp", matrix_type(x)))
  iter_function(iterate_matrix(x@matrix), iterate_matrix(x@mask), x@invert)
})

setMethod("short_description", "MatrixMask", function(x) {
  c(
    short_description(x@matrix),
    sprintf(
      "Mask entries according to matrix %s %s",
      class(x@mask),
      ifelse(x@invert, "(inverted)", "")
    )
  )
})

#' Mask matrix entries to zero
#' Set matrix entries to zero given a mask matrix of the 
#' same dimensions. Normally, non-zero values in the mask
#' will set the matrix entry to zero. If inverted, zero
#' values in the mask matrix will set the matrix entry to zero.
#' @param mat Data matrix (IterableMatrix)
#' @param mask Mask matrix (IterableMatrix or dgCMatrix)
#' @keywords internal
mask_matrix <- function(mat, mask, invert=FALSE) {
  assert_is(mat, "IterableMatrix")
  assert_is(invert, "logical")
  assert_is(mask, c("IterableMatrix", "dgCMatrix"))
  assert_true(nrow(mat) == nrow(mask) && ncol(mat) == ncol(mask))
  if (is(mask, "dgCMatrix")) {
    if (mat@transpose)
      mask <- t(as(t(mask), "IterableMatrix"))
    else
      mask <- as(mask, "IterableMatrix")
  }
  
  if (mat@transpose != mask@transpose) stop("Cannot mask matrices with different internal transpose states.\nPlease use transpose_storage_order().")
  mask <- convert_matrix_type(mask, "uint32_t")

  wrapMatrix("MatrixMask",
    mat,
    mask = mask,
    invert = invert
  )
}

setClass("MatrixRankTransform",
  contains = "IterableMatrix",
  slots = c(
    matrix = "IterableMatrix"
  ),
  prototype = list(
    matrix = NULL
  )
)
setMethod("matrix_type", signature(x = "MatrixRankTransform"), function(x) "double")
setMethod("iterate_matrix", "MatrixRankTransform", function(x) {
  iter_function <- get(sprintf("iterate_matrix_rank_%s_cpp", matrix_type(x@matrix)))
  iter_function(iterate_matrix(x@matrix))
})

setMethod("short_description", "MatrixRankTransform", function(x) {
  c(
    short_description(x@matrix),
    sprintf(
      "Rank transform each matrix %s",
      ifelse(x@transpose, "row", "col")
    )
  )
})

#' Rank-transform a matrix
#' 
#' Rank the values within each row/col of a matrix, and output
#' the rank values as a new matrix. Rank values are offset such
#' that the rank of a 0 value is 0, and ties are handled by
#' averaging ranks.
#'
#' Note that efficient rank calculation depends on the storage order
#' of a matrix, so it may be necessary to call transpose_storage_order()
#'
#' @param mat Data matrix (IterableMatrix)
#' @param axis Axis to rank values within. "col" to rank values within each column,
#'     and "row" to rank values within each row.
#' @keywords internal
rank_transform <- function(mat, axis) {
  assert_is(mat, "IterableMatrix")
  assert_is_character(axis)
  assert_len(axis, 1)
  assert_true(axis %in% c("row", "col"))
  assert_true(storage_order(mat) == axis)

  wrapMatrix("MatrixRankTransform", mat)
}

# Row sums and row means

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate rowSums
#' @return * `rowSums()`: vector of row sums
#' @examples
#' #######################################################################
#' ## rowSums() example
#' #######################################################################
#' rowSums(mat)
#' 
#' 
setMethod("rowSums", signature(x = "IterableMatrix"), function(x) {
  iter <- iterate_matrix(convert_matrix_type(x, "double"))
  if (x@transpose) {
    res <- col_sums_double_cpp(iter)
  } else {
    res <- row_sums_double_cpp(iter)
  }
  names(res) <- rownames(x)
  res
})

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate colSums
#' @return * `colSums()`: vector of col sums
#' @examples
#' #######################################################################
#' ## colSums() example
#' #######################################################################
#' colSums(mat)
#' 
#' 
setMethod("colSums", signature(x = "IterableMatrix"), function(x) {
  iter <- iterate_matrix(convert_matrix_type(x, "double"))
  if (x@transpose) {
    res <- row_sums_double_cpp(iter)
  } else {
    res <- col_sums_double_cpp(iter)
  }
  names(res) <- colnames(x)
  res
})

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate rowMeans
#' @return * `rowMeans()`: vector of row means
#' @examples
#' #######################################################################
#' ## rowMeans() example
#' #######################################################################
#' rowMeans(mat)
#' 
#' 
setMethod("rowMeans", signature(x = "IterableMatrix"), function(x) rowSums(x) / ncol(x))

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate colMeans
#' @return * `colMeans()`: vector of col means
#' @examples
#' #######################################################################
#' ## colMeans() example
#' #######################################################################
#' colMeans(mat)
#' 
#' 
setMethod("colMeans", signature(x = "IterableMatrix"), function(x) colSums(x) / nrow(x))


# Strategy for rowVars and colVars:
#  - matrixStats::rowVars and matrixStats::colVars are not generic, and don't accept BPCells objects
#  - The bioconductor MatrixGenerics package would be inconvenient to add as a hard dependency because it's not in CRAN
#  - BPCells registers S3 methods BPCells::rowVars and BPCells::colVars which work with BPCells objects, and will
#    fall back to matrixStats or MatrixGenerics implementations if available (otherwise erroring for other matrix inputs)
#  - If MatrixGenerics in installed, BPCells will also register as a generic with it
#  - In summary, BPCells::rowVars and BPCells::colVars will work on all inputs, and so will MatrixGenerics::rowVars and
#    MatrixGenerics::colVars. matrixStats::rowVars and matrixStats::colVars will only work on base R matrix objects.

#' @describeIn IterableMatrix-methods Calculate colVars (replacement for `matrixStats::colVars()`)
#' @return * `colVars()`: vector of col variance
#' @examples
#' #######################################################################
#' ## colVars() example
#' #######################################################################
#' colVars(mat)
#' 
#' 
#' @export
colVars <- function(x, rows = NULL, cols = NULL, na.rm = FALSE, center = NULL, ..., useNames = TRUE) UseMethod("colVars")
#' @export
colVars.default <- function(x, rows = NULL, cols = NULL, na.rm = FALSE, center = NULL, ..., useNames = TRUE) {
  if (requireNamespace("MatrixGenerics", quietly = TRUE)) {
    MatrixGenerics::colVars(x, rows=rows, cols=cols, na.rm=na.rm, center=center, ..., useNames=useNames)
  } else if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colVars(x, rows=rows, cols=cols, na.rm=na.rm, center=center, ..., useNames=useNames)
  } else {
    stop("Can't run colVars on a non-BPCells object unless MatrixGenerics or matrixStats are installed.")
  }
}
#' @export
colVars.IterableMatrix <- function(x, rows = NULL, cols = NULL, na.rm = FALSE, center = NULL, ..., useNames = TRUE) {
  if (!is.null(rows) || !is.null(cols) || !isFALSE(na.rm) || !is.null(center) || !isTRUE(useNames)) {
    stop("colVars(IterableMatrix) doesn't support extra arguments rows, cols, na.rm, center, or useNames")
  }
  matrix_stats(x, col_stats="variance")$col_stats["variance",]
}
rlang::on_load({
  if (requireNamespace("MatrixGenerics", quietly=TRUE)) {
    setMethod(MatrixGenerics::colVars, "IterableMatrix", colVars.IterableMatrix)
  }
})

#' @describeIn IterableMatrix-methods Calculate rowVars (replacement for `matrixStats::rowVars()`)
#' @return * `rowVars()`: vector of row variance
#' @export
rowVars <- function(x, rows = NULL, cols = NULL, na.rm = FALSE, center = NULL, ..., useNames = TRUE) UseMethod("rowVars")
#' @export
rowVars.default <- function(x, rows = NULL, cols = NULL, na.rm = FALSE, center = NULL, ..., useNames = TRUE) {
  if (requireNamespace("MatrixGenerics", quietly = TRUE)) {
    MatrixGenerics::rowVars(x, rows=rows, cols=cols, na.rm=na.rm, center=center, ..., useNames=useNames)
  } else if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowVars(x, rows=rows, cols=cols, na.rm=na.rm, center=center, ..., useNames=useNames)
  } else {
    stop("Can't run rowVars on a non-BPCells object unless MatrixGenerics or matrixStats are installed.")
  }
}
#' @export
rowVars.IterableMatrix <- function(x, rows = NULL, cols = NULL, na.rm = FALSE, center = NULL, ..., useNames = TRUE) {
  if (!is.null(rows) || !is.null(cols) || !isFALSE(na.rm) || !is.null(center) || !isTRUE(useNames)) {
    stop("rowVars(IterableMatrix) doesn't support extra arguments rows, cols, na.rm, center, or useNames")
  }
  matrix_stats(x, row_stats="variance")$row_stats["variance",]
}
rlang::on_load({
  if (requireNamespace("MatrixGenerics", quietly=TRUE)) {
    setMethod(MatrixGenerics::rowVars, "IterableMatrix", rowVars.IterableMatrix)
  }
})

#' Get the max of each row in an iterable matrix
#' @param x IterableMatrix object/dgCMatrix object
#' @return * `rowMaxs()`: vector of maxes for every row
#' @describeIn IterableMatrix-methods Calculate rowMaxs (replacement for `matrixStats::rowMaxs()`)
#' @examples
#' #######################################################################
#' ## rowMaxs() example
#' #######################################################################
#' rowMaxs(mat)
#' 
#' 
#' @export
rowMaxs <- function(x, rows = NULL, cols = NULL, na.rm = FALSE, ...) UseMethod("rowMaxs")
#' @export
rowMaxs.default <- function(x, rows = NULL, cols = NULL, na.rm = FALSE, ...) {
  if (requireNamespace("MatrixGenerics", quietly = TRUE)) {
    MatrixGenerics::rowMaxs(x, rows = rows, cols = cols, na.rm = na.rm, ...)
  } else if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowMaxs(x, rows = rows, cols = cols, na.rm = na.rm, ...)
  }
  else {
    stop("Can't run rowMaxs on a non-BPCells object unless MatrixGenerics or matrixStats are installed.")
  }
}
#' @export
rowMaxs.IterableMatrix <- function(x, rows = NULL, cols = NULL, na.rm = FALSE, ...) {
  if(!is.null(rows) || !is.null(cols) || !isFALSE(na.rm)) {
    stop("rowMaxs(IterableMatrix) doesn't support extra arguments rows, cols, or na.rm")
  }
  iter <- iterate_matrix(convert_matrix_type(x, "double"))
  if(x@transpose == TRUE) {
    res <- matrix_max_per_col_cpp(iter)
  } else {
    res <- matrix_max_per_row_cpp(iter)
  }
  names(res) <- rownames(x)
  return(res)
}
rlang::on_load({
  if (requireNamespace("MatrixGenerics", quietly=TRUE)) {
    setMethod(MatrixGenerics::rowMaxs, "IterableMatrix", rowMaxs.IterableMatrix)
  }
})

#' Get the max of each col in an interable matrix
#' @param x IterableMatrix/dgCMatrix object
#' @return * `colMaxs()`: vector of column maxes
#' @describeIn IterableMatrix-methods Calculate colMax (replacement for `matrixStats::colMax()`)
#' @examples
#' #######################################################################
#' ## colMaxs() example
#' #######################################################################
#' colMaxs(mat)
#' 
#' 
#' @export
colMaxs <- function(x, rows = NULL, cols = NULL, na.rm = FALSE, ...) UseMethod("colMaxs")
#' @export
colMaxs.default <- function(x, rows = NULL, cols = NULL, na.rm = FALSE, ...) {
  if (requireNamespace("MatrixGenerics", quietly = TRUE)) {
    MatrixGenerics::colMaxs(x, rows = rows, cols = cols, na.rm = na.rm, ...)
  } else if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colMaxs(x, rows = rows, cols = cols, na.rm = na.rm, ...)
  }
  else {
    stop("Can't run colMaxs on a non-BPCells object unless MatrixGenerics or matrixStats are installed.")
  }
}
#' @export
colMaxs.IterableMatrix <- function(x, rows = NULL, cols = NULL, na.rm = FALSE, ...) {
  iter <- iterate_matrix(convert_matrix_type(x, "double"))
  if(x@transpose == TRUE) {
    res <- matrix_max_per_row_cpp(iter)
  } else {
    res <- matrix_max_per_col_cpp(iter)
  }
  names(res) <- colnames(x)
  return(res)
}
rlang::on_load({
  if (requireNamespace("MatrixGenerics", quietly=TRUE)) {
    setMethod(MatrixGenerics::colMaxs, "IterableMatrix", colMaxs.IterableMatrix)
  }
})

# Index subsetting
setClass("MatrixSubset",
  contains = "IterableMatrix",
  slots = c(
    matrix = "IterableMatrix",
    row_selection = "integer",
    col_selection = "integer",
    zero_dims = "logical"
  ),
  prototype = list(
    row_selection = integer(0),
    col_selection = integer(0),
    zero_dims = c(FALSE, FALSE)
  )
)
setMethod("matrix_type", signature(x = "MatrixSubset"), function(x) matrix_type(x@matrix))

# Helper function to convert logical/character indexing into numeric indexing
selection_index <- function(selection, dim_len, dimnames) {
  if (!rlang::is_missing(selection)) {
    indices <- seq_len(dim_len)
    names(indices) <- dimnames
    if (!is.logical(selection)) assert_distinct(selection, n = 2)
    res <- vctrs::vec_slice(indices, selection)
    names(res) <- NULL
    return(res)
  } else {
    return(rlang::missing_arg())
  }
}
# Helper function for fixing dims and dimnames after selection.
# i and j must be integer indexes or missing
selection_fix_dims <- function(x, i, j) {
  if (!rlang::is_missing(i)) {
    x@dim[1] <- length(i)
    x@dimnames[1] <- list(rownames(x)[i])
  }
  if (!rlang::is_missing(j)) {
    x@dim[2] <- length(j)
    x@dimnames[2] <- list(colnames(x)[j])
  }
  x
}

# Helper function to split a selection into two components:
# First a subset, followed by a reorder.
# Set subset or reorder to rlang::missing_arg() if they are an identity transformation
# When we push down subset operations, we'll just push down the subset, and not the reorder
# Use instead of `selection_index()`. Convert these results into the equivalent from `selection_index()` by
# calling `unsplit_selection()`
split_selection_index <- function(selection, dim_len, dimnames) {
  selection <- selection_index(selection, dim_len, dimnames)
  if (rlang::is_missing(selection)) {
    return(list(subset=rlang::missing_arg(), reorder=rlang::missing_arg()))
  }
  res <- list(subset = sort(selection), reorder=rank(selection, ties.method="first"))
  if (length(res$subset) == dim_len) res$subset <- rlang::missing_arg()
  if (all(res$reorder == seq_along(res$reorder))) res$reorder <- rlang::missing_arg()
  return(res)
}

# Reverse split_selection
unsplit_selection <- function(selection) {
  if (rlang::is_missing(selection$subset)) return(rlang::maybe_missing(selection$reorder))
  if (rlang::is_missing(selection$reorder)) return(selection$subset)
  selection$subset[selection$reorder]
}

setMethod("[", "IterableMatrix", function(x, i, j, ...) {
  if (missing(x)) stop("x is missing in matrix selection")
  if (rlang::is_missing(i) && rlang::is_missing(j)) {
    return(x)
  }

  ret <- wrapMatrix("MatrixSubset", x)
  i <- selection_index(i, nrow(x), rownames(x))
  j <- selection_index(j, ncol(x), colnames(x))
 
  # Check for a trivial selection that doesn't introduce any subset
  if ((rlang::is_missing(i) || isTRUE(all.equal(i, seq_len(nrow(x))))) && 
      (rlang::is_missing(j) || isTRUE(all.equal(j, seq_len(ncol(x)))))) {
    return(x)
  }
  ret <- selection_fix_dims(ret, rlang::maybe_missing(i), rlang::maybe_missing(j))

  if (x@transpose) {
    tmp <- rlang::maybe_missing(i)
    i <- rlang::maybe_missing(j)
    j <- rlang::maybe_missing(tmp)
  }
  if (!rlang::is_missing(i)) {
    ret@row_selection <- i
    ret@zero_dims[1] <- length(i) == 0
  }
  if (!rlang::is_missing(j)) {
    ret@col_selection <- j
    ret@zero_dims[2] <- length(j) == 0
  }
  ret
})

# Simulate assigning to a subset of the matrix.
# We concatenate the un-modified matrix subsets with the new values,
# then reorder rows/columns appropriately
setMethod("[<-", "IterableMatrix", function(x, i, j, ..., value) {
  # Do type conversions if needed
  if (is.matrix(value)) value <- as(value, "dgCMatrix")
  if (is(value, "dgCMatrix")) {
    if (x@transpose) {
      value <- t(as(t(value), "IterableMatrix"))
    } else {
      value <- as(value, "IterableMatrix")
    }
    if (matrix_type(value) != matrix_type(x)) {
      rlang::warn(c(
        "Converting input matrix type to match destination",
        sprintf("input type: %s", matrix_type(value)),
        sprintf("destination type: %s", matrix_type(x))
      ))
      value <- convert_matrix_type(value, matrix_type(x))
    }
  }

  rownames(value) <- if (rlang::is_missing(i)) rownames(x) else rownames(x)[i]
  colnames(value) <- if (rlang::is_missing(j)) colnames(x) else colnames(x)[j]
  
  if (!rlang::is_missing(i)) {
    i <- selection_index(i, nrow(x), rownames(x))
    ni <- if (length(i) > 0) seq_len(nrow(x))[-i] else seq_len(nrow(x))
    
    x_i <- x[i,]
    x_ni <- x[ni,]
    
    if (rlang::is_missing(j)) {
      if (any(dim(x_i) != dim(value))) {
        stop("Mismatched dimensions in assignment to subset")
      }
      x_i <- value
    } else {
      x_i[,j] <- value
    }
    rownames(x_i) <- rownames(x)[i]
    x <- rbind(x_i, x_ni)[order(c(i, ni)),]
  } else if(!rlang::is_missing(j)) {
    j <- selection_index(j, ncol(x), colnames(x))
    nj <- if (length(j) > 0) seq_len(ncol(x))[-j] else seq_len(ncol(x))
    
    x_j <- x[,j]
    x_nj <- x[,nj]
    if (any(dim(x_j) != dim(value))) {
      stop("Mismatched dimensions in assignment to subset")
    }
    x_j <- value
    colnames(x_j) <- colnames(x)[j]
    x <- cbind(x_j, x_nj)[,order(c(j, nj))]
  } else {
    if (any(dim(x) != dim(value))) {
      stop("Mismatched dimensions in assignment to subset")
    }
    x <- value
  }
  
  return(x)
})

setMethod("[", "MatrixSubset", function(x, i, j, ...) {
  if (missing(x)) stop("x is missing in matrix selection")
  
  if (x@transpose) {
    return(t(t(x)[rlang::maybe_missing(j), rlang::maybe_missing(i)]))
  }

  i <- selection_index(i, nrow(x), rownames(x))
  j <- selection_index(j, ncol(x), colnames(x))
  x <- selection_fix_dims(x, rlang::maybe_missing(i), rlang::maybe_missing(j))

  if (!rlang::is_missing(i)) {
    if (length(x@row_selection) == 0 && !x@zero_dims[1]) {
      x@row_selection <- i
    } else {
      x@row_selection <- x@row_selection[i]
    }
    x@zero_dims[1] <- length(i) == 0
  }
  if (!rlang::is_missing(j)) {
    if (length(x@col_selection) == 0 && !x@zero_dims[2]) {
      x@col_selection <- j
    } else {
      x@col_selection <- x@col_selection[j]
    }
    x@zero_dims[2] <- length(j) == 0
  }
  # Apply the consolidated selection to the inner matrix in case there is some
  # new subsetting that should be pushed down further
  i <- if(length(x@row_selection) == 0 && !x@zero_dims[1]) rlang::missing_arg() else x@row_selection
  j <- if(length(x@col_selection) == 0 && !x@zero_dims[2]) rlang::missing_arg() else x@col_selection
  res <- x@matrix[rlang::maybe_missing(i), rlang::maybe_missing(j)]
  # If we are the final operation (signified by having non-null dimnames), make sure
  # those dimnames are preserved as the inner matrix becomes the new outer matrix.
  # If we are an inner operation (signified by having null dimnames), then setting
  # things here will cause C++ to no longer find the correct dimnames. (See issue #97)
  if (!is.null(rownames(x))) rownames(res) <- rownames(x)
  if (!is.null(colnames(x))) colnames(res) <- colnames(x)
  res
})


setMethod("iterate_matrix", "MatrixSubset", function(x) {
  assert_true(matrix_type(x) %in% c("uint32_t", "float", "double"))
  
  iter_row_function <- get(sprintf("iterate_matrix_row_select_%s_cpp", matrix_type(x)))
  iter_col_function <- get(sprintf("iterate_matrix_col_select_%s_cpp", matrix_type(x)))

  ret <- iterate_matrix(x@matrix)

  if (length(x@row_selection) != 0 || x@zero_dims[1]) ret <- iter_row_function(ret, x@row_selection - 1L)
  if (length(x@col_selection) != 0 || x@zero_dims[2]) ret <- iter_col_function(ret, x@col_selection - 1L)
  
  ret
})

setMethod("short_description", "MatrixSubset", function(x) {
  if (x@transpose) {
    rows <- if (x@zero_dims[2]) "none" else x@col_selection
    cols <- if (x@zero_dims[1]) "none" else x@row_selection
  } else {
    rows <- if (x@zero_dims[1]) "none" else x@row_selection
    cols <- if (x@zero_dims[2]) "none" else x@col_selection
  }
  c(
    short_description(x@matrix),
    sprintf(
      "Select rows: %s and cols: %s",
      pretty_print_vector(rows, empty = "all"),
      pretty_print_vector(cols, empty = "all")
    )
  )
})

# Renaming dims
setClass("RenameDims",
  contains = "IterableMatrix",
  slots = c(
    matrix = "IterableMatrix"
  )
)
setMethod("matrix_type", "RenameDims", function(x) matrix_type(x@matrix))
setMethod("iterate_matrix", "RenameDims", function(x) {
  if (x@transpose) {
    return(iterate_matrix(t(x)))
  }
  assert_true(matrix_type(x) %in% c("uint32_t", "float", "double"))
  
  iter_function <- get(sprintf("iterate_matrix_rename_dims_%s_cpp", matrix_type(x)))
  row_names <- if (is.null(rownames(x))) character(0) else rownames(x)
  col_names <- if (is.null(colnames(x))) character(0) else colnames(x)

  iter_function(iterate_matrix(x@matrix), row_names, col_names, is.null(rownames(x)), is.null(colnames(x)))
})

setMethod("[", "RenameDims", function(x, i, j, ...) {
  if (missing(x)) stop("x is missing in matrix selection")

  i <- selection_index(i, nrow(x), rownames(x))
  j <- selection_index(j, ncol(x), colnames(x))
  x <- selection_fix_dims(x, rlang::maybe_missing(i), rlang::maybe_missing(j))

  x@matrix <- x@matrix[rlang::maybe_missing(i), rlang::maybe_missing(j)]

  if (x@transpose) {
    tmp <- rlang::maybe_missing(i)
    i <- rlang::maybe_missing(j)
    j <- rlang::maybe_missing(tmp)
  }
  x
})
setMethod("short_description", "RenameDims", function(x) {
  c(
    short_description(x@matrix),
    sprintf("Reset dimnames")
  )
})
setMethod("dimnames<-", signature(x = "IterableMatrix", value = "list"), function(x, value) {
  if (identical(dimnames(x), value)) return(x)
  d <- dim(x)
  has_error <- FALSE
  if (!is.list(value) || length(value) != 2) has_error <- TRUE
  if (!is.null(value[[1]]) && length(value[[1]]) != d[1]) has_error <- TRUE
  if (!is.null(value[[2]]) && length(value[[2]]) != d[2]) has_error <- TRUE
  if (has_error) stop("Invalid dimnames supplied")

  if (!is(x, "RenameDims")) {
    x <- wrapMatrix("RenameDims", x)
  }
  if (!is.null(value[[1]])) {
    x@dimnames[[1]] <- as.character(value[[1]])
  } else {
    x@dimnames[1] <- list(NULL)
  }
  if (!is.null(value[[2]])) {
    x@dimnames[[2]] <- as.character(value[[2]])
  } else {
    x@dimnames[2] <- list(NULL)
  }
  x
})
setMethod("dimnames<-", signature(x = "IterableMatrix", value = "NULL"), function(x, value) {
  if (identical(dimnames(x), value)) return(x)
  if (!is(x, "RenameDims")) {
    x <- wrapMatrix("RenameDims", x)
  }
  x@dimnames <- list(NULL, NULL)
  x
})

# Concatenating matrices by row or by column

#' Helper function for rbind/cbind concatenating dimnames
#' @keywords internal
concat_dimnames <- function(x, y, len_x, len_y, warning_prefix, dim_type) {
  if (is.null(x) && is.null(y)) {
    return(NULL)
  }
  if (!is.null(x) && !is.null(y)) {
    return(c(x, y))
  }
  warning(sprintf(
    "%s: %s names present on some but not all matrices. Setting missing names to \"\"",
    warning_prefix, dim_type
  ), call. = FALSE)
  if (is.null(x)) x <- rep_len("", len_x)
  if (is.null(y)) y <- rep_len("", len_y)
  c(x, y)
}

#' Helper function for rbind/cbind merging dimnames
#' @keywords internal
merge_dimnames <- function(x, y, warning_prefix, dim_type) {
  if (!is.null(x) && !is.null(y) && !all(x == y)) {
    warning(sprintf("%s: %s names are mismatched. Setting names to match first matrix", warning_prefix, dim_type), call. = FALSE)
  }
  if (!is.null(x)) {
    return(x)
  }
  if (!is.null(y)) {
    return(y)
  }
  return(NULL)
}

setClass("RowBindMatrices",
  contains = "IterableMatrix",
  slots = c(
    matrix_list = "list",
    threads = "integer"
  ),
  prototype = list(
    matrix_list = list(),
    threads = 0L
  )
)
setMethod("matrix_type", signature(x = "RowBindMatrices"), function(x) matrix_type(x@matrix_list[[1]]))

setMethod("iterate_matrix", "RowBindMatrices", function(x) {
  iter_function <- get(sprintf("iterate_matrix_row_bind_%s_cpp", matrix_type(x)))
  iterators <- lapply(x@matrix_list, iterate_matrix)
  iter_function(iterators, x@threads)
})

setMethod("matrix_inputs", "RowBindMatrices", function(x) x@matrix_list)
setMethod("matrix_inputs<-", "RowBindMatrices", function(x, ..., value) {
  assert_is(value, "list")
  assert_len(value, length(x@matrix_list))
  for (v in value) assert_is(v, "IterableMatrix")
  x@matrix_list <- value
  return(x)
})

setMethod("short_description", "RowBindMatrices", function(x) {
  sprintf(
    "Concatenate %s of %d matrix objects with classes%s (threads=%d)",
    ifelse(x@transpose, "cols", "rows"),
    length(x@matrix_list),
    pretty_print_vector(vapply(x@matrix_list, class, character(1)), prefix = ": ", max_len = 3),
    x@threads
  )
})

setMethod("rbind2", signature(x = "IterableMatrix", y = "IterableMatrix"), function(x, y, ...) {
  if (x@transpose != y@transpose) stop("Cannot merge matrices with different internal transpose states.\nPlease use transpose_storage_order().")
  if (matrix_type(x) != matrix_type(y)) {
    # error out if matrix type x or y are not "double", "float", or "uint32_t"
    if(!(matrix_type(x) %in% c("uint32_t", "float", "double"))) {
      stop("rbind2(): Cannot merge matrices with different unspported data types. Please use convert_matrix_type().")
    }
    if(!(matrix_type(y) %in% c("uint32_t", "float", "double"))) {
      stop("rbind2(): Cannot merge matrices with different unspported data types. Please use convert_matrix_type().")
    }
    # upcast if mismatching types
    rlang::warn(
      sprintf(
        paste0("rbind2(): Mismatching matrix types (%s vs. %s). Upcasting to double"),
        matrix_type(x), matrix_type(y)
      )
    )
    x <- convert_matrix_type(x, "double")
    y <- convert_matrix_type(y, "double")
  }
  if (x@transpose) {
    return(t(cbind2(t(x), t(y))))
  }

  if (ncol(x) != ncol(y)) stop("Error in rbind: matrices must have equal number of columns")
  if (nrow(x) == 0) return(y)
  if (nrow(y) == 0) return(x)

  # Handle dimnames
  col_names <- merge_dimnames(colnames(x), colnames(y), "rbind", "column")
  row_names <- concat_dimnames(rownames(x), rownames(y), nrow(x), nrow(y), "rbind", "row")
  if (matrix_is_transform(x) && !is(x, "RenameDims")) x@dimnames <- list(NULL, NULL)
  if (matrix_is_transform(y) && !is(y, "RenameDims")) y@dimnames <- list(NULL, NULL)

  matrix_list <- list()
  if (is(x, "RowBindMatrices")) {
    matrix_list <- c(matrix_list, x@matrix_list)
  } else {
    matrix_list <- c(matrix_list, x)
  }
  if (is(y, "RowBindMatrices")) {
    matrix_list <- c(matrix_list, y@matrix_list)
  } else {
    matrix_list <- c(matrix_list, y)
  }

  new("RowBindMatrices", matrix_list = matrix_list, dim = c(nrow(x) + nrow(y), ncol(x)), dimnames = list(row_names, col_names), transpose = FALSE)
})

# Helper methods for unusual rbind2 calls
setMethod("rbind2", signature(x = "IterableMatrix", y = "missing"), function(x, y, ...) x)
setMethod("rbind2", signature(x = "IterableMatrix", y = "dgCMatrix"), function(x, y, ...) rbind2(x, as(y, "IterableMatrix")))
setMethod("rbind2", signature(x = "dgCMatrix", y = "IterableMatrix"), function(x, y, ...) rbind2(as(x, "IterableMatrix"), y))

#' Set matrix op thread count
#'
#' Set number of threads to use for sparse-dense multiply and matrix_stats.
#'
#' Only valid for concatenated matrices
#' @param mat IterableMatrix, product of rbind or cbind
#' @param threads Number of threads to use for execution
#' @keywords internal
set_threads <- function(mat, threads=0L) {
  assert_is(mat, c("RowBindMatrices", "ColBindMatrices"))
  assert_is_wholenumber(threads)
  assert_true(threads >= 0)
  mat@threads <- as.integer(threads)
  mat
}

setClass("ColBindMatrices",
  contains = "IterableMatrix",
  slots = c(
    matrix_list = "list",
    threads = "integer"
  ),
  prototype = list(
    matrix_list = list(),
    threads = 0L
  )
)
setMethod("matrix_type", signature(x = "ColBindMatrices"), function(x) matrix_type(x@matrix_list[[1]]))

setMethod("iterate_matrix", "ColBindMatrices", function(x) {
  iter_function <- get(sprintf("iterate_matrix_col_bind_%s_cpp", matrix_type(x)))
  iterators <- lapply(x@matrix_list, iterate_matrix)
  iter_function(iterators, x@threads)
})

setMethod("matrix_inputs", "ColBindMatrices", function(x) x@matrix_list)
setMethod("matrix_inputs<-", "ColBindMatrices", function(x, ..., value) {
  assert_is(value, "list")
  assert_len(value, length(x@matrix_list))
  for (v in value) assert_is(v, "IterableMatrix")
  x@matrix_list <- value
  return(x)
})


setMethod("short_description", "ColBindMatrices", function(x) {
  sprintf(
    "Concatenate %s of %d matrix objects with classes%s (threads=%d)",
    ifelse(x@transpose, "rows", "cols"),
    length(x@matrix_list),
    pretty_print_vector(vapply(x@matrix_list, class, character(1)), prefix = ": ", max_len = 3),
    x@threads
  )
})

setMethod("cbind2", signature(x = "IterableMatrix", y = "IterableMatrix"), function(x, y, ...) {
  if (x@transpose != y@transpose) stop("Cannot merge matrices with different internal transpose states.\nPlease use transpose_storage_order().")
  if (matrix_type(x) != matrix_type(y)) {
    # error out if matrix type x or y are not "double", "float", or "uint32_t"
    if(!(matrix_type(x) %in% c("uint32_t", "float", "double"))) {
      stop("cbind2(): Cannot merge matrices with different unsupported data types. Please use convert_matrix_type().")
    }
    if(!(matrix_type(y) %in% c("uint32_t", "float", "double"))) {
      stop("cbind2(): Cannot merge matrices with different unsupported data types. Please use convert_matrix_type().")
    }
    # upcast if mismatching types
    rlang::warn(
      sprintf(
        paste0("cbind2(): Mismatching matrix types (%s vs. %s). Upcasting to double"),
        matrix_type(x), matrix_type(y)
      )
    )
    x <- convert_matrix_type(x, "double")
    y <- convert_matrix_type(y, "double")
  }
  if (x@transpose) {
    return(t(rbind2(t(x), t(y))))
  }

  if (nrow(x) != nrow(y)) stop("Error in cbind: matrices must have equal number of columns")
  if (ncol(x) == 0) return(y)
  if (ncol(y) == 0) return(x)
  # Handle dimnames
  row_names <- merge_dimnames(rownames(x), rownames(y), "cbind", "row")
  col_names <- concat_dimnames(colnames(x), colnames(y), ncol(x), ncol(y), "cbind", "column")
  if (matrix_is_transform(x) && !is(x, "RenameDims")) x@dimnames <- list(NULL, NULL)
  if (matrix_is_transform(y) && !is(y, "RenameDims")) y@dimnames <- list(NULL, NULL)

  matrix_list <- list()
  if (is(x, "ColBindMatrices")) {
    matrix_list <- c(matrix_list, x@matrix_list)
  } else {
    matrix_list <- c(matrix_list, x)
  }
  if (is(y, "ColBindMatrices")) {
    matrix_list <- c(matrix_list, y@matrix_list)
  } else {
    matrix_list <- c(matrix_list, y)
  }

  new("ColBindMatrices", matrix_list = matrix_list, dim = c(nrow(x), ncol(x) + ncol(y)), dimnames = list(row_names, col_names), transpose = FALSE)
})

# Helper methods for unusual cbind2 calls
setMethod("cbind2", signature(x = "IterableMatrix", y = "missing"), function(x, y, ...) x)
setMethod("cbind2", signature(x = "IterableMatrix", y = "dgCMatrix"), function(x, y, ...) cbind2(x, as(y, "IterableMatrix")))
setMethod("cbind2", signature(x = "dgCMatrix", y = "IterableMatrix"), function(x, y, ...) cbind2(as(x, "IterableMatrix"), y))

# Row bind needs specialization because there's not a default row-seek operation
setMethod("[", "RowBindMatrices", function(x, i, j, ...) {
  if (missing(x)) stop("x is missing in matrix selection")
  # Handle transpose via recursive call
  if (x@transpose) {
    return(t(t(x)[rlang::maybe_missing(j), rlang::maybe_missing(i)]))
  }

  i <- split_selection_index(i, nrow(x), rownames(x))
  j <- split_selection_index(j, ncol(x), colnames(x))


  # If we're just reordering rows/cols, do a standard matrix selection
  if (rlang::is_missing(i$subset) && rlang::is_missing(j$subset)) {
    return(callNextMethod(x, unsplit_selection(i), unsplit_selection(j)))
  }

  # if the length of our row selection is 0, do a standard matrix selection
  if (!rlang::is_missing(i$subset) && length(i$subset) == 0) {
    return(callNextMethod(x, unsplit_selection(i), unsplit_selection(j)))
  }

  x <- selection_fix_dims(x, rlang::maybe_missing(i$subset), rlang::maybe_missing(j$subset))

  # Calculate helper variables to extract and transform the relevant parts of i$subset for each sub-matrix
  if (!rlang::is_missing(i$subset)) {
    rows <- vapply(x@matrix_list, nrow, integer(1))
    local_i_offset <- cumsum(c(0, rows))
    # Find the range of indices in i$subset that correspond to each matrix in x@matrix_list
    local_i_range <- findInterval(local_i_offset, i$subset)
  }

  new_mats <- list()
  for (k in seq_along(x@matrix_list)) {
    mat <- x@matrix_list[[k]]
    if (!rlang::is_missing(i$subset)) {
      if (local_i_range[k] == local_i_range[k+1]) {
        local_i <- integer(0)
      } else {
        local_i <- i$subset[(local_i_range[k]+1):(local_i_range[k+1])] - local_i_offset[k]
      }
      mat <- mat[local_i,]
    }
    if (!rlang::is_missing(j$subset)) {
      # Only pass through the subset operation to a lower-level, not the shuffle
      mat <- mat[,j$subset]
    }
    if (nrow(mat) > 0) {
      new_mats <- c(new_mats, mat)
    }
  }
  if (length(new_mats) > 1) {
    x@matrix_list <- new_mats
  } else if(length(new_mats) == 1) {
    # Only set dimnames here if we know non-null ones, otherwise it might
    # prevent name propagation at the C++ level. (See "[" for MatrixSubset)
    if (!is.null(rownames(x))) rownames(new_mats[[1]]) <- rownames(x)
    if (!is.null(colnames(x))) colnames(new_mats[[1]]) <- colnames(x)
    x <- new_mats[[1]]
  } else {
    stop("Subset RowBindMatrix error: got 0-length matrix_list after subsetting (please report this BPCells bug)")
  }
  if (!rlang::is_missing(i$reorder)) {
    x <- x[i$reorder,]
  }
  if (!rlang::is_missing(j$reorder)) {
    x <- x[,j$reorder]
  }
  return(x)
})

setMethod("[", "ColBindMatrices", function(x, i, j, ...) {
  if (missing(x)) stop("x is missing in matrix selection")
  # Handle transpose via recursive call
  if (x@transpose) {
    return(t(t(x)[rlang::maybe_missing(j), rlang::maybe_missing(i)]))
  }

  i <- split_selection_index(i, nrow(x), rownames(x))
  j <- split_selection_index(j, ncol(x), colnames(x))


  # If we're just reordering rows/cols, do a standard matrix selection
  if (rlang::is_missing(i$subset) && rlang::is_missing(j$subset)) {
    return(callNextMethod(x, unsplit_selection(i), unsplit_selection(j)))
  }

  # if the length of our col selection is 0, do a standard matrix selection
  if (!rlang::is_missing(j$subset) && length(j$subset) == 0) {
    return(callNextMethod(x, unsplit_selection(i), unsplit_selection(j)))
  }

  x <- selection_fix_dims(x, rlang::maybe_missing(i$subset), rlang::maybe_missing(j$subset))

  # Calculate helper variables to extract and transform the relevant parts of j$subset for each sub-matrix
  if (!rlang::is_missing(j$subset)) {
    cols <- vapply(x@matrix_list, ncol, integer(1))
    local_j_offset <- cumsum(c(0, cols))
    # Find the range of indices in j$subset that correspond to each matrix in x@matrix_list
    local_j_range <- findInterval(local_j_offset, j$subset)
  }

  new_mats <- list()
  for (k in seq_along(x@matrix_list)) {
    mat <- x@matrix_list[[k]]
    if (!rlang::is_missing(j$subset)) {
      if (local_j_range[k] == local_j_range[k+1]) {
        local_j <- integer(0)
      } else {
        local_j <- j$subset[(local_j_range[k]+1):(local_j_range[k+1])] - local_j_offset[k]
      }
      mat <- mat[,local_j]
    }
    if (!rlang::is_missing(i$subset)) {
      # Only pass through the subset operation to a lower-level, not the shuffle
      mat <- mat[i$subset,]
    }
    if (ncol(mat) > 0) {
      new_mats <- c(new_mats, mat)
    }
  }
  if (length(new_mats) > 1) {
    x@matrix_list <- new_mats
  } else if(length(new_mats) == 1) {
    # Only set dimnames here if we know non-null ones, otherwise it might
    # prevent name propagation at the C++ level. (See "[" for MatrixSubset)
    if (!is.null(rownames(x))) rownames(new_mats[[1]]) <- rownames(x)
    if (!is.null(colnames(x))) colnames(new_mats[[1]]) <- colnames(x)
    x <- new_mats[[1]]
  } else {
    stop("Subset ColBindMatrix error: got 0-length matrix_list after subsetting (please report this BPCells bug)")
  }
  if (!rlang::is_missing(i$reorder)) {
    x <- x[i$reorder,]
  }
  if (!rlang::is_missing(j$reorder)) {
    x <- x[,j$reorder]
  }
  return(x)
})
#' Prepare a matrix for multi-threaded operation
#' 
#' Transforms a matrix such that `matrix_stats` or matrix multiplies with
#' a vector/dense matrix will be evaluated in parallel. This only speeds up
#' those specific operations, not reading or writing the matrix in general.
#' The parallelism is not guaranteed to work if additional operations are
#' applied after the parallel split.
#'
#' @param mat IterableMatrix
#' @param threads Number of execution threads
#' @param chunks Number of chunks to use (>= threads)
#' @return IterableMatrix which will perform certain operations in parallel
#' @keywords internal
parallel_split <- function(mat, threads, chunks=threads) {
  assert_is(mat, "IterableMatrix")
  assert_is_wholenumber(threads)
  assert_is_wholenumber(chunks)
  assert_true(chunks >= threads)

  if (threads <= 1L) {
    return(mat)
  }

  if (mat@transpose) {
    return(t(parallel_split(t(mat), threads, chunks)))
  }

  start_col <- 1
  mats <- list()
  for (i in seq_len(chunks)) {
    col_count <- (ncol(mat) - start_col + 1) %/% (chunks - i + 1)
    indices <- seq(start_col, start_col + col_count - 1)
    start_col <- start_col + col_count
    subset <- mat[,indices]
    mats <- c(mats, list(mat[,indices]))
  }
  ret <- do.call(cbind, mats)
  ret@threads <- as.integer(threads)
  return(ret)
}

# Packed integer matrix
setClass("PackedMatrixMemBase",
  contains = "IterableMatrix",
  slots = c(
    # Leave out val storage since it's datatype-dependent
    index_data = "integer",
    index_starts = "integer",
    index_idx = "integer",
    index_idx_offsets = "numeric",
    idxptr = "numeric",
    version = "character"
  ),
  prototype = list(
    # Leave out val storage since it's datatype-dependent
    index_data = integer(0),
    index_starts = integer(0),
    index_idx = integer(0),
    index_idx_offsets = numeric(0),
    idxptr = numeric(0),
    version = character(0)
  )
)
setMethod("short_description", "PackedMatrixMemBase", function(x) {
  "Load compressed matrix from memory"
})
setMethod("matrix_inputs", "PackedMatrixMemBase", function(x) list())

setClass("PackedMatrixMem_uint32_t",
  contains = "PackedMatrixMemBase",
  slots = c(
    val_data = "integer",
    val_idx = "integer",
    val_idx_offsets = "numeric"
  ),
  prototype = list(
    val_data = integer(0),
    val_idx = integer(0),
    val_idx_offsets = numeric(0)
  )
)
setMethod("matrix_type", "PackedMatrixMem_uint32_t", function(x) "uint32_t")
setMethod("iterate_matrix", "PackedMatrixMem_uint32_t", function(x) {
  if (x@transpose) x <- t(x)
  x@dimnames <- denormalize_dimnames(x@dimnames)
  iterate_packed_matrix_mem_uint32_t_cpp(x, x@dimnames[[1]], x@dimnames[[2]], nrow(x))
})

setClass("PackedMatrixMem_float",
  contains = "PackedMatrixMemBase",
  slots = c(val = "integer"),
  prototype = list(val = integer(0))
)
setMethod("matrix_type", "PackedMatrixMem_float", function(x) "float")
setMethod("iterate_matrix", "PackedMatrixMem_float", function(x) {
  if (x@transpose) x <- t(x)
  x@dimnames <- denormalize_dimnames(x@dimnames)
  iterate_packed_matrix_mem_float_cpp(x, x@dimnames[[1]], x@dimnames[[2]], nrow(x))
})

setClass("PackedMatrixMem_double",
  contains = "PackedMatrixMemBase",
  slots = c(val = "numeric"),
  prototype = list(val = numeric(0))
)
setMethod("matrix_type", "PackedMatrixMem_double", function(x) "double")
setMethod("iterate_matrix", "PackedMatrixMem_double", function(x) {
  if (x@transpose) x <- t(x)
  x@dimnames <- denormalize_dimnames(x@dimnames)
  iterate_packed_matrix_mem_double_cpp(x, x@dimnames[[1]], x@dimnames[[2]], nrow(x))
})

setClass("UnpackedMatrixMemBase",
  contains = "IterableMatrix",
  slots = c(
    # Leave out val storage since it's data-type dependent
    index = "integer",
    idxptr = "numeric",
    version = "character"
  ),
  prototype = list(
    index = integer(0),
    idxptr = numeric(0),
    version = character(0)
  )
)
setMethod("short_description", "UnpackedMatrixMemBase", function(x) {
  "Load uncompressed matrix from memory"
})
setMethod("matrix_inputs", "UnpackedMatrixMemBase", function(x) list())

setClass("UnpackedMatrixMem_uint32_t",
  contains = "UnpackedMatrixMemBase",
  slots = c(val = "integer"),
  prototype = list(val = integer())
)
setMethod("matrix_type", "UnpackedMatrixMem_uint32_t", function(x) "uint32_t")
setMethod("iterate_matrix", "UnpackedMatrixMem_uint32_t", function(x) {
  if (x@transpose) x <- t(x)
  x@dimnames <- denormalize_dimnames(x@dimnames)
  iterate_unpacked_matrix_mem_uint32_t_cpp(x, x@dimnames[[1]], x@dimnames[[2]], nrow(x))
})

setClass("UnpackedMatrixMem_float",
  contains = "UnpackedMatrixMemBase",
  slots = c(val = "integer"),
  prototype = list(val = integer(0))
)
setMethod("matrix_type", "UnpackedMatrixMem_float", function(x) "float")
setMethod("iterate_matrix", "UnpackedMatrixMem_float", function(x) {
  if (x@transpose) x <- t(x)
  x@dimnames <- denormalize_dimnames(x@dimnames)
  iterate_unpacked_matrix_mem_float_cpp(x, x@dimnames[[1]], x@dimnames[[2]], nrow(x))
})

setClass("UnpackedMatrixMem_double",
  contains = "UnpackedMatrixMemBase",
  slots = c(val = "numeric"),
  prototype = list(val = numeric(0))
)
setMethod("matrix_type", "UnpackedMatrixMem_double", function(x) "double")
setMethod("iterate_matrix", "UnpackedMatrixMem_double", function(x) {
  if (x@transpose) x <- t(x)
  x@dimnames <- denormalize_dimnames(x@dimnames)
  iterate_unpacked_matrix_mem_double_cpp(x, x@dimnames[[1]], x@dimnames[[2]], nrow(x))
})


#' Transpose the storage order for a matrix
#' @param matrix Input matrix
#' @param outdir Directory to store the output
#' @param tmpdir Temporary directory to use for intermediate storage
#' @param load_bytes The minimum contiguous load size during the merge sort passes
#' @param sort_bytes The amount of memory to allocate for re-sorting chunks of entries
#' @details This re-sorts the entries of a matrix to change the storage order
#' from row-major to col-major. For large matrices, this can be slow -- around 2
#' minutes to transpose a 500k cell RNA-seq matrix The default load_bytes (4MiB)
#' and sort_bytes (1GiB) parameters allow ~85GB of data to be sorted with two
#' passes through the data, and ~7.3TB of data to be sorted in three passes
#' through the data.
#' @return MatrixDir object with a copy of the input matrix, but the storage order flipped
#' @examples
#' mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
#' rownames(mat) <- paste0("gene", seq_len(10))
#' colnames(mat) <- paste0("cell", seq_len(5))
#' mat <- mat %>% as("dgCMatrix") %>% as("IterableMatrix")
#' mat
#' 
#' ## A regular transpose operation switches a user's rows and cols 
#' t(mat)
#' 
#' ## Running `transpose_storage_order()` instead changes whether the storage is in row-major or col-major,
#' ## but does not switch the rows and cols
#' transpose_storage_order(mat)
#' @export
transpose_storage_order <- function(matrix, outdir = tempfile("transpose"), tmpdir = tempdir(), load_bytes = 4194304L, sort_bytes = 1073741824L) {
  assert_true(matrix_type(matrix) %in% c("uint32_t", "float", "double"))

  write_function <- get(sprintf("write_matrix_transpose_%s_cpp", matrix_type(matrix)))

  outdir <- normalizePath(outdir, mustWork = FALSE)
  tmpdir <- normalizePath(tmpdir, mustWork = FALSE)
  tmpdir <- tempfile("transpose_tmp", tmpdir = tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE, expand = FALSE))

  it <- iterate_matrix(matrix)
  write_function(it, outdir, tmpdir, load_bytes, sort_bytes, !matrix@transpose)

  open_matrix_dir(outdir)
}

#' Read/write sparse matrices
#'
#' BPCells matrices are stored in sparse format, meaning only the non-zero entries
#' are stored. Matrices can store integer counts data or decimal numbers (float or double).
#' See details for more information.
#'
#' ### Storage locations
#' Matrices can be stored in a directory on disk, in memory, or in an HDF5 file.
#' Saving in a directory on disk is a good default for local analysis, as it provides
#' the best I/O performance and lowest memory usage. The HDF5 format
#' allows saving within existing hdf5 files to group data together, and the in
#' memory format provides the fastest performance in the event memory usage is
#' unimportant.
#'
#' ### Bitpacking Compression
#' For typical RNA counts matrices holding integer counts, this bitpacking
#' compression will result in 6-8x less space than an R dgCMatrix, and 4-6x
#' smaller than a scipy csc_matrix. The compression will be more effective when
#' the count values in the matrix are small, and when the rows of the matrix are
#' sorted by rowMeans. In tests on RNA-seq data optimal ordering could save up
#' to 40% of storage space. On non-integer data only the row indices are
#' compressed, not the values themselves so space savings will be smaller.
#'
#' For non-integer data matrices, bitpacking compression is much less effective,
#' as it can only be applied to the indexes of each entry but not the values.
#' There will still be some space savings, but far less than for counts matrices.
#'
#' @param matrix Input matrix, either IterableMatrix or dgCMatrix
#' @param compress Whether or not to compress the data.
#' @return BPCells matrix object
#' @examples
#' ## Create temporary directory to keep demo matrix
#' data_dir <- file.path(tempdir(), "mat")
#' if (dir.exists(data_dir)) unlink(data_dir, recursive = TRUE)
#' dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
#' 
#' mat <- get_demo_mat()
#' mat
#' 
#' #######################################################################
#' ## write_matrix_memory() example
#' #######################################################################
#' mat_memory <- write_matrix_memory(mat)
#' mat_memory
#' 
#' 
#' @rdname matrix_io
#' @export
write_matrix_memory <- function(mat, compress = TRUE) {
  assert_is(mat, c("IterableMatrix", "dgCMatrix"))
  if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")

  assert_true(matrix_type(mat) %in% c("uint32_t", "float", "double"))
  if (compress && matrix_type(mat) != "uint32_t") {
    rlang::inform(c(
      "Warning: Matrix compression performs poorly with non-integers.",
      "Consider calling convert_matrix_type if a compressed integer matrix is intended."
    ), .frequency = "regularly", .frequency_id = "matrix_compress_non_integer")
  }

  write_function <- get(sprintf("write_%s_matrix_mem_%s_cpp", ifelse(compress, "packed", "unpacked"), matrix_type(mat)))
  class <- sprintf("%sMatrixMem_%s", ifelse(compress, "Packed", "Unpacked"), matrix_type(mat))

  it <- iterate_matrix(mat)
  m <- write_function(it, mat@transpose)
  m[["dimnames"]] <- normalized_dimnames(m$row_names, m$col_names)
  m$dim <- m$shape
  m$transpose <- m$storage_order == "row"
  m$row_names <- NULL
  m$col_names <- NULL
  m$shape <- NULL
  m$storage_order <- NULL
  res <- do.call(new, c(class, m))
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
setMethod("matrix_inputs", "MatrixDir", function(x) list())

setMethod("iterate_matrix", "MatrixDir", function(x) {
  if (x@transpose) x <- t(x)
  x@dimnames <- denormalize_dimnames(x@dimnames)

  iter_function <- get(sprintf("iterate_%s_matrix_file_%s_cpp", ifelse(x@compressed, "packed", "unpacked"), matrix_type(x)))

  iter_function(x@dir, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]], nrow(x))
})

setMethod("short_description", "MatrixDir", function(x) {
  sprintf(
    "Load %s matrix from directory %s",
    if (x@compressed) "compressed" else "uncompressed",
    x@dir
  )
})


#' @rdname matrix_io
#' @param dir Directory to save the data into
#' @param buffer_size For performance tuning only. The number of items to be buffered
#' in memory before calling writes to disk.
#' @param overwrite If `TRUE`, write to a temp dir then overwrite existing data. Alternatively,
#'   pass a temp path as a string to customize the temp dir location.
#' @examples
#' #######################################################################
#' ## write_matrix_dir() example
#' #######################################################################
#' mat %>% write_matrix_dir(
#'  file.path(data_dir, "demo_mat"),
#'  overwrite = TRUE
#' )
#' 
#' 
#' @export
write_matrix_dir <- function(mat, dir, compress = TRUE, buffer_size = 8192L, overwrite = FALSE) {
  assert_is(mat, c("IterableMatrix", "dgCMatrix"))
  if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")

  assert_is(dir, "character")
  assert_is(compress, "logical")
  assert_is(buffer_size, "integer")
  assert_is(overwrite, c("logical", "character"))
  if (is(overwrite, "character")) {
    assert_true(dir.exists(overwrite))
    overwrite_path <- tempfile("overwrite", tmpdir=overwrite)
    overwrite <- TRUE
  } else if (overwrite) {
    overwrite_path <- tempfile("overwrite")
  }


  assert_true(matrix_type(mat) %in% c("uint32_t", "float", "double"))
  if (compress && matrix_type(mat) != "uint32_t") {
    rlang::inform(c(
      "Warning: Matrix compression performs poorly with non-integers.",
      "Consider calling convert_matrix_type if a compressed integer matrix is intended."
    ), .frequency = "regularly", .frequency_id = "matrix_compress_non_integer")
  }

  dir <- normalizePath(dir, mustWork = FALSE)
  
  did_tmp_copy <- FALSE
  if (overwrite && dir.exists(dir)) {
    mat <- write_matrix_dir(mat, overwrite_path, compress, buffer_size)
    did_tmp_copy <- TRUE
  }

  it <- iterate_matrix(mat)

  write_function <- get(sprintf("write_%s_matrix_file_%s_cpp", ifelse(compress, "packed", "unpacked"), matrix_type(mat)))
  write_function(it, dir, buffer_size, overwrite, mat@transpose)

  if (did_tmp_copy) {
    unlink(overwrite_path, recursive=TRUE)
  }
  open_matrix_dir(dir, buffer_size)
}

#' @rdname matrix_io
#' @examples
#' #######################################################################
#' ## open_matrix_dir() example
#' #######################################################################
#' mat <- open_matrix_dir(
#'  file.path(data_dir, "demo_mat")
#' )
#' mat
#' 
#' 
#' @export
open_matrix_dir <- function(dir, buffer_size = 8192L) {
  assert_is_file(dir)
  assert_is(buffer_size, "integer")

  dir <- normalizePath(dir, mustWork = FALSE)
  info <- dims_matrix_file_cpp(dir, buffer_size)
  new("MatrixDir",
    dir = dir, dim = info$dims, compressed = info$compressed, buffer_size = buffer_size,
    dimnames = normalized_dimnames(info$row_names, info$col_names), type = info$type,
    transpose = info$transpose
  )
}

setClass("EXPERIMENTAL_MatrixDirCompressedCol",
  contains = "IterableMatrix",
  slots = c(
    dir = "character",
    buffer_size = "integer"
  ),
  prototype = list(
    dir = character(0),
    buffer_size = integer(0)
  )
)
setMethod("matrix_type", "EXPERIMENTAL_MatrixDirCompressedCol", function(x) "uint32_t")
setMethod("matrix_inputs", "EXPERIMENTAL_MatrixDirCompressedCol", function(x) list())

setMethod("iterate_matrix", "EXPERIMENTAL_MatrixDirCompressedCol", function(x) {
  if (x@transpose) x <- t(x)

  EXPERIMENTAL_iterate_packed_sparse_column_matrix_file_uint32_t_cpp(
    x@dir, x@buffer_size
  )
})

setMethod("short_description", "EXPERIMENTAL_MatrixDirCompressedCol", function(x) {
  sprintf(
    "Load EXPERIMENTAL sparse column format matrix from directory %s",
    x@dir
  )
})


#' Write to experimental sparse-column format integer matrix
#'
#' The experimental sparse-column format is designed to handle storage of matrices
#' with many columns of all-zero, and less than 2^32-1 non-zero entries.
#'
#' @param dir Directory to save the data into
#' @param overwrite If `TRUE`, write to a temp dir then overwrite existing data. Alternatively,
#'   pass a temp path as a string to customize the temp dir location.
#' @keywords internal
EXPERIMENTAL_write_matrix_dir <- function(mat, dir, buffer_size = 8192L, overwrite = FALSE) {
  assert_is(mat, c("IterableMatrix", "dgCMatrix"))
  if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")

  assert_is(dir, "character")
  assert_is(buffer_size, "integer")
  assert_is(overwrite, c("logical", "character"))
  if (is(overwrite, "character")) {
    assert_true(dir.exists(overwrite))
    overwrite_path <- tempfile("overwrite", tmpdir=overwrite)
    overwrite <- TRUE
  } else if (overwrite) {
    overwrite_path <- tempfile("overwrite")
  }

  assert_true(matrix_type(mat) == "uint32_t")

  dir <- normalizePath(dir, mustWork = FALSE)
  
  did_tmp_copy <- FALSE
  if (overwrite && dir.exists(dir)) {
    mat <- EXPERIMENTAL_write_matrix_dir(mat, overwrite_path, buffer_size)
    did_tmp_copy <- TRUE
  }

  it <- iterate_matrix(mat)

  EXPERIMENTAL_write_packed_sparse_column_matrix_file_uint32_t_cpp(
    it, dir, buffer_size, overwrite, mat@transpose
  )

  if (did_tmp_copy) {
    unlink(overwrite_path, recursive=TRUE)
  }
  EXPERIMENTAL_open_matrix_dir(dir, buffer_size)
}

#' Open experimental sparse-column format integer matrix
#'
#' The experimental sparse-column format is designed to handle storage of matrices
#' with many columns of all-zero, and less than 2^32-1 non-zero entries.
#' @param dir Directory to load data from
#' @param buffer_size For performance tuning only. The number of items to be buffered
#' in memory before calling writes to disk.
#' @keywords internal
EXPERIMENTAL_open_matrix_dir <- function(dir, buffer_size = 8192L) {
  assert_is_file(dir)
  assert_is(buffer_size, "integer")

  dir <- normalizePath(dir, mustWork = FALSE)
  info <- EXPERIMENTAL_dims_packed_sparse_column_matrix_file_cpp(dir, buffer_size)
  new("EXPERIMENTAL_MatrixDirCompressedCol",
    dir = dir, dim = info$dims, buffer_size = buffer_size,
    dimnames = normalized_dimnames(info$row_names, info$col_names),
    transpose = info$transpose
  )
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
setMethod("matrix_inputs", "MatrixH5", function(x) list())

setMethod("iterate_matrix", "MatrixH5", function(x) {
  if (x@transpose) x <- t(x)
  x@dimnames <- denormalize_dimnames(x@dimnames)

  iter_function <- get(sprintf("iterate_%s_matrix_hdf5_%s_cpp", ifelse(x@compressed, "packed", "unpacked"), matrix_type(x)))

  iter_function(x@path, x@group, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]], nrow(x))
})

setMethod("short_description", "MatrixH5", function(x) {
  sprintf(
    "Load %s matrix in hdf5 file %s, group %s",
    if (x@compressed) "compressed" else "uncompressed",
    x@path,
    x@group
  )
})

#' @rdname matrix_io
#' @inheritParams write_fragments_hdf5
#' @examples
#' #######################################################################
#' ## write_matrix_hdf5() example
#' #######################################################################
#' mat %>% write_matrix_hdf5(path = file.path(data_dir, "demo_mat.h5"), group = "mat")
#' 
#' 
#' @export
write_matrix_hdf5 <- function(
    mat, 
    path, 
    group, 
    compress = TRUE, 
    buffer_size = 8192L, 
    chunk_size = 1024L, 
    overwrite = FALSE,
    gzip_level = 0L
) {
  assert_is(mat, c("IterableMatrix", "dgCMatrix"))
  if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")

  assert_is(path, "character")
  assert_is(group, "character")
  assert_is(compress, "logical")
  assert_is(buffer_size, "integer")
  assert_is(chunk_size, "integer")
  assert_is(gzip_level, "integer")
  assert_is(overwrite, c("logical", "character"))
  if (is(overwrite, "character")) {
    assert_true(dir.exists(overwrite))
    overwrite_path <- tempfile("overwrite", tmpdir=overwrite)
    overwrite <- TRUE
  } else if (overwrite) {
    overwrite_path <- tempfile("overwrite")
  }

  if (gzip_level != 0L && compress) {
     rlang::inform(c(
        "Warning: Mixing gzip compression (gzip_level > 0) with bitpacking compression (compress=TRUE) may be slower than bitpacking compression alone, with little space savings"
     ))
  }

  assert_true(matrix_type(mat) %in% c("uint32_t", "float", "double"))
  if (compress && matrix_type(mat) != "uint32_t") {
    rlang::inform(c(
      "Warning: Matrix compression performs poorly with non-integers.",
      "Consider calling convert_matrix_type if a compressed integer matrix is intended."
    ), .frequency = "regularly", .frequency_id = "matrix_compress_non_integer")
  }

  path <- normalizePath(path, mustWork = FALSE)
  did_tmp_copy <- FALSE
  if (overwrite && hdf5_group_exists_cpp(path, group)) {
    rlang::inform(c(
      "Warning: Overwriting an hdf5 dataset does not free old storage"
    ), .frequency = "regularly", .frequency_id = "hdf5_overwrite")
    did_tmp_copy <- TRUE
    mat <- write_matrix_dir(mat, overwrite_path, compress, buffer_size)
  }

  it <- iterate_matrix(mat)

  write_function <- get(sprintf("write_%s_matrix_hdf5_%s_cpp", ifelse(compress, "packed", "unpacked"), matrix_type(mat)))
  write_function(it, path, group, buffer_size, chunk_size, overwrite, mat@transpose, gzip_level)

  if (did_tmp_copy) {
    unlink(overwrite_path, recursive=TRUE)
  }
  open_matrix_hdf5(path, group, buffer_size)
}

#' @rdname matrix_io
#' @inheritParams open_fragments_hdf5
#' @examples
#' #######################################################################
#' ## open_matrix_hdf5() example
#' #######################################################################
#' mat_hdf5 <- open_matrix_hdf5(
#'  file.path(data_dir, "demo_mat.h5"),
#'  group = 'mat'
#' )
#' mat_hdf5
#' 
#' 
#' @export
open_matrix_hdf5 <- function(path, group, buffer_size = 16384L) {
  assert_is_file(path)
  assert_is(group, "character")
  assert_is(buffer_size, "integer")

  path <- normalizePath(path, mustWork = FALSE)
  info <- dims_matrix_hdf5_cpp(path, group, buffer_size)
  new("MatrixH5",
    path = path, group = group, dim = info$dims, compressed = info$compressed, buffer_size = buffer_size,
    dimnames = normalized_dimnames(info$row_names, info$col_names), type = info$type,
    transpose = info$transpose
  )
}

setClass("10xMatrixH5",
  contains = "IterableMatrix",
  slots = c(
    path = "character",
    group = "character",
    type = "character",
    buffer_size = "integer"
  ),
  prototype = list(
    path = character(0),
    group = "matrix",
    type = "uint32_t",
    buffer_size = integer(0)
  )
)
setMethod("matrix_type", "10xMatrixH5", function(x) x@type)
setMethod("matrix_inputs", "10xMatrixH5", function(x) list())
setMethod("iterate_matrix", "10xMatrixH5", function(x) {
  if (x@transpose) x <- t(x)
  x@dimnames <- denormalize_dimnames(x@dimnames)
  iterate_matrix_10x_hdf5_cpp(x@path, x@group, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]])
})
setMethod("short_description", "10xMatrixH5", function(x) {
  sprintf("10x HDF5 feature matrix in file %s", x@path)
})

#' Read/write a 10x feature matrix
#'
#' @inheritParams open_matrix_hdf5
#' @param feature_type Optional selection of feature types to include in output matrix.
#'    For multiome data, the options are "Gene Expression" and "Peaks". This option is
#'    only compatible with files from cellranger 3.0 and newer.
#' @return BPCells matrix object
#' @details The 10x format makes use of gzip compression for the matrix data,
#' which can slow down read performance. Consider writing into another format
#' if the read performance is important to you.
#' @examples
#' ## Download example matrices from pbmc 500 dataset and save in temp directory
#' data_dir <- file.path(tempdir(), "mat_10x")
#' dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
#' url_base <- "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_X/"
#' mat_file <- "500_PBMC_3p_LT_Chromium_X_filtered_feature_bc_matrix.h5"
#' rna_url <- paste0(url_base, mat_file)
#' if (!file.exists(file.path(data_dir, mat_file))) {
#'  download.file(rna_url, file.path(data_dir, mat_file), mode="wb")
#' }
#' 
#' #######################################################################
#' ## open_matrix_10x_hdf5() example
#' #######################################################################
#' mat <- open_matrix_10x_hdf5(
#'  file.path(data_dir, mat_file)
#' )
#' mat
#' 
#' 
#' @export
open_matrix_10x_hdf5 <- function(path, feature_type = NULL, buffer_size = 16384L) {
  assert_is_file(path)
  assert_is(buffer_size, "integer")
  if (!is.null(feature_type)) assert_is(feature_type, "character")
  path <- normalizePath(path, mustWork = FALSE)
  all_groups <- hdf5_group_objnames_cpp(path, "/")
  if ("matrix" %in% all_groups) {
    # If we detect something that looks like a new-style 10x matrix, 
    # ignore any other top-level groups in the file
    all_groups <- "matrix"
  }
  mats <- list()
  for (group in all_groups) {
    # Note: to find old-style 10x files for testing purposes, try these links from 10x:
    # Dataset description: https://www.10xgenomics.com/datasets/100-1-1-mixture-of-fresh-frozen-human-hek-293-t-and-mouse-nih-3-t-3-cells-2-standard-2-1-0
    # HDF5 file url: https://cf.10xgenomics.com/samples/cell-exp/2.1.0/hgmm_100/hgmm_100_raw_gene_bc_matrices_h5.h5
    info <- dims_matrix_10x_hdf5_cpp(path, group, buffer_size)
    res <- new("10xMatrixH5",
               path = path, group = group, type = info$type, dim = info$dims, buffer_size = buffer_size,
               dimnames = normalized_dimnames(info$row_names, info$col_names)
    )
    if (!is.null(feature_type)) {
      valid_features <- read_hdf5_string_cpp(path, file.path(group, "features/feature_type"), buffer_size) %in% feature_type # nolint
      res <- res[valid_features, ]
    }
    mats[[group]] <- res
  }
  if (length(mats) == 1) {
    return(mats[[1]])
  }
  rlang::inform("The input HDF5 is detected as an old-style multi-genome file, returning a list of matrices")
  return(mats)
}

#' @rdname open_matrix_10x_hdf5
#' @inheritParams write_matrix_hdf5
#' @param mat IterableMatrix
#' @param barcodes Vector of names for the cells
#' @param feature_ids Vector of IDs for the features
#' @param feature_names Vector of names for the features
#' @param feature_types String or vector of feature types
#' @param feature_metadata Named list of additional metadata vectors
#' to store for each feature
#' @param gzip_level Gzip compression level. Default is 0 (no compression)
#' @param type Data type of the output matrix. Default is `uint32_t` to match a 
#' matrix of 10x UMI counts. Non-integer data types include `float` and 
#' `double`. If `auto`, will use the data type of `mat`.
#' 
#' @details Input matrices must be in column-major storage order,
#' and if the rownames and colnames are not set, names must be
#' provided for the relevant metadata parameters. Some of the
#' metadata parameters are not read by default in BPCells, but
#' it is possible to export them for use with other tools.
#' @examples
#' #######################################################################
#' ## write_matrix_10x_hdf5() example
#' #######################################################################
#' mat <- write_matrix_10x_hdf5(
#'  mat,
#'  file.path(data_dir, paste0("new", mat_file))
#' )
#' mat
#' 
#' 
#' @export
write_matrix_10x_hdf5 <- function(
    mat,
    path,
    barcodes = colnames(mat),
    feature_ids = rownames(mat),
    feature_names = rownames(mat),
    feature_types = "Gene Expression",
    feature_metadata = list(),
    buffer_size = 16384L,
    chunk_size = 1024L,
    gzip_level = 0L,
    type = c("uint32_t", "double", "float", "auto")
) {
  type <- match.arg(type)
  assert_is(mat, "IterableMatrix")
  assert_is(path, "character")
  if (mat@transpose) {
    stop(
      "Matrix must have column-major storage order.\n", 
      "Call t() or transpose_storage_order() first."
    )
  }
  if (type == "auto") {
    type <- matrix_type(mat)
  }
  if (matrix_type(mat) != type) {
    warning(
      "Converting from ", matrix_type(mat), " to ", type, 
      " matrix for output to 10x format"
    )
    mat <- convert_matrix_type(mat, type)
  }
  assert_is(barcodes, "character")
  assert_len(barcodes, ncol(mat))
  assert_is(feature_ids, "character")
  assert_len(feature_ids, nrow(mat))
  assert_is(feature_names, "character")
  assert_len(feature_names, nrow(mat))
  assert_is(feature_types, "character")
  if (!(length(feature_types) %in% c(1, nrow(mat)))) {
    stop("feature_types must have length 1 or nrow(mat)")
  }
  if (length(feature_types) == 1) {
    feature_types <- rep_len(feature_types, nrow(mat))
  }
  assert_is(feature_metadata, "list")
  if (length(feature_metadata) != 0) {
    assert_not_null(names(feature_metadata))
    for (name in names(feature_metadata)) {
      assert_len(feature_metadata[[name]], nrow(mat))
      assert_is(feature_metadata[[name]], "character")
    }
  } else {
    names(feature_metadata) <- character(0)
  }
  assert_is(buffer_size, "integer")
  assert_is(chunk_size, "integer")
  assert_is(gzip_level, "integer")
  
  path <- normalizePath(path, mustWork = FALSE)
  it <- iterate_matrix(mat)
  write_matrix_10x_hdf5_cpp(
    it,
    path,
    type = matrix_type(x = mat),
    barcodes,
    feature_ids,
    feature_names,
    feature_types,
    feature_metadata,
    buffer_size,
    chunk_size,
    gzip_level
  )
  open_matrix_10x_hdf5(path, buffer_size = buffer_size)
}

setClass("AnnDataMatrixH5",
  contains = "IterableMatrix",
  slots = c(
    path = "character",
    group = "character",
    type = "character",
    buffer_size = "integer"
  ),
  prototype = list(
    path = character(0),
    group = "matrix",
    type = character(0),
    buffer_size = integer(0)
  )
)
setMethod("matrix_type", "AnnDataMatrixH5", function(x) x@type)
setMethod("matrix_inputs", "AnnDataMatrixH5", function(x) list())
setMethod("iterate_matrix", "AnnDataMatrixH5", function(x) {
  if (x@transpose) x <- t(x)
  x@dimnames <- denormalize_dimnames(x@dimnames)
  iterate_matrix_anndata_hdf5_cpp(x@path, x@group, x@type, x@buffer_size, x@dimnames[[1]], x@dimnames[[2]])
})
setMethod("short_description", "AnnDataMatrixH5", function(x) {
  sprintf(
    "AnnData HDF5 matrix in file %s, group %s",
    x@path, x@group
  )
})

#' Read/write AnnData matrix
#'
#' @description
#' Read or write a matrix from an anndata hdf5 file. These functions will
#' automatically transpose matrices when converting to/from the AnnData
#' format. This is because the AnnData convention stores cells as rows, whereas the R
#' convention stores cells as columns. If this behavior is undesired, call `t()`
#' manually on the matrix inputs and outputs of these functions. 
#'
#' Most users writing to AnnData files should default to `write_matrix_anndata_hdf5()` rather
#' than the dense variant (see details for more information).
#'
#' @inheritParams open_matrix_hdf5
#' @return AnnDataMatrixH5 object, with cells as the columns.
#' @details 
#'   **Efficiency considerations**: Reading from a dense AnnData matrix will generally be slower
#'   than sparse for single cell datasets, so it is recommended to re-write any dense AnnData
#'   inputs to a sparse format early in processing.
#'
#'   `write_matrix_anndata_hdf5()` should be used by default, as it always writes in the more efficient sparse format.
#'   `write_matrix_anndata_hdf5_dense()` writes in the AnnData dense format, and can be used for smaller matrices 
#'   when efficiency and file size are less of a concern than increased portability (e.g. writing to `obsm` or `varm` matrices).
#'   See the [AnnData docs](https://anndata.readthedocs.io/en/latest/fileformat-prose.html#dense-arrays) for format details.
#'
#'   **Dimension names:** Dimnames are inferred from `obs/_index` or `var/_index` based on length matching.
#'   This helps to infer dimnames for `obsp`,` varm`, etc. If the number of `len(obs) == len(var)`,
#'   dimname inference will be disabled.
#' 
#'   **Signed integers:** When int32 and int64 matrices are read, they will be converted to float and double matrices respectively.
#'   This is because BPCells only supports unsigned integer matrices, and signed integer matrices would have their negative values
#'   misinterpreted as zeros.
#' @examples
#' ## Create temporary directory to keep demo matrix
#' data_dir <- file.path(tempdir(), "mat_anndata")
#' if (dir.exists(data_dir)) unlink(data_dir, recursive = TRUE)
#' dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
#' mat <- get_demo_mat()
#' 
#' 
#' #######################################################################
#' ## write_matrix_anndata_hdf5() example
#' #######################################################################
#' mat <- write_matrix_anndata_hdf5(
#'  mat,
#'  file.path(data_dir, paste0("new_demo_mat.h5"))
#' )
#' mat
#' 
#' 
#' #######################################################################
#' ## open_matrix_anndata_hdf5() example
#' #######################################################################
#' mat <- open_matrix_anndata_hdf5(
#'  file.path(data_dir, paste0("new_demo_mat.h5"))
#' )
#' mat
#' 
#' 
#' @export
open_matrix_anndata_hdf5 <- function(path, group = "X", buffer_size = 16384L) {
  assert_is_file(path)
  assert_is(buffer_size, "integer")

  path <- normalizePath(path, mustWork = FALSE)
  info <- dims_matrix_anndata_hdf5_cpp(path, group, buffer_size)
  res <- new("AnnDataMatrixH5",
    path = path, dim = info$dims, buffer_size = buffer_size,
    group = group, dimnames = normalized_dimnames(info$row_names, info$col_names),
    transpose = info$transpose,
    type = info$type
  )

  return(t(res))
}

#' @rdname open_matrix_anndata_hdf5
#' @inheritParams open_matrix_anndata_hdf5
#' @inheritParams write_matrix_hdf5
#' @param gzip_level Gzip compression level. Default is 0 (no compression)
#' @export
write_matrix_anndata_hdf5 <- function(mat, path, group = "X", buffer_size = 16384L, chunk_size = 1024L, gzip_level = 0L) {
  assert_is(mat, "IterableMatrix")
  assert_is(path, "character")
  mat <- t(mat)
  write_matrix_anndata_hdf5_cpp(
    iterate_matrix(mat),
    path,
    group,
    matrix_type(mat),
    mat@transpose,
    buffer_size,
    chunk_size,
    gzip_level
  )
  open_matrix_anndata_hdf5(path, group, buffer_size)
}

#' @rdname open_matrix_anndata_hdf5
#' @inheritParams write_matrix_anndata_hdf5
#' 
#' @param dataset The dataset within the hdf5 file to write the matrix to. Used for `write_matrix_anndata_hdf5_dense`
#' @examples
#' #######################################################################
#' ## write_matrix_anndata_hdf5_dense() example
#' #######################################################################
#' mat <- write_matrix_anndata_hdf5_dense(
#'  mat,
#'  file.path(data_dir, paste0("new_demo_mat_dense.h5"))
#' )
#' mat
#' 
#' 
#' @export
write_matrix_anndata_hdf5_dense <- function(mat, path, dataset = "X", buffer_size = 16384L, chunk_size = 1024L, gzip_level = 0L) {
  assert_is(mat, "IterableMatrix")
  assert_is(path, "character")
  mat <- t(mat)
  write_matrix_anndata_hdf5_dense_cpp(
    iterate_matrix(mat),
    path,
    dataset,
    matrix_type(mat),
    mat@transpose,
    chunk_size,
    gzip_level
  )
  open_matrix_anndata_hdf5(path, dataset, buffer_size)
}

#' Import MatrixMarket files
#'
#' Read a sparse matrix from a MatrixMarket file. This is a text-based format used by 
#' 10x, Parse, and others to store sparse matrices. 
#' Format details on the [NIST website](https://math.nist.gov/MatrixMarket/formats.html).
#' @inheritParams transpose_storage_order
#' @param mtx_path Path of mtx or mtx.gz file
#' @param row_names Character vector of row names
#' @param col_names Character vector of col names
#' @param row_major If true, store the matrix in row-major orientation
#' @return MatrixDir object with the imported matrix
#' @details Import MatrixMarket mtx files to the BPCells format. This implementation ensures
#'   fixed memory usage even for very large inputs by doing on-disk sorts. It will be
#'   much slower than hdf5 inputs, so only use MatrixMarket format when absolutely necessary.
#'
#'   As a rough speed estimate, importing the 17GB Parse 
#'   [1M PBMC](https://www.parsebiosciences.com/datasets/pbmc/single-cell-rna-sequencing-of-1-million-human-cells-in-a-single-experiment)
#'   `DGE_1M_PBMC.mtx` file takes about 4 minutes and 1.3GB of RAM, producing a compressed output matrix of 1.5GB. `mtx.gz`
#'   files will be slower to import due to gzip decompression.
#'   
#'   When importing from 10x mtx files, the row and column names can be read automatically
#'   using the `import_matrix_market_10x()` convenience function.
#' @export
import_matrix_market <- function(
  mtx_path, outdir = tempfile("matrix_market"), row_names = NULL, col_names = NULL, row_major = FALSE,
  tmpdir = tempdir(), load_bytes = 4194304L, sort_bytes = 1073741824L
) {
  if (!is.null(row_names)) {
    row_names <- as.character(row_names)
  } else {
    row_names <- character(0)
  }
  if (!is.null(col_names)) {
    col_names <- as.character(col_names)
  } else {
    col_names <- character(0)
  }

  mtx_path <- normalizePath(mtx_path, mustWork = FALSE)

  outdir <- normalizePath(outdir, mustWork = FALSE)
  tmpdir <- normalizePath(tmpdir, mustWork = FALSE)
  tmpdir <- tempfile("matrix_market_tmp", tmpdir = tmpdir)
  #on.exit(unlink(tmpdir, recursive = TRUE, expand = FALSE))


  import_matrix_market_cpp(mtx_path, row_names, col_names, outdir, tmpdir, load_bytes, sort_bytes, row_major)
  m <- open_matrix_dir(outdir)
  return(m)
}

#' @rdname import_matrix_market
#' @param mtx_dir Directory holding matrix.mtx.gz, barcodes.tsv.gz, and features.tsv.gz
#' @param feature_type String or vector of feature types to include. (cellranger 3.0 and newer)
#' @export
import_matrix_market_10x <- function(
  mtx_dir, outdir = tempfile("matrix_market"), feature_type=NULL, row_major = FALSE, 
  tmpdir = tempdir(), load_bytes = 4194304L, sort_bytes = 1073741824L
) {
  mtx_dir <- normalizePath(mtx_dir, mustWork=TRUE)
  assert_is_file(file.path(mtx_dir, c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")), multiple_ok=TRUE)
  col_types <- readr::cols(.default=readr::col_character())
  col_names <- readr::read_tsv(file.path(mtx_dir, "barcodes.tsv.gz"), col_names=FALSE, col_types=col_types, progress=FALSE)[[1]]
  features <- readr::read_tsv(file.path(mtx_dir, "features.tsv.gz"), col_names=FALSE, col_types=col_types, progress=FALSE)
  row_names <- features[[1]]
  if (!is.null(feature_type)) {
    assert_true(ncol(features) >= 3)
  }
  
  m <- import_matrix_market(
    file.path(mtx_dir, "matrix.mtx.gz"), 
    row_names = row_names, 
    col_names = col_names,
    row_major = row_major,
    outdir = outdir,
    tmpdir = tmpdir,
    load_bytes = load_bytes,
    sort_bytes = sort_bytes
  )

  if (!is.null(feature_type)) {
    valid_features <- features[[3]] %in% feature_type # nolint
    m <- m[valid_features, ]
  }

  m
}

# Overlap matrix from fragments
setClass("PeakMatrix",
  contains = "IterableMatrix",
  slots = c(
    fragments = "IterableFragments",
    chr_id = "integer",
    start = "integer",
    end = "integer",
    chr_levels = "character",
    mode = "character"
  ),
  prototype = list(
    fragments = NULL,
    chr_id = integer(0),
    start = integer(0),
    end = integer(0),
    chr_levels = character(0),
    mode = "insertions"
  )
)
setMethod("matrix_type", "PeakMatrix", function(x) "uint32_t")
setMethod("matrix_inputs", "PeakMatrix", function(x) list())

#' Calculate ranges x cells overlap matrix
#' @param fragments Input fragments object. Must have cell names and chromosome names defined
#' @param ranges `r document_granges("Peaks/ranges to overlap,")`
#' @param mode Mode for counting peak overlaps. (See "value" section for more details)
#' @inheritParams convert_to_fragments
#' @param explicit_peak_names Boolean for whether to add rownames to the output matrix in format e.g
#'  chr1:500-1000, where start and end coords are given in a 0-based coordinate system.
#'  Note that either way, peak names will be written when the matrix is saved.
#' @note When calculating the matrix directly from a fragments tsv, it's necessary to first call `select_chromosomes()` in order to
#'     provide the ordering of chromosomes to expect while reading the tsv.
#' @return Iterable matrix object with dimension ranges x cells. When saved,
#'   the column names of the output matrix will be in the format chr1:500-1000,
#'   where start and end coords are given in a 0-based coordinate system.
#'
#' **`mode` options**
#'
#' - `"insertions"`: Start and end coordinates are separately overlapped with each peak
#' - `"fragments"`: Like `"insertions"`, but each fragment can contribute at most 1 count
#'    to each peak, even if both the start and end coordinates overlap
#' - `"overlaps"`: Like `"fragments"`, but an overlap is also counted if the fragment fully
#'    spans the peak even if neither the start or end falls within the peak
#' @examples
#' ## Prep demo data
#' frags <- get_demo_frags(subset = FALSE)
#' chrom_sizes <- read_ucsc_chrom_sizes(file.path(tempdir(), "references"), genome="hg38")
#' blacklist <- read_encode_blacklist(file.path(tempdir(), "references"), genome="hg38")
#' frags_filter_blacklist <- frags %>% select_regions(blacklist, invert_selection = TRUE)
#' peaks <- call_peaks_tile(
#'   frags_filter_blacklist, 
#'   chrom_sizes,
#'   effective_genome_size = 2.8e9
#' )
#' top_peaks <- head(peaks, 5000)
#' top_peaks <- top_peaks[order_ranges(top_peaks, chrNames(frags)),]
#' 
#' ## Get peak matrix
#' peak_matrix(frags_filter_blacklist, top_peaks, mode="insertions")
#' @export
peak_matrix <- function(fragments, ranges, mode = c("insertions", "fragments", "overlaps"), zero_based_coords = !is(ranges, "GRanges"), explicit_peak_names = TRUE) {
  assert_is(fragments, "IterableFragments")
  ranges <- normalize_ranges(ranges, zero_based_coords = zero_based_coords)

  assert_is(zero_based_coords, "logical")

  assert_not_null(cellNames(fragments))
  assert_not_na(cellNames(fragments))
  assert_not_null(chrNames(fragments))
  assert_not_na(chrNames(fragments))

  mode <- match.arg(mode)

  peak_order <- order_ranges(ranges, chrNames(fragments))
  if (!all(peak_order == seq_along(peak_order))) {
    rlang::warn("Peaks given out of order. Reordering with order_ranges()")
    ranges$chr <- ranges$chr[peak_order]
    ranges$start <- ranges$start[peak_order]
    ranges$end <- ranges$end[peak_order]
  }

  chr_id <- as.integer(factor(as.character(ranges$chr), chrNames(fragments))) - 1L

  res <- new("PeakMatrix", fragments = fragments, chr_id = chr_id, start = ranges$start, end = ranges$end, mode = mode)
  res@chr_levels <- chrNames(fragments)
  res@dim <- c(length(cellNames(fragments)), length(res@chr_id))
  res@dimnames[[1]] <- cellNames(fragments)
  if (explicit_peak_names) {
    res@dimnames[[2]] <- stringr::str_c(res@chr_levels[chr_id + 1], ":", ranges$start, "-", ranges$end)
  }
  return(t(res))
}

setMethod("iterate_matrix", "PeakMatrix", function(x) {
  it <- iterate_fragments(x@fragments)
  iterate_peak_matrix_cpp(it, x@chr_id, x@start, x@end, x@chr_levels, x@mode)
})
setMethod("short_description", "PeakMatrix", function(x) {
  # Subset strings first to avoid a very slow string concatenation process
  indices <- c(head(seq_along(x@chr_id), 3), tail(seq_along(x@chr_id), 1))
  labels <- paste0(x@chr_levels[1 + x@chr_id[indices]], ":", x@start[indices] + 1, "-", x@end[indices])
  c(
    short_description(x@fragments),
    sprintf(
      "Calculate %d peaks over %d ranges%s", ncol(x), length(x@chr_id),
      pretty_print_vector(labels, prefix = ": ", max_len = 2)
    )
  )
})

setMethod("[", "PeakMatrix", function(x, i, j, ...) {
    if (missing(x)) stop("x is missing in matrix selection")
  # Handle transpose via recursive call
  if (x@transpose) {
    return(t(t(x)[rlang::maybe_missing(j), rlang::maybe_missing(i)]))
  }

  i <- selection_index(i, nrow(x), rownames(x))
  j <- split_selection_index(j, ncol(x), colnames(x))


  x <- selection_fix_dims(x, rlang::maybe_missing(i), rlang::maybe_missing(j$subset))

  if (!rlang::is_missing(i)) {
    x@fragments <- select_cells(x@fragments, i)
  }
  if (!rlang::is_missing(j$subset)) {
    x@chr_id <- x@chr_id[j$subset]
    x@start <- x@start[j$subset]
    x@end <- x@end[j$subset]
  }
  if (!rlang::is_missing(j$reorder)) {
    x <- callNextMethod(x,,j$reorder)
  }

  x
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
    chr_levels = "character",
    mode = "character"
  ),
  prototype = list(
    fragments = NULL,
    chr_id = integer(0),
    start = integer(0),
    end = integer(0),
    tile_width = integer(0),
    chr_levels = character(0),
    mode = character(0)
  )
)
setMethod("matrix_type", "TileMatrix", function(x) "uint32_t")
setMethod("matrix_inputs", "TileMatrix", function(x) list())


#' Calculate ranges x cells tile overlap matrix
#' @param fragments Input fragments object
#' @param ranges `r document_granges("Tiled regions", extras=c("tile_width"="Size of each tile in this region in basepairs"))`  
#'  
#'  Must be non-overlapping and sorted by
#'  (chr, start), with chromosomes ordered according to the chromosome names of `fragments`
#' @inheritParams convert_to_fragments
#' @param explicit_tile_names Boolean for whether to add rownames to the output matrix in format e.g
#'  chr1:500-1000, where start and end coords are given in a 0-based coordinate system. For
#'  whole-genome Tile matrices the names will take ~5 seconds to generate and take up 400MB of memory.
#'  Note that either way, tile names will be written when the matrix is saved.
#' @param mode Mode for counting tile overlaps. (See "value" section for more detail)
#' @note When calculating the matrix directly from a fragments tsv, it's necessary to first call `select_chromosomes()` in order to
#'     provide the ordering of chromosomes to expect while reading the tsv.
#' @return Iterable matrix object with dimension ranges x cells. When saved,
#'   the column names will be in the format chr1:500-1000,
#'   where start and end coords are given in a 0-based coordinate system.
#'
#' **`mode` options**
#'
#' - `"insertions"`: Start and end coordinates are separately overlapped with each tile
#' - `"fragments"`: Like `"insertions"`, but each fragment can contribute at most 1 count
#'    to each tile, even if both the start and end coordinates overlap
#' @examples
#' ## Prep demo data
#' frags <- get_demo_frags(subset = FALSE)
#' chrom_sizes <- read_ucsc_chrom_sizes(file.path(tempdir(), "references"), genome="hg38")
#' blacklist <- read_encode_blacklist(file.path(tempdir(), "references"), genome="hg38")
#' frags_filter_blacklist <- frags %>% select_regions(blacklist, invert_selection = TRUE)
#' ranges <- tibble::tibble(
#'   chr = "chr4",
#'   start = 0,
#'   end = "190214555", 
#'   tile_width = 200
#' )
#' 
#' 
#' ## Get tile matrix
#' tile_matrix(frags_filter_blacklist, ranges)
#' @export
tile_matrix <- function(fragments, ranges, mode = c("insertions", "fragments"), zero_based_coords = !is(ranges, "GRanges"), explicit_tile_names = FALSE) {
  assert_is(fragments, "IterableFragments")
  ranges <- normalize_ranges(ranges, metadata_cols = "tile_width", zero_based_coords = zero_based_coords)

  assert_is(zero_based_coords, "logical")

  assert_not_null(cellNames(fragments))
  assert_not_na(cellNames(fragments))
  assert_not_null(chrNames(fragments))
  assert_not_na(chrNames(fragments))

  mode <- match.arg(mode)

  ranges$tile_width <- as.integer(ranges$tile_width)
  assert_true(all(as.character(ranges$chr) %in% chrNames(fragments)))

  tile_order <- order_ranges(ranges, chrNames(fragments))
  if (!all(tile_order == seq_along(tile_order))) {
    rlang::warn("Tiles given out of order. Reordering with order_ranges()")
    ranges$chr <- ranges$chr[tile_order]
    ranges$start <- ranges$start[tile_order]
    ranges$end <- ranges$end[tile_order]
    ranges$tile_width <- ranges$tile_width[tile_order]
  }
  chr_id <- as.integer(factor(as.character(ranges$chr), chrNames(fragments))) - 1L

  # Check to make sure tiles are non-overlapping
  pair_1 <- seq_len(length(chr_id) - 1)
  pair_2 <- pair_1 + 1
  if (any(chr_id[pair_1] == chr_id[pair_2] & ranges$end[pair_1] >= ranges$start[pair_2])) {
    rlang::abort("Tile regions must be non-overlapping")
  }

  # Construct tile matrix
  res <- new("TileMatrix", fragments = fragments, chr_id = chr_id, start = ranges$start, end = ranges$end, tile_width = ranges$tile_width, mode=mode)
  res@chr_levels <- chrNames(fragments)

  tiles <- as.integer(ceiling((ranges$end - ranges$start) / ranges$tile_width))
  res@dim <- c(length(cellNames(fragments)), sum(tiles))
  res@dimnames[[1]] <- cellNames(fragments)
  if (explicit_tile_names) {
    res@dimnames[[2]] <- get_tile_names_cpp(chr_id, ranges$start, ranges$end, ranges$tile_width, res@chr_levels)
  }
  return(t(res))
}

#' Get ranges corresponding to selected tiles of a tile matrix
#' @keywords internal
tile_ranges <- function(tile_matrix, selection) {
  # Handle manually transposed tile_matrix objects
  if (!tile_matrix@transpose) {
    tile_matrix <- t(tile_matrix)
  }
  indices <- seq_len(nrow(tile_matrix))
  selection <- vctrs::vec_slice(indices, selection) - 1
  ranges <- get_tile_ranges_cpp(
    tile_matrix@chr_id, tile_matrix@start, tile_matrix@end,
    tile_matrix@tile_width, tile_matrix@chr_levels,
    selection
  )
  ranges$chr <- factor(tile_matrix@chr_levels[ranges$chr + 1], levels = tile_matrix@chr_levels)
  return(tibble::as_tibble(ranges))
}

setMethod("iterate_matrix", "TileMatrix", function(x) {
  it <- iterate_fragments(x@fragments)
  iterate_tile_matrix_cpp(it, x@chr_id, x@start, x@end, x@tile_width, x@chr_levels, x@mode)
})

setMethod("short_description", "TileMatrix", function(x) {
  # Subset strings first to avoid a very slow string concatenation process
  indices <- c(head(seq_along(x@chr_id), 3), tail(seq_along(x@chr_id), 1))
  labels <- paste0(x@chr_levels[1 + x@chr_id[indices]], ":", x@start[indices] + 1, "-", x@end[indices], " (", x@tile_width[indices], "bp)")
  tiles <- sum(as.integer(ceiling((x@end - x@start) / x@tile_width)))
  c(
    short_description(x@fragments),
    sprintf(
      "Calculate %d tiles over %d ranges%s", tiles, length(x@chr_id),
      pretty_print_vector(labels, prefix = ": ", max_len = 2)
    )
  )
})

setMethod("[", "TileMatrix", function(x, i, j, ...) {
  if (missing(x)) stop("x is missing in matrix selection")

  # Handle transpose via recursive call
  if (x@transpose) {
    return(t(t(x)[rlang::maybe_missing(j), rlang::maybe_missing(i)]))
  }
  
  i <- selection_index(i, nrow(x), rownames(x))
  j <- split_selection_index(j, ncol(x), colnames(x))

  x_orig <- x
  x <- selection_fix_dims(x, rlang::maybe_missing(i), rlang::maybe_missing(j$subset))

  if (!rlang::is_missing(i)) {
    x@fragments <- select_cells(x@fragments, i)
  }
  if (!rlang::is_missing(j$subset)) {
    if (mean(diff(j$subset) == 1) > .9) {
      # If 90% contiguous, then use a slice
      new_tiles <- subset_tiles_cpp(
        x@chr_id, x@start, x@end, x@tile_width, x@chr_levels,
        j$subset - 1
      )
      x@chr_id <- new_tiles$chr_id
      x@start <- new_tiles$start
      x@end <- new_tiles$end
      x@tile_width <- new_tiles$tile_width
    } else {
      # Otherwise, convert to a peak matrix
      peaks <- tile_ranges(x_orig, j$subset)
      x <- t(peak_matrix(x@fragments, peaks))
    }
  }
  if (!rlang::is_missing(j$reorder)) {
    x <- callNextMethod(x,,j$reorder)
  }
  x
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
setMethod("matrix_type", signature(x = "ConvertMatrixType"), function(x) x@type)
setMethod("iterate_matrix", "ConvertMatrixType", function(x) {
  iter_function <- get(sprintf("convert_matrix_%s_%s_cpp", matrix_type(x@matrix), matrix_type(x)))
  it <- iterate_matrix(x@matrix)
  iter_function(it)
})
setMethod("short_description", "ConvertMatrixType", function(x) {
  c(
    short_description(x@matrix),
    sprintf("Convert type from %s to %s", matrix_type(x@matrix), matrix_type(x))
  )
})

setMethod("[", "ConvertMatrixType", function(x, i, j, ...) {
  if (missing(x)) stop("x is missing in matrix selection")

  i <- selection_index(i, nrow(x), rownames(x))
  j <- selection_index(j, ncol(x), colnames(x))
  x <- selection_fix_dims(x, rlang::maybe_missing(i), rlang::maybe_missing(j))

  x@matrix <- x@matrix[rlang::maybe_missing(i), rlang::maybe_missing(j)]
  x
})

#' Convert the type of a matrix
#' @param matrix IterableMatrix object input
#' @param type One of uint32_t (unsigned 32-bit integer), float (32-bit real number),
#'   or double (64-bit real number)
#' @return IterableMatrix object
#' @examples
#' mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
#' rownames(mat) <- paste0("gene", seq_len(10))
#' colnames(mat) <- paste0("cell", seq_len(5))
#' mat <- mat %>% as("dgCMatrix") %>% as("IterableMatrix")
#' mat
#' convert_matrix_type(mat, "float")
#' @export
convert_matrix_type <- function(matrix, type = c("uint32_t", "double", "float")) {
  assert_is(matrix, c("dgCMatrix", "IterableMatrix"))
  type <- match.arg(type)
  if (is(matrix, "dgCMatrix")) {
    matrix <- as(matrix, "IterableMatrix")
  }
  if (matrix_type(matrix) == type) {
    return(matrix)
  } else if (is(matrix, "ConvertMatrixType")) {
    if (matrix_type(matrix@matrix) == type) {
      ret <- matrix@matrix
      # Restore dimnames that would have been cleared with wrapMatrix
      dimnames(ret) <- dimnames(matrix) 
      return(ret)
    } else {
      matrix@type <- type
      return(matrix)
    }
  } else if (matrix@transpose) {
    return(t(convert_matrix_type(t(matrix), type)))
  }
  wrapMatrix("ConvertMatrixType", matrix, type = type)
}

# Conversions with dgCMatrix

#' Convert between BPCells matrix and R objects.
#'
#' BPCells matrices can be interconverted with Matrix package 
#' dgCMatrix sparse matrices, as well as base R
#' dense matrices (though this may result in high memory usage for large matrices)
#'
#' @usage
#' # Convert to R from BPCells
#' as(bpcells_mat, "dgCMatrix") # Sparse matrix conversion
#' as.matrix(bpcells_mat) # Dense matrix conversion
#' 
#' # Convert to BPCells from R
#' as(dgc_mat, "IterableMatrix")
#' @examples
#' mat <- get_demo_mat()[1:2, 1:2]
#' mat
#' 
#' 
#' #######################################################################
#' ## as(bpcells_mat, "dgCMatrix") example
#' #######################################################################
#' mat_dgc <- as(mat, "dgCMatrix")
#' mat_dgc
#' 
#' 
#' ## as.matrix(bpcells_mat) example
#' as.matrix(mat)
#' 
#' ## Alternatively, can also use function as()
#' as(mat, "matrix")
#' 
#' #######################################################################
#' ## as(dgc_mat, "IterableMatrix") example
#' #######################################################################
#' as(mat_dgc, "IterableMatrix")
#' 
#' 
#' @name matrix_R_conversion
NULL

setClass("Iterable_dgCMatrix_wrapper",
  contains = "IterableMatrix",
  slots = c(
    mat = "dgCMatrix"
  ),
  prototype = list(
    mat = NULL
  )
)
setMethod("matrix_type", signature(x = "Iterable_dgCMatrix_wrapper"), function(x) "double")
setMethod("matrix_inputs", "Iterable_dgCMatrix_wrapper", function(x) list())

setAs("dgCMatrix", "IterableMatrix", function(from) {
  new("Iterable_dgCMatrix_wrapper", dim = dim(from), dimnames = dimnames(from), transpose = FALSE, mat = from)
})
setMethod("iterate_matrix", "Iterable_dgCMatrix_wrapper", function(x) {
  if (x@transpose) x <- t(x)
  x@dimnames <- denormalize_dimnames(x@dimnames)
  iterate_csparse_matrix_cpp(x@mat, x@dimnames[[1]], x@dimnames[[2]])
})
setMethod("short_description", "Iterable_dgCMatrix_wrapper", function(x) {
  "Load dgCMatrix from memory"
})

setMethod("iterate_matrix", "dgCMatrix", function(x) {
  iterate_csparse_matrix_cpp(x)
})

setAs("IterableMatrix", "dgCMatrix", function(from) {
  res <- write_matrix_memory(convert_matrix_type(from, "double"), compress=FALSE)
  if (length(res@index) >= 2^31-1) {
    rlang::abort(c(
      "Error converting IterableMatrix to dgCMatrix",
      "dgCMatrix objects cannot hold more than 2^31 non-zero entries",
      sprintf("Input matrix has %0.f entries", length(res@index))
    ))
  }
  if (from@transpose) {
    res <- Matrix::sparseMatrix(
      j = res@index,
      p = res@idxptr,
      x = res@val,
      index1=FALSE,
      dims = dim(res),
      dimnames = dimnames(res)
    )
    return(as(res, "dgCMatrix"))
  } else {
    return(Matrix::sparseMatrix(
      i = res@index,
      p = res@idxptr,
      x = res@val,
      index1=FALSE,
      dims = dim(res),
      dimnames = dimnames(res)
    ))
  }
})

# Add conversion to and from base R dense matrices
setAs("IterableMatrix", "matrix", function(from) {
  rlang::inform(c(
      "Warning: Converting to a dense matrix may use excessive memory"
    ), .frequency = "regularly", .frequency_id = "matrix_dense_conversion")
  # `mat` will always be numeric mode
  mat <- as.matrix(as(from, "dgCMatrix"))
  # to keep the original mode, we transform it when necessary
  if (matrix_type(from) == "uint32_t") {
      mat <- matrix_to_integer(mat)
  }
  mat
})

setAs("matrix", "IterableMatrix", function(from) mat <- as(as(from, "dgCMatrix"), "IterableMatrix"))

matrix_to_integer <- function(matrix) { # a numeric matrix
    if (is.integer(matrix)) return(matrix) # styler: off
    if (all(matrix <= .Machine$integer.max)) {
        storage.mode(matrix) <- "integer"
    } else {
      warning(
        "Using `double` mode since some values exceed `.Machine$integer.max`"
      )
    }
    matrix
}

#' @exportS3Method base::as.matrix
as.matrix.IterableMatrix <- function(x, ...) as(x, "matrix")

#' @export
setMethod("as.matrix", signature(x = "IterableMatrix"), function(x, ...) as(x, "matrix"))

#' Calculate matrix stats
#' @param matrix Input matrix object
#' @param row_stats Which row statistics to compute
#' @param col_stats Which col statistics to compute
#' @param threads Number of threads to use during execution
#' @return List of row_stats: matrix of n_stats x n_rows,
#'          col_stats: matrix of n_stats x n_cols
#' @details The statistics will be calculated in a single pass over the matrix,
#' so this method is desirable to use for efficiency purposes compared to
#' the more standard rowMeans or colMeans if multiple statistics are needed.
#' The stats are ordered by complexity: nonzero, mean, then variance. All
#' less complex stats are calculated in the process of calculating a more complicated stat.
#' So to calculate mean and variance simultaneously, just ask for variance,
#' which will compute mean and nonzero counts as a side-effect
#' @examples
#' mat <- matrix(rpois(100, lambda = 5), nrow = 10)
#' rownames(mat) <- paste0("gene", 1:10)
#' colnames(mat) <- paste0("cell", 1:10)
#' mat <- mat %>% as("dgCMatrix") %>% as("IterableMatrix")
#' 
#' ## By default, no row or column stats are calculated
#' res_none <- matrix_stats(mat)
#' res_none
#' 
#' ## Request row variance (automatically computes mean and nonzero too)
#' res_row_var <- matrix_stats(mat, row_stats = "variance")
#' res_row_var
#' 
#' ## Request both row variance and column variance
#' res_both_var <- matrix_stats(
#'   mat = mat,
#'   row_stats = "variance",
#'   col_stats = "mean"
#' )
#' res_both_var
#' @export
matrix_stats <- function(matrix,
                         row_stats = c("none", "nonzero", "mean", "variance"),
                         col_stats = c("none", "nonzero", "mean", "variance"),
                         threads = 0L
                         ) {
  if (!is(matrix, "IterableMatrix")) {
    if (canCoerce(matrix, "IterableMatrix")) {
      matrix <- as(matrix, "IterableMatrix")
    } else {
      rlang::abort("Input matrix cannot be converted to an IterableMatrix object")
    }
  }
  assert_is_wholenumber(threads)

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

  it <- matrix %>%
    convert_matrix_type("double") %>%
    parallel_split(threads, threads*4) %>%
    iterate_matrix()
  res <- matrix_stats_cpp(it, row_stats_number, col_stats_number)
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


#' @export
svds <- function (A, k, nu = k, nv = k, opts = list(), ...) UseMethod("svds")

# RSpectra exports svds as an S3 Generic, but is a suggested dependency
# With this approach, IterableMatrix objects will work with RSpectra::svds,
# And we will additionally create + register BPCells::svds which will 
# default to calling RSpectra::svds if possible at run-time.
# This way, users should be fine calling either RSpectra::svds or BPCells::svds in the event both are installed
rlang::on_load({
  vctrs::s3_register("RSpectra::svds", "IterableMatrix")
})

#' @export
svds.default <- function(A, k, nu = k, nv = k, opts = list(), ...) {
  if (requireNamespace("RSpectra", quietly = TRUE)) {
    RSpectra::svds(A=A, k=k, nu=nu, nv=nv, opts=opts, ...)
  } else {
    stop("Can't run svds on a non-BPCells object unless RSpectra is installed.")
  }
}

#' @export
svds.IterableMatrix <- function(A, k, nu = k, nv = k, opts = list(), threads=0, ...) {
  assert_is_wholenumber(threads)
  assert_is_wholenumber(k)
  assert_true(k == nu)
  assert_true(k == nv)
  assert_true(k < min(nrow(A), ncol(A)))
  if (min(nrow(A), ncol(A)) < 3) stop("Must have at least 3 rows and cols")

  assert_true(!("center" %in% names(opts)))
  assert_true(!("scale" %in% names(opts)))

  solver_params <- list(
    ncv = min(min(nrow(A), ncol(A)), max(2*k+1, 20)),
    tol = 1e-5,
    maxitr = 1000
  )
  solver_params[names(opts)] <- opts

  it <- A %>%
    convert_matrix_type("double") %>%
    parallel_split(threads, threads*4) %>%
    iterate_matrix()
  
  svds_cpp(
    it, 
    k, 
    solver_params[["ncv"]],
    solver_params[["maxitr"]],
    solver_params[["tol"]],
    A@transpose
  )
}


#' Calculate the MD5 checksum of an IterableMatrix
#'
#' @description
#'   Calculate the MD5 checksum of an IterableMatrix and return the checksum in
#'   hexidecimal format.
#' @details
#'   `checksum()` converts the non-zero elements of the sparse input matrix to double
#'   precision, concatenates each element value with the element row and column index words,
#'   and uses these 16-byte blocks along with the matrix dimensions and row and column
#'   names to calculate the checksum. The checksum value depends on the storage order so
#'   column- and row-order matrices with the same element values give different checksum
#'   values. `checksum()` uses element and index values in little-endian CPU storage order.
#'   It converts to little-endian order on big-endian architecture although this has not
#'   been tested.
#' @param matrix IterableMatrix object
#' @return MD5 checksum string in hexidecimal format.
#' @examples
#' library(Matrix)
#' library(BPCells)
#' m1 <- matrix(seq(1,12), nrow=3)
#' m2 <- as(m1, 'dgCMatrix')
#' m3 <- as(m2, 'IterableMatrix')
#' checksum(m3)
#' @export
checksum <- function(matrix) {
    assert_is(matrix, "IterableMatrix")

    iter <- iterate_matrix(BPCells:::convert_matrix_type(matrix, "double"))
    checksum_double_cpp(iter)
}

#' Apply a function to summarize rows/cols
#'
#' Apply a custom R function to each row/col of a BPCells matrix. This
#' will run slower than the builtin C++-backed functions, but will 
#' keep most of the memory benefits from disk-backed operations.
#' @param mat IterableMatrix object
#' @param fun `function(val, row, col)` that takes in a row/col of values and returns a summary output. Argument details:
#'
#'   1. `val` - Vector length (# non-zero values) with the value for each non-zero matrix entry
#'   2. `row` - one-based row index (`apply_by_col`: vector length (# non-zero values), `apply_by_row`: single integer)
#'   3. `col` - one-based col index (`apply_by_col`: single integer, `apply_by_row`: vector length (# non-zero values))
#'   4. `...` - Optional additional arguments (should not be named row, col, or val)
#' @param ... Optional additional arguments passed to `fun`
#' @return **apply_by_row** - A list of length `nrow(matrix)` with the results returned by `fun()` on each row
#' @details These functions require row-major matrix storage for apply_by_row and col-major storage for apply_by_col,
#' so matrices stored in the wrong order may neeed a re-ordered copy created using `transpose_storage_order()` first.
#' This is required to be able to keep memory-usage low and allow calculating the result with a single streaming pass of the
#' input matrix.
#' 
#' If vector/matrix outputs are desired instead of lists, calling `unlist(x)` or `do.call(cbind, x)` or `do.call(rbind, x)`
#' can convert the list output.
#'
#' @seealso For an interface more similar to `base::apply`, see the [BPCellsArray](https://github.com/Yunuuuu/BPCellsArray/)
#' project. For calculating colMeans on a sparse single cell RNA matrix it is about 8x slower than `apply_by_col`, due to the
#' `base::apply` interface not being sparsity-aware. (See [pull request #104](https://github.com/bnprks/BPCells/pull/104) for benchmarking.)
#' @examples
#' mat <- matrix(rbinom(40, 1, 0.5) * sample.int(5, 40, replace = TRUE), nrow = 4)
#' rownames(mat) <- paste0("gene", 1:4)
#' mat
#' 
#' mat <- mat %>% as("dgCMatrix") %>% as("IterableMatrix")
#' 
#' #######################################################################
#' ## apply_by_row() example
#' #######################################################################
#' ## Get mean of every row
#' 
#' ## expect an error in the case that col-major matrix is passed
#' apply_by_row(mat, function(val, row, col) {sum(val) / nrow(mat)}) %>% 
#'  unlist()
#' 
#' ## Need to transpose matrix to make sure it is in row-order
#' mat_row_order <- transpose_storage_order(mat)
#' 
#' ## works as expected for row major
#' apply_by_row(mat_row_order, 
#'  function(val, row, col) sum(val) / ncol(mat_row_order)
#' ) %>% unlist()
#' 
#' # Also analogous to running rowMeans() without names
#' rowMeans(mat)
#' 
#' 
#' @export
apply_by_row <- function(mat, fun, ...) {
  assert_is(mat, "IterableMatrix")
  if (storage_order(mat) != "row") {
    rlang::abort("Cannot call apply_by_row on a col-major matrix. Please call transpose_storage_order() first")
  }
  if (length(list(...)) > 0) {
    f <- function(val, row, col) {fun(val, row, col, ...)}
  } else {
    f <- fun
  }
  apply_matrix_double_cpp(iterate_matrix(convert_matrix_type(mat, "double")), f, TRUE)
}

#' @return **apply_by_col** - A list of length `ncol(matrix)` with the results returned by `fun()` on each row
#' @rdname apply_by_row
#' @examples
#' #######################################################################
#' ## apply_by_col() example
#' #######################################################################
#' ## Get argmax of every col
#' apply_by_col(mat, 
#'  function(val, row, col) if (length(val) > 0) row[which.max(val)] else 1L
#' ) %>% unlist()
#' 
#' 
#' @export
apply_by_col <- function(mat, fun, ...) {
  if (storage_order(mat) != "col") {
    rlang::abort("Cannot call apply_by_col on a row-major matrix. Please call transpose_storage_order() first")
  }
  if (length(list(...)) > 0) {
    f <- function(val, row, col) {fun(val, row, col, ...)}
  } else {
    f <- fun
  }
  apply_matrix_double_cpp(iterate_matrix(convert_matrix_type(mat, "double")), f, FALSE)
}

