
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
setMethod("dimnames<-", signature(x = "IterableMatrix", value = "list"), function(x, value) {
  d <- dim(x)
  has_error <- FALSE
  if (!is.list(value) || length(value) != 2) has_error <- TRUE
  if (!is.null(value[[1]]) && length(value[[1]]) != d[1]) has_error <- TRUE
  if (!is.null(value[[2]]) && length(value[[2]]) != d[2]) has_error <- TRUE
  if (has_error) stop("Invalid dimnames supplied")

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
  x@dimnames <- list(NULL, NULL)
  x
})

#' Construct an S4 matrix object wrapping another matrix object
#'
#' Helps to avoid duplicate storage of dimnames
#' @keywords internal
wrapMatrix <- function(class, m, ...) {
  dimnames <- dimnames(m)
  if (matrix_is_transform(m)) m@dimnames <- list(NULL, NULL)
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
setGeneric("matrix_type", function(x) standardGeneric("matrix_type"))

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

setMethod("short_description", "IterableMatrix", function(x) {
  character(0)
})


#' @describeIn IterableMatrix-methods Display an IterableMatrix
#' @param object IterableMatrix object
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
    return(t(dense_multiply_left_cpp(x@xptr, t(y))))
  } else {
    return(dense_multiply_right_cpp(x@xptr, y))
  }
})

setMethod("%*%", signature(x = "matrix", y = "LinearOperator"), function(x, y) {
  if (y@transpose) {
    return(t(dense_multiply_right_cpp(y@xptr, t(x))))
  } else {
    return(dense_multiply_left_cpp(y@xptr, x))
  }
})

setMethod("%*%", signature(x = "LinearOperator", y = "numeric"), function(x, y) {
  if (x@transpose) {
    return(vec_multiply_left_cpp(x@xptr, y))
  } else {
    return(vec_multiply_right_cpp(x@xptr, y))
  }
})

setMethod("%*%", signature(x = "numeric", y = "LinearOperator"), function(x, y) {
  if (y@transpose) {
    return(vec_multiply_right_cpp(y@xptr, x))
  } else {
    return(vec_multiply_left_cpp(y@xptr, x))
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
  if (x@transpose != y@transpose) stop("Cannot multiply matrices with different interal transpose states.\nPlease use transpose_storage_order().")
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

  i <- selection_index(i, nrow(x), rownames(x))
  j <- selection_index(j, ncol(x), colnames(x))
  x <- selection_fix_dims(x, rlang::maybe_missing(i), rlang::maybe_missing(j))

  # Selection will be a no-op if i or j is missing
  x@left <- x@left[rlang::maybe_missing(i), ]
  x@right <- x@right[, rlang::maybe_missing(j)]

  x
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
  
  if (mat@transpose != mask@transpose) stop("Cannot mask matrices with different interal transpose states.\nPlease use transpose_storage_order().")
  mask <- convert_matrix_type(mask, "uint32_t")

  wrapMatrix("MatrixMask",
    mat,
    mask = mask,
    invert = invert
  )
}


# Row sums and row means

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate rowSums
#' @return * `rowSums()`: vector of row sums
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
setMethod("rowMeans", signature(x = "IterableMatrix"), function(x) rowSums(x) / ncol(x))

#' @param x IterableMatrix object
#' @describeIn IterableMatrix-methods Calculate colMeans
#' @return * `colMeans()`: vector of col means
setMethod("colMeans", signature(x = "IterableMatrix"), function(x) colSums(x) / nrow(x))


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
setMethod("matrix_type", signature(x = "MatrixSubset"), function(x) matrix_type(x@matrix))

# Helper function to convert logical/character indexing into numeric indexing
selection_index <- function(selection, dim_len, dimnames) {
  if (!rlang::is_missing(selection)) {
    indices <- seq_len(dim_len)
    names(indices) <- dimnames
    if (!is.logical(selection)) assert_distinct(selection, n = 2)
    return(vctrs::vec_slice(indices, selection))
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
setMethod("[", "IterableMatrix", function(x, i, j, ...) {
  if (missing(x)) stop("x is missing in matrix selection")
  if (rlang::is_missing(i) && rlang::is_missing(j)) {
    return(x)
  }

  ret <- wrapMatrix("MatrixSubset", x)
  i <- selection_index(i, nrow(x), rownames(x))
  j <- selection_index(j, ncol(x), colnames(x))
  ret <- selection_fix_dims(ret, rlang::maybe_missing(i), rlang::maybe_missing(j))

  if (x@transpose) {
    tmp <- rlang::maybe_missing(i)
    i <- rlang::maybe_missing(j)
    j <- rlang::maybe_missing(tmp)
  }
  if (!rlang::is_missing(i)) ret@row_selection <- i
  if (!rlang::is_missing(j)) ret@col_selection <- j
  ret
})

setMethod("[", "MatrixSubset", function(x, i, j, ...) {
  if (missing(x)) stop("x is missing in matrix selection")
  i <- selection_index(i, nrow(x), rownames(x))
  j <- selection_index(j, ncol(x), colnames(x))
  x <- selection_fix_dims(x, rlang::maybe_missing(i), rlang::maybe_missing(j))

  if (x@transpose) {
    tmp <- rlang::maybe_missing(i)
    i <- rlang::maybe_missing(j)
    j <- rlang::maybe_missing(tmp)
  }
  if (!rlang::is_missing(i)) {
    if (length(x@row_selection) == 0) {
      x@row_selection <- i
    } else {
      x@row_selection <- x@row_selection[i]
    }
  }
  if (!rlang::is_missing(j)) {
    if (length(x@col_selection) == 0) {
      x@col_selection <- j
    } else {
      x@col_selection <- x@col_selection[j]
    }
  }
  x
})


setMethod("iterate_matrix", "MatrixSubset", function(x) {
  assert_true(matrix_type(x) %in% c("uint32_t", "float", "double"))
  
  iter_row_function <- get(sprintf("iterate_matrix_row_select_%s_cpp", matrix_type(x)))
  iter_col_function <- get(sprintf("iterate_matrix_col_select_%s_cpp", matrix_type(x)))

  ret <- iterate_matrix(x@matrix)

  if (length(x@row_selection) != 0) ret <- iter_row_function(ret, x@row_selection - 1L)
  if (length(x@col_selection) != 0) ret <- iter_col_function(ret, x@col_selection - 1L)
  
  ret
})

setMethod("short_description", "MatrixSubset", function(x) {
  if (x@transpose) {
    rows <- x@col_selection
    cols <- x@row_selection
  } else {
    rows <- x@row_selection
    cols <- x@col_selection
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
    "%s: %s names presenent on some but not all matrices. Setting missing names to \"\"",
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
    matrix_list = "list"
  ),
  prototype = list(
    matrix_list = list()
  )
)
setMethod("matrix_type", signature(x = "RowBindMatrices"), function(x) matrix_type(x@matrix_list[[1]]))

setMethod("iterate_matrix", "RowBindMatrices", function(x) {
  iter_function <- get(sprintf("iterate_matrix_row_bind_%s_cpp", matrix_type(x)))
  iterators <- lapply(x@matrix_list, iterate_matrix)
  iter_function(iterators)
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
    "Concatenate rows of %d matrix objects with classes%s",
    length(x@matrix_list),
    pretty_print_vector(vapply(x@matrix_list, class, character(1)), prefix = ": ", max_len = 3)
  )
})

setMethod("rbind2", signature(x = "IterableMatrix", y = "IterableMatrix"), function(x, y, ...) {
  if (x@transpose != y@transpose) stop("Cannot merge matrices with different interal transpose states.\nPlease use transpose_storage_order().")
  if (x@transpose) {
    return(t(cbind2(t(x), t(y))))
  }

  if (ncol(x) != ncol(y)) stop("Error in rbind: matrices must have equal number of columns")
  # Handle dimnames
  col_names <- merge_dimnames(colnames(x), colnames(y), "rbind", "column")
  row_names <- concat_dimnames(rownames(x), rownames(y), nrow(x), nrow(y), "rbind", "row")
  if (matrix_is_transform(x)) x@dimnames <- list(NULL, NULL)
  if (matrix_is_transform(y)) y@dimnames <- list(NULL, NULL)

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

setClass("ColBindMatrices",
  contains = "IterableMatrix",
  slots = c(
    matrix_list = "list"
  ),
  prototype = list(
    matrix_list = list()
  )
)
setMethod("matrix_type", signature(x = "ColBindMatrices"), function(x) matrix_type(x@matrix_list[[1]]))

setMethod("iterate_matrix", "ColBindMatrices", function(x) {
  iter_function <- get(sprintf("iterate_matrix_col_bind_%s_cpp", matrix_type(x)))
  iterators <- lapply(x@matrix_list, iterate_matrix)
  iter_function(iterators)
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
    "Concatenate columns of %d matrix objects with classes%s",
    length(x@matrix_list),
    pretty_print_vector(vapply(x@matrix_list, class, character(1)), prefix = ": ", max_len = 3)
  )
})

setMethod("cbind2", signature(x = "IterableMatrix", y = "IterableMatrix"), function(x, y, ...) {
  if (x@transpose != y@transpose) stop("Cannot merge matrices with different interal transpose states.\nPlease use transpose_storage_order().")
  if (x@transpose) {
    return(t(rbind2(t(x), t(y))))
  }

  if (nrow(x) != nrow(y)) stop("Error in cbind: matrices must have equal number of columns")
  # Handle dimnames
  row_names <- merge_dimnames(rownames(x), rownames(y), "cbind", "row")
  col_names <- concat_dimnames(colnames(x), colnames(y), ncol(x), ncol(y), "cbind", "column")
  if (matrix_is_transform(x)) x@dimnames <- list(NULL, NULL)
  if (matrix_is_transform(y)) y@dimnames <- list(NULL, NULL)

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

# Packed integer matrix
setClass("PackedMatrixMemBase",
  contains = "IterableMatrix",
  slots = c(
    # Leave out val storage since it's datatype-dependent
    index_data = "integer",
    index_starts = "integer",
    index_idx = "integer",
    idxptr = "integer",
    version = "character"
  ),
  prototype = list(
    # Leave out val storage since it's datatype-dependent
    index_data = integer(0),
    index_starts = integer(0),
    index_idx = integer(0),
    idxptr = integer(0),
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
    val_idx = "integer"
  ),
  prototype = list(
    val_data = integer(0),
    val_idx = integer(0)
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
    idxptr = "integer",
    version = "character"
  ),
  prototype = list(
    index = integer(0),
    idxptr = integer(0),
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
#' @export
transpose_storage_order <- function(matrix, outdir = tempfile("transpose"), tmpdir = tempdir(), load_bytes = 4194304L, sort_bytes = 1073741824L) {
  assert_true(matrix_type(matrix) %in% c("uint32_t", "float", "double"))

  write_function <- get(sprintf("write_matrix_transpose_%s_cpp", matrix_type(matrix)))

  outdir <- normalizePath(outdir, mustWork = FALSE)
  tmpdir <- normalizePath(tmpdir, mustWork = FALSE)
  tmpdir <- tempfile("transpose_tmp", tmpdir = tmpdir)

  it <- iterate_matrix(matrix)
  write_function(it, outdir, tmpdir, load_bytes, sort_bytes, !matrix@transpose)

  unlink(tmpdir, recursive = TRUE, expand = FALSE)
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
#' @param overwrite If `TRUE`, overwrite any pre-existing data
#' @export
write_matrix_dir <- function(mat, dir, compress = TRUE, buffer_size = 8192L, overwrite = FALSE) {
  assert_is(mat, c("IterableMatrix", "dgCMatrix"))
  if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")

  assert_is(dir, "character")
  assert_is(compress, "logical")
  assert_is(buffer_size, "integer")
  assert_is(overwrite, "logical")

  assert_true(matrix_type(mat) %in% c("uint32_t", "float", "double"))
  if (compress && matrix_type(mat) != "uint32_t") {
    rlang::inform(c(
      "Warning: Matrix compression performs poorly with non-integers.",
      "Consider calling convert_matrix_type if a compressed integer matrix is intended."
    ), .frequency = "regularly", .frequency_id = "matrix_compress_non_integer")
  }

  dir <- normalizePath(dir, mustWork = FALSE)
  it <- iterate_matrix(mat)

  write_function <- get(sprintf("write_%s_matrix_file_%s_cpp", ifelse(compress, "packed", "unpacked"), matrix_type(mat)))
  write_function(it, dir, buffer_size, overwrite, mat@transpose)

  open_matrix_dir(dir, buffer_size)
}

#' @rdname matrix_io
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
#' @export
write_matrix_hdf5 <- function(mat, path, group, compress = TRUE, buffer_size = 8192L, chunk_size = 1024L, overwrite=FALSE) {
  assert_is(mat, c("IterableMatrix", "dgCMatrix"))
  if (is(mat, "dgCMatrix")) mat <- as(mat, "IterableMatrix")

  assert_is(path, "character")
  assert_is(group, "character")
  assert_is(compress, "logical")
  assert_is(buffer_size, "integer")
  assert_is(chunk_size, "integer")
  assert_is(overwrite, "logical")

  assert_true(matrix_type(mat) %in% c("uint32_t", "float", "double"))
  if (compress && matrix_type(mat) != "uint32_t") {
    rlang::inform(c(
      "Warning: Matrix compression performs poorly with non-integers.",
      "Consider calling convert_matrix_type if a compressed integer matrix is intended."
    ), .frequency = "regularly", .frequency_id = "matrix_compress_non_integer")
  }

  path <- normalizePath(path, mustWork = FALSE)
  if (overwrite && hdf5_group_exists_cpp(path, group)) {
    rlang::inform(c(
      "Warning: Overwriting an hdf5 dataset does not free old storage"
    ), .frequency = "regularly", .frequency_id = "hdf5_overwrite")
  }

  it <- iterate_matrix(mat)

  write_function <- get(sprintf("write_%s_matrix_hdf5_%s_cpp", ifelse(compress, "packed", "unpacked"), matrix_type(mat)))
  write_function(it, path, group, buffer_size, chunk_size, overwrite, mat@transpose)

  open_matrix_hdf5(path, group, buffer_size)
}

#' @rdname matrix_io
#' @inheritParams open_fragments_hdf5
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
    buffer_size = "integer"
  ),
  prototype = list(
    path = character(0),
    buffer_size = integer(0)
  )
)
setMethod("matrix_type", "10xMatrixH5", function(x) "uint32_t")
setMethod("matrix_inputs", "10xMatrixH5", function(x) list())
setMethod("iterate_matrix", "10xMatrixH5", function(x) {
  iterate_matrix_10x_hdf5_cpp(x@path, x@buffer_size)
})
setMethod("short_description", "10xMatrixH5", function(x) {
  sprintf("10x HDF5 feature matrix in file %s", x@path)
})

#' Read/write a 10x feature matrix
#'
#' @inheritParams open_matrix_hdf5
#' @param feature_type Optional selection of feature types to include in output matrix.
#'    For multiome data, the options are "Gene Expression" and "Peaks". This option is
#'    only compatible
#' @return BPCells matrix object
#' @details The 10x format makes use of gzip compression for the matrix data,
#' which can slow down read performance. Consider writing into another format
#' if the read performance is important to you.
#' @export
open_matrix_10x_hdf5 <- function(path, feature_type = NULL, buffer_size = 16384L) {
  assert_is_file(path)
  assert_is(buffer_size, "integer")
  if (!is.null(feature_type)) assert_is(feature_type, "character")

  path <- normalizePath(path, mustWork = FALSE)
  info <- dims_matrix_10x_hdf5_cpp(path, buffer_size)
  res <- new("10xMatrixH5",
    path = path, dim = info$dims, buffer_size = buffer_size,
    dimnames = normalized_dimnames(info$row_names, info$col_names)
  )

  if (!is.null(feature_type)) {
    valid_features <- read_hdf5_string_cpp(path, "matrix/features/feature_type", buffer_size) %in% feature_type # nolint
    res <- res[valid_features, ]
  }

  return(res)
}

#' @rdname open_matrix_10x_hdf5
#' @inheritParams write_matrix_hdf5
#' @param mat IterableMatrix
#' @param barcodes Vector of names for the cells
#' @param feature_ids Vector of IDs for the features
#' @param feature_names Vector of names for the features
#' @param feature_type String or vector of feature types
#' @param feature_metadata Named list of additional metadata vectors
#' to store for each feature
#' @details Input matrices must be in column-major storage order,
#' and if the rownames and colnames are not set, names must be
#' provided for the relevant metadata parameters. Some of the
#' metadata parameters are not read by default in BPCells, but
#' it is possible to export them for use with other tools.
#' @export
write_matrix_10x_hdf5 <- function(mat,
                                  path,
                                  barcodes = colnames(mat),
                                  feature_ids = rownames(mat),
                                  feature_names = rownames(mat),
                                  feature_types = "Gene Expression",
                                  feature_metadata = list(),
                                  buffer_size = 16384L,
                                  chunk_size = 1024L) {
  assert_is(mat, "IterableMatrix")
  assert_is(path, "character")
  if (mat@transpose) {
    stop("Matrix must have column-major storage order.\nCall t() or transpose_storage_order() first.")
  }
  if (matrix_type(mat) != "uint32_t") {
    warning("Converting to integer matrix for output to 10x format")
    mat <- convert_matrix_type(mat, "uint32_t")
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
  }
  assert_is(buffer_size, "integer")
  assert_is(chunk_size, "integer")

  path <- normalizePath(path, mustWork = FALSE)
  it <- iterate_matrix(mat)
  write_matrix_10x_hdf5_cpp(
    it,
    path,
    barcodes,
    feature_ids,
    feature_names,
    feature_types,
    feature_metadata,
    buffer_size,
    chunk_size
  )
  open_matrix_10x_hdf5(path, buffer_size)
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
setMethod("matrix_inputs", "AnnDataMatrixH5", function(x) list())
setMethod("iterate_matrix", "AnnDataMatrixH5", function(x) {
  iterate_matrix_anndata_hdf5_cpp(x@path, x@group, x@buffer_size)
})
setMethod("short_description", "AnnDataMatrixH5", function(x) {
  sprintf(
    "AnnData HDF5 matrix in file %s, group %s",
    x@path, x@group
  )
})

#' Read AnnData matrix
#'
#' Read a sparse integer matrix from an anndata matrix in an hdf5 file.
#' @inheritParams open_matrix_hdf5
#' @return AnnDataMatrixH5 object, with cells as the columns.
#' @details Since AnnData stores RNA matrices as cells x genes, whereas BPCells
#'   stores RNA matrices as genes x cells, the returned matrix will be transposed
#'   relative to the native AnnData matrix.
#' @export
open_matrix_anndata_hdf5 <- function(path, group = "X", buffer_size = 16384L) {
  assert_is_file(path)
  assert_is(buffer_size, "integer")

  path <- normalizePath(path, mustWork = FALSE)
  info <- dims_matrix_anndata_hdf5_cpp(path, group, buffer_size)
  res <- new("AnnDataMatrixH5",
    path = path, dim = info$dims, buffer_size = buffer_size,
    group = group, dimnames = normalized_dimnames(info$row_names, info$col_names),
    transpose = info$transpose
  )

  # We do the reverse of what transpose says, because anndata files are usually
  # stored row-major with cells as rows, whereas BPCells will work best with
  # col-major with cells as cols.
  if (info[["transpose"]]) {
    return(t(res))
  } else {
    return(res)
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
#' @param ranges GRanges object with the ranges to overlap, or list/data frame with columns chr, start, & end.
#' @param mode Mode for counting peak overlaps. (See "value" section for more details)
#' @inheritParams convert_to_fragments
#' @param explicit_peak_names Boolean for whether to add rownames to the output matrix in format e.g
#'  chr1:500-1000, where start and end coords are given in a 0-based coordinate system.
#'  Note that either way, peak names will be written when the matrix is saved.
#' @note When calculating the matrix directly from a fragments tsv, it's necessary to first call `select_chromosomes` in order to
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
setMethod("matrix_inputs", "TileMatrix", function(x) list())


#' Calculate ranges x cells tile overlap matrix
#' @param fragments Input fragments object
#' @param ranges GRanges object with the ranges to overlap including a metadata column tile_width,
#'  or a list/data frame with columns chr, start, end, and tile_width. Must be non-overlapping and sorted by
#'  (chr, start), with chromosomes ordered according to the chromosome names of `fragments`
#' @inheritParams convert_to_fragments
#' @param explicit_tile_names Boolean for whether to add rownames to the output matrix in format e.g
#'  chr1:500-1000, where start and end coords are given in a 0-based coordinate system. For
#'  whole-genome Tile matrices the names will take ~5 seconds to generate and take up 400MB of memory.
#'  Note that either way, tile names will be written when the matrix is saved.
#' @note When calculating the matrix directly from a fragments tsv, it's necessary to first call `select_chromosomes` in order to
#'     provide the ordering of chromosomes to expect while reading the tsv.
#' @return Iterable matrix object with dimension ranges x cells. When saved,
#'   the column names will be in the format chr1:500-1000,
#'   where start and end coords are given in a 0-based coordinate system.
#' @export
tile_matrix <- function(fragments, ranges, zero_based_coords = !is(ranges, "GRanges"), explicit_tile_names = FALSE) {
  assert_is(fragments, "IterableFragments")
  ranges <- normalize_ranges(ranges, metadata_cols = "tile_width", zero_based_coords = zero_based_coords)

  assert_is(zero_based_coords, "logical")

  assert_not_null(cellNames(fragments))
  assert_not_na(cellNames(fragments))
  assert_not_null(chrNames(fragments))
  assert_not_na(chrNames(fragments))

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
  res <- new("TileMatrix", fragments = fragments, chr_id = chr_id, start = ranges$start, end = ranges$end, tile_width = ranges$tile_width)
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
tile_ranges <- function(tile_matrix, selection) {
  # Handle manually transposed tile_matrix objects
  if (!tile_matrix@transpose) {
    tile_matrix <- tile_matrix@transpose
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
  iterate_tile_matrix_cpp(it, x@chr_id, x@start, x@end, x@tile_width, x@chr_levels)
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
#' @export
convert_matrix_type <- function(matrix, type = c("uint32_t", "double", "float")) {
  assert_is(matrix, c("dgCMatrix", "IterableMatrix"))
  type <- match.arg(type)
  if (is(matrix, "dgCMatrix")) {
    matrix <- as(matrix, "IterableMatrix")
  }
  if (matrix_type(matrix) == type) {
    return(matrix)
  } else if (matrix@transpose) {
    return(t(convert_matrix_type(t(matrix), type)))
  }
  wrapMatrix("ConvertMatrixType", matrix, type = type)
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
  iter <- iterate_matrix(convert_matrix_type(from, "double"))
  res <- build_csparse_matrix_double_cpp(iter)
  if (from@transpose) {
    res <- t(res)
  }
  res@Dimnames <- from@dimnames
  return(res)
})

# Add conversion to base R dense matrices
setAs("IterableMatrix", "matrix", function(from) {
  rlang::inform(c(
      "Warning: Converting to a dense matrix may use excessive memory"
    ), .frequency = "regularly", .frequency_id = "matrix_dense_conversion")
  as(from, "dgCMatrix") %>% as.matrix()  
})

#' @exportS3Method base::as.matrix
as.matrix.IterableMatrix <- function(x, ...) as(x, "matrix")

#' @export
setMethod("as.matrix", signature(x = "IterableMatrix"), function(x, ...) as(x, "matrix"))

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
matrix_stats <- function(matrix,
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
