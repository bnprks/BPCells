

setClass("TransformedMatrix",
  contains = "IterableMatrix",
  slots = c(
    matrix = "IterableMatrix",
    row_params = "matrix",
    col_params = "matrix",
    global_params = "numeric"
  ),
  prototype = list(
    row_params = matrix(0, 0, 0),
    col_params = matrix(0, 0, 0),
    global_params = numeric(0)
  )
)
setMethod("matrix_type", "TransformedMatrix", function(x) "double")

# Subsetting on TransformedMatrix objects
setMethod("[", "TransformedMatrix", function(x, i, j, ...) {
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

  # Subset the row/col params
  if (!rlang::is_missing(i) && ncol(x@row_params) != 0) {
    x@row_params <- x@row_params[, i, drop = FALSE]
  }
  if (!rlang::is_missing(j) && ncol(x@col_params) != 0) {
    x@col_params <- x@col_params[, j, drop = FALSE]
  }
  x
})


# log1p method support. The SIMD method may be ever so slightly less precise,
# but it can be substantially faster depending on the CPU SIMD features
# (Should still provide 32-bit float accuracy)
setClass("TransformLog1p", contains = "TransformedMatrix")
setMethod("iterate_matrix", "TransformLog1p", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_log1psimd_cpp(ptr(it)), it)
})
setMethod("short_description", "TransformLog1p", function(x) {
  c(
    short_description(x@matrix),
    "Transform log1p"
  )
})
setMethod("log1p", "IterableMatrix", function(x) {
  wrapMatrix("TransformLog1p", convert_matrix_type(x, "double"))
})

setClass("TransformLog1pSlow", contains = "TransformedMatrix")
setMethod("iterate_matrix", "TransformLog1pSlow", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_log1p_cpp(ptr(it)), it)
})
setMethod("short_description", "TransformLog1pSlow", function(x) {
  c(
    short_description(x@matrix),
    "Transform log1p (non-SIMD implementation)"
  )
})
log1p_slow <- function(x) {
  wrapMatrix("TransformLog1pSlow", convert_matrix_type(x, "double"))
}

setClass("TransformMin", contains = "TransformedMatrix")
setMethod("iterate_matrix", "TransformMin", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_min_cpp(ptr(it), x@global_params[1]), it)
})
setMethod("short_description", "TransformMin", function(x) {
  c(
    short_description(x@matrix),
    sprintf("Transform min(x, %d)", x@global_params[1])
  )
})
#' Take the elementwise min with a constant
#' @param mat IterableMatrix
#' @param val Single positive numeric value
#' @return IterableMatrix
#' @description Take the elementwise minimum with a positive value.
#' This has the effect of capping the maximum value in the matrix
#' @export
min_scalar <- function(mat, val) {
  assert_is(mat, "IterableMatrix")
  assert_is(val, "numeric")
  assert_greater_than_zero(val)
  assert_len(val, 1)

  res <- wrapMatrix("TransformMin", convert_matrix_type(mat, "double"))
  res@global_params <- val
  res
}

# Scaling + shifting support (Scale first, then shift)
#
# To provide greatest ease of use, this class will try its hardest to compose operations with itself
# such that multiple +,-,*,/ operations can be coalesced into as few transforms as possible.
#
# Slot active_transforms: A 3x2 logical matrix for whether certain transformations active -
#   rows are Row,Col,Global(scalar), cols are Scale,Shift. e.g. element (0,0) says if there is
#   an active RowScale operation
setClass("TransformScaleShift",
  contains = "TransformedMatrix",
  slots = c(
    active_transforms = "matrix"
  ),
  prototype = list(
    active_transforms = matrix(FALSE, nrow = 3, ncol = 2, dimnames = list(c("row", "col", "global"), c("scale", "shift")))
  )
)
setMethod("iterate_matrix", "TransformScaleShift", function(x) {
  res <- iterate_matrix(x@matrix)
  if (any(x@active_transforms[, "scale"])) {
    scale_row <- matrix(0, 0, 0)
    scale_col <- matrix(0, 0, 0)
    if (x@active_transforms["row", "scale"]) scale_row <- x@row_params[1, , drop = FALSE]
    if (x@active_transforms["col", "scale"]) scale_col <- x@col_params[1, , drop = FALSE]
    if (x@active_transforms["global", "scale"]) {
      if (x@active_transforms["row", "scale"]) {
        scale_row <- scale_row * x@global_params[1]
      } else if (x@active_transforms["col", "scale"]) {
        scale_col <- scale_col * x@global_params[1]
      } else {
        scale_row <- matrix(x@global_params[1], nrow = 1, ncol = nrow(x))
      }
    }
    res <- wrapMat_double(iterate_matrix_scale_cpp(ptr(res), scale_row, scale_col), res)
  }
  if (any(x@active_transforms[, "shift"])) {
    shift_row <- matrix(0, 0, 0)
    shift_col <- matrix(0, 0, 0)
    if (x@active_transforms["row", "shift"]) shift_row <- x@row_params[2, , drop = FALSE]
    if (x@active_transforms["col", "shift"]) shift_col <- x@col_params[2, , drop = FALSE]
    if (x@active_transforms["global", "shift"]) {
      if (x@active_transforms["row", "shift"]) {
        shift_row <- shift_row + x@global_params[2]
      } else if (x@active_transforms["col", "shift"]) {
        shift_col <- shift_col + x@global_params[2]
      } else {
        shift_row <- matrix(x@global_params[2], nrow = 1, ncol = nrow(x))
      }
    }
    if (nrow(shift_row) != 0) res <- wrapMat_double(iterate_matrix_row_shift_cpp(ptr(res), shift_row), res)
    if (nrow(shift_col) != 0) res <- wrapMat_double(iterate_matrix_col_shift_cpp(ptr(res), shift_col), res)
  }
  res
})
setMethod("short_description", "TransformScaleShift", function(x) {
  # Return multiple lines, one for each transform active
  res <- short_description(x@matrix)

  # Subset the row + col params matrices for faster pretty printing of
  # large parameter sets
  print_entries <- 3
  if (ncol(x@row_params) > print_entries + 1) {
    x@row_params <- x@row_params[,c(1:print_entries, ncol(x@row_params))]
  }
  if (ncol(x@col_params) > print_entries + 1) {
    x@col_params <- x@col_params[,c(1:print_entries, ncol(x@col_params))]
  }

  # Handle scale transforms
  if (x@active_transforms["global", "scale"]) {
    res <- c(res, sprintf("Scale by %.3g", x@global_params[1]))
  }
  if (x@active_transforms["row", "scale"]) {
    res <- c(res, sprintf(
      "Scale %s by %s",
      if (x@transpose) "columns" else "rows",
      pretty_print_vector(sprintf("%.3g", x@row_params[1, ]), max_len = 3)
    ))
  }
  if (x@active_transforms["col", "scale"]) {
    res <- c(res, sprintf(
      "Scale %s by %s",
      if (x@transpose) "rows" else "columns",
      pretty_print_vector(sprintf("%.3g", x@col_params[1, ]), max_len = 3)
    ))
  }

  # Handle shift transforms
  if (x@active_transforms["global", "shift"]) {
    res <- c(res, sprintf("Shift by %.3g", x@global_params[2]))
  }
  if (x@active_transforms["row", "shift"]) {
    res <- c(res, sprintf(
      "Shift %s by %s",
      if (x@transpose) "columns" else "rows",
      pretty_print_vector(sprintf("%.3g", x@row_params[2, ]), max_len = 3)
    ))
  }
  if (x@active_transforms["col", "shift"]) {
    res <- c(res, sprintf(
      "Shift %s by %s",
      if (x@transpose) "rows" else "columns",
      pretty_print_vector(sprintf("%.3g", x@col_params[2, ]), max_len = 3)
    ))
  }
  res
})

# Basic dispatch for scaling/shifting (Create TransformScaleShift and then apply function to it)
setMethod("*", signature(e1 = "IterableMatrix", e2 = "numeric"), function(e1, e2) {
  e1 <- wrapMatrix("TransformScaleShift", convert_matrix_type(e1, "double"))
  e1 * e2
})
setMethod("*", signature(e1 = "numeric", e2 = "IterableMatrix"), function(e1, e2) {
  e2 <- wrapMatrix("TransformScaleShift", convert_matrix_type(e2, "double"))
  e2 * e1
})
setMethod("+", signature(e1 = "IterableMatrix", e2 = "numeric"), function(e1, e2) {
  e1 <- wrapMatrix("TransformScaleShift", convert_matrix_type(e1, "double"))
  e1 + e2
})
setMethod("+", signature(e1 = "numeric", e2 = "IterableMatrix"), function(e1, e2) {
  e2 <- wrapMatrix("TransformScaleShift", convert_matrix_type(e2, "double"))
  e2 + e1
})
# Note: we skip numeric / IterableMatrix as it would result in a lot of infinities for dividing by 0.
setMethod("/", signature(e1 = "IterableMatrix", e2 = "numeric"), function(e1, e2) {
  e1 * (1 / e2)
})
setMethod("-", signature(e1 = "IterableMatrix", e2 = "numeric"), function(e1, e2) {
  e1 + (-e2)
})
setMethod("-", signature(e1 = "numeric", e2 = "IterableMatrix"), function(e1, e2) {
  e2 * -1 + e1
})

# Full dispatch for scaling/shifting
setMethod("*", signature(e1 = "TransformScaleShift", e2 = "numeric"), function(e1, e2) {
  # Convenience renaming - x is matrix, y is vector/scalar
  x <- e1
  y <- e2
  # 1. Error checking - dimensions match
  # Note that since x@dim is swapped upon transpose, the relevant dimension is
  # always nrow(x) regardless of whether x@transpose is true
  assert_true(length(y) == 1 || length(y) == nrow(x))

  # 2. Handle multiplying by a scalar
  if (length(y) == 1) {
    # Initialize global_params if necessary
    if (!any(x@active_transforms["global", ])) x@global_params <- c(1, 0)
    x@global_params[1] <- y * x@global_params[1]
    x@global_params[2] <- x@global_params[2] * y
    x@active_transforms["global", "scale"] <- TRUE
    return(x)
  }
  # 3. Check the transform: if we're trying to scale row/col after having set a shift on col/row, then make new layer
  if (x@transpose && x@active_transforms["row", "shift"]) {
    res <- wrapMatrix("TransformScaleShift", x)
    res@active_transforms["col", "scale"] <- TRUE
    res@col_params <- matrix(c(1, 0), nrow = 2, ncol = nrow(x)) # Note since x is transposed ncol has the underlying row count
    res@col_params[1, ] <- y
    return(res)
  }
  if (!x@transpose && x@active_transforms["col", "shift"]) {
    res <- wrapMatrix("TransformScaleShift", x)
    res@active_transforms["col", "scale"] <- TRUE
    res@row_params <- matrix(c(1, 0), nrow = 2, ncol = nrow(x))
    res@row_params[1, ] <- y
    return(res)
  }
  # 4. Otherwise, update existing parameters (multiply the shift+scale appropriately)
  if (x@transpose) {
    # Initialize col_params if necessary
    if (!any(x@active_transforms["col", ])) x@col_params <- matrix(c(1, 0), nrow = 2, ncol = nrow(x))

    # Update scale
    x@col_params[1, ] <- x@col_params[1, ] * y
    x@active_transforms["col", "scale"] <- TRUE

    # Update shift
    x@col_params[2, ] <- x@col_params[2, ] * y
  } else {
    # Initialize row_params if necessary
    if (!any(x@active_transforms["row", ])) x@row_params <- matrix(c(1, 0), nrow = 2, ncol = nrow(x))

    # Update scale
    x@row_params[1, ] <- x@row_params[1, ] * y
    x@active_transforms["row", "scale"] <- TRUE

    # Update shift
    x@row_params[2, ] <- x@row_params[2, ] * y
  }
  return(x)
})
setMethod("+", signature(e1 = "TransformScaleShift", e2 = "numeric"), function(e1, e2) {
  # Convenience renaming - x is matrix, y is vector/scalar
  x <- e1
  y <- e2
  # 1. Error checking - dimensions match
  # Note that since x@dim is swapped upon transpose, the relevant dimension is
  # always nrow(x) regardless of whether x@transpose is true
  assert_true(length(y) == 1 || length(y) == nrow(x))

  # 2. Handle shifting by a scalar
  if (length(y) == 1) {
    # Initialize global_params if necessary
    if (!any(x@active_transforms["global", ])) x@global_params <- c(1, 0)
    x@global_params[2] <- y + x@global_params[2]
    x@active_transforms["global", "shift"] <- TRUE
    return(x)
  }

  # 3. Otherwise, update existing parameters
  if (x@transpose) {
    # Initialize col_params if necessary
    if (!any(x@active_transforms["col", ])) x@col_params <- matrix(c(1, 0), nrow = 2, ncol = nrow(x))

    # Update shift
    x@col_params[2, ] <- x@col_params[2, ] + y
    x@active_transforms["col", "shift"] <- TRUE
  } else {
    # Initialize row_params if necessary
    if (!any(x@active_transforms["row", ])) x@row_params <- matrix(c(1, 0), nrow = 2, ncol = nrow(x))

    # Update scale
    x@row_params[2, ] <- x@row_params[2, ] + y
    x@active_transforms["row", "shift"] <- TRUE
  }
  return(x)
})
# Just take advantage of commutative property to only implement half
setMethod("*", signature(e1 = "numeric", e2 = "TransformScaleShift"), function(e1, e2) {
  e2 * e1
})
setMethod("+", signature(e1 = "numeric", e2 = "TransformScaleShift"), function(e1, e2) {
  e2 + e1
})



#' Broadcasting vector arithmetic
#'
#' Convenience functions for adding or multiplying
#' each row / column of a mtarix by a number.
#'
#' @rdname mat_norm
#'
#' @param mat Matrix-like object
#' @param vec Numeric vector
#' @return Matrix-like object
add_rows <- function(mat, vec) {
  assert_is(mat, c("dgCMatrix", "IterableMatrix", "matrix"))
  assert_is_numeric(vec)
  assert_true(length(vec) == nrow(mat))
  mat + vec
}
#' @rdname mat_norm
#' @export
add_cols <- function(mat, vec) {
  assert_is(mat, c("dgCMatrix", "IterableMatrix", "matrix"))
  assert_is_numeric(vec)
  assert_true(length(vec) == ncol(mat))
  t(t(mat) + vec)
}
#' @rdname mat_norm
#' @export
multiply_rows <- function(mat, vec) {
  assert_is(mat, c("dgCMatrix", "IterableMatrix", "matrix"))
  assert_is_numeric(vec)
  assert_true(length(vec) == nrow(mat))
  mat * vec
}
#' @rdname mat_norm
#' @export
multiply_cols <- function(mat, vec) {
  assert_is(mat, c("dgCMatrix", "IterableMatrix", "matrix"))
  assert_is_numeric(vec)
  assert_true(length(vec) == ncol(mat))
  t(t(mat) * vec)
}
