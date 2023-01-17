

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

#################
# Log1p and Expm1
#################

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


setClass("TransformExpm1", contains = "TransformedMatrix")
setMethod("iterate_matrix", "TransformExpm1", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_expm1simd_cpp(ptr(it)), it)
})
setMethod("short_description", "TransformExpm1", function(x) {
  c(
    short_description(x@matrix),
    "Transform expm1"
  )
})
setMethod("expm1", "IterableMatrix", function(x) {
  wrapMatrix("TransformExpm1", convert_matrix_type(x, "double"))
})

setClass("TransformExpm1Slow", contains = "TransformedMatrix")
setMethod("iterate_matrix", "TransformExpm1Slow", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_expm1_cpp(ptr(it)), it)
})
setMethod("short_description", "TransformExpm1Slow", function(x) {
  c(
    short_description(x@matrix),
    "Transform expm1 (non-SIMD implementation)"
  )
})
expm1_slow <- function(x) {
  wrapMatrix("TransformExpm1Slow", convert_matrix_type(x, "double"))
}


#################
# Pow and Square
#################

setClass("TransformSquare", contains="TransformedMatrix")
setMethod("iterate_matrix", "TransformSquare", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_square_cpp(ptr(it)), it)
})
setMethod("short_description", "TransformSquare", function(x) {
  c(
    short_description(x@matrix),
    "Square elements"
  )
})

setClass("TransformPow", contains="TransformedMatrix")
setMethod("iterate_matrix", "TransformPow", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_powsimd_cpp(ptr(it), x@global_params[1]), it)
})
setMethod("short_description", "TransformPow", function(x) {
  c(
    short_description(x@matrix),
    sprintf("Raise elements to the power of %.2g", x@global_params[1])
  )
})

setMethod("^", signature(e1 = "IterableMatrix", e2 = "numeric"), function(e1, e2) {
  assert_len(e2, 1)
  assert_true(e2 != 0)
  if (e2 == 2) {
    wrapMatrix("TransformSquare", convert_matrix_type(e1, "double"))
  } else {
    wrapMatrix("TransformPow", convert_matrix_type(e1, "double"), global_params=e2)
  }
})

setClass("TransformPowSlow", contains="TransformedMatrix")
setMethod("iterate_matrix", "TransformPowSlow", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_pow_cpp(ptr(it), x@global_params[1]), it)
})
setMethod("short_description", "TransformPowSlow", function(x) {
  c(
    short_description(x@matrix),
    sprintf("Raise elements to the power of %.2g (non-SIMD implementation)", x@global_params[1])
  )
})

pow_slow <- function(x, exponent) {
  wrapMatrix("TransformPowSlow", convert_matrix_type(x, "double"), global_params=exponent)
}

#################
# Min
#################

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

#' Elementwise minimum
#'
#' Take the minimum value of a matrix with a per-row, per-col, or global
#' constant. This constant must be >0 to preserve sparsity of the matrix.
#' This has the effect of capping the maximum value in the matrix.
#'
#' @param mat IterableMatrix
#' @param val Single positive numeric value
#' @return IterableMatrix
#' @description **min_scalar**: Take minumum with a global constant
#' @rdname min_elementwise
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

setClass("TransformMinByRow", contains = "TransformedMatrix")
setMethod("iterate_matrix", "TransformMinByRow", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_min_by_row_cpp(ptr(it), x@row_params), it)
})
setMethod("short_description", "TransformMinByRow", function(x) {
  # Subset the row + col params matrices for faster pretty printing of
  # large parameter sets
  print_entries <- 3
  if (ncol(x@row_params) > print_entries + 1) {
    x@row_params <- x@row_params[,c(1:print_entries, ncol(x@row_params)),drop=FALSE]
  }
  
  c(
    short_description(x@matrix),
    sprintf("Transform min by row: %s", pretty_print_vector(sprintf("%.3g", x@row_params[1, ]), max_len = 3))
  )
})

#' @rdname min_elementwise
#' @description **min_by_row**: Take the minimum with a per-row constant
#' @export
min_by_row <- function(mat, vals) {
  assert_is(mat, "IterableMatrix")
  assert_is(vals, "numeric")
  assert_greater_than_zero(vals)
  assert_len(vals, nrow(mat))

  wrapMatrix("TransformMinByRow", convert_matrix_type(mat, "double"), row_params=matrix(vals, nrow=1))
}

setClass("TransformMinByCol", contains = "TransformedMatrix")
setMethod("iterate_matrix", "TransformMinByCol", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_min_by_col_cpp(ptr(it), x@col_params), it)
})
setMethod("short_description", "TransformMinByCol", function(x) {
  # Subset the row + col params matrices for faster pretty printing of
  # large parameter sets
  print_entries <- 3
  if (ncol(x@col_params) > print_entries + 1) {
    x@col_params <- x@col_params[,c(1:print_entries, ncol(x@col_params)), drop=FALSE]
  }
  
  c(
    short_description(x@matrix),
    sprintf("Transform min by col: %s", pretty_print_vector(sprintf("%.3g", x@col_params[1, ]), max_len = 3))
  )
})

#' @rdname min_elementwise
#' @description **min_by_col**: Take the minimum with a per-col constant
#' @export
min_by_col <- function(mat, vals) {
  assert_is(mat, "IterableMatrix")
  assert_is(vals, "numeric")
  assert_greater_than_zero(vals)
  assert_len(vals, ncol(mat))

  wrapMatrix("TransformMinByCol", convert_matrix_type(mat, "double"), col_params=matrix(vals, nrow=1))
}


#################
# SCTransform
#################

setClass("SCTransformPearson", contains = "TransformedMatrix")
setMethod("iterate_matrix", "SCTransformPearson", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_sctransform_pearson_simd_cpp(ptr(it), x@row_params, x@col_params, x@global_params), it)
})
setMethod("short_description", "SCTransformPearson", function(x) {
  c(
    short_description(x@matrix),
    "SCTransform (Pearson Residuals)"
  )
})

setClass("SCTransformPearsonTranspose", contains = "TransformedMatrix")
setMethod("iterate_matrix", "SCTransformPearsonTranspose", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_sctransform_pearson_transpose_simd_cpp(ptr(it), x@row_params, x@col_params, x@global_params), it)
})
setMethod("short_description", "SCTransformPearsonTranspose", function(x) {
  c(
    short_description(x@matrix),
    "SCTransform (Pearson Residuals with transposed orientation)"
  )
})


setClass("SCTransformPearsonSlow", contains = "TransformedMatrix")
setMethod("iterate_matrix", "SCTransformPearsonSlow", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_sctransform_pearson_cpp(ptr(it), x@row_params, x@col_params, x@global_params), it)
})
setMethod("short_description", "SCTransformPearsonSlow", function(x) {
  c(
    short_description(x@matrix),
    "SCTransform slow (Pearson Residuals)"
  )
})

setClass("SCTransformPearsonTransposeSlow", contains = "TransformedMatrix")
setMethod("iterate_matrix", "SCTransformPearsonTransposeSlow", function(x) {
  it <- iterate_matrix(x@matrix)
  wrapMat_double(iterate_matrix_sctransform_pearson_transpose_cpp(ptr(it), x@row_params, x@col_params, x@global_params), it)
})
setMethod("short_description", "SCTransformPearsonTransposeSlow", function(x) {
  c(
    short_description(x@matrix),
    "SCTransform slow (Pearson Residuals with transposed orientation)"
  )
})

#' SCTransform Pearson Residuals
#'
#' Calculate pearson residuals of a negative binomial sctransform model.
#' Normalized values are calculated as `(X - mu) / sqrt(mu + mu^2/theta)`.
#' mu is calculated as `cell_read_counts * gene_beta`.
#' 
#' The parameterization used is somewhat simplified compared to the original 
#' SCTransform paper, in particular it uses a linear-scale rather than
#' log-scale to represent the cell_read_counts and gene_beta variables. It also
#' does not support the addition of arbitrary cell metadata (e.g. batch) to add to the 
#' negative binomial regression.
#'
#' @param mat IterableMatrix (raw counts)
#' @param gene_theta Vector of per-gene thetas (overdispersion values)
#' @param gene_beta Vector of per-gene betas (expression level values)
#' @param cell_read_counts Vector of total reads per (umi count for RNA)
#' @param min_var Minimum value for clipping variance
#' @param clip_range Length 2 vector of min and max clipping range
#' @param columns_are_cells Whether the columns of the matrix correspond to cells (default) or genes
#' @param slow If TRUE, use a 10x slower but more precise implementation (default FALSE)
#' @return IterableMatrix
#' @export
sctransform_pearson <- function(mat, gene_theta, gene_beta, cell_read_counts, min_var = -Inf, clip_range = c(-10, 10), columns_are_cells=TRUE, slow=FALSE) {
  assert_is(mat, "IterableMatrix")
  assert_is_numeric(gene_theta)
  assert_is_numeric(gene_beta)
  assert_is_numeric(cell_read_counts)
  assert_is_numeric(min_var)
  assert_is_numeric(clip_range)
  assert_is(columns_are_cells, "logical")

  # Check dimensions
  if (columns_are_cells) {
    assert_true(length(gene_theta) == 1 || length(gene_theta) == nrow(mat))
    assert_len(gene_beta, nrow(mat))
    assert_len(cell_read_counts, ncol(mat))
  } else {
    assert_true(length(gene_theta) == 1 || length(gene_theta) == ncol(mat))
    assert_len(gene_beta, ncol(mat))
    assert_len(cell_read_counts, nrow(mat))
  }
  if (!is.null(clip_range)) {
    assert_len(clip_range, 2)
  } else {
    clip_range <- c(-Inf, Inf)
  }
  if (!is.null(min_var)) {
    assert_len(min_var, 1)
  } else {
    min_var <- -Inf
  }

  # Re-scale gene_beta and cell_read_counts to be similar magnitudes
  ratio <- exp(0.5*(mean(log(gene_beta)) - mean(log(cell_read_counts))))
  gene_beta <- gene_beta / ratio
  cell_read_counts <- cell_read_counts * ratio

  # Determine which implementation to use in the backend
  matrix_class <- sprintf(
    "SCTransformPearson%s%s", 
    ifelse(mat@transpose == columns_are_cells, "Transpose", ""), 
    ifelse(slow, "Slow", "")
  )
  
  col_params <- matrix(cell_read_counts, nrow=1)
  row_params <- rbind(1/gene_theta, t(gene_beta))

  if (mat@transpose == columns_are_cells) {
    # Transposed orientation
    tmp <- row_params
    row_params <- col_params
    col_params <- tmp
    rm(tmp)
  }
  wrapMatrix(
    matrix_class, 
    convert_matrix_type(mat, "double"),
    row_params = row_params,
    col_params = col_params,
    global_params = c(1/sqrt(min_var), clip_range)
  )
}



#################
# Scale + Shift
#################

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



#################
# Arithmetic helpers
#################

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
