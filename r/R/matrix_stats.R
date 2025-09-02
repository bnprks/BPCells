# Copyright 2024 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

#' Find the nth quantile value(s) of each row in a matrix. Only supports transposed matrices.
#' @param x IterableMatrix object or a matrix-like object.
#' @param probs (Numeric) Quantile value(s) to be computed, between 0 and 1.
#' @param type (Integer) between 4 and 9 selecting which quantile algorithm to use, detailed in `matrixStats::rowQuantiles()`
#' @return  - `rowQuantiles():` If `length(probs) == 1`, return a numeric with number of entries equal to the number of rows in the matrix.
#' Else, return a Matrix of quantile values, with cols representing each quantile, and each row representing a row in the input matrix.
#' @describeIn IterableMatrix-methods Calculate rowQuantiles (replacement for `matrixStats::rowQuantiles`)
#' @usage rowQuantiles(
#'   x,
#'   rows = NULL,
#'   cols = NULL,
#'   probs = seq(from = 0, to = 1, by = 0.25),
#'   na.rm = FALSE,
#'   type = 7L,
#'   digits = 7L,
#'   ...,
#'   useNames = TRUE,
#'   drop = TRUE
#' )
#' @examples
#' #######################################################################
#' ## rowQuantiles() example
#' #######################################################################
#' rowQuantiles(transpose_storage_order(mat))
#' 
#' 
#' @export
rowQuantiles <- function(x, rows = NULL, cols = NULL,
                         probs = seq(from = 0, to = 1, by = 0.25),
                         na.rm = FALSE, type = 7L, digits = 7L, ...,
                         useNames = TRUE, drop = TRUE) {
  UseMethod("rowQuantiles")
}
#' @export
rowQuantiles.default <- function(x, rows = NULL, cols = NULL,
                                 probs = seq(from = 0, to = 1, by = 0.25),
                                 na.rm = FALSE, type = 7L, digits = 7L, ...,
                                 useNames = TRUE, drop = TRUE) {
  if (requireNamespace("MatrixGenerics", quietly = TRUE)) {
    MatrixGenerics::rowQuantiles(x, rows = rows, cols = cols, probs = probs, na.rm = na.rm, type = type, digits = digits, ..., useNames = useNames, drop = drop)
  } else if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowQuantiles(x, rows = rows, cols = cols, probs = probs, na.rm = na.rm, type = type, digits = digits, ..., useNames = useNames, drop = drop)
  } else {
    rlang::abort("Cannot run rowQuantiles on a non-BPCells object unless MatrixGenerics or matrixStats is installed.")
  }
}
#' @export
rowQuantiles.IterableMatrix <- function(x, rows = NULL, cols = NULL,
                                        probs = seq(from = 0, to = 1, by = 0.25),
                                        na.rm = FALSE, type = 7L, digits = 7L, ...,
                                        useNames = TRUE, drop = TRUE) {
  if (!is.null(rows) || !is.null(cols) || isTRUE(na.rm) || isFALSE(useNames) || isFALSE(drop)) {
    rlang::abort("rowQuantiles(IterableMatrix) doesn't support extra arguments rows, cols, na.rm, useNames, or drop.")
  }
  if (!x@transpose) {
    rlang::abort("rowQuantiles(IterableMatrix) only supports row-major matrices. Please call transpose_storage_order() first, then run rowQuantiles.")
  }
  assert_is(x, "IterableMatrix")
  assert_is(probs, "numeric")
  assert_true(all(probs >= 0 & probs <= 1))
  assert_is(type, c("integer", "numeric"))
  assert_true(length(type) == 1)
  assert_true(type <= 9 || type >= 4)
  type <- as.integer(type)
  # Convert type to alpha and beta values
  alpha_beta <- switch(as.character(type),
    `4` = c(0, 1),
    `5` = c(0.5, 0.5),
    `6` = c(0, 0),
    `7` = c(1, 1),
    `8` = c(1 / 3, 1 / 3),
    `9` = c(3 / 8, 3 / 8)
  )
  alpha <- alpha_beta[1]
  beta <- alpha_beta[2]
  iter <- iterate_matrix(convert_matrix_type(x, "double"))
  res <- matrix_quantile_per_col_cpp(iter, probs, alpha = alpha, beta = beta)
  rownames(res) <- rownames(x)
  if (length(probs) == 1) return(res[,1])
  # `quantile()` from base R returns rownames as percentages, so we follow that convention
  colnames(res) <- paste0(format(100 * probs, trim = TRUE, digits = digits), "%")
  return(res)
}


#' Find the nth quantile value(s) of each column in a matrix. Only supports non-transposed matrices.
#' @return - `colQuantiles():` If `length(probs) == 1`, return a numeric with number of entries equal to the number of columns in the matrix. 
#' Else, return a Matrix of quantile values, with cols representing each quantile, and each row representing a col in the input matrix.
#' @describeIn IterableMatrix-methods Calculate colQuantiles (replacement for `matrixStats::colQuantiles`)
#' @inheritParams rowQuantiles
#' @usage colQuantiles(
#'   x,
#'   rows = NULL,
#'   cols = NULL,
#'   probs = seq(from = 0, to = 1, by = 0.25),
#'   na.rm = FALSE,
#'   type = 7L,
#'   digits = 7L,
#'   ...,
#'   useNames = TRUE,
#'   drop = TRUE
#' )
#' @examples
#' #######################################################################
#' ## colQuantiles() example
#' #######################################################################
#' colQuantiles(mat)
#' 
#' 
#' @export
colQuantiles <- function(x, rows = NULL, cols = NULL,
                         probs = seq(from = 0, to = 1, by = 0.25),
                         na.rm = FALSE, type = 7L, digits = 7L, ...,
                         useNames = TRUE, drop = TRUE) {
  UseMethod("colQuantiles")
}
#' @export
colQuantiles.default <- function(x, rows = NULL, cols = NULL, 
                                 probs = seq(from = 0, to = 1, by = 0.25),
                                 na.rm = FALSE, type = 7L, digits = 7L, ...,
                                 useNames = TRUE, drop = TRUE) {
  if (requireNamespace("MatrixGenerics", quietly = TRUE)) {
    MatrixGenerics::colQuantiles(x, rows = rows, cols = cols, probs = probs, na.rm = na.rm, type = type, digits = digits, ..., useNames = useNames, drop = drop)
  } else if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colQuantiles(x, rows = rows, cols = cols, probs = probs, na.rm = na.rm, type = type, digits = digits, ..., useNames = useNames, drop = drop)
  } else {
    rlang::abort("Cannot run colQuantiles on a non-BPCells object unless MatrixGenerics or matrixStats is installed.")
  }
}
#' @export
colQuantiles.IterableMatrix <- function(x, rows = NULL, cols = NULL, 
                                        probs = seq(from = 0, to = 1, by = 0.25), 
                                        na.rm = FALSE, type = 7L, digits = 7L, ...,
                                        useNames = TRUE, drop = TRUE) {
  if (!is.null(rows) || !is.null(cols) || isTRUE(na.rm) || isFALSE(useNames) || isFALSE(drop)) {
    rlang::abort("colQuantiles(IterableMatrix) doesn't support extra arguments rows, cols, na.rm, useNames, or drop.")
  }
  assert_is(x, "IterableMatrix")
  if (x@transpose) {
    rlang::abort("colQuantiles(IterableMatrix) does not support row-major matrices.\nPlease call transpose_storage_order() first.")
  }
  assert_is(probs, "numeric")
  assert_true(all(probs >= 0 & probs <= 1))
  assert_is(type, c("integer", "numeric"))
  assert_true(length(type) == 1)
  assert_true(type <= 9 || type >= 4)
  type <- as.integer(type)
  # Convert type to alpha and beta values
  alpha_beta <- switch(as.character(type),
    `4` = c(0, 1),
    `5` = c(0.5, 0.5),
    `6` = c(0, 0),
    `7` = c(1, 1),
    `8` = c(1 / 3, 1 / 3),
    `9` = c(3 / 8, 3 / 8)
  )
  alpha <- alpha_beta[1]
  beta <- alpha_beta[2]
  iter <- iterate_matrix(convert_matrix_type(x, "double"))
  res <- matrix_quantile_per_col_cpp(iter, probs, alpha = alpha, beta = beta)
  rownames(res) <- colnames(x)
  if (length(probs) == 1) return(res[,1])
  # `quantile()` from base R returns rownames as percentages, so we follow that convention
  colnames(res) <- paste0(format(100 * probs, trim = TRUE, digits = digits), "%")
  return(res)
}

rlang::on_load({
  if (requireNamespace("MatrixGenerics", quietly = TRUE)) {
    setMethod(MatrixGenerics::colQuantiles, "IterableMatrix", colQuantiles.IterableMatrix)
    setMethod(MatrixGenerics::rowQuantiles, "IterableMatrix", rowQuantiles.IterableMatrix)
  }
})