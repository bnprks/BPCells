# This file is distributed under the MPL-2.0 license <https://www.mozilla.org/en-US/MPL/2.0/>
# It is adapted from the RSpectra package, which is distributed under the MPL-2.0 license. 
# The RSpectra package is copyright Yixuan Qiu 2016 
# Original source code available here: https://github.com/yixuan/RSpectra
# Modifications are copyright 2024 BPCells contributors
# 
# SPDX-License-Identifier: MPL-2.0

#' Calculate svds
#'
#' Use the C++ Spectra solver (same as RSpectra package), in order to
#' compute the largest k values and corresponding singular vectors.
#' Empirically, memory usage is much lower than using `irlba::irlba()`, likely
#' due to avoiding R garbage creation while solving due to the pure-C++ solver.
#' This documentation is a slightly-edited version of the `RSpectra::svds()` 
#' documentation.
#'
#' @param A The matrix whose truncated SVD is to be computed.
#' @param k Number of singular values requested.
#' @param nu Number of left singular vectors to be computed. This must be between 0 and 'k'. (Must be equal to 'k' for BPCells IterableMatrix) 
#' @param nu Number of right singular vectors to be computed. This must be between 0 and 'k'. (Must be equal to 'k' for BPCells IterableMatrix) 
#' @param opts Control parameters related to computing algorithm. See *Details* below
#' @param threads Control threads to use calculating mat-vec producs (BPCells specific)
#' @return A list with the following components:
##' \item{d}{A vector of the computed singular values.}
##' \item{u}{An \code{m} by \code{nu} matrix whose columns contain
##'          the left singular vectors. If \code{nu == 0}, \code{NULL}
##'          will be returned.}
##' \item{v}{An \code{n} by \code{nv} matrix whose columns contain
##'          the right singular vectors. If \code{nv == 0}, \code{NULL}
##'          will be returned.}
##' \item{nconv}{Number of converged singular values.}
##' \item{niter}{Number of iterations used.}
##' \item{nops}{Number of matrix-vector multiplications used.}
#' @details
#' When RSpectra is installed, this function will just add a method to
#' `RSpectra::svds()` for the `IterableMatrix` class.
#'
#' The \code{opts} argument is a list that can supply any of the
#' following parameters:
#'
#' \describe{
#' \item{\code{ncv}}{Number of Lanzcos basis vectors to use. More vectors
#'                   will result in faster convergence, but with greater
#'                   memory use. \code{ncv} must be satisfy
#'                   \eqn{k < ncv \le p}{k < ncv <= p} where
#'                   \code{p = min(m, n)}.
#'                   Default is \code{min(p, max(2*k+1, 20))}.}
#' \item{\code{tol}}{Precision parameter. Default is 1e-10.}
#' \item{\code{maxitr}}{Maximum number of iterations. Default is 1000.}
#' \item{\code{center}}{Either a logical value (\code{TRUE}/\code{FALSE}), or a numeric
#'                      vector of length \eqn{n}. If a vector \eqn{c} is supplied, then
#'                      SVD is computed on the matrix \eqn{A - 1c'}{A - 1 * c'},
#'                      in an implicit way without actually forming this matrix.
#'                      \code{center = TRUE} has the same effect as
#'                      \code{center = colMeans(A)}. Default is \code{FALSE}. Ignored in BPCells}
#' \item{\code{scale}}{Either a logical value (\code{TRUE}/\code{FALSE}), or a numeric
#'                     vector of length \eqn{n}. If a vector \eqn{s} is supplied, then
#'                     SVD is computed on the matrix \eqn{(A - 1c')S}{(A - 1 * c')S},
#'                     where \eqn{c} is the centering vector and \eqn{S = diag(1/s)}.
#'                     If \code{scale = TRUE}, then the vector \eqn{s} is computed as
#'                     the column norm of \eqn{A - 1c'}{A - 1 * c'}.
#'                     Default is \code{FALSE}. Ignored in BPCells}
#' }
#' @references Qiu Y, Mei J (2022). _RSpectra: Solvers for Large-Scale Eigenvalue and SVD Problems_. R package version 0.16-1, <https://CRAN.R-project.org/package=RSpectra>.
#' @usage svds(A, k, nu = k, nv = k, opts = list(), threads=0L, ...)
#' @examples
#' mat <- matrix(rnorm(500), nrow = 50, ncol = 10)
#' rownames(mat) <- paste0("gene", seq_len(50))
#' colnames(mat) <- paste0("cell", seq_len(10))
#' mat <- mat %>% as("dgCMatrix") %>% as("IterableMatrix")
#' 
#' svd_res <- svds(mat, k = 5)
#' 
#' names(svd_res)
#' 
#' svd_res$d
#' 
#' dim(svd_res$u)
#' 
#' dim(svd_res$v)
#' # Can also pass in values directly into RSpectra::svds
#' svd_res <- svds(mat, k = 5, opts=c(maxitr = 500))
#' @name svds
NULL