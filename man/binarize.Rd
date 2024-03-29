% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transforms.R
\name{binarize}
\alias{binarize}
\title{Convert matrix elements to zeros and ones}
\usage{
binarize(mat, threshold = 0, strict_inequality = TRUE)
}
\arguments{
\item{mat}{IterableMatrix}

\item{threshold}{A numeric value that determines whether the elements
of x are set to zero or one.}

\item{strict_inequality}{A logical value determining whether the
comparison to the threshold is >= (strict_inequality=FALSE)
or > (strict_inequality=TRUE).}
}
\value{
binarized IterableMatrix object
}
\description{
Binarize compares the matrix element values to the
threshold value and sets the output elements to either zero
or one. By default, element values greater than the threshold
are set to one; otherwise, set to zero. When strict_inequality
is set to FALSE, element values greater than or equal to the
threshold are set to one. As an alternative, the \code{<}, \code{<=}, \code{>},
and \code{>=} operators are also supported.
}
