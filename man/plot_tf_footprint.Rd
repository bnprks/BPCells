% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_tf_footprint}
\alias{plot_tf_footprint}
\title{Plot TF footprint}
\usage{
plot_tf_footprint(
  fragments,
  motif_positions,
  cell_groups = rlang::rep_along(cellNames(fragments), "all"),
  flank = 250L,
  smooth = 0L,
  zero_based_coords = !is(genes, "GRanges"),
  colors = discrete_palette("stallion"),
  return_data = FALSE,
  apply_styling = TRUE
)
}
\arguments{
\item{fragments}{IterableFragments object}

\item{motif_positions}{Coordinate ranges for motifs (must include strand) and
have constant width}

\item{cell_groups}{Character or factor assigning a group to each cell, in order of
\code{cellNames(fragments)}}

\item{flank}{Number of flanking basepairs to include on either side of the motif}

\item{smooth}{(optional) Sparse matrix of dimensions cells x cells with cell-cell distance
weights for smoothing.}

\item{zero_based_coords}{If true, coordinates start and 0 and the end coordinate is not included in the range.
If false, coordinates start at 1 and the end coordinate is included in the range}

\item{return_data}{If true, return data from just before plotting rather than a plot.}

\item{apply_styling}{If false, return a plot without pretty styling applied}
}
\description{
Plot the footprinting around TF motif sites
}
\seealso{
\code{footprint()}, \code{plot_tss_profile()}
}
