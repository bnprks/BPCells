% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/singlecell_utils.R
\name{marker_features}
\alias{marker_features}
\title{Test for marker features}
\usage{
marker_features(mat, groups, method = "wilcoxon")
}
\arguments{
\item{mat}{IterableMatrix object of dimensions features x cells}

\item{groups}{Character/factor vector of cell groups/clusters. Length #cells}

\item{method}{Test method to use. Current options are:
\itemize{
\item \code{wilcoxon}: Wilconxon rank-sum test a.k.a Mann-Whitney U test
}}
}
\value{
tibble with the following columns:
\itemize{
\item \strong{foreground}: Group ID used for the foreground
\item \strong{background}: Group ID used for the background (or NA if comparing to rest of cells)
\item \strong{feature}: ID of the feature
\item \strong{p_val_raw}: Unadjusted p-value for differential test
\item \strong{foreground_mean}: Average value in the foreground group
\item \strong{background_mean}: Average value in the background group
}
}
\description{
Given a features x cells matrix, perform one-vs-all differential
tests to find markers.
}
\details{
Tips for using the values from this function:
\itemize{
\item Use \code{dplyr::mutate()} to add columns for e.g. adjusted p-value and log fold change.
\item Use \code{dplyr::filter()} to get only differential genes above some given threshold
\item To get adjusted p-values, use R \code{p.adjust()}, recommended method is "BH"
\item To get log2 fold change: if your input matrix was already log-transformed,
calculate \code{(foreground_mean - background_mean)/log(2)}. If your input
matrix was not log-transformed, calculate \code{log2(forground_mean/background_mean)}
}
}
