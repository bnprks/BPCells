% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_tss_scatter}
\alias{plot_tss_scatter}
\title{TSS Enrichment vs. Fragment Counts plot}
\usage{
plot_tss_scatter(
  atac_qc,
  min_frags = NULL,
  min_tss = NULL,
  bins = 100,
  apply_styling = TRUE
)
}
\arguments{
\item{atac_qc}{Tibble as returned by \code{qc_scATAC()}. Must have columns \code{nFrags} and \code{TSSEnrichment}}

\item{min_frags}{Minimum fragment count cutoff}

\item{min_tss}{Minimum TSS Enrichment cutoff}

\item{bins}{Number of bins for density calculation}

\item{apply_styling}{If false, return a plot without pretty styling applied}
}
\description{
Density scatter plot with log10(fragment_count) on the x-axis and TSS
enrichment on the y-axis. This plot is most useful to select which cell barcodes
in an experiment correspond to high-quality cells
}
