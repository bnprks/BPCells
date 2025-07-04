% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneScores.R
\name{gene_score_tiles_archr}
\alias{gene_score_tiles_archr}
\title{Calculate gene-tile distances for ArchR gene activities}
\usage{
gene_score_tiles_archr(
  genes,
  chromosome_sizes = NULL,
  tile_width = 500,
  addArchRBug = FALSE
)
}
\arguments{
\item{genes}{Gene coordinates given as GRanges, data.frame, or list. See \code{help("genomic-ranges-like")} for details on format and coordinate systems. Required attributes:
\itemize{
\item \code{chr}, \code{start}, \code{end}: genomic position
\item \code{strand}: +/- or TRUE/FALSE for positive or negative strand
}}

\item{chromosome_sizes}{(optional) Size of chromosomes as a \link{genomic-ranges} object}

\item{tile_width}{Size of tiles to consider}

\item{addArchRBug}{Replicate ArchR bug in handling nested genes}
}
\value{
Tibble with one range per tile, with additional metadata
columns gene_idx (row index of the gene this tile corresponds to) and
distance.

Distance is a signed distance calculated such that if the tile has a smaller
start coordinate than the gene and the gene is on the + strand, distance will
be negative. The distance of adjacent but non-overlapping regions is 1bp, counting
up from there.
}
\description{
ArchR-style gene activity scores are based on a weighted sum of each tile
according to the signed distance from the tile to a gene body. This function
calculates the signed distances according to ArchR's default parameters.
}
\details{
ArchR's tile distance algorithm works as follows
\enumerate{
\item Genes are extended 5kb upstream
\item Genes are linked to any tiles 1kb-100kb upstream + downstream, but tiles
beyond a neighboring gene are not considered
}
}
\examples{
## Prep data
directory  <- file.path(tempdir(), "references")
genes <- read_gencode_genes(
    directory,
    release = "42",
    annotation_set = "basic",    
)


## Get gene scores by tile
gene_score_tiles_archr(
    genes
)
}
