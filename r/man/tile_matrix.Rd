% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix.R
\name{tile_matrix}
\alias{tile_matrix}
\title{Calculate ranges x cells tile overlap matrix}
\usage{
tile_matrix(
  fragments,
  ranges,
  mode = c("insertions", "fragments"),
  zero_based_coords = !is(ranges, "GRanges"),
  explicit_tile_names = FALSE
)
}
\arguments{
\item{fragments}{Input fragments object}

\item{ranges}{Tiled regions given as GRanges, data.frame, or list. See \code{help("genomic-ranges-like")} for details on format and coordinate systems. Required attributes:
\itemize{
\item \code{chr}, \code{start}, \code{end}: genomic position
\item \code{tile_width}: Size of each tile in this region in basepairs
}

Must be non-overlapping and sorted by
(chr, start), with chromosomes ordered according to the chromosome names of \code{fragments}}

\item{mode}{Mode for counting tile overlaps. (See "value" section for more detail)}

\item{zero_based_coords}{Whether to convert the ranges from a 1-based end-inclusive
coordinate system to a 0-based end-exclusive coordinate system. Defaults to true
for GRanges and false for other formats
(see this \href{https://web.archive.org/web/20210920203703/http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/}{archived UCSC blogpost})}

\item{explicit_tile_names}{Boolean for whether to add rownames to the output matrix in format e.g
chr1:500-1000, where start and end coords are given in a 0-based coordinate system. For
whole-genome Tile matrices the names will take ~5 seconds to generate and take up 400MB of memory.
Note that either way, tile names will be written when the matrix is saved.}
}
\value{
Iterable matrix object with dimension ranges x cells. When saved,
the column names will be in the format chr1:500-1000,
where start and end coords are given in a 0-based coordinate system.

\strong{\code{mode} options}
\itemize{
\item \code{"insertions"}: Start and end coordinates are separately overlapped with each tile
\item \code{"fragments"}: Like \code{"insertions"}, but each fragment can contribute at most 1 count
to each tile, even if both the start and end coordinates overlap
}
}
\description{
Calculate ranges x cells tile overlap matrix
}
\note{
When calculating the matrix directly from a fragments tsv, it's necessary to first call \code{select_chromosomes()} in order to
provide the ordering of chromosomes to expect while reading the tsv.
}
\examples{
## Prep demo data
frags <- get_demo_frags(subset = FALSE)
chrom_sizes <- read_ucsc_chrom_sizes(file.path(tempdir(), "references"), genome="hg38")
blacklist <- read_encode_blacklist(file.path(tempdir(), "references"), genome="hg38")
frags_filter_blacklist <- frags \%>\% select_regions(blacklist, invert_selection = TRUE)
ranges <- tibble::tibble(
  chr = "chr4",
  start = 0,
  end = "190214555", 
  tile_width = 200
)


## Get tile matrix
tile_matrix(frags_filter_blacklist, ranges)
}
