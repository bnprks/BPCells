#' @useDynLib BPCells, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods is new as
#' @importFrom magrittr %>%
NULL

#' @importMethodsFrom Matrix t 
#' @export 
Matrix::t
#' @importMethodsFrom Matrix rowSums
#' @export
Matrix::rowSums
#' @importMethodsFrom Matrix colSums
#' @export
Matrix::colSums
#' @importMethodsFrom Matrix rowMeans
#' @export
Matrix::rowMeans
#' @importMethodsFrom Matrix colMeans
#' @export
Matrix::colMeans



#' Genomic range formats
#'
#' @name genomic-ranges-like
#' @description
#' BPCells accepts a flexible set of genomic ranges-like objects as input, either
#' GRanges, data.frame, lists, or character vectors. These objects must specify chromosome, start,
#' and end coordinates along with optional metadata about each range.

#'  With the exception of `GenomicRanges::GRanges` objects, BPCells assumes all
#'  objects use a zero-based, end-exclusive coordinate system (see below for details).
#'
#' ## Valid Range-like objects
#' BPCells can interpret the following types as ranges:
#' - `list()`, `data.frame()`, with columns:
#'      - `chr`: Character or factor of chromosome names
#'      - `start`: Start coordinates (0-based)
#'      - `end`: End coordinates (exclusive)
#'      - (optional) `strand`: `"+"`/`"-"` or `TRUE`/`FALSE` for pos/neg strand
#'      - (optional) Additional metadata as named `list` entries or `data.frame` columns
#' - `GenomicRanges::GRanges`
#'      - `start(x)` is interpreted as a 1-based start coordinate
#'      - `end(x)` is interpreted as an inclusive end coordinate
#'      - `strand(x)`: `"*"` entries are interpeted as postive strand
#'      - (optional) `mcols(x)` holds additional metadata
#' - `character`
#'      - Given in format `"chr1:1000-2000"` or `"chr1:1,000-2,000"`
#'      - Uses 0-based, end-exclusive coordinate system
#'      - Cannot be used for ranges where additional metadata is required
#'
#' ## Range coordinate systems
#' There are two main conventions for the coordinate systems:
#'
#' ***One-based, end-inclusive ranges***
#' - The first base of a chromosome is numbered 1
#' - The last base in a range is equal to the end coordinate
#' - e.g. 1-5 describes the first 5 bases of the chromosome
#' - Used in formats such as SAM, GTF
#' - In BPCells, used when reading or writing `GenomicRanges::GRanges` objects
#'
#' ***Zero-based, end-exclusive ranges***
#' - The first base of a chromosome is numbered 0
#' - The last base in a range is one less than the end coordinate
#' - e.g. 0-5 describes the first 5 bases of the chromosome
#' - Used in formats such as BAM, BED
#' - In BPCells, used for all other range objects
#'
NULL
