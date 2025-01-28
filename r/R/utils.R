# Copyright 2021 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

# Helper-function to get random seed. Adapted from withr:::get_seed
get_seed <- function() {
  if (exists(".Random.seed", globalenv(), mode = "integer", inherits = FALSE)) {
    return(get(".Random.seed", globalenv(), mode = "integer", inherits = FALSE))
  } else {
    return(NULL)
  }
}
# Helper-function to set random seed. Adapted from withr:::restore_seed / withr:::rm_seed
restore_seed <- function(seed) {
  if (is.null(seed)) {
    rm(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    assign(".Random.seed", seed, envir = globalenv(), inherits = FALSE)
  }
}

# Helper-function for pkgdown documentation about genomic ranges inputs
document_granges <- function(
  intro_noun="Genomic regions", 
  position="`chr`, `start`, `end`: genomic position",
  strand=NULL,
  extras=NULL
) {
  if (!is.null(strand) && strand == "default") {
    strand <- "`strand`: +/- or TRUE/FALSE for positive or negative strand"
  }
  if (!is.null(extras)) {
    extras <- sprintf("`%s`: %s", names(extras), as.character(extras))
  }
  bullets <- paste0(sprintf("  - %s", c(position, strand, extras)), collapse="\n")
  sprintf(paste0(
    "%s given as GRanges, data.frame, or list. ",
    "See `help(\"genomic-ranges-like\")` for details on format and coordinate systems. ",
    "Required attributes:\n\n",
    "%s"
  ), intro_noun, bullets)
}

# Function which prints a message using shell echo.
# Useful for printing messages from inside mclapply when running in Rstudio.
log_progress <- function(msg, add_timestamp = TRUE){
  if (add_timestamp) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S "), msg)
  }
  if (.Platform$GUI == "RStudio") {
    system(sprintf('echo "%s"', paste0(msg, collapse="")))
  } else {
    message(msg)
  }
}

#' Prepare a test matrix and test fragments for BPCells.
#'
#' @param directory (character) Where the input/output data should be stored.  If NULL, a temporary directory is created.
#' @param mat_name (character) Name of the RNA matrix file. If NULL, the matrix is named "test_mat."
#' @param frags_name (character) Name of the ATAC fragments file. If NULL, the fragments are named "test_frags".
#' @return (list) A list with the RNA matrix under the name "mat", and the ATAC fragments under the name "frags".
#' @details
#' This function downloads the 10x Genomics PBMC 3k dataset, and filters the fragments and matrix to cells with at least 1000 reads.
#' Following, both fragments and the matrix is subset to only genes and insertions on chromosomes 4 and 11.
#' The RNA matrix is 1 MB and the fragments are 12.5 MB, after BPCells compression.
#' @keywords internal
prepare_test_data <- function(directory = NULL, mat_name = NULL, frags_name = NULL) {
    if (is.null(directory)) {
        directory <- file.path(tempdir())
        dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    }
    if (is.null(mat_name)) {
        mat_name <- "test_rna"
    }
    if (is.null(frags_name)) {
        frags_name <- "test_frags"
    }
    url_base <- "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/"
    rna_raw_url <- paste0(url_base, "pbmc_granulocyte_sorted_3k_raw_feature_bc_matrix.h5")
    atac_raw_url <- paste0(url_base, "pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz")
    options(timeout=300)
    if (!file.exists(file.path(directory,"pbmc_3k_10x.h5"))) {
      download.file(rna_raw_url, file.path(directory, "pbmc_3k_10x.h5"), mode="wb")
    }
    if (!file.exists(file.path(directory,"pbmc_3k_10x.fragments.tsv.gz"))) {
      download.file(atac_raw_url, file.path(directory, "pbmc_3k_10x.fragments.tsv.gz"), mode="wb")
    }
    if (!file.exists(file.path(directory,"pbmc_3k_rna_raw"))) {
      mat_raw <- open_matrix_10x_hdf5(file.path(directory, "pbmc_3k_10x.h5"), feature_type="Gene Expression") %>% 
        write_matrix_dir(file.path(directory,"pbmc_3k_rna_raw"))
    } else {
      mat_raw <- open_matrix_dir(file.path(directory,"pbmc_3k_rna_raw"))
    }
    # Check if we already ran import
    if (!file.exists(file.path(directory,"pbmc_3k_frags"))) {
      frags_raw <- open_fragments_10x(file.path(directory,"pbmc_3k_10x.fragments.tsv.gz")) %>%
          write_fragments_dir(file.path(directory,"pbmc_3k_frags"))
    } else {
      frags_raw <- open_fragments_dir(file.path(directory,"pbmc_3k_frags"))
    }
    # for atac transcripts
    transcripts <- read_gencode_transcripts(
      file.path(directory,"references_transcripts"), 
      release="42", 
      transcript_choice="MANE_Select",
      annotation_set = "basic", 
      features="transcript"
    )
    # for RNA genes
    genes_test <- read_gencode_genes(
        file.path(directory,"./reference_genes"),
        release = "42",
        annotation_set = "basic",
    )
    # Filter to only cells that have at least 1000 reads on the RNA side
    # and only genes/fragments that exist on chr 4 and 11
    filtered_cells <- colnames(mat_raw)[reads_per_cell >= 1e3]
    filtered_genes  <- genes[genes$chr %in% c("chr4", "chr11"),]$gene_id
    # remove version numbers
    filtered_genes <- gsub("\\..*", "", filtered_genes)
    mat <- mat_raw[which(rownames(mat_raw) %in% filtered_genes), pass_rna]
    frags <- select_cells(frags_raw, pass_rna) %>% select_chromosomes(c("chr4", "chr11"))
    mat <- write_matrix_dir(mat, file.path(directory, mat_name))
    frags <- write_fragments_dir(frags, file.path(directory, frags_name))
    return(list(mat = mat, frags = frags))
}