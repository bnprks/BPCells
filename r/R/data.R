# Copyright 2022 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

#' Prepare a demo matrix and demo fragments for BPCells.
#'
#' @param directory (character) Where the input/output data should be stored.  If NULL, a temporary directory is created.
#' @param mat_name (character) Name of the RNA matrix file. If NULL, the matrix is named "demo_mat."
#' @param frags_name (character) Name of the ATAC fragments file. If NULL, the fragments are named "demo_frags".
#' @param timeout (numeric) Timeout for downloading files in seconds.
#' @param remove_input_data (logical) Whether to remove the downloaded non-procesed matrix, frags, gencode transcripts, and gencode genes
#' after processing.
#' @return (list) A list with the RNA matrix under the name "mat", and the ATAC fragments under the name "frags".
#' @details
#' This function downloads the 10x Genomics PBMC 3k dataset, and filters the fragments and matrix to cells with at least 1000 reads.
#' Following, both fragments and the matrix is subset to only genes and insertions on chromosomes 4 and 11.
#' The RNA matrix is 1 MB and the fragments are 12.5 MB, after BPCells compression.
#' @keywords internal
prepare_demo_data <- function(directory = NULL, mat_name = NULL, frags_name = NULL, timeout = 300, remove_input_data = TRUE) {
    if (is.null(directory)) {
        directory <- file.path(tempdir())
        dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    }
    if (is.null(mat_name)) {
        mat_name <- "demo_mat"
    }
    if (is.null(frags_name)) {
        frags_name <- "demo_frags"
    }
    url_base <- "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/"
    rna_raw_url <- paste0(url_base, "pbmc_granulocyte_sorted_3k_raw_feature_bc_matrix.h5")
    atac_raw_url <- paste0(url_base, "pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz")
    ensure_downloaded(file.path(directory, "pbmc_3k_10x.h5"), rna_raw_url, timeout = timeout)
    ensure_downloaded(file.path(directory, "pbmc_3k_10x.fragments.tsv.gz"), atac_raw_url, timeout = timeout)
    if (!file.exists(file.path(directory,"pbmc_3k_rna_raw"))) {
      mat_raw <- open_matrix_10x_hdf5(file.path(directory, "pbmc_3k_10x.h5"), feature_type="Gene Expression") %>% 
        write_matrix_dir(file.path(directory, "pbmc_3k_rna_raw"))
    } else {
      mat_raw <- open_matrix_dir(file.path(directory, "pbmc_3k_rna_raw"))
    }
    # Check if we already ran import
    if (!file.exists(file.path(directory, "pbmc_3k_frags"))) {
      frags_raw <- open_fragments_10x(file.path(directory, "pbmc_3k_10x.fragments.tsv.gz")) %>%
          write_fragments_dir(file.path(directory, "pbmc_3k_frags"))
    } else {
      frags_raw <- open_fragments_dir(file.path(directory, "pbmc_3k_frags"))
    }
    # for atac transcripts
    transcripts <- read_gencode_transcripts(
      file.path(directory, "references_transcripts"), 
      release = "42", 
      transcript_choice = "MANE_Select",
      annotation_set = "basic", 
      features = "transcript"
    )
    # for RNA genes
    genes_demo <- read_gencode_genes(
        file.path(directory, "./reference_genes"),
        release = "42",
        annotation_set = "basic",
    )
    # Filter to only cells that have at least 1000 reads on the RNA side
    # and only genes/fragments that exist on chr 4 and 11
    reads_per_cell <- colSums(mat_raw)
    filtered_cells <- colnames(mat_raw)[reads_per_cell >= 1e3]
    filtered_genes  <- genes_demo[genes_demo$chr %in% c("chr4", "chr11"),]$gene_id
    # remove version numbers
    filtered_genes <- gsub("\\..*", "", filtered_genes)
    mat <- mat_raw[which(rownames(mat_raw) %in% filtered_genes), filtered_cells]
    frags <- select_cells(frags_raw, filtered_cells) %>% select_chromosomes(c("chr4", "chr11"))
    mat <- write_matrix_dir(mat, file.path(directory, mat_name), overwrite = TRUE)
    frags <- write_fragments_dir(frags, file.path(directory, frags_name), overwrite = TRUE)
    if (remove_input_data) {
        unlink(file.path(directory, "pbmc_3k_10x.h5"))
        unlink(file.path(directory, "pbmc_3k_10x.fragments.tsv.gz"))
        unlink(file.path(directory, "pbmc_3k_rna_raw"))
        unlink(file.path(directory, "pbmc_3k_frags"))
    }
    return(list(mat = mat, frags = frags))
}

#' Retrieve BPCells demo data
#' 
#' The demo dataset is a subset of the  10x Genomics PBMC 3k dataset, and filters both the matrix and the fragments
#' to cells with at least 1000 reads.  Both the matrix and the fragments are subset to only genes on chromosomes 4 and 11.  
#' @rdname demo_data
#' @return 
#' - `get_demo_mat()`: (IterableMatrix) A `(features x cells)` matrix of shape `(1984 x 2724)`.
#' @details 
#' The first time either `get_demo_mat()` are ran `get_demo_frags()`, the demo data is downloaded and stored in the BPCells data directory 
#' (under `file.path(tools::R_user_dir("BPcells", which="data"), "demo_data")`).  Subsequent calls to this function will use the previously downloaded matrix.
#' The preperation of this matrix can be reproduced by running the internal function `prepare_demo_data()`.  
#' 
#' In the case that demo data is not pre-downloaded and demo data download fails, `prepare_demo_data()` will be run,
#'  which manually builds the demo dataset from the 10x Genomics PBMC 3k dataset.
#' 
#' Both the matrix from `get_demo_mat()` and the fragments from `get_demo_frags()` 
#' may be removed by running `remove_demo_data()`.
#' 
#' - `get_demo_mat()`: Retrieve a 1 MB demo matrix, representing a subset of the 10X Genomics PBMC 3k dataset.  
#' @export
get_demo_mat <- function() {
    # Use the data directory for BPCells
    data_dir <- file.path(tools::R_user_dir("BPCells", which = "data"), "demo_data")
    if (!dir.exists(data_dir)) {
        dir.create(data_dir, recursive = TRUE)
    }
    if (!dir.exists(file.path(data_dir, "demo_mat"))) {
        url <- "https://pub-c4e56988ff67429e9856ffa33aecb0c1.r2.dev/demo_mat.tar.gz"
        download.file(url, file.path(data_dir, "demo_mat.tar.gz"))
        # Check if file download failed
        if (!file.exists(file.path(data_dir, "demo_mat.tar.gz"))) {
            prepare_demo_data(data_dir)
        } else {
            untar(file.path(data_dir, "demo_mat.tar.gz"), exdir=data_dir)
            file.remove(file.path(data_dir, "demo_mat.tar.gz"))
        }
    }
    return(open_matrix_dir(file.path(data_dir, "demo_mat")))
}

#' @rdname demo_data
#' @return
#' - `get_demo_frags()`: (IterableFragments) A Fragments object with 2724 cells and fragments on chromosomes 4 and 11.
#' @details
#' - `get_demo_frags()`: Retrieve a 12.5 MB demo fragments object, representing a subset of the 10X Genomics PBMC 3k dataset.
#' @export 
get_demo_frags <- function() {
    data_dir <- file.path(tools::R_user_dir("BPCells", which = "data"), "demo_data")
    if (!dir.exists(data_dir)) {
        dir.create(data_dir, recursive = TRUE)
    }
    if (!dir.exists(file.path(data_dir, "demo_frags"))) {
        url <- "https://pub-c4e56988ff67429e9856ffa33aecb0c1.r2.dev/demo_frags.tar.gz"
        download.file(url, file.path(data_dir, "demo_frags.tar.gz"))
        if (!file.exists(file.path(data_dir, "demo_frags.tar.gz"))) {
            prepare_demo_data(data_dir)
        } else {
            untar(file.path(data_dir, "demo_frags.tar.gz"), exdir = data_dir)
            file.remove(file.path(data_dir, "demo_frags.tar.gz"))
        }
    }
    return(open_fragments_dir(file.path(data_dir, "demo_frags")))
}

#' @rdname demo_data
#' @return
#' - `remove_demo_data()`: NULL
#' @details 
#'  - `remove_demo_data()`: Remove the demo data from the BPCells data directory.
#' @export
remove_demo_data <- function() {
    data_dir <- file.path(tools::R_user_dir("BPCells", which = "data"), "demo_data")
    if (dir.exists(data_dir)) {
        unlink(data_dir, recursive = TRUE)
    }
}

#' Gene Symbol Mapping data
#'
#' Mapping of the canonical gene symbols corresponding to each
#' unambiguous alias, previous symbol, ensembl ID, or entrez ID. 
#'
#' @format **human_gene_mapping**
#'
#' A named character vector. Names are aliases or IDs and values
#' are the corresponding canonical gene symbol
#'
#' @details See the source code in `data-raw/human_gene_mapping.R` and
#' `data-raw/mouse_gene_mapping.R` for exactly how these mappings were made.
#' @source **human_gene_mapping**
#'
#' <http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt>
#'
#' @rdname gene_mapping
"human_gene_mapping"

#' @rdname gene_mapping
#' @format **mouse_gene_mapping**
#'
#' A named character vector. Names are aliases or IDs and values
#' are the corresponding canonical gene symbol
#'
#' @source **mouse_gene_mapping**
#'
#' <http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt>
#' <http://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt>
"mouse_gene_mapping"


