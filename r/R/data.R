# Copyright 2022 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

#' Create a small demo matrix and fragment object.
#' 
#' Downloads a 
#' [10x Genomics dataset](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0), 
#'  consisting of 3k cells then performs optional QC and subsetting.  Holds subsetted objects in disk,
#' and returns a list with both the matrix and fragments. 
#' @param directory (character) The directory where all the input/output data will be stored.
#' Downloaded intermediates will be stored in subdir `intermediates`.
#' If `NULL`, a temporary directory is created.
#' @param filter_qc (bool) Whether to filter both the RNA and ATAC data using QC information.
#' @param subset (bool) Whether to subset to only genes/insertions on chromosome 4 and 11.
#' @param timeout (numeric) Timeout for downloading files in seconds.
#' @return (list) A list with the RNA matrix under the name `mat`, and the ATAC fragments under the name `frags`.
#' @details
#' This function downloads the 10x Genomics PBMC 3k dataset. 
#' Filtering using QC information on the fragments and matrix provides cells with at least 1000 reads, 1000 frags, and a minimum tss enrichment of 10.
#' Subsetting provides only genes and insertions on chromosomes 4 and 11.
#' The name of the matrix and fragments folders are `demo_mat` and `demo_frags` respectively.  
#' Additionally, choosing to qc filter appends a `_filtered`, and choosing to subset data appends a `_subsetted` to the name.
#' @keywords internal
prepare_demo_data <- function(directory = NULL, filter_qc = TRUE, subset = TRUE, timeout = 300) {
    if (is.null(directory)) {
      directory <- file.path(tempdir())
    }
    intermediate_dir <- file.path(directory, "intermediates")
    dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(intermediate_dir, recursive = TRUE))

    mat_name <- "demo_mat"
    frags_name <- "demo_frags"

    # Download matrix/frags if not done previously, and open
    url_base <- "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/"
    # Recreate mat if mat is malformed
    mat <- NULL
    frags <- NULL
    tryCatch({
        mat <- open_matrix_dir(file.path(directory, mat_name))
    }, error = function(e) {
        rna_raw_url <- paste0(url_base, "pbmc_granulocyte_sorted_3k_raw_feature_bc_matrix.h5")
        ensure_downloaded(file.path(intermediate_dir, "pbmc_3k_10x.h5"), rna_raw_url, timeout = timeout)
        mat <<- open_matrix_10x_hdf5(file.path(intermediate_dir, "pbmc_3k_10x.h5"), feature_type="Gene Expression") %>% 
            write_matrix_dir(file.path(directory, mat_name), overwrite = TRUE)
    })
    # Recreate frags if frags are malformed
    tryCatch({
        frags <- open_fragments_dir(file.path(directory, frags_name))
    }, error = function(e) {
        atac_raw_url <- paste0(url_base, "pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz")
        ensure_downloaded(file.path(intermediate_dir, "pbmc_3k_10x.fragments.tsv.gz"), atac_raw_url, timeout = timeout)
        frags <<- open_fragments_10x(file.path(intermediate_dir, "pbmc_3k_10x.fragments.tsv.gz")) %>%
            write_fragments_dir(file.path(directory, frags_name), overwrite = TRUE)
    })
    if (filter_qc) {
        # Download annotations for transcripts
        transcripts <- read_gencode_transcripts(
            intermediate_dir, 
            release = "42", 
            transcript_choice = "MANE_Select",
            annotation_set = "basic", 
            features = "transcript"
        )
        blacklist <- read_encode_blacklist(intermediate_dir, genome="hg38")
        atac_qc <- qc_scATAC(frags, transcripts, blacklist)
        # Filter to only cells that have at least 1000 reads on the RNA side
        # a minimum of 1000 frag reads, and greater than 10 tss enrichment
        pass_atac <- atac_qc %>%
            dplyr::filter(nFrags > 1000, TSSEnrichment > 10) %>%
            dplyr::pull(cellName)
        pass_rna <- colnames(mat)[colSums(mat) > 1000]
        filtered_cells <- intersect(pass_atac, pass_rna)
        frags <- select_cells(frags, filtered_cells)
        mat <- mat[, filtered_cells]
    }
    if (subset) {
        # Subset to only genes/fragments that exist on chr4 and 11
        genes_demo <- read_gencode_genes(
            intermediate_dir,
            release = "42",
            annotation_set = "basic",
        )
        filtered_genes  <- genes_demo[genes_demo$chr %in% c("chr4", "chr11"),]$gene_id
        # remove version numbers
        filtered_genes <- gsub("\\..*", "", filtered_genes)
        mat <- mat[which(rownames(mat) %in% filtered_genes), ]
        frags <- frags %>% select_chromosomes(c("chr4", "chr11"))
    }

    # Rename mat and frags depending on state of filtering and subsetting
    if (filter_qc) {
        mat_name <- paste0(mat_name, "_filtered")
        frags_name <- paste0(frags_name, "_filtered")
    }
    if (subset) {
        mat_name <- paste0(mat_name, "_subsetted")
        frags_name <- paste0(frags_name, "_subsetted")
    }
    # Write changes to directory
    if (filter_qc || subset) {
        mat <- write_matrix_dir(mat, file.path(directory, mat_name), overwrite = TRUE)
        frags <- write_fragments_dir(frags, file.path(directory, frags_name), overwrite = TRUE)
    }
    return(list(mat = mat, frags = frags))
}

#' Retrieve BPCells demo data
#' 
#' `r lifecycle::badge("experimental")` \cr Functions to download matrices and fragments derived from a 
#' [10X Genomics PBMC 3k dataset](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0),
#' with options to filter with common qc metrics, and to subset genes and fragments to only chromosome 4 and 11.
#' @rdname demo_data
#' @param filter_qc (bool) Whether to filter both the RNA and ATAC data using qc metrics (described in `details`).
#' @param subset (bool) Whether to subset to only genes/insertions on chromosome 4 and 11.
#' @return 
#' - `get_demo_mat()`: (IterableMatrix) A `(features x cells)` matrix.
#' @details
#' These data functions are experimental. 
#' The interface, as well as the demo dataset itself will likely undergo changes in the near future.
#' 
#' **Data Processing**:
#' 
#' The first time either `get_demo_mat()`, or `get_demo_frags()`, are ran 
#' demo data is downloaded and stored in the BPCells data directory  (under `file.path(tools::R_user_dir("BPcells", which="data"), "demo_data")`).
#' 
#' Subsequent calls to this function will use the previously downloaded matrix/fragments, given that the same combination of filtering and 
#' subsetting has been performed previously.
#' 
#' The preparation of this matrix can be reproduced by running the internal function `prepare_demo_data()` with `directory` set to the BPCells data directory.
#' 
#' In the case that demo data is not pre-downloaded and demo data download fails, `prepare_demo_data()` will act
#' as a fallback.
#' 
#' Both the matrix from `get_demo_mat()` and the fragments from `get_demo_frags()`
#' may be removed by running `remove_demo_data()`.
#' 
#' Filtering using QC information on the fragments and matrix object chooses cells with at least 1000 reads, 1000 frags, and a minimum tss enrichment of 10.
#' Subsetting provides only genes and insertions on chromosomes 4 and 11.
#' 
#' **Dimensions**:
#' \tabular{lll}{
#'   \strong{Condition} \tab \strong{RNA matrix (features x cells)} \tab \strong{Fragments (chromosomes x cells)} \cr
#'   Raw                \tab (36601 x 650165)                       \tab (39 x 462264)                            \cr
#'   Filter             \tab (36601 x 2600)                         \tab (39 x 2600)                              \cr
#'   Subset             \tab (3582 x 650165)                        \tab (2 x 462264)                             \cr
#'   Filter + Subset    \tab (3582 x 2600)                          \tab (2 x 2600)                               \cr
#' }
#' **Data size**:
#' \tabular{lll}{
#'   \strong{Condition}          \tab \strong{RNA matrix (MB)} \tab \strong{Fragments (MB)} \cr
#'   Raw                         \tab 31.9                     \tab 200                     \cr
#'   Filter                      \tab 9.4                      \tab 137                     \cr
#'   Subset                      \tab 18.3                     \tab 25.6                    \cr
#'   Filter + Subset             \tab 1.2                      \tab 12.3                    \cr
#' }
#' 
#' **Function Description**:
#' 
#' - `get_demo_mat()`: Retrieve a demo `IterableMatrix` object representing the 10X Genomics PBMC 3k dataset.
#' @examples
#' #######################################################################
#' ## get_demo_mat() example
#' #######################################################################
#' get_demo_mat()
#' 
#' 
#' @export
get_demo_mat <- function(filter_qc = TRUE, subset = TRUE) {
    # Use the data directory for BPCells
    data_dir <- file.path(tools::R_user_dir("BPCells", which = "data"), "demo_data")
    if (!dir.exists(data_dir)) {
        dir.create(data_dir, recursive = TRUE)
    }
    mat_name = "demo_mat"
    if (filter_qc) mat_name <- paste0(mat_name, "_filtered")
    if (subset) mat_name <- paste0(mat_name, "_subsetted")
    if (!dir.exists(file.path(data_dir, mat_name))) {
        url <- paste0("https://pub-c4e56988ff67429e9856ffa33aecb0c1.r2.dev/", mat_name, ".tar.gz")
        suppressWarnings(try(download.file(url, file.path(data_dir, paste0(mat_name, ".tar.gz"))), silent = TRUE))
        # Check if file download failed
        if (!file.exists(file.path(data_dir, paste0(mat_name, ".tar.gz")))) {
            prepare_demo_data(data_dir, filter_qc = filter_qc, subset = subset)
        } else {
            untar(file.path(data_dir, paste0(mat_name, ".tar.gz")), exdir=data_dir)
            file.remove(file.path(data_dir, paste0(mat_name, ".tar.gz")))
        }
    }
    return(open_matrix_dir(file.path(data_dir, mat_name)))
}

#' @rdname demo_data
#' @return
#' - `get_demo_frags()`: (IterableFragments) A Fragments object.
#' @details
#' - `get_demo_frags()`: Retrieve a demo `IterableFragments` object representing the 10X Genomics PBMC 3k dataset.
#' @examples
#' #######################################################################
#' ## get_demo_frags() example
#' #######################################################################
#' get_demo_frags()
#' 
#' 
#' @export 
get_demo_frags <- function(filter_qc = TRUE, subset = TRUE) {
    data_dir <- file.path(tools::R_user_dir("BPCells", which = "data"), "demo_data")
    if (!dir.exists(data_dir)) {
        dir.create(data_dir, recursive = TRUE)
    }
    frags_name <- "demo_frags"
    if (filter_qc) frags_name <- paste0(frags_name, "_filtered")
    if (subset) frags_name <- paste0(frags_name, "_subsetted")
    if (!dir.exists(file.path(data_dir, frags_name))) {
        url <- paste0("https://pub-c4e56988ff67429e9856ffa33aecb0c1.r2.dev/", frags_name, ".tar.gz")
        suppressWarnings(try(download.file(url, file.path(data_dir, paste0(frags_name, ".tar.gz"))), silent = TRUE))
        if (!file.exists(file.path(data_dir, paste0(frags_name, ".tar.gz")))) {
            prepare_demo_data(data_dir)
        } else {
            untar(file.path(data_dir, paste0(frags_name, ".tar.gz")), exdir = data_dir)
            file.remove(file.path(data_dir, paste0(frags_name, ".tar.gz")))
        }
    }
    return(open_fragments_dir(file.path(data_dir, frags_name)))
}

#' @rdname demo_data
#' @return
#' - `remove_demo_data()`: `NULL`
#' @details 
#' - `remove_demo_data()`: Remove the demo data from the BPCells data directory.
#' @examples
#' #######################################################################
#' ## remove_demo_data() example
#' #######################################################################
#' remove_demo_data()
#' 
#' 
#' ## Demo data folder is now empty
#' data_dir <- file.path(tools::R_user_dir("BPCells", which = "data"), "demo_data")
#' list.files(data_dir)
#' 
#' 
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
#' @examples
#' #######################################################################
#' ## human_gene_mapping
#' head(human_gene_mapping)
#' #######################################################################
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
#' @examples
#' #######################################################################
#' ## mouse_gene_mapping
#' head(mouse_gene_mapping)
"mouse_gene_mapping"


