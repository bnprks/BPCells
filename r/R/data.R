# Copyright 2022 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

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
