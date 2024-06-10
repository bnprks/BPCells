# Copyright 2023 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

# TODO in progress
suppressPackageStartupMessages({
  library(ArchR)
  library(BPCells)
  library(tidyverse)
})

archr_dir <- "/oak/stanford/groups/wjg/bparks/BPCells/04_data/test_real_data/archr_snyder_panc"
bpcells_dir <- "/oak/stanford/groups/wjg/bparks/BPCells/04_data/test_real_data/bpcells_snyder_panc"
min_tss <- 5 # Minimum TSS cutoff

proj <- loadArchRProject(file.path(archr_dir, "project"))
ArrowFiles <- getArrowFiles(proj)

#' Test matrix equality across arrow files
#' @param ArrowFiles Vector of arrow file path names
#' @param ArchRMatrix Name of ArchR Matrix to load from arrows
#' @param bpcells_mat BPCells matrix of orientation feature x cells
#' @param bpcells_transform If not NULL, a function of 1 parameter that transforms a 
#'   dgCMatrix from BPCells prior to comparison with ArchR
test_matrix_equality <- function(ArrowFiles, ArchRMatrix, bpcells_mat, bpcells_transform=NULL) {
  cell_names <- NULL
  for (i in seq_along(ArrowFiles)) {
    cat(sprintf("Testing matrix %s on file %s\n", ArchRMatrix, ArrowFiles[i]))
    mat_archr <- ArchR::getMatrixFromArrow(ArrowFiles[i], useMatrix=ArchRMatrix, logFile = NULL)
    cell_names <- c(cell_names, colnames(mat_archr))
    mat_bpcells <- bpcells_mat[, colnames(mat_archr)] %>% as("dgCMatrix")
    if (!is.null(bpcells_transform)) {
      mat_bpcells <- bpcells_transform(mat_bpcells)
    }
    mat_archr <- assay(mat_archr)
    if (is.null(rownames(mat_archr))) {
      rownames(mat_bpcells) <- NULL
    }
    stopifnot(all.equal(mat_bpcells, mat_archr))
  }
  stopifnot(length(cell_names) == ncol(bpcells_mat))
  stopifnot(length(setdiff(cell_names, colnames(bpcells_mat))) == 0)
}

##########################
# QC Stats Equality
##########################
cat("Checking QC stats\n")
archr_qc_path <- Sys.glob(file.path(archr_dir, "/QualityControl/*/*-Pre-Filter-Metadata.rds"))
qc_archr <- lapply(archr_qc_path, readRDS) %>%
  do.call(rbind, .) %>%
  tibble::as_tibble()

qc_bpcells <- readr::read_tsv(file.path(bpcells_dir, "cell_qc.tsv.gz"))
stopifnot(length(setdiff(qc_archr$cellNames, qc_bpcells$cellName)) == 0)
rownames(qc_bpcells) <- qc_bpcells$cellName
qc_bpcells_compat <- qc_bpcells[qc_archr$cellNames, ] %>%
  dplyr::mutate(TSSEnrichment = round(TSSEnrichment, 3)) %>%
  dplyr::rename(
    cellNames = cellName,
    nMonoFrags = subNucleosomal,
    nDiFrags = monoNucleosomal,
    nMultiFrags = multiNucleosomal
  ) %>%
  dplyr::mutate(
    BlacklistRatio = ReadsInBlacklist / (2 * nFrags),
    PromoterRatio = ReadsInPromoter / (2 * nFrags),
    Keep = 0 + (TSSEnrichment >= min_tss)
  ) %>%
  dplyr::select(colnames(qc_archr))

stopifnot(all.equal(qc_bpcells_compat, qc_archr))

#################################
# Tile matrix equality
#################################
cat("Checking tile matrix\n")
# Get the tile labels for ArchR's matrix
chrom_sizes <- getGenomeAnnotation(proj)$chromSizes %>%
    {set_names(end(.), seqnames(.))}
tile_metadata <- ArchR:::.getFeatureDF(ArrowFiles, subGroup="TileMatrix")
tile_names_archr <- stringr::str_c(
    as.character(tile_metadata$seqnames), ":", 
    as.integer(tile_metadata$start), "-", 
    as.integer(pmin(tile_metadata$start + 500L, chrom_sizes[as.character(tile_metadata$seqnames)]))
)

# Get which tiles are blacklisted. Note that ArchR has an incorrect blacklist calculation in its TileMatrix code
buggy_blacklist <- getGenomeAnnotation(proj)$blacklist
width(buggy_blacklist) <- ((width(buggy_blacklist)-1)%/% 250) * 250 + 1
start(buggy_blacklist) <- start(buggy_blacklist) + 1
end(buggy_blacklist) <- end(buggy_blacklist) + 1
blacklist_tiles <- unique(BPCells:::range_overlaps(tile_names_archr, buggy_blacklist)$from)
blacklist_mask <- ifelse(seq_along(tile_names_archr) %in% blacklist_tiles, 0, 1)

tile_mat_bpcells <- open_matrix_dir(file.path(bpcells_dir, "tile_mat")) %>%
  .[tile_names_archr,] %>%
  multiply_rows(blacklist_mask)

test_matrix_equality(ArrowFiles, "TileMatrix", tile_mat_bpcells, bpcells_transform=Matrix::drop0)

#################################
# Peak matrix equality
#################################
cat("Checking peak matrix\n")

peak_metadata <- ArchR:::.getFeatureDF(ArrowFiles, subGroup="PeakMatrix")
peak_names_archr <- sprintf(
  "%s:%d-%d", 
  peak_metadata$seqnames, 
  peak_metadata$start-1, 
  peak_metadata$end
)
#TODO: Test this actually works after deduplicating peak set
peak_mat_bpcells <- open_matrix_dir(file.path(bpcells_dir, "peak_mat")) %>%
  .[peak_names_archr,]
test_matrix_equality(ArrowFiles, "PeakMatrix", peak_mat_bpcells)

##########################
# Footprint Counts Equality
##########################
cat("Checking footprints\n")
foot_archr_raw <- readRDS(file.path(archr_dir, "project/footprints.Rds"))
foot_archr_mat <- SummarizedExperiment::assay(foot_archr_raw)[SummarizedExperiment::rowData(foot_archr_raw)$type == "footprint", ]

foot_archr <- tibble::as_tibble(foot_archr_mat) %>%
  dplyr::mutate(position = -2000:2000) %>%
  tidyr::pivot_longer(!position, names_to = "group", values_to = "count")

foot_bpcells <- readr::read_tsv(file.path(bpcells_dir, "footprint.tsv.gz"))

foot <- dplyr::inner_join(foot_archr, foot_bpcells, by = c("position", "group"), suffix = c(".archr", ".bpcells"))

stopifnot(all.equal(foot$count.archr, foot$count.bpcells))
stopifnot(nrow(foot) == nrow(foot_bpcells))
stopifnot(nrow(foot_bpcells) == nrow(foot_archr))


#################################
# Gene activity score equality
#################################
cat("Checking gene activity scores\n")
gene_metadata <- ArchR:::.getFeatureDF(ArrowFiles, subGroup="GeneScoreMatrix")

gene_mat_bpcells <- open_matrix_dir(file.path(bpcells_dir, "gene_score_mat")) %>%
  .[gene_metadata$name,]

transform_gene_mat <- function(mat) {
  mat@x <- round(mat@x, 3)
  Matrix::drop0(mat)
}


test_matrix_equality(ArrowFiles, "GeneScoreMatrix", gene_mat_bpcells, bpcells_transform=transform_gene_mat)