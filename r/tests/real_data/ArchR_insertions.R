# Copyright 2024 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

devtools::load_all("/mnt/c/Users/Immanuel/PycharmProjects/BPCells/r")
devtools::load_all("/mnt/c/Users/Immanuel/PycharmProjects/ArchR")

# Set up temp dir in case it's not already set
create_temp_dir <- function(dir = NULL) {
  if (is.null(dir)) {
    dir <- file.path(tempdir(), "lsi_test")
    if (dir.exists(dir)) unlink(dir, recursive = TRUE)
    dir.create(dir)
  }
  return(dir)
}

fix_granges_syntax_for_archr <- function(gr) {
    mcols(gr)$RG <- gsub("PBSmall#", "", mcols(gr)$RG)
    names(mcols(gr))[names(mcols(gr)) == "RG"] <- "cell"
    mcols(gr)$barcode <- 1
    gr <- tibble::tibble(as.data.frame(gr))
    # remove width and strand columns
    gr <- gr %>% dplyr::select(!c(width, strand)) %>% dplyr::rename(chr = seqnames)
    return(gr)
}

#' Subset an archr project into two, and create a new archr project with each subset
#' as a sample.
#' @param final_dir Where to store the final archr project with two samples. Chr vector of length 1.
#' @param int_dir Where to put intermediate files for generation of final archr project. Chr vector of length 1.
create_subsetted_test_project <- function(final_dir, int_dir) {
    setwd(int_dir)
    files <- list.files(int_dir, full.names = TRUE)
    unlink(files, recursive=TRUE)
    proj <- getTestProject()
    # remove previous files from archr project
    proj <- saveArchRProject(proj, outputDirectory = file.path(int_dir, "ArchRIntermediates"))
    subset <- sample(x=c(TRUE, FALSE), size = length(proj$Sample), replace = TRUE)
    
    proj_1 <- proj[subset,]
    proj_2 <- proj[!subset,]
    proj_granges <- getFragmentsFromProject(proj)
    proj_1_granges <- getFragmentsFromProject(proj, cellNames = getCellNames(proj_1))
    proj_2_granges <- getFragmentsFromProject(proj_2, cellNames = getCellNames(proj_2))
    proj_all <- fix_granges_syntax_for_archr(proj_granges$PBSmall)
    readr::write_delim(fix_granges_syntax_for_archr(proj_1_granges$PBSmall), paste0(int_dir, "/ArchRIntermediates/frag_1.tsv"),
                       delim = "\t", col_names = FALSE)
    readr::write_delim(fix_granges_syntax_for_archr(proj_2_granges$PBSmall), paste0(int_dir, "/ArchRIntermediates/frag_2.tsv"), 
                       delim = "\t", col_names = FALSE)
    readr::write_delim(proj_all, paste0(int_dir, "/ArchRIntermediates/frag_all.tsv"), delim = "\t", col_names = FALSE)
    
    reformatFragmentFiles(
        fragmentFiles = c(paste0(int_dir, "/ArchRIntermediates/frag_1.tsv"),
                          paste0(int_dir, "/ArchRIntermediates/frag_2.tsv"),
                          paste0(int_dir, "/ArchRIntermediates/frag_all.tsv")),
        checkChrPrefix = FALSE
    )
    createArrowFiles(
        inputFiles = c(paste0(int_dir, "/ArchRIntermediates/frag_1.tsv.gz"), paste0(int_dir, "/ArchRIntermediates/frag_2.tsv.gz")),
        sampleNames = c("frag_1", "frag_2"),
        minFrags = 10,
        force=TRUE
    )
    createArrowFiles(
        inputFiles = c(paste0(int_dir, "/ArchRIntermediates/frag_all.tsv.gz")), sampleNames = c("frag_all"),
        minFrags = 10,
        force=TRUE
    )
    final_dir <- file.path(final_dir, "ArchRInsertionTestProject")
    if (dir.exists(final_dir)) unlink(final_dir, recursive = TRUE)
    dir.create(final_dir)
    proj_subsetted <- ArchRProject(
        ArrowFiles = c(paste0(int_dir, "/frag_1.arrow"), paste0(int_dir, "/frag_2.arrow")),
        outputDirectory = final_dir,
        copyArrows = TRUE,
    )
    return(proj_subsetted)
}

#' Do dim reduction with iterative LSI, add clusters, and provide insertions by group to
#' an archr project.
#' @param proj An archr project.
give_insertions_to_archr_project <- function(proj) {
    # add iterative lsi for dim reduction
    proj<- addIterativeLSI(
        proj,
        useMatrix = "TileMatrix",
        name = "IterativeLSI",
        iterations = 2,
        clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10),
        varFeatures = 4471,
        dimsToUse = 1:30
    )
    # add clusters
    proj <- addClusters(
        proj,
        reducedDims = "IterativeLSI",
        method= "Seurat",
        name = "Clusters", 
        resolution = 0.8
    )
    addGroupCoverages(
      proj,
      groupBy = "Clusters",
      useLabels = TRUE,
      minCells = 10,
      minReplicates=2,
      maxReplicates=2
    )
    return(proj)
}

#' Turn the insertions rle file created by archr into a bedfile read by macs.
#' @param proj An archr project.
create_archr_insertion_beds <- function(proj) {
  proj_path <- getOutputDirectory(proj)
  coverage_path <- file.path(proj_path, "GroupCoverages", "Clusters")
  out_path <- file.path(proj_path, "InsertionBeds")
  if (dir.exists(out_path)) unlink(out_path, recursive = TRUE)
  dir.create(out_path)
  for (file_name in list.files(coverage_path))
    .writeCoverageToBed(
      coverageFile = file.path(coverage_path, file_name),
      out = file.path(out_path, (gsub(".h5", ".bed", file_name)))
    )
}

#' Use the non-subsetted archr fragment file and use the cluster assignemnts from proj (including sample assignments).
#' Then form insertion bed files with those cluster assignments.
#' @param proj An archr project.
#' @param final_dir Where to store the final BPcells insertions, held in a sub directory called BPCellsInsertions. Chr vector of length 1.
#' @param int_dir Where the intermediate archr files were stored. Chr vector of length 1.
write_bpcells_insertion_beds <- function(proj, final_dir, int_dir) {
  frags <- readr::read_tsv(paste0(int_dir, "/ArchRIntermediates/frag_all.tsv.gz"), col_names = c("chr", "start", "end", "cell_id", "barcode"))
  frags_bpcells <- convert_to_fragments(frags)
  getCellNames(proj)
  mapping <- tibble::tibble(
    archr_names = getCellNames(proj),
    sample = gsub("#.*", "", getCellNames(proj)),
    cell = gsub(".*#", "", getCellNames(proj)),
    cluster = proj$Clusters,
    group = paste0(cluster, "._.", sample),
  ) %>% arrange(match(cell, cellNames(frags_bpcells)))
  dir.create(file.path(final_dir, "BPCellsInsertions"))
  bpcells_out_path <- paste0(final_dir, "/BPCellsInsertions/", levels(as.factor(mapping$group)), ".insertions.coverage.bed")
  names(bpcells_out_path) <- levels(as.factor(mapping$group))
  write_insertion_bed(
    fragments = frags_bpcells,
    cell_groups = mapping$group,
    path = bpcells_out_path, insertion_mode = "both"
  )
}

#' Test to ensure that the insertions grouped by pseudobulked and written by ArchR and BPCells are the same 
#' given the same fragments file. 
#' @details
#' 
#' Steps: 
#'  - Use the test archr project of 127 cells, and split into two.  This is done to prevent replicate creation during pseudobulk creation.
#'  - Run iterative LSI, cluster, and create insertion beds from archr.
#'  - Use the same non-subsetted fragment file from archr, as well as the cluster/sample information to determine pseudobulks.  
#'    Use this information for pseudobulk assignment in BPCells.
#'  - Compare insertion beds between BPCells and ArchR.
#'  
compare_bpcells_archr_insertion_writing <- function() {
  set.seed(123)
  final_dir <- "/oak/stanford/groups/wjg/iabdi/"
  final_dir <- file.path(final_dir, "insertion_test")
  if (dir.exists(final_dir)) unlink(final_dir, recursive = TRUE)
  dir.create(final_dir)
  int_dir <- file.path(tempdir(), "insertion_test")
  dir.create(int_dir)
  proj_subsetted <- create_subsetted_test_project(final_dir, int_dir)
  proj_subsetted <- give_insertions_to_archr_project(proj_subsetted)
  create_archr_insertion_beds(proj_subsetted)
  write_bpcells_insertion_beds(proj_subsetted, final_dir, int_dir)
  for (file_name in  list.files(file.path(final_dir, "BPCellsInsertions"))) {
    bp_insertion <- readr::read_tsv(file.path(final_dir, "BPCellsInsertions", file_name), col_names = FALSE)
    archr_insertion <- readr::read_tsv(file.path(final_dir, "ArchRInsertionTestProject", "InsertionBeds", file_name), col_names = FALSE) %>%
      dplyr::arrange(X1, X2, X3)
    testthat::expect_equal(bp_insertion, archr_insertion)
  }
}
compare_bpcells_archr_insertion_writing()


