# Copyright 2023 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(GenomicRanges)
  library(ArchR)
})

args <- c(
  "/oak/stanford/groups/wjg/bparks/BPCells/01_raw_data/ENCODE/samples/snyder_*_PANC/fragments.tsv.gz",
  "hg38",
  "2000",
  "5",
  "/oak/stanford/groups/wjg/bparks/BPCells/04_data/test_real_data/ArchR_snyder_panc",
  "/oak/stanford/groups/wjg/bparks/BPCells/01_raw_data/peaksets/ENCFF305BGH.shuffle.bed",
  "8"
)
get_abspath <- function(path) {
  file.path(normalizePath(dirname(path)), basename(path))
}

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7)

input_glob <- args[1] # Glob pattern for input fragments.tsv.gz files
genome <- args[2] # Name of genome
min_frags <- as.numeric(args[3]) # Minimum fragments cutoff
min_tss <- as.numeric(args[4]) # Minimum TSS cutoff
output_dir <- get_abspath(args[5]) # Output directory
peak_path <- get_abspath(args[6]) # Input peak set
threads <- as.numeric(args[7])

output_arrow_dir <- file.path(output_dir, "arrows")
output_qc <- file.path(output_dir, "QualityControl")
output_log <- file.path(output_dir, "logs")
output_project <- file.path(output_dir, "project")

addArchRGenome(genome)
addArchRThreads(threads)

# Read peaks here to fail fast
peak_set <- read_tsv(peak_path, col_names = FALSE, col_types = NULL) %>%
  transmute(chr = X1, start = X2 + 1, end = X3) %>%
  filter(chr %in% seqnames(getGenomeAnnotation()$chromSizes)) %>%
  distinct() %>%
  GenomicRanges::makeGRangesFromDataFrame()

dir.create(output_arrow_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_log, recursive = TRUE, showWarnings = FALSE)
# Write the arguments passed while running the code
write_lines(args, file.path(output_dir, "arguments.txt"))
# Write the current script as well
cmdArgs <- commandArgs(trailingOnly = FALSE)
scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]

# Make sure these come out as absolute paths since we're going to change directories
# during actual Arrow creation
input_paths <- normalizePath(Sys.glob(input_glob))
sample_names <- basename(dirname(input_paths))
output_arrow_path <- file.path(normalizePath(output_arrow_dir), str_c(sample_names, ".arrow"))

# Change to an actual tempdir so that we can have a more reliable temporary filesystem on OAK
original_dir <- getwd()
tmp <- tempdir()
setwd(tmp)

arrow_path <- createArrowFiles(
  inputFiles = input_paths,
  sampleName = sample_names,
  outputNames = str_remove(output_arrow_path, fixed(".arrow")),
  excludeChr = c("chrM"),
  QCDir = output_qc,
  minTSS = min_tss,
  minFrags = min_frags,
  maxFrags = Inf,
  addTileMat = TRUE,
  TileMatParams = list(binarize = FALSE),
  addGeneScoreMat = TRUE,
  logFile = createLogFile("createArrowFiles", output_log)
)

setwd(original_dir)

proj <- ArchRProject(
  ArrowFiles = arrow_path,
  outputDirectory = output_project,
  copyArrows = FALSE,
  showLogo = FALSE,
  threads = 1
)

proj <- saveArchRProject(proj, output_project, logFile = createLogFile("saveProject_1", output_log))

# Add peak matrix
proj <- addPeakSet(proj, peak_set)

proj <- addPeakMatrix(
  ArchRProj = proj,
  ceiling = 1e9,
  binarize = FALSE,
  verbose = TRUE,
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addPeakMatrix", output_log)
)

proj <- saveArchRProject(proj, output_project, logFile = createLogFile("saveProject_2", output_log))

# Add dummy groups based on first base of cell barcode
barcode_group <- stringr::str_match(proj$cellNames, "#(.)")[, 2]
proj <- addCellColData(proj, barcode_group, "barcode_group", cells = proj$cellNames)

# Trying my darndest to prevent weird pseudoreplicate creation
cellGroups <- addGroupCoverages(
  proj,
  returnGroups = TRUE,
  groupBy = "barcode_group",
  minCells = 0,
  maxCells = 1e9,
  maxFragments = 1e10,
  maxRep = 1e9,
  logFile = createLogFile("addGroupCoverages", output_log)
)
saveRDS(cellGroups, file.path(output_project, "barcode_group_membership.Rds"), compress = FALSE)
proj <- addGroupCoverages(
  proj,
  groupBy = "barcode_group",
  minCells = 0,
  maxCells = 1e9,
  maxFragments = 1e10,
  maxRep = 1e9,
  logFile = createLogFile("addGroupCoverages", output_log)
)

proj <- saveArchRProject(proj, output_project, logFile = createLogFile("saveProject_3", output_log))

footprints <- getFootprints(
  proj,
  positions = GRangesList(TSS = getTSS(proj)),
  groupBy = "barcode_group",
  flank = 2000
)

saveRDS(footprints, file.path(output_project, "footprints.Rds"), compress = FALSE)
