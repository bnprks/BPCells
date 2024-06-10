# Copyright 2023 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

suppressPackageStartupMessages({
  library(BPCells)
  library(magrittr)
})



args <- c(
  "/oak/stanford/groups/wjg/bparks/BPCells/01_raw_data/ENCODE/samples/snyder_*_PANC/fragments.tsv.gz",
  "hg38",
  "2000",
  "5",
  "/oak/stanford/groups/wjg/bparks/BPCells/04_data/test_real_data/bpcells_snyder_panc",
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

# Read peaks here to fail fast
peak_set <- BPCells::read_bed(peak_path) %>%
  dplyr::distinct()

# Get ArchR annotation data
genome_annotation <- getExportedValue("ArchR", stringr::str_c("genomeAnno", stringr::str_to_title(genome)))
gene_annotation <- getExportedValue("ArchR", stringr::str_c("geneAnno", stringr::str_to_title(genome)))

standard_chr <- GenomicRanges::seqnames(genome_annotation$chromSizes) %>% as.character()

#################################
# Import + merge fragments
#################################
cat(sprintf("Importing and merging fragments: %s\n", Sys.time()))
input_paths <- normalizePath(Sys.glob(input_glob))
sample_names <- basename(dirname(input_paths))
temp_frag_paths <- lapply(sample_names, tempfile)


fragments_list <- parallel::mclapply(seq_along(input_paths), mc.cores = threads, function(i) {
  open_fragments_10x(input_paths[i]) %>%
    prefix_cell_names(paste0(sample_names[i], "#")) %>%
    select_chromosomes(standard_chr) %>%
    write_fragments_dir(temp_frag_paths[[i]])
})
fragments <- do.call(c, fragments_list) %>% write_fragments_dir(file.path(output_dir, "fragments"))
for (p in temp_frag_paths) {
  unlink(p, recursive = TRUE)
}

#################################
# QC data + filter
#################################
cat(sprintf("Computing QC data and filtering: %s\n", Sys.time()))
qc_res <- fragments %>%
  shift_fragments(shift_end = -1) %>% # Match ArchR bug in fragment coords
  qc_scATAC(gene_annotation$TSS, genome_annotation$blacklist)
qc_res_promoter <- fragments %>%
  shift_fragments(shift_end = -1) %>% # Match ArchR bug in fragment coords
  qc_scATAC(gene_annotation$genes, genome_annotation$blacklist)
qc_res$ReadsInPromoter <- qc_res_promoter$ReadsInPromoter
readr::write_tsv(qc_res, file.path(output_dir, "cell_qc.tsv.gz"))

keeper_cells <- dplyr::filter(qc_res, TSSEnrichment >= min_tss, nFrags >= min_frags) %>%
  dplyr::pull(cellName)

fragments <- select_cells(fragments, keeper_cells) %>%
  write_fragments_dir(file.path(output_dir, "fragments_filtered"))

#################################
# Calculate tile matrix
#################################
cat(sprintf("Calculating tile matrix: %s\n", Sys.time()))
tile_coords <- genome_annotation$chromSizes %>%
  tibble::as_tibble() %>%
  dplyr::transmute(chr = seqnames, start = 0, end = end, tile_width = 500)

tile_mat <- fragments %>%
  shift_fragments(shift_start = 1) %>%
  # Mimic a quirk of ArchR's tile alignment
  tile_matrix(tile_coords) %>%
  write_matrix_dir(file.path(output_dir, "tile_mat"))

#################################
# Calculate peak matrix
#################################
cat(sprintf("Calculating peak matrix: %s\n", Sys.time()))
ordered_peaks <- BPCells:::order_ranges(peak_set, chrNames(fragments))
peak_matrix <- fragments %>%
  shift_fragments(shift_end = -1) %>%
  # Since ArchR incorrectly sets the end-coordinate of 10x fragments
  peak_matrix(peak_set[ordered_peaks, ]) %>%
  write_matrix_dir(file.path(output_dir, "peak_mat"))

#################################
# Calculate footprints
#################################
cat(sprintf("Calculating footprints: %s\n", Sys.time()))
groups <- stringr::str_c(
  stringr::str_match(cellNames(fragments), "#(.)")[, 2],
  "._.",
  stringr::str_split_fixed(cellNames(fragments), "#", 2)[, 1]
)
footprint_data <- fragments %>%
  shift_fragments(shift_end = -1) %>%
  footprint(
    gene_annotation$TSS,
    cell_groups = groups,
    flank = 2000
  )
readr::write_tsv(footprint_data, file.path(output_dir, "footprint.tsv.gz"))

#################################
# Calculate gene activity score
#################################
cat(sprintf("Calculating gene activity matrix: %s\n", Sys.time()))

# Shift gene and blacklist coordinates to account for ArchR's shift of tile coordinates
gene_regions <- gene_annotation$genes
gene_regions <- gene_regions[as.character(GenomicRanges::seqnames(gene_regions)) %in% standard_chr]
gene_regions <- gene_regions[!is.na(gene_regions$symbol)]
GenomicRanges::start(gene_regions) <- GenomicRanges::start(gene_regions) + 1
GenomicRanges::end(gene_regions) <- GenomicRanges::end(gene_regions) + 1

blacklist <- genome_annotation$blacklist
GenomicRanges::start(blacklist) <- GenomicRanges::start(blacklist) + 1
GenomicRanges::end(blacklist) <- GenomicRanges::end(blacklist) + 1

gene_score_matrix <- fragments %>%
  shift_fragments(shift_start = 1) %>% # Mimic a quirk of ArchR's tile alinment
  gene_score_archr(
    genes=gene_regions, chromosome_sizes = tile_coords, tile_width = 500, blacklist = blacklist,
    gene_name_column="symbol", addArchRBug = TRUE
  )
gene_score_matrix <- write_matrix_dir(gene_score_matrix, file.path(output_dir, "gene_score_mat"))


cat(sprintf("Done: %s\n", Sys.time()))