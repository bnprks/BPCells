# library(BPCells)
suppressPackageStartupMessages({
  library(testthat)
  library(patchwork)
  library(tidyverse)
})

devtools::load_all("~/Dropbox/greenleaf/playground/fragment_io/BPCells/")

source_dir <- "/Volumes/BP Archive 21/PackedInsertionsData/"
# fragments_file <- file.path(source_dir, "01_raw_data/atac_pbmc_10k_2.0/fragments.tsv.gz")
# system.time({
#   fragments2 <- open_10x_fragments(fragments_file) %>% write_packed_fragments()
# })
# all.equal(fragments, fragments2)
fragments <- readRDS(file.path(source_dir, "atac_pbmc_10_2.0.rds"))

time_op <- function(expr) {
  gc()
  mem_before <- system(sprintf("ps -o rss= %d", Sys.getpid()), intern = TRUE) %>% as.integer()
  res <- system.time(expr)
  mem_after <- system(sprintf("ps -o rss= %d", Sys.getpid()), intern = TRUE) %>% as.integer()
  gc()
  as.list(c(res, mem=mem_after-mem_before))
}

results_list <- list()
for (flank in c(0, 2, 10, 20)) {
#for (flank in c(0, 1)) {
  tss <- ArchR::geneAnnoHg38$TSS %>% unique() %>%
    GenomicRanges::resize(2*flank + 1, fix="center") %>%
    GenomicRanges::slidingWindows(1) %>%
    unlist()
  for (version in 1:6) {
    cat(sprintf("Running flank %d, version %d\n", flank, version))
    res <- time_op({
      overlapMatrix(fragments, tss, version=version) %>% matrixStats("none", "none")
    })
    results_list <- c(results_list, list(c(res, version=as.character(version), flank=flank)))
  }
  tss_ranges <- ArchR::geneAnnoHg38$TSS %>% unique() %>%
    GenomicRanges::resize(2*flank + 1, fix="center")
  cat(sprintf("Running flank %d, version tile\n", flank))
  res <- time_op({
    tileMatrix(fragments, tss_ranges, tile_width = 1) %>% matrixStats("none", "none")
  })
  res$version <- "tile"
  res$flank <- flank
  results_list <- c(results_list, list(res))
}

results_list <- lapply(results_list, function(x) {
  x <- as.list(x); x[["version"]] <- as.character(x[["version"]]); 
  x
})
time_plot <- bind_rows(results_list) %>% 
  ggplot(aes(flank, elapsed, color=as.character(version))) +
  geom_point() +
  labs(y="Time (seconds)", color="version")

mem_plot <- bind_rows(results_list) %>% 
  ggplot(aes(flank, mem/1e6, color=as.character(version))) +
  geom_point() +
  labs(y="Memory (GB)", color="version")

plot <- time_plot + mem_plot + 
  patchwork::plot_layout(guides="collect") +
  patchwork::plot_annotation(title="Scaling TSS profile by flanking basepairs")

ggsave("timing_plot.png", plot, width=7, height=3)
write_csv(bind_rows((results_list)), "timing_data.csv")
flank <- 1
tss <- ArchR::geneAnnoHg38$TSS %>% unique() %>%
  GenomicRanges::resize(2*flank + 1, fix="center") %>%
  GenomicRanges::slidingWindows(1) %>%
  unlist()

tss_ranges <- ArchR::geneAnnoHg38$TSS %>% unique() %>%
  GenomicRanges::resize(2*flank + 1, fix="center")

cat("Checking correctness\n")

m1 <- overlapMatrix(fragments, tss[1:3000,], version=1) %>% as("dgCMatrix") 
m2 <- overlapMatrix(fragments, tss[1:3000,], version=2) %>% as("dgCMatrix") 
m3 <- overlapMatrix(fragments, tss[1:3000,], version=3) %>% as("dgCMatrix") 
m4 <- overlapMatrix(fragments, tss[1:3000,], version=4) %>% as("dgCMatrix") 
m5 <- tileMatrix(fragments, tss_ranges[1:1000,], tile_width = 1) %>% as("dgCMatrix")
stopifnot(all.equal(m1, m2))
stopifnot(all.equal(m1, m3))
stopifnot(all.equal(m1, m4))
stopifnot(all.equal(m1, m5))


