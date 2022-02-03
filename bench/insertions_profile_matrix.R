# library(BPCells)
suppressPackageStartupMessages({
  library(testthat)
  library(patchwork)
  library(tidyverse)
  library(GenomicRanges)
})

devtools::load_all("~/Dropbox/greenleaf/playground/fragment_io/BPCells/")

source_dir <- "~/Downloads/PackedInsertionsData/"
# fragments_file <- file.path(source_dir, "01_raw_data/atac_pbmc_10k_2.0/fragments.tsv.gz")
# system.time({
#   fragments2 <- open_10x_fragments(fragments_file) %>% write_packed_fragments()
# })
# all.equal(fragments, fragments2)
fragments <- readRDS(file.path(source_dir, "atac_pbmc_10_2.0.rds"))
peaks <- readr::read_tsv(file.path(source_dir, "01_raw_data/atac_pbmc_10k_2.0/peaks.bed"),
                        comment = "#", col_names=c("chr", "start", "end"), col_types="cii")

read_10x_h5 <- function(path, feature_type="Peaks") {
  h5 <- hdf5r::H5File$new(path, mode="r")
  
  valid_feature_types <- unique(h5[["matrix/features/feature_type"]][])
  if (!all(feature_type %in% valid_feature_types))
    stop(sprintf("Feature type must be one of: %s", paste0(valid_feature_types, collapse=", ")))
  
  mat <- Matrix::sparseMatrix(
    i=h5[["matrix/indices"]][]+1, x=h5[["matrix/data"]][], p=h5[["matrix/indptr"]][], dims = h5[["matrix/shape"]][]
  )
  rownames(mat) <- h5[["matrix/features/id"]][]
  colnames(mat) <- h5[["matrix/barcodes"]][]
  mat <- mat[h5[["matrix/features/feature_type"]][] %in% feature_type,]
  h5$close()
  mat
}
time_op <- function(expr) {
  gc()
  mem_before <- system(sprintf("ps -o rss= %d", Sys.getpid()), intern = TRUE) %>% as.integer()
  res <- system.time(expr)
  mem_after <- system(sprintf("ps -o rss= %d", Sys.getpid()), intern = TRUE) %>% as.integer()
  gc()
  as.list(c(res, mem=mem_after-mem_before))
}

# Test peak matrix creation for correctness
correct_peaks <- read_10x_h5(
  file.path(source_dir, "01_raw_data/atac_pbmc_10k_2.0/raw_peak_bc_matrix.h5"),
  "Peaks"
) %>% t()

for (version in 1:5) {
  cat(sprintf("Checking correctness of peaks version %d\n", version))
  mat <- fragments %>% 
    shift_fragments(shift_end=1) %>%
    select_cells(rownames(correct_peaks)) %>%
    overlapMatrix(as.list(peaks), convert_to_0_based_coords = FALSE,
                version = version) %>% as("dgCMatrix")
  stopifnot(all.equal(mat, correct_peaks, check.attributes = FALSE))
}
mat <- fragments %>% 
  shift_fragments(shift_end=1) %>% 
  select_cells(rownames(correct_peaks)) %>%
  tileMatrix(as.list(peaks), peaks$end - peaks$start, convert_to_0_based_coords = FALSE) %>%
  as("dgCMatrix")
stopifnot(all.equal(mat, correct_peaks, check.attributes = FALSE))

mat <- fragments %>% 
  shift_fragments(shift_end=1) %>% 
  select_cells(rownames(correct_peaks)) %>%
  tileMatrix(as.list(peaks), peaks$end - peaks$start, convert_to_0_based_coords = FALSE,
             version = 2) %>%
  as("dgCMatrix")
stopifnot(all.equal(mat, correct_peaks, check.attributes = FALSE))

rm(mat, correct_peaks)

peak_timing <- NULL
for (version in 1:5) {
  cat(sprintf("Checking performance of peaks version %d\n", version))
  timing <- time_op({
    overlapMatrix(fragments, as.list(peaks), convert_to_0_based_coords = FALSE, 
                  version=version) %>% matrixStats()
  })
  timing$version <- as.character(version)
  peak_timing <- bind_rows(peak_timing, timing)
}
timing <- time_op({
  tileMatrix(fragments, as.list(peaks), peaks$end - peaks$start, convert_to_0_based_coords = FALSE) %>%
    matrixStats()
})
timing$version <- "tile"
peak_timing <- bind_rows(peak_timing, timing)

timing <- time_op({
  tileMatrix(fragments, as.list(peaks), peaks$end - peaks$start, convert_to_0_based_coords = FALSE,
             version=2) %>%
    matrixStats()
})
timing$version <- "tile2"
peak_timing <- bind_rows(peak_timing, timing)

peak_timing_plot <- ggplot(peak_timing, aes(version, elapsed)) +
  geom_point() +
  ylim(0, NA) +
  labs(y="time (seconds)")
peak_mem_plot <- ggplot(peak_timing, aes(version, mem/1e3)) +
  geom_point() +
  ylim(0, NA) +
  labs(y="memory (MB)")

(peak_timing_plot + peak_mem_plot) +
  patchwork::plot_annotation(title="Peak matrix creation benchmarking")


# Test 500bp tile matrix for correctness
chr_sizes <- read_tsv(file.path(source_dir, "01_raw_data/hg38.chrom.sizes"),
                      col_names=c("chr", "size"), col_types="ci") %>%
  filter(chr %in% chrNames(fragments))
chrs <- GenomicRanges::GRanges(chr_sizes$chr, IRanges::IRanges(1, chr_sizes$size))
chr_tiles <- GenomicRanges::slidingWindows(chrs, width=500, step=500) %>% unlist()

# Test correctness on just chr1 + chr10
ref_ans <- overlapMatrix(fragments, chr_tiles[seqnames(chr_tiles) %in% c("chr1", "chr10")], 
                         version = 1) %>%
  as("dgCMatrix")

for (version in 2:5) {
  cat(sprintf("Testing correctness on version %d\n", version))
  alt_ans <- overlapMatrix(fragments, chr_tiles[seqnames(chr_tiles) %in% c("chr1", "chr10")], 
                           version = version) %>%
    as("dgCMatrix")
  stopifnot(all.equal(ref_ans, alt_ans))
}
alt_ans <- fragments %>% 
  #shift_fragments(shift_start = 0, shift_end = -1) %>%
  tileMatrix(chrs[seqnames(chrs) %in% c("chr1", "chr10")], tile_width = 500) %>%
  as("dgCMatrix")
stopifnot(all.equal(ref_ans, alt_ans))

alt_ans <- fragments %>% 
  #shift_fragments(shift_start = 0, shift_end = -1) %>%
  tileMatrix(chrs[seqnames(chrs) %in% c("chr1", "chr10")], tile_width = 500, version=2) %>%
  as("dgCMatrix")
stopifnot(all.equal(ref_ans, alt_ans))

rm(ref_ans, alt_ans)


tile_timing <- NULL
for (version in 1:5) {
  cat(sprintf("Checking performance of tiles version %d\n", version))
  timing <- time_op({
    overlapMatrix(fragments, chr_tiles, 
                  version=version) %>% matrixStats()
  })
  timing$version <- as.character(version)
  tile_timing <- bind_rows(tile_timing, timing)
}
timing <- time_op({
  tileMatrix(fragments, chrs, tile_width=500, convert_to_0_based_coords = FALSE) %>%
    matrixStats()
})
timing$version <- "tile"
tile_timing <- bind_rows(tile_timing, timing)

timing <- time_op({
  tileMatrix(fragments, chrs, tile_width=500, convert_to_0_based_coords = FALSE,
             version=2) %>%
    matrixStats()
})
timing$version <- "tile2"
tile_timing <- bind_rows(tile_timing, timing)


tile_timing_plot <- ggplot(tile_timing, aes(version, elapsed)) +
  geom_point() +
  ylim(0, NA) +
  labs(y="time (seconds)")
tile_mem_plot <- ggplot(tile_timing, aes(version, mem/1e3)) +
  geom_point() +
  ylim(0, NA) +
  labs(y="memory (MB)")

(tile_timing_plot + tile_mem_plot) +
  patchwork::plot_annotation(title="Tile matrix creation benchmarking")




# Test timing + correctness of insertions iterator
timing1 <- time_op({ res1 <- scan_fragments_modulo_cpp(iterate_fragments(fragments)) })
timing2 <- time_op({ res2 <- scan_insertions_cpp(iterate_fragments(fragments)) })
timing3 <- time_op({ res3 <- scan_insertions2_cpp(iterate_fragments(fragments)) })
stopifnot(res1[4] == res2[3])
dplyr::bind_rows(timing1, timing2, timing3)

system.time({
start_profiler("/Users/ben/Dropbox/greenleaf/playground/fragment_io/insertions2_profile.prof")
scan_insertions2_cpp(iterate_fragments(fragments))
stop_profiler()
})

start_profiler("insertions2_profile.prof")
scan_insertions2_cpp(iterate_fragments(fragments))
stop_profiler()

start_profiler("insertions1_profile.prof")
scan_insertions_cpp(iterate_fragments(fragments))
stop_profiler()

time_op({overlapMatrix(fragments, chr_tiles, version = 1) %>% as("dgCMatrix")}) %>%
  as_tibble()
x <- time_op({
  overlapMatrix(fragments, as.list(peaks), convert_to_0_based_coords = FALSE) %>%
    matrixStats("none", "none")
})

flank <- 2000
tss_ranges <- ArchR::geneAnnoHg38$TSS %>% unique() %>%
  GenomicRanges::resize(2*flank + 1, fix="center")

time_op({
  tileMatrix(fragments, sample(tss_ranges, 5), tile_width = 1) %>% matrixStats("none", "none")
})

peaks_to_use <- 1:100000
system.time({
  res1 <- overlapMatrix(fragments, as.list(peaks[peaks_to_use,]), convert_to_0_based_coords = FALSE) %>%
    matrixStats("none", "none")
})

#devtools::load_all("~/Dropbox/greenleaf/playground/fragment_io/BPCells/")
system.time({
  res2 <- overlapMatrix(fragments, as.list(peaks[peaks_to_use,]), convert_to_0_based_coords = FALSE,
                        version = 2) %>%
    matrixStats("none", "none")
})

system.time({
  res3 <- overlapMatrix(fragments, as.list(peaks[peaks_to_use,]), convert_to_0_based_coords = FALSE,
                        version = 3) %>%
    matrixStats("none", "none")
})





time_op({
  junk <- overlapMatrix(fragments, tss, version=4) %>% matrixStats("none", "none")
})

results_list <- list()
for (flank in c(0, 2)) { #, 10, 20)) {
  tss <- ArchR::geneAnnoHg38$TSS %>% unique() %>%
    GenomicRanges::resize(2*flank + 1, fix="center") %>%
    GenomicRanges::slidingWindows(1) %>%
    unlist()
  for (version in 1:4) {
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
  #results_list_old <- results_list
results_list <- lapply(results_list, function(x) {
  x <- as.list(x); x[["version"]] <- as.character(x[["version"]]); 
  #x$user.self <- as.num
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

time_plot + mem_plot + 
  patchwork::plot_layout(guides="collect") +
  patchwork::plot_annotation(title="Scaling TSS profile by flanking basepairs")




time_op({
  overlapMatrix(fragments, tss, version=1) %>% matrixStats("none", "none")
})
system.time({
  overlapMatrix(fragments, tss, version=1) %>% matrixStats("none", "none")
})

gc()
mem_before <- system(sprintf("ps -o rss= %d", Sys.getpid()), intern = TRUE) %>% as.integer()
system.time({
  overlapMatrix(fragments, tss, version=2) %>% matrixStats("none", "none")
})
mem_after <- system(sprintf("ps -o rss= %d", Sys.getpid()), intern = TRUE) %>% as.integer()
gc()
print(mem_after - mem_before)

gc()
mem_before <- system(sprintf("ps -o rss= %d", Sys.getpid()), intern = TRUE) %>% as.integer()
system.time({
  overlapMatrix(fragments, tss, version=3) %>% matrixStats("none", "none")
})
mem_after <- system(sprintf("ps -o rss= %d", Sys.getpid()), intern = TRUE) %>% as.integer()
gc()


flank <- 1
tss <- ArchR::geneAnnoHg38$TSS %>% unique() %>%
  GenomicRanges::resize(2*flank + 1, fix="center") %>%
  GenomicRanges::slidingWindows(1) %>%
  unlist()

tss_ranges <- ArchR::geneAnnoHg38$TSS %>% unique() %>%
  GenomicRanges::resize(2*flank + 1, fix="center")


m1 <- overlapMatrix(fragments, tss[1:3000,], version=1) %>% as("dgCMatrix") 
m2 <- overlapMatrix(fragments, tss[1:3000,], version=2) %>% as("dgCMatrix") 
m3 <- overlapMatrix(fragments, tss[1:3000,], version=3) %>% as("dgCMatrix") 
m4 <- overlapMatrix(fragments, tss[1:3000,], version=4) %>% as("dgCMatrix") 
m5 <- overlapMatrix(fragments, tss[1:3000,], version=5) %>% as("dgCMatrix") 
m6 <- tileMatrix(fragments, tss_ranges[1:1000,], tile_width = 1) %>% as("dgCMatrix")
all.equal(m1, m2)
all.equal(m1, m3)
all.equal(m1, m4)
all.equal(m1, m5)
all.equal(m1, m6)



groups_matrix[groups_indices] <- cell_normalization_factors[groups_indices[,2]]

# Get a matrix of cell_groups x all positions
motif_bp <- unlist(GenomicRanges::slidingWindows(motif_positions, 1))
mat_all <- groups_matrix %*% overlapMatrix(fragments, motif_bp)


mat1 <- overlapMatrix(fragments, as.list(peaks[1:1000,]), convert_to_0_based_coords = FALSE) %>%
  as("dgCMatrix")
mat2 <- overlapMatrix(fragments, as.list(peaks[1:1000,]), convert_to_0_based_coords = FALSE,
                      version=2) %>%
  as("dgCMatrix")
mat3 <- overlapMatrix(fragments, as.list(peaks[1:1000,]), convert_to_0_based_coords = FALSE,
                      version=3) %>%
  as("dgCMatrix")
all.equal(mat1, mat2)
all.equal(mat1, mat3)

x <- fragments %>% select_chromosomes("chr1") %>% as("GRanges")

# Got next coord: 629935, last_ouptut: 634019
sub <- subsetByOverlaps(x, GRanges("chr1", ranges=IRanges(start=628000, end=635000)))

sub[start(sub) %in% c(629935, 634019)]

generate_sparse_matrix <- function(nrow, ncol, fraction_nonzero=0.5, max_val=10) {
  m <- matrix(rbinom(nrow*ncol, 1, fraction_nonzero)*sample.int(max_val, nrow*ncol, replace=TRUE), nrow=nrow)
  as(m, "dgCMatrix")
}
generate_dense_matrix <- function(nrow, ncol) {
  m <- matrix(runif(nrow*ncol), nrow=nrow)
}

to_matrix <- function(x) {
  x <- as.matrix(x)
  attr(x, "dimnames") <- NULL
  x
}
to_vector <- function(x) {
  x <- as.numeric(x)
  attributes(x) <- NULL
  x
}

tf_idf_manual <- function(x, scaleTo=1e4) {
  x <- t( t(x/rowMeans(x)) /colSums(x) )
  log(1 + scaleTo * x)
}


withr::local_seed(195123)
t <- Matrix::t

withr::local_seed(195123)
t <- Matrix::t

m1 <- generate_sparse_matrix(100, 200)
m2 <- t(m1)

i1 <- as(m1, "IterableMatrix")
i2 <- t(i1)

i1_tf <- normalize_TFIDF(i1)
i2_tf <- normalize_TFIDF(i2, scaleTo=5)

m1_tf <- tf_idf_manual(as.matrix(m1))
m2_tf <- tf_idf_manual(as.matrix(m2), scaleTo=5)

b <- generate_dense_matrix(200, 3)
b2 <- generate_dense_matrix(100, 3)

expect_equal(i1_tf %*% b, m1_tf %*% b)
expect_equal(t(b2) %*% i1_tf, t(b2) %*% m1_tf)

expect_equal(t(b) %*% i2_tf, t(b) %*% m2_tf)
expect_equal(i2_tf %*% b2, m2_tf %*% b2) 



x <- t(b) %*% i2_tf
y <- t(b) %*% m2_tf

I <- matrix(0, 5, 5)
diag(I) <- 1

m1 <- generate_sparse_matrix(20, 5)
m2 <- t(m1)

i1 <- as(m1, "IterableMatrix")
i2 <- t(i1)

i1_tf <- normalize_TFIDF(i1)
i2_tf <- normalize_TFIDF(i2, scaleTo=5)

m1_tf <- tf_idf_manual(as.matrix(m1))
m2_tf <- tf_idf_manual(as.matrix(m2), scaleTo=5)


all.equal(m1_tf, i1_tf %*% I)
all.equal(m2_tf, I %*% i2_tf)


flank <- 0
tss <- ArchR::geneAnnoHg38$TSS %>% unique() %>%
  GenomicRanges::resize(2*flank + 1, fix="center") %>%
  slidingWindows(1) %>%
  unlist()

tssY <- tss[seqnames(tss) == "chrY"]

chrY <- fragments %>% select_chromosomes("chrY") %>% as("GRanges")
subsetByOverlaps(chrY, as("chrY:26411100-26411200", "GRanges"))
