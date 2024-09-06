devtools::load_all("/mnt/c/Users/Immanuel/PycharmProjects/BPCells/r")
library(Rgperftools)
# Force devtools to compile without debug flags for optimal performance profiling
library(pkgbuild)
flags <- pkgbuild::compiler_flags(debug=FALSE)
new_compiler_flags <- function(debug=FALSE) {flags}
assignInNamespace('compiler_flags', new_compiler_flags, 'pkgbuild')

dir <- "/tmp/RtmpG0kkA9/file1f2a70a57b47"
dir <- withr::local_tempdir()
# These parameters should make it virtually certain we have
# counts of 0, 1, and >1 represented
chr1 <- tibble::tibble(
  chr = "chr1",
  start = sort(sample.int(1000, 1000, replace=TRUE)),
  end = start + sample.int(150, 1000, replace=TRUE),
  cell_id = sample(LETTERS, 1000, replace=TRUE)
)
chr2 <- tibble::tibble(
  chr = "chr2",
  start = sort(sample.int(1000, 1000, replace=TRUE)),
  end = start + sample.int(150, 1000, replace=TRUE),
  cell_id = sample(LETTERS, 1000, replace=TRUE)
)
frags <- convert_to_fragments(dplyr::bind_rows(chr1, chr2))
cell_groups <- dplyr::if_else(cellNames(frags) %in% c("A", "E", "I", "O", "U"), "vowel", "consonant")
frag_table <- dplyr::bind_rows(chr1, chr2) %>%
  dplyr::mutate(cell_group = dplyr::if_else(cell_id %in% c("A", "E", "I", "O", "U"), "vowel", "consonant"))
coverage_start <- frag_table %>%
  dplyr::mutate(end=start+1L) %>%
  dplyr::group_by(chr, start, end, cell_group) %>%
  dplyr::summarize(value = dplyr::n(), .groups="drop")

coverage_end <- frag_table %>%
  dplyr::mutate(start=end-1L) %>%
  dplyr::group_by(chr, start, end, cell_group) %>%
  dplyr::summarize(value = dplyr::n(), .groups="drop")

coverage <- dplyr::bind_rows(coverage_start, coverage_end) %>%
  dplyr::group_by(chr, start, end, cell_group) %>%
  dplyr::summarize(value = sum(value), .groups="drop")
answers <- list(
  "start_only" = coverage_start,
  "end_only" = coverage_end,
  "both" = coverage
)
# Test start + end insertions
for (mode in c("start_only", "end_only", "both")) {
  write_insertion_bedgraph(
    frags,
    c("vowel"=file.path(dir, "vowel.bg"), "consonant"=file.path(dir, "consonant.bg.gz")),
    cell_groups,
    mode
  )
  # Check for gzip encoding
  expect_false(readBin(file.path(dir, "vowel.bg"), "integer", size=2, endian="big") == 0x1f8b)
  expect_true(readBin(file.path(dir, "consonant.bg.gz"), "integer", size=2, endian="big") == 0x1f8b)
  result_vowel <- readr::read_tsv(file.path(dir, "vowel.bg"), col_names=c("chr", "start", "end", "value"), col_types="ciii")
  result_consonant <- readr::read_tsv(file.path(dir, "consonant.bg.gz"), col_names=c("chr", "start", "end", "value"), col_types="ciii")
}



data_dir <- file.path("/mnt/c/Users/Immanuel/PycharmProjects/BPCells/r/data/pbmc-3k")
setwd(data_dir)
# try with using pbmc data
frags <- open_fragments_10x("/mnt/c/Users/Immanuel/PycharmProjects/BPCells/r/data/pbmc-3k/pbmc_3k_10x.fragments.tsv.gz") %>%
  write_fragments_memory()
plot_fragment_length(frags)
genes <- read_gencode_transcripts(
  "./references", 
  release="42", 
  transcript_choice="MANE_Select",
  annotation_set = "basic", 
  features="transcript" # Make sure to set this so we don't get exons as well
)
# Check if we already ran import
if (!file.exists("pbmc_3k_rna_raw")) {
  mat_raw <- open_matrix_10x_hdf5("pbmc_3k_10x.h5", feature_type="Gene Expression") %>% 
    write_matrix_dir("pbmc_3k_rna_raw")
} else {
  mat_raw <- open_matrix_dir("pbmc_3k_rna_raw")
}
blacklist <- read_encode_blacklist("./references", genome="hg38")
atac_qc <- qc_scATAC(frags, genes, blacklist)
pass_atac <- atac_qc %>%
    dplyr::filter(nFrags > 1000, TSSEnrichment > 10) %>%
    dplyr::pull(cellName)
pass_rna <- colnames(mat_raw)[Matrix::colSums(mat_raw) > 1e3]
keeper_cells <- intersect(pass_atac, pass_rna)
frags <- frags %>% select_cells(keeper_cells)
keeper_genes <- Matrix::rowSums(mat_raw) > 3
mat <- mat_raw[keeper_genes,keeper_cells]
# Normalize by reads-per-cell
mat <- multiply_cols(mat, 1/Matrix::colSums(mat))

# Log normalization
mat <- log1p(mat * 10000) # Log normalization

stats <- matrix_stats(mat, row_stats="variance")
variable_genes <- order(stats$row_stats["variance",], decreasing=TRUE) %>% 
  head(1000) %>% 
  sort()

mat_norm <- mat[variable_genes,]
mat_norm <- mat_norm %>% write_matrix_memory(compress=FALSE)
gene_means <- stats$row_stats["mean",variable_genes]
gene_vars <- stats$row_stats["variance", variable_genes]
mat_norm <- (mat_norm - gene_means) / gene_vars
svd <- BPCells::svds(mat_norm, k=50)
# Alternate option: irlba::irlba(mat_norm, nv=50)
pca <- multiply_cols(svd$v, svd$d)

cat(sprintf("PCA dimensions: %s\n", toString(dim(pca))))
pca[1:4,1:3]
clusts <- knn_hnsw(pca, ef=500) %>% # Find approximate nearest neighbors
  knn_to_snn_graph() %>% # Convert to a SNN graph
  cluster_graph_louvain() # Perform graph-based clustering
umap <- uwot::umap(pca)
plot_embedding(clusts, umap)
cluster_annotations <- c(
  "1" = "T",
  "2" = "CD8-T",
  "3" = "B",
  "4" = "T",
  "5" = "NK",
  "6" = "Mono",
  "7" = "Mono",
  "8" = "Mono",
  "9" = "T",
  "10" = "DC",
  "11" = "Mono",
  "12" = "DC"
)
cell_types <- cluster_annotations[clusts]
cell_types_factor <- as.factor(cell_types)
levels(cell_types_factor)
paths <- file.path(data_dir, levels(cell_types_factor))
names(paths) <- levels(cell_types_factor)
# Save the frags object
write_fragments_dir(frags, file.path(data_dir, "frags"))
#save paths and cell tyupes
saveRDS(paths, file.path(data_dir, "paths.rds"))
saveRDS(cell_types, file.path(data_dir, "cell_types.rds"))



#################################################################
########## Downstream analysis of the data ######################
#################################################################

data_dir <- file.path("/mnt/c/Users/Immanuel/PycharmProjects/BPCells/r/data/pbmc-3k")
devtools::load_all("/mnt/c/Users/Immanuel/PycharmProjects/BPCells/r")
setwd(data_dir)
paths <- readRDS(file.path(data_dir, "paths.rds"))
paths_gz <- paste0(paths, ".bed.gz")
names(paths_gz) <- names(paths)
paths_bed <- paste0(paths, ".bed")
names(paths_bed) <- names(paths)
cell_types <- readRDS(file.path(data_dir, "cell_types.rds"))
frags <- open_fragments_dir(file.path(data_dir, "frags"))
profile_insertion_bed_writing <- function(frags, cell_types, paths, name, compress) {
  # add time
  time.start <- Sys.time()
  #start_profiler(name)
  call_macs_peaks(
    fragments=frags,
    macs_executable="macs3",
    cell_groups = cell_types,
    path = data_dir,
    effective_genome_size = 2.7e9,
    insertion_mode = "both",
    step = "all",
    threads = 4,
    compress_inputs = compress
  )
  #write_insertion_bed(frags, paths, cell_types, "both", threads = 1, verbose = TRUE)
  #stop_profiler()
  #time.end <- Sys.time()
  #time.taken <- time.end - time.start
  #print(time.taken)
}
profile_insertion_bed_writing(frags, cell_types, paths_bed, "/oak/stanford/groups/wjg/iabdi/profiling/insertion_bed_without_gz.out", compress = FALSE)
profile_insertion_bed_writing(frags, cell_types, paths_gz, "/oak/stanford/groups/wjg/iabdi/profiling/insertion_bed_with_gz.out", compress = TRUE)

# append .bed.gz to paths
paths <- paste0(paths, ".bed.gz")
call_macs_peaks(
  fragments=frags,
  cell_names = cell_types,
  genome_size = 2.7e9,
  path = data_dir,
  insertion_mode = "start_only",
  quiet = TRUE,
  step = "prep-inputs",
  threads = 4
)

call_macs_peaks(
  genome_size = 2.7e9,
  path = dir,
  step = "run-macs",
)

macs_output <- call_macs_peaks(
  genome_size = 2.7e9,
  path = dir,
  step = "read-outputs",
)

data_macs <- call_macs_peaks(
  fragments=frags,
  cell_names = cell_types,
  genome_size = 2.7e9,
  path = data_dir,
  insertion_mode = "start_only",
  step = "all",
  threads = 8
)
