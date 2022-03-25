devtools::load_all()
library(tidyverse)

open_fragments_10x("~/Downloads/pbmc_500/atac_pbmc_500_nextgem_fragments.tsv.gz") %>%
  write_fragments_dir("~/Downloads/pbmc_500/fragments")


f <- open_fragments_dir("pbmc_500/fragments") %>%
  shift_fragments(shift_end=1) %>%
  write_fragments_memory()
m <- open_matrix_10x_hdf5("pbmc_500/atac_pbmc_500_nextgem_raw_peak_bc_matrix.h5") 
mi <- t(m) %>% as("dgCMatrix") %>% .[cellNames(f),] %>% as("IterableMatrix") %>% as("mat_uint32_t")
peaks <- str_split_fixed(rownames(m), "-|:", n=3) %>%
  {tibble(chr=.[,1], start=as.integer(.[,2]), end=as.integer(.[,3]))}

m2 <- peakMatrix(f, peaks, zero_based_coords = TRUE)


matrix_identical_uint32_t_cpp(iterate_matrix(mi), iterate_matrix(m2))

start_profiler("profile_1.out")
s2 <- colSums(m2)
stop_profiler()



open_fragments_10x("pbmc_10k/atac_pbmc_10k_nextgem_fragments.tsv.gz") %>%
  write_fragments_dir("pbmc_10k/fragments")
f <- open_fragments_dir("pbmc_10k/fragments")

m <- open_matrix_10x_hdf5("pbmc_10k/atac_pbmc_10k_nextgem_raw_peak_bc_matrix.h5") 
mi <- t(m) %>% as("dgCMatrix") %>% .[cellNames(f),] %>% as("IterableMatrix") %>% as("mat_uint32_t")
peaks <- str_split_fixed(rownames(m), "-|:", n=3) %>%
  {tibble(chr=.[,1], start=as.integer(.[,2]), end=as.integer(.[,3]))}

m2 <- peakMatrix(f, peaks, zero_based_coords = TRUE)

matrix_identical_uint32_t_cpp(iterate_matrix(mi), iterate_matrix(m2))

start_profiler("profile_4.out")
s2 <- colSums(m2)
stop_profiler()

