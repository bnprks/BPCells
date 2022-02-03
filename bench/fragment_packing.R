devtools::load_all("~/Dropbox/greenleaf/playground/fragment_io/BPCells/")

input <- open_10x_fragments("~/Downloads/pbmc_unsorted_3k_raw_feature_bc_matrix.tar.gz") %>%
    write_fragments_memory()