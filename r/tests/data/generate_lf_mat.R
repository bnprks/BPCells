# Simple test script to generate the lf and crlf matrix files
# The mini_mat_crlf must be generated on Windows with an 
# older BPCells version (pre-CRLF fix from pull #257).
# Version 0.3.0 is sufficiently old
library(BPCells)

mat <- matrix(
  c(1, 2, 0, 3, 0,
    0, 0, 2, 2, 1,
    0, 0, 2, 0, 2),
  nrow=5, ncol=3
)
rownames(mat) <- paste0("row", seq_len(nrow(mat)))
colnames(mat) <- paste0("col", seq_len(ncol(mat)))

mat <- as(mat, "dgCMatrix") |>
  as("IterableMatrix")

write_matrix_dir(mat, "mini_mat_lf")