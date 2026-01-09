# Copyright 2026 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

# Pulls transcripts from gencode, filters to only chr4 and saves as package data.
# Also runs qc_scATAC on the demo fragments filtered to chr4 and saves as package data.

library(dplyr)
library(BPCells)


# Get transcripts
transcripts <- read_gencode_transcripts(
  file.path(tempdir(), "references"), release = "42",
  annotation_set = "basic",
  features = "transcript"
)
transcripts_filtered_example_chr_4 <- transcripts %>% dplyr::filter(chr %in% "chr4")



frags <- get_demo_frags() %>% select_chromosomes("chr4")
blacklist <- read_encode_blacklist(file.path(tempdir(), "references"), genome="hg38")
qc_results_filtered_example_chr_4 <- qc_scATAC(frags, transcripts_filtered_example_chr_4, blacklist)
readr::write_delim(
  transcripts_filtered_example_chr_4, file.path("./inst/extdata/transcripts_filtered_example_chr_4.tsv.gz"), 
  delim = "\t"
)
readr::write_delim(
  qc_results_filtered_example_chr_4, file.path("./inst/extdata/qc_results_filtered_example_chr_4.tsv.gz"), 
  delim = "\t"
)

