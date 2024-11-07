# Copyright 2024 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

# Set up k-v pairs for other tests
config_csv <- read.csv("config.csv")
config <- as.list(config_csv)$value
names(config) <- as.list(config_csv)$key

# Set up temp dir in case it's not already set
create_temp_dir <- function(dir = NULL) {
  if (is.null(dir)) {
    dir <- file.path(tempdir(), "lsi_test")
    if (dir.exists(dir)) unlink(dir, recursive = TRUE)
    dir.create(dir)
  }
  return(dir)
}