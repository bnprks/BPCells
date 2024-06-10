# Copyright 2022 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

# Pull gene lists from MGI, and make a named vector that maps non-canonical gene names/
# symbols to their canonical names.
# Since this involves more file stitching than for human genes, there's also
# more checks that the data is self-consistent

library(dplyr)

# A full list of non-withdrawn gene IDs
mgi_list <- readr::read_tsv(
  "http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt"
) %>%
  select(mgi_id = "MGI Accession ID", symbol = "Marker Symbol", status = "Status", name = "Marker Name", type = "Marker Type", alias = "Marker Synonyms (pipe-separated)")


# Contains some withdrawn IDs, and a superset of the current mgi_list
mgi_entrez <- readr::read_tsv(
  "http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt",
  col_names = c("mgi_id", "symbol", "status", "name", "cM_pos", "chr", "type", "secondary_ids", "entrez_id", "alias", "feature_type", "start", "end", "strand", "biotype")
) %>%
  select(mgi_id, symbol, status, name, secondary_ids, entrez_id, alias) %>%
  filter(status == "O")

# Check that symbol names match between the two sources
stopifnot(all.equal(sort(mgi_list$symbol), sort(mgi_entrez$symbol)))

# Check that given aliases match
mgi_list_mapping <- select(mgi_list, symbol, alt = alias) %>%
  tidyr::separate_rows(alt, sep = "\\|") %>%
  filter(!is.na(alt)) %>%
  arrange(symbol, alt)

mgi_entrez_mapping <- select(mgi_entrez, symbol, alt = alias) %>%
  tidyr::separate_rows(alt, sep = "\\|") %>%
  filter(!is.na(alt)) %>%
  arrange(symbol, alt)
stopifnot(all.equal(mgi_list_mapping, mgi_entrez_mapping))

mgi_ensembl <- readr::read_tsv(
  "http://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt",
  col_names = c("mgi_id", "symbol", "name", "cM_pos", "chr", "ensembl_id", "ensembl_transcript", "ensembl_prot", "feature_type", "start", "end", "strand", "biotypes")
) %>% select(mgi_id, symbol, name, ensembl_id)

# Check that ensembl symbols are a subset of the full list
stopifnot(length(setdiff(mgi_ensembl$symbol, mgi_entrez$symbol)) == 0)

mouse_gene_mapping <- bind_rows(
  alias = mgi_entrez_mapping,
  entrez = mgi_entrez %>% transmute(symbol, alt = as.character(entrez_id)),
  ensembl = mgi_ensembl %>% select(symbol, alt = ensembl_id),
  symbol = transmute(mgi_entrez_mapping, symbol, alt = symbol),
  .id = "source"
) %>%
  filter(!is.na(alt)) %>%
  dplyr::group_by(alt) %>%
  dplyr::filter(dplyr::n_distinct(symbol) == 1 | source == "symbol") %>%
  dplyr::ungroup() %>%
  dplyr::distinct(symbol, alt) %>%
  dplyr::arrange(alt) %>%
  dplyr::pull(symbol, alt)

usethis::use_data(mouse_gene_mapping, overwrite = TRUE, compress="xz")
