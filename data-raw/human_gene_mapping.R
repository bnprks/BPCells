library(magrittr)

# Pull data from HGNC, and make a named vector that maps non-canonical gene names/
# symbols to their canonical names. Only unambiguous mappings will be stored
hgnc <- readr::read_tsv(
  "http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt",
  col_types=readr::cols(.default=readr::col_character())
)

human_gene_mapping <- dplyr::bind_rows(
  alias = dplyr::select(hgnc, symbol, alt = alias_symbol) %>% tidyr::separate_rows(alt, sep = "\\|"),
  prev = dplyr::select(hgnc, symbol, alt = prev_symbol) %>% tidyr::separate_rows(alt, sep = "\\|"),
  entrez = dplyr::select(hgnc, symbol, alt = entrez_id) %>% dplyr::mutate(alt = as.character(alt)),
  ensembl = dplyr::select(hgnc, symbol, alt = ensembl_gene_id),
  symbol = dplyr::transmute(hgnc, symbol, alt = symbol),
  .id = "source"
) %>%
  dplyr::filter(!is.na(alt), alt != "") %>%
  dplyr::group_by(alt) %>%
  dplyr::filter(dplyr::n_distinct(symbol) == 1 | source == "symbol") %>%
  dplyr::ungroup() %>%
  dplyr::distinct(symbol, alt) %>%
  dplyr::arrange(alt) %>%
  dplyr::pull(symbol, alt)

usethis::use_data(human_gene_mapping, overwrite = TRUE, compress="xz")
