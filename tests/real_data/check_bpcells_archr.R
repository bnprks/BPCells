# TODO in progress

archr_dir <- "/oak/stanford/groups/wjg/bparks/BPCells/04_data/test_real_data/ArchR_snyder_panc"
bpcells_dir <- "/oak/stanford/groups/wjg/bparks/BPCells/04_data/test_real_data/bpcells_snyder_panc"
min_tss <- 5 # Minimum TSS cutoff



##########################
# QC Stats Equality
##########################
archr_qc_path <- Sys.glob(file.path(archr_dir, "/QualityControl/*/*-Pre-Filter-Metadata.rds"))
qc_archr <- lapply(archr_qc_path, readRDS) %>%
    do.call(rbind, .) %>%
    tibble::as_tibble()

qc_bpcells <- readr::read_tsv(file.path(bpcells_dir, 'cell_qc.tsv.gz'))
stopifnot(length(setdiff(qc_archr$cellNames, qc_bpcells$cellName)) == 0)
rownames(qc_bpcells) <- qc_bpcells$cellName
qc_bpcells_compat <- qc_bpcells[qc_archr$cellNames,] %>%
    dplyr::mutate(TSSEnrichment = round(TSSEnrichment, 3)) %>%
    dplyr::rename(
        cellNames=cellName, 
        nMonoFrags=subNucleosomal,
        nDiFrags=monoNucleosomal,
        nMultiFrags=multiNucleosomal
    ) %>%
    dplyr::mutate(
        BlacklistRatio=ReadsInBlacklist/(2*nFrags),
        PromoterRatio=ReadsInPromoter/(2*nFrags),
        Keep=0 + (TSSEnrichment >= min_tss)
    ) %>% 
    dplyr::select(colnames(qc_archr))

stopifnot(all.equal(qc_bpcells_compat, qc_archr))


##########################
# Footprint Counts Equality
##########################

foot_archr_raw <- readRDS(file.path(archr_dir, "project/footprints.Rds"))
foot_archr_mat <- SummarizedExperiment::assay(foot_archr_raw)[SummarizedExperiment::rowData(foot_archr_raw)$type == "footprint",]

foot_archr <- tibble::as_tibble(foot_archr_mat) %>%
    dplyr::mutate(pos=-2000:2000) %>%
    tidyr::pivot_longer(!pos, names_to="group", values_to="value")

foot_bpcells <- readr::read_tsv(file.path(bpcells_dir, 'footprint.tsv.gz'))

foot <- dplyr::inner_join(foot_archr, foot_bpcells, by=c("pos", "group"), suffix=c(".archr", ".bpcells"))

stopifnot(all.equal(foot$value.archr, foot$value.bpcells))