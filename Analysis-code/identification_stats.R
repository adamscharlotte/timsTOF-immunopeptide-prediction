# ssh cadams@10.152.135.57

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
# pool <- "TUM_HLA2_88"

# base_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/" # nolint
base_path <- "/Users/adams/Projects/300K/2022-library-run/"
meta_path <- paste(base_path, "Metadata/full-pool-sequence.txt", sep = "")
meta_qc_path <- paste(base_path, "Metadata/qc-peptides.txt", sep = "")
mapped_precursor_path <- paste(base_path, "Annotation/precursor-mapped/", pool, ".csv", sep = "")  # nolint
txt_output_path <- paste(base_path, "Metadata/identification-summary.txt", sep = "") # nolint

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

tbl_obs_sequence <- fread(mapped_precursor_path) %>%
    as_tibble() %>%
    dplyr::rename(obs_sequence = Sequence) %>%
    rename_with(toupper) %>%
    spaceless() %>%
    select(OBS_SEQUENCE, PROTEINS, SCAN_NUMBER) %>%
    distinct()

tbl_pool_all <- fread(meta_path) %>%
    as_tibble() %>%
    filter(pool_name == pool) %>%
    dplyr::rename(PROTEINS = pool_name) %>%
    rename_with(toupper) %>%
    select(SEQUENCE, PROTEINS) %>%
    distinct()

tbl_qc_peptides <- fread(meta_qc_path) %>%
    as_tibble() %>%
    rename_with(toupper) %>%
    select(SEQUENCE, PROTEINS)

tbl_peptides <- rbind(tbl_pool_all, tbl_qc_peptides)

tbl_mapped_sequences <- tbl_obs_sequence %>%
    merge(tbl_peptides) %>%
    as_tibble()

tbl_filtered_sequences <- tbl_mapped_sequences %>%
    filter(endsWith(SEQUENCE, OBS_SEQUENCE)) %>%
    filter(str_length(OBS_SEQUENCE) <= str_length(SEQUENCE)) %>%
    group_by(OBS_SEQUENCE) %>%
    filter(str_length(SEQUENCE) == min(str_length(SEQUENCE))) %>%
    ungroup()

tbl_full_length <- tbl_filtered_sequences %>%
    filter(SEQUENCE == OBS_SEQUENCE) %>%
    select(-SCAN_NUMBER) %>%
    distinct()

tbl_n_trunc <- tbl_filtered_sequences %>%
    filter(endsWith(SEQUENCE, OBS_SEQUENCE)) %>%
    filter(str_length(OBS_SEQUENCE) < str_length(SEQUENCE)) %>%
    select(-SCAN_NUMBER) %>%
    distinct()

tbl_wrong_identification <- tbl_obs_sequence %>%
    filter(!OBS_SEQUENCE %in% tbl_filtered_sequences$OBS_SEQUENCE)

tbl_c_trunc <- tbl_mapped_sequences %>%
    filter(startsWith(SEQUENCE, OBS_SEQUENCE)) %>%
    filter(str_length(OBS_SEQUENCE) < str_length(SEQUENCE)) %>%
    group_by(OBS_SEQUENCE) %>%
    filter(str_length(SEQUENCE) == min(str_length(SEQUENCE))) %>%
    ungroup() %>%
    select(-SCAN_NUMBER) %>%
    distinct()

tbl_internal <- tbl_mapped_sequences %>%
    filter(!startsWith(SEQUENCE, OBS_SEQUENCE)) %>%
    filter(!endsWith(SEQUENCE, OBS_SEQUENCE)) %>%
    filter(str_detect(SEQUENCE, OBS_SEQUENCE)) %>%
    group_by(OBS_SEQUENCE) %>%
    filter(str_length(SEQUENCE) == min(str_length(SEQUENCE))) %>%
    ungroup() %>%
    select(-SCAN_NUMBER) %>%
    distinct()

tbl_wrong_identification %>% count(PROTEINS)

tbl_info <- tibble(pool_name = pool, pool_peptides = nrow(tbl_pool_all),
    total_psm = nrow(tbl_obs_sequence),
    pool_correct = nrow(filter(tbl_filtered_sequences, PROTEINS == pool)),
    pool_fl = nrow(filter(tbl_full_length, PROTEINS == pool)),
    pool_trunc = nrow(filter(tbl_n_trunc, PROTEINS == pool)),
    pool_wrong = nrow(filter(tbl_wrong_identification, PROTEINS == pool)),
    pool_cterm = nrow(filter(tbl_c_trunc, PROTEINS == pool)),
    pool_internal = nrow(filter(tbl_internal, PROTEINS == pool)),
    qc_peptides = nrow(tbl_qc_peptides),
    qc_correct = nrow(filter(tbl_filtered_sequences, !PROTEINS == pool)),
    qc_fl = nrow(filter(tbl_full_length, !PROTEINS == pool)),
    qc_trunc = nrow(filter(tbl_n_trunc, !PROTEINS == pool)),
    qc_wrong = nrow(filter(tbl_wrong_identification, !PROTEINS == pool)),
    qc_cterm = nrow(filter(tbl_c_trunc, !PROTEINS == pool)),
    qc_internal = nrow(filter(tbl_internal, !PROTEINS == pool)),
    rt_trunc = nrow(filter(tbl_n_trunc, PROTEINS == "QC_RT_QC_Peptide"))
)

fwrite(tbl_info, txt_output_path, append = TRUE)