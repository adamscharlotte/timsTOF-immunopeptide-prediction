# ssh cadams@10.152.135.57

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
# pool <- "TUM_HLA_16"

# base_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/" # nolint
base_path <- "/Users/adams/Projects/300K/2022-library-run/"
meta_path <- paste(base_path, "Metadata/full-pool-sequence.txt", sep = "")
meta_qc_path <- paste(base_path, "Metadata/qc-peptides.txt", sep = "")
txt_output_path <- paste(base_path, "Metadata/info-sum.txt", sep = "")
mapped_precursor_path <- paste(base_path, "Annotation/precursor-mapped/", pool, ".csv", sep = "")  # nolint

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

tbl_obs_sequence <- fread(mapped_precursor_path) %>%
    as_tibble() %>%
    dplyr::rename(obs_sequence = Sequence) %>%
    rename_with(toupper) %>%
    spaceless() %>%
    select(OBS_SEQUENCE) %>%
    distinct()

tbl_pool_all <- fread(meta_path) %>%
    as_tibble() %>%
    # filter(pool_name == pool) %>%
    mutate(trunc_sequence = str_sub(Sequence, start = -7)) %>%
    rename_with(toupper) %>%
    select(POOL_NAME, SEQUENCE, TRUNC_SEQUENCE) %>%
    distinct()

tbl_pool_trunc <- tbl_pool_all %>%
    add_count(POOL_NAME, TRUNC_SEQUENCE) %>%
    filter(n == 1) %>%
    select(-n)

tbl_pool_full <- tbl_pool_all %>%
    add_count(POOL_NAME, TRUNC_SEQUENCE) %>%
    filter(n > 1) %>%
    select(-c(n, TRUNC_SEQUENCE)) %>%
    mutate(TRUNC_SEQUENCE = POOL_NAME)

tbl_pool_peptides <- rbind(tbl_pool_trunc, tbl_pool_full)

tbl_qc_peptides <- fread(meta_qc_path) %>%
    as_tibble() %>%
    rename_with(toupper) %>%
    select(SEQUENCE, TRUNC_SEQUENCE)

tbl_peptides <- rbind(tbl_pool_peptides, tbl_qc_peptides)


tbl_peptides %>%
    count(TRUNC_SEQUENCE) %>%
    arrange((desc(n)))

tbl_mapped_sequences <- tbl_obs_sequence %>%
    merge(tbl_peptides) %>%
    as_tibble()

tbl_filtered_sequences <- tbl_mapped_sequences %>%
    filter(stringr::str_ends(OBS_SEQUENCE, TRUNC_SEQUENCE))

tbl_filtered_sequences %>%
    filter(!OBS_SEQUENCE == SEQUENCE) %>%
    print(n = 150)

tbl_filtered_sequences %>%
    count(OBS_SEQUENCE) %>%
    arrange((desc(n)))

tbl_filtered_sequences %>%
    filter(OBS_SEQUENCE == c("YELTQPPSV", "YVLTQPPSV")) %>%
    distinct()

tbl_annotation <- tbl_mapped_precursor %>%
    filter(obs_sequence %in% tbl_filtered_precursors$obs_sequence &
        Precursor %in% tbl_filtered_precursors$Precursor) %>%
    filter(Score >= 70) %>%
    rename_with(toupper) %>%
    group_by(PRECURSOR) %>%
    mutate(combined_INTENSITIES = paste0(INTENSITIES, collapse = ";")) %>%
    mutate(combined_MZ = paste0(MZ, collapse = ";")) %>%
    mutate(RETENTION_TIME = median(RETENTION_TIME)) %>%
    ungroup() %>%
    select(RAW_FILE, SCAN_NUMBER, MODIFIED_SEQUENCE, CHARGE, FRAGMENTATION,
    MASS_ANALYZER, SCAN_EVENT_NUMBER, MASS, SCORE, REVERSE, RETENTION_TIME,
    combined_MZ, combined_INTENSITIES, COLLISION_ENERGY, PRECURSOR) %>%
    distinct()

output_path <- paste(base_path, "Annotation/precursor-consensus/un-annotated/", pool, ".csv", sep = "") # nolint
fwrite(tbl_annotation, output_path)