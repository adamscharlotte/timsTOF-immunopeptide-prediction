# ssh cadams@10.152.135.57

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
# pool <- "TUM_lysn_33"

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
    select(OBS_SEQUENCE, PROTEINS) %>%
    distinct()

tbl_pool_all <- fread(meta_path) %>%
    as_tibble() %>%
    filter(pool_name == pool) %>%
    mutate(trunc_sequence = str_sub(Sequence, start = -7)) %>%
    dplyr::rename(PROTEINS = pool_name) %>%
    rename_with(toupper) %>%
    select(SEQUENCE, TRUNC_SEQUENCE, PROTEINS) %>%
    distinct()

# Remove all sequences that have the same truncated sequence
tbl_pool_trunc <- tbl_pool_all %>%
    add_count(TRUNC_SEQUENCE) %>%
    filter(n == 1) %>%
    select(-n)

tbl_pool_all %>%
    add_count(TRUNC_SEQUENCE) %>%
    filter(!n == 1) %>%
    select(-n) %>%
    arrange(TRUNC_SEQUENCE)

# Use the full length sequence as the truncated
tbl_pool_full <- tbl_pool_all %>%
    add_count(TRUNC_SEQUENCE) %>%
    filter(n > 1) %>%
    select(-c(n, TRUNC_SEQUENCE)) %>%
    mutate(TRUNC_SEQUENCE = SEQUENCE) %>%
    distinct()

tbl_pool_peptides <- rbind(tbl_pool_trunc, tbl_pool_full)

tbl_qc_peptides <- fread(meta_qc_path) %>%
    as_tibble() %>%
    rename_with(toupper) %>%
    select(SEQUENCE, TRUNC_SEQUENCE, PROTEINS)

tbl_peptides <- rbind(tbl_pool_peptides, tbl_qc_peptides)

tbl_mapped_sequences <- tbl_obs_sequence %>%
    merge(tbl_peptides) %>%
    as_tibble()

tbl_sequences_full <- tbl_mapped_sequences %>%
    filter(OBS_SEQUENCE == SEQUENCE) %>%
    select(OBS_SEQUENCE, SEQUENCE, PROTEINS)

tbl_sequences_trunc <- tbl_mapped_sequences %>%
    filter(!OBS_SEQUENCE %in% tbl_sequences_full$OBS_SEQUENCE) %>%
    filter(stringr::str_ends(OBS_SEQUENCE, TRUNC_SEQUENCE)) %>%
    select(OBS_SEQUENCE, SEQUENCE, PROTEINS)

tbl_sequences_trunc %>% distinct() %>% count(PROTEINS)

tbl_filtered_sequences <- rbind(tbl_sequences_full, tbl_sequences_trunc)

tbl_filtered_sequences %>%
    add_count(OBS_SEQUENCE) %>%
    filter(n > 1) %>%
    arrange(OBS_SEQUENCE)

output_path <- paste(base_path, "Annotation/full-length-map/", pool, ".csv", sep = "") # nolint
fwrite(tbl_filtered_sequences, output_path)