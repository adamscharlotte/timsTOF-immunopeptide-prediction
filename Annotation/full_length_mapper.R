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
    filter(str_length(OBS_SEQUENCE) <= str_length(SEQUENCE))

output_path <- paste(base_path, "Annotation/full-length-map/", pool, ".csv", sep = "") # nolint
fwrite(tbl_filtered_sequences, output_path)