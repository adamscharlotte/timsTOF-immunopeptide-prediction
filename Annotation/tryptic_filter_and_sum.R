library(tidyverse)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
pool <- "TUM_first_pool_2"

base_path <- "/Users/adams/Projects/300K/2022-library-run/"
meta_path <- paste(base_path, "Metadata/tryptic-pool-sequence.txt", sep = "")
meta_qc_path <- paste(base_path, "Metadata/qc-peptides.txt", sep = "")
mapped_precursor_path <- paste(base_path, "Annotation/precursor-mapped/", pool, ".csv", sep = "")  # nolint

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

tbl_mapped_precursor <- fread(mapped_precursor_path) %>%
    as_tibble() %>%
    dplyr::rename(obs_sequence = Sequence) %>%
    spaceless()

# -----------------------------------------------------------------------------

tbl_pool_peptides <- fread(meta_path) %>%
    as_tibble() %>%
    filter(pool_name == pool) %>%
    dplyr::rename(Proteins = pool_name) %>%
    distinct() %>%
    select(Proteins, Sequence)

tbl_qc_peptides <- fread(meta_qc_path) %>%
    as_tibble() %>%
    select(Proteins, Sequence)

tbl_peptides <- rbind(tbl_pool_peptides, tbl_qc_peptides)

tbl_psms <- tbl_mapped_precursor %>%
    merge(tbl_peptides, by = "Proteins") %>%
    as_tibble()

tbl_mapped_precursor %>%
    filter(obs_sequence %in% tbl_peptides$Sequence)

tbl_filtered_psms <- tbl_psms %>%
    filter(stringr::str_ends(obs_sequence, trunc_sequence)) %>%
    select(Proteins, Scan_number, frame, Precursor, obs_sequence) %>%
    distinct()

tbl_filtered_precursors <- tbl_filtered_psms %>%
    select(Precursor, frame, obs_sequence) %>%
    distinct() %>%
    add_count(Precursor, frame) %>%
    filter(n == 1)

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