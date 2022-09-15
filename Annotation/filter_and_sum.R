# ssh cadams@10.152.135.57

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(version = '3.12')
# BiocManager::install("MSnbase")

library(tidyverse)
library(data.table)
# library(MSnbase)
# library(timsr)

args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
# pool <- "TUM_HLA_16"

base_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/" # nolint
# base_path <- "/Users/adams/Projects/300K/2022-library-run/"
meta_path <- paste(base_path, "Metadata/full-pool-sequence.txt", sep = "")
meta_qc_path <- paste(base_path, "Metadata/qc-peptides.txt", sep = "")
txt_output_path <- paste(base_path, "Metadata/info-sum.txt", sep = "")
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
    mutate(trunc_sequence = str_sub(Sequence, start = -7)) %>%
    dplyr::rename(Proteins = pool_name) %>%
    distinct()

tbl_qc_peptides <- fread(meta_qc_path) %>% as_tibble()

tbl_peptides <- rbind(tbl_pool_peptides, tbl_qc_peptides)

tbl_psms <- tbl_mapped_precursor %>%
    merge(tbl_peptides, by = "Proteins") %>%
    as_tibble()

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
    select(RAW_FILE, SCAN_NUMBER, MODIFIED_SEQUENCE, CHARGE, FRAGMENTATION,
    MASS_ANALYZER, SCAN_EVENT_NUMBER, MASS, SCORE, REVERSE, RETENTION_TIME,
    MZ, INTENSITIES, COLLISION_ENERGY, PRECURSOR, FRAME)

output_path <- paste(base_path, "Annotation/precursor-consensus/un-annotated/", pool, ".csv", sep = "") # nolint
fwrite(tbl_annotation, output_path)

colnames(tbl_annotation)

tbl_annotation_list <- tbl_annotation %>%
    mutate(MZ_list = as.list(strsplit(MZ, ";"))) %>%
    mutate(INTENSITIES_list = as.list(strsplit(INTENSITIES, ";"))) %>%
    select(MZ_list, INTENSITIES_list, MODIFIED_SEQUENCE, PRECURSOR) %>%
    distinct()

tbl_precursor_mz <- tbl_annotation_list %>%
    filter(PRECURSOR == 19980) %>%
    mutate(spectrum_name = as.character(PRECURSOR))




tbl_precursor_mz %>%
    mutate(spectrum = new("Spectrum1", mz = as.numeric(unlist(MZ_list)), intensity = as.numeric(unlist(INTENSITIES_list))))

sp1 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.01),
    intensity = ints1)
tbl_precursor_mz %>%
    mutate(MZ_num = as.numeric(unlist(MZ_list)))
class(tbl_precursor_mz$MZ_list)

x <- list(as.character(1:3))
x

# -----------------------------------------------------------------------------

tbl_qc_psms <- tbl_mapped_precursor %>%
    # filter(!Proteins == pool) %>%
    filter(stringr::str_starts(Proteins, "QC")) %>%
    distinct()

tbl_meta <- fread(meta_path) %>% as_tibble()
tbl_pool <- tbl_meta %>%
    filter(pool_name == pool) %>%
    mutate(trunc_sequence = str_sub(Sequence, start = -7)) %>%
    dplyr::rename(Sequence_pool = Sequence) %>%
    distinct()

tbl_qc_peptides <- fread(meta_qc_path) %>% as_tibble()

tbl_filtered_qc_frames <- tbl_mapped_precursor %>%
    filter(stringr::str_starts(Proteins, "QC")) %>%
    filter(Sequence %in% tbl_qc_peptides$Sequence) %>%
    select(Scan_number, frame, Precursor, Sequence) %>%
    distinct()

tbl_pool_psms <- tbl_mapped_precursor %>%
    mutate(pool_name = pool) %>%
    filter(Proteins == pool) %>%
    merge(tbl_pool, by = "pool_name") %>%
    as_tibble()

tbl_filtered_pool_frames <- tbl_pool_psms %>%
    filter(stringr::str_ends(Sequence, trunc_sequence)) %>%
    select(Scan_number, frame, Precursor, Sequence) %>%
    distinct()

tbl_psms <- tbl_mapped_precursor %>% select(Scan_number, Sequence) %>% distinct()
tbl_pool_psms <- tbl_pool_psms %>% select(Scan_number, Sequence) %>% distinct()
tbl_qc_psms <- tbl_qc_psms %>% select(Scan_number, Sequence) %>% distinct()
tbl_filtered_pool_psms <- tbl_filtered_pool_frames %>% select(Scan_number, Sequence) %>% distinct()
tbl_filtered_qc_psms <- tbl_filtered_qc_frames %>% select(Scan_number, Sequence) %>% distinct()
tbl_unique_pool_peptides <- tbl_filtered_pool_frames %>% select(Sequence) %>% distinct()
tbl_full_length <- tbl_filtered_pool_frames %>% filter(Sequence %in% tbl_pool$Sequence_pool) %>%
    select(Sequence) %>% distinct
tbl_filtered_qc_precursors <- tbl_filtered_qc_frames %>%
    select(Precursor, frame, Sequence) %>%
    distinct() %>%
    add_count(Precursor, frame) %>%
    filter(n == 1)

tbl_filtered_pool_precursors <- tbl_filtered_pool_frames %>%
    select(Precursor, frame, Sequence) %>%
    distinct() %>%
    add_count(Precursor, frame) %>%
    filter(n == 1)

tbl_filtered_precursors <- rbind(tbl_filtered_qc_precursors, tbl_filtered_pool_precursors)

# -----------------------------------------------------------------------------

tbl_info <- tibble(pool_name = pool, PSMs = nrow(tbl_psms),
    max_length = max(tbl_pool_psms$Length),
    pool_PSMs = nrow(tbl_pool_psms), qc_PSMs = nrow(tbl_qc_psms),
    filtered_pool_psms = nrow(tbl_filtered_pool_psms), 
    filtered_qc_psms = nrow(tbl_filtered_qc_psms), 
    pool_unique_peptides = nrow(tbl_unique_pool_peptides),
    full_length = nrow(tbl_full_length),
    # truncated = nrow(tbl_truncated),
    precursor = nrow(tbl_filtered_precursors),
    frame = nrow(tbl_filtered_pool_frames),
    # wrong_identification = nrow(tbl_wrong),
    fasta_size = nrow(tbl_pool),
    percentage = full_length / fasta_size * 100,
    min_score = min(tbl_mapped_precursor$Score),
    max_score = max(tbl_mapped_precursor$Charge),
    total_annotation = nrow(tbl_annotation)
)

fwrite(tbl_info, txt_output_path, append = TRUE)

# -----------------------------------------------------------------------------