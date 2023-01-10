library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
# pool <- "TUM_first_pool_5"
# pool <- "TUM_aspn_11"

# base_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/" # nolint
base_path <- "/Users/adams/Projects/300K/2022-library-run/"
meta_path <- paste(base_path, "Metadata/full-pool-sequence.txt", sep = "")
meta_qc_path <- paste(base_path, "Metadata/qc-peptides.txt", sep = "")
# mapped_precursor_path <- paste(base_path, "Annotation/extract-psm/", pool, ".csv", sep = "")  # nolint
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
    distinct()

tbl_qc_peptides <- fread(meta_qc_path) %>%
    as_tibble() %>%
    select(Proteins, Sequence)

tbl_peptides <- rbind(tbl_pool_peptides, tbl_qc_peptides)

# tbl_mapped_precursor %>%
#     select(frame, Scan_number, Precursor, Proteins, obs_sequence) %>%
#     distinct()

tbl_psms <- tbl_mapped_precursor %>%
    select(frame, Scan_number, Precursor, Proteins, obs_sequence) %>%
    distinct() %>%
    merge(tbl_peptides, by = "Proteins") %>%
    as_tibble()

tbl_filtered_psms <- tbl_psms %>%
    filter(endsWith(obs_sequence, Sequence)) %>%
    filter(str_length(obs_sequence) <= str_length(Sequence)) %>%
    group_by(obs_sequence) %>%
    filter(str_length(Sequence) == min(str_length(Sequence))) %>%
    ungroup()

# tbl_filtered_psms %>%
#     distinct() %>%
#     select() %>%
#     distinct()

# tbl_filtered_psms %>%
#     select(Scan_number, Precursor, obs_sequence, frame, Protein_group_IDs) %>%
#     distinct()

# tbl_mapped_precursor %>% summarise_all(n_distinct) %>% as.data.frame()

if (!nrow(tbl_filtered_psms %>%
    group_by(Scan_number, Precursor, frame) %>%
    filter(n() == 1) %>%
    ungroup()
    ) ==
    nrow(tbl_filtered_psms)) warning(
    "Some scans are present multiple times. This can cause duplicates.")

tbl_annotation <- tbl_mapped_precursor %>%
    select(Proteins, Scan_number, frame, Precursor, Raw_file,
        Modified_sequence, Charge, Fragmentation, Mass_analyzer,
        Scan_event_number, Mass, Score, Reverse, retention_time,
        collision_energy, mz, intensities, obs_sequence) %>%
    merge(tbl_filtered_psms,
        by = c("Proteins", "Scan_number", "frame", "Precursor",
        "obs_sequence")) %>%
    as_tibble() %>%
    filter(Score >= 70) %>%
    rename_with(toupper) %>%
    # Group by scan numbers instead of precursor
    group_by(SCAN_NUMBER) %>%
    mutate(median_CE = median(COLLISION_ENERGY)) %>%
    mutate(combined_INTENSITIES = paste0(INTENSITIES, collapse = ";")) %>%
    mutate(combined_MZ = paste0(MZ, collapse = ";")) %>%
    mutate(median_RETENTION_TIME = median(RETENTION_TIME)) %>%
    ungroup() %>%
    select(RAW_FILE, SCAN_NUMBER, MODIFIED_SEQUENCE, CHARGE, FRAGMENTATION,
        MASS_ANALYZER, SCAN_EVENT_NUMBER, MASS, SCORE, REVERSE,
        median_RETENTION_TIME, combined_MZ, combined_INTENSITIES, median_CE,
        OBS_SEQUENCE, SEQUENCE) %>%
    distinct()

output_path <- paste(base_path, "Annotation/total-scan-consensus/filtered/", pool, ".csv", sep = "") # nolint
fwrite(tbl_annotation, output_path)