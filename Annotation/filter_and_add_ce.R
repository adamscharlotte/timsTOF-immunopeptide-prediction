# ssh cadams@10.152.135.57

library(tidyverse)
library(data.table)
library(timsr)

args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
plate <- args[2]
# pool <- "TUM_aspn_1"
# plate <- "20220621_AspNlysN"

base_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/" # nolint
meta_path <- paste(base_path, "Metadata/full-pool-sequence.txt", sep = "")
# meta_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/full-pool-sequence.txt"
mapped_precursor_path <- paste(base_path, "Annotation/precursor-mapped/", pool, ".csv", sep = "")
# mapped_precursor_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/precursor-mapped/TUM_aspn_1.csv"
pasef_path <- paste(base_path, "Searches/", plate, "/", pool, "/combined/txt/pasefMsmsScans.txt", sep = "") # nolint
# pasef_path <- "/Users/adams/Projects/300K/2022-library-run/Searches/txt/pasefMsmsScans.txt"
meta_qc_path <- paste(base_path, "Metadata/qc-peptides.txt", sep = "")

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

tbl_mapped_precursor <- fread(mapped_precursor_path) %>% as_tibble() %>% spaceless()
tbl_pasef <- fread(pasef_path) %>%
    as_tibble() %>%
    select(Precursor, CollisionEnergy) %>%
    distinct()

tbl_qc_peptides <- fread(meta_qc_path) %>% as_tibble()

tbl_pool <- fread(meta_path) %>%
    as_tibble() %>%
    filter(pool_name == pool) %>%
    mutate(trunc_sequence = str_sub(Sequence, start = -7)) %>%
    rename(Sequence_pool = Sequence) %>%
    distinct()

tbl_qc_psms <- tbl_mapped_precursor %>% filter(stringr::str_starts(Proteins, "QC"))

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

tbl_mapped_ce <- tbl_mapped_precursor %>%
    select(-collision_energy) %>%
    merge(tbl_pasef) %>%
    as_tibble() %>%
    distinct()

tbl_annotation <- tbl_mapped_ce %>%
    filter(Sequence %in% tbl_filtered_precursors$Sequence &
        Precursor %in% tbl_filtered_precursors$Precursor) %>%
    filter(Score >= 70) %>%
    rename_with(toupper) %>%
    rename(COLLISION_ENERGY = COLLISIONENERGY) %>%
    select(RAW_FILE, SCAN_NUMBER, MODIFIED_SEQUENCE, CHARGE, FRAGMENTATION,
    MASS_ANALYZER, SCAN_EVENT_NUMBER, MASS, SCORE, REVERSE, RETENTION_TIME,
    MZ, INTENSITIES, COLLISION_ENERGY, PRECURSOR, FRAME)


output_path <- paste(base_path, "Annotation/full-truncated-qc/un-annotated-ce/", pool, ".csv", sep = "")
fwrite(tbl_annotation, output_path)