# ssh cadams@10.152.135.57

library(tidyverse)
library(data.table)
library(timsr)

# -----------------------------------------------------------------------------
path_to_bruker_dll <- "/media/c2m/5_user_files/cadams/libtimsdata.so"
setup_bruker_so(path_to_bruker_dll)

accept_Bruker_EULA_and_on_Windows_or_Linux <- TRUE

if (accept_Bruker_EULA_and_on_Windows_or_Linux) {
    path_to_bruker_dll <- "/media/c2m/5_user_files/cadams/libtimsdata.so"
    setup_bruker_so(path_to_bruker_dll)
    all_columns <- c("frame", "scan", "tof", "intensity", "mz",
        "inv_ion_mobility", "retention_time")
} else {
    all_columns <- c("frame", "scan", "tof", "intensity", "retention_time")
}
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
# pool <- "TUM_HLA_16"

base_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/" # nolint
meta_path <- paste(base_path, "Metadata/full-meta-map.txt", sep = "")
# meta_qc_path <- paste(base_path, "Metadata/qc-peptides.txt", sep = "")

meta_pool_path <- paste(base_path, "Metadata/full-pool-sequence.txt", sep = "")
tbl_meta <- fread(meta_path) %>% as_tibble()
plate <- tbl_meta %>% filter(pool_name == pool) %>% pull(plate)
raw_folder <- tbl_meta %>% filter(pool_name == pool) %>% pull(folder_names)

raw_path <- paste(base_path, "raw-library-run/", plate, "/", raw_folder, sep = "") # nolint
txt_path <- paste(base_path, "Searches/", plate, "/", pool, "/combined/txt/", sep = "") # nolint
setwd(txt_path)

tbl_msms <- fread("msms.txt") %>% as_tibble()

tbl_precursors <- fread("accumulatedMsmsScans.txt") %>%
    as_tibble %>%
    filter(`Scan number` %in% tbl_msms$`Scan number`) %>%
    mutate(Precursor = `PASEF precursor IDs`) %>%
    # Split the precursor ID and create a new row for each
    mutate(Precursor = strsplit(as.character(Precursor), ";")) %>%
    unnest(Precursor) %>%
    distinct()

scan_precursor_map <- tbl_precursors %>%
    select(`Scan number`, Precursor) %>%
    distinct()

tbl_pasef <- fread("pasefMsmsScans.txt") %>%
    as_tibble %>%
    filter(Precursor %in% tbl_precursors$Precursor) %>%
    mutate(frame = Frame) %>%
    mutate(collision_energy = CollisionEnergy) %>%
    select(Precursor, frame, ScanNumBegin, ScanNumEnd, collision_energy) %>%
    distinct()

frame_min <- min(tbl_pasef$frame)
frame_max <- max(tbl_pasef$frame)

df_tims <- TimsR(raw_path) # get data handle

df_raw <- df_tims[frame_min : frame_max, all_columns]

tbl_raw <- df_raw %>%
    as_tibble() %>%
    distinct()

tbl_raw_pasef <- merge(tbl_raw, tbl_pasef, by = "frame") %>% as_tibble()

tbl_raw_group <- tbl_raw_pasef %>%
    filter(ScanNumBegin <= scan & scan <= ScanNumEnd) %>%
    distinct() %>%
    group_by(Precursor, frame) %>%
    arrange(mz) %>%
    mutate(intensities = paste(intensity, collapse = ";")) %>%
    mutate(mz = paste(mz, collapse = ";")) %>%
    ungroup() %>%
    select(Precursor, frame, intensities, mz, retention_time,
        collision_energy, inv_ion_mobility) %>%
    distinct()

tbl_msms_precursor <- merge(tbl_msms, scan_precursor_map) %>% as_tibble()

tbl_msms_raw <- tbl_msms_precursor %>%
    select(-c(`Retention time`, Intensities)) %>%
    merge(tbl_raw_group, by = "Precursor") %>%
    as_tibble()

output_path <- paste(base_path, "Annotation/precursor-mapped/", pool, ".csv", sep = "")
fwrite(tbl_msms_raw, output_path)