# ssh cadams@10.152.135.57

library(tidyverse)
library(data.table)
library(timsr)

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

args <- commandArgs(trailingOnly = TRUE)
d_path <- args[1]
txt_path <- args[2] # Should contain accumulatedMsmsScans.txt, msms.txt and pasefMsmsScans.txt
path_to_bruker_dll <- args[3]

# d_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/d-folder/211113_SS_malignant_HNSCC_Tue39L243_20%_DDA_Rep1.d" # nolint
# txt_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/reresults/test-d" # nolint
# path_to_bruker_dll <- "/media/c2m/5_user_files/cadams/libtimsdata.so"

# -----------------------------------------------------------------------------
setup_bruker_so(path_to_bruker_dll)

accept_Bruker_EULA_and_on_Windows_or_Linux <- TRUE

if (accept_Bruker_EULA_and_on_Windows_or_Linux) {
    path_to_bruker_dll <- path_to_bruker_dll
    setup_bruker_so(path_to_bruker_dll)
    all_columns <- c("frame", "scan", "tof", "intensity", "mz",
        "inv_ion_mobility", "retention_time")
} else {
    all_columns <- c("frame", "scan", "tof", "intensity", "retention_time")
}
# -----------------------------------------------------------------------------

tbl_msms <- fread(paste(txt_path, "/msms.txt", sep = "")) %>%
    as_tibble() %>%
    rename_with(toupper) %>%
    spaceless()

tbl_precursors <- fread(paste(txt_path,
        "/accumulatedMsmsScans.txt", sep = "")) %>%
    as_tibble() %>%
    rename_with(toupper) %>%
    spaceless() %>%
    filter(SCAN_NUMBER %in% tbl_msms$SCAN_NUMBER) %>%
    mutate(PRECURSOR = PASEF_PRECURSOR_IDS) %>%
    # Split the precursor ID and create a new row for each
    mutate(PRECURSOR = strsplit(as.character(PRECURSOR), ";")) %>%
    unnest(PRECURSOR) %>%
    distinct()

scan_precursor_map <- tbl_precursors %>%
    select(SCAN_NUMBER, PRECURSOR) %>%
    distinct()

tbl_pasef <- fread(paste(txt_path, "/pasefMsmsScans.txt", sep = "")) %>%
    as_tibble() %>%
    rename_with(toupper) %>%
    spaceless() %>%
    filter(PRECURSOR %in% tbl_precursors$PRECURSOR) %>%
    mutate(COLLISION_ENERGY = COLLISIONENERGY) %>%
    select(PRECURSOR, FRAME, SCANNUMBEGIN, SCANNUMEND, COLLISION_ENERGY) %>%
    distinct()

df_tims <- TimsR(d_path) # get data handle
df_raw <- df_tims[tbl_pasef$FRAME, all_columns]

tbl_raw <- df_raw %>%
    as_tibble() %>%
    rename_with(toupper) %>%
    spaceless() %>%
    distinct()

dt_raw <- data.table(tbl_raw)
dt_pasef <- data.table(tbl_pasef)

tbl_raw_pasef <- merge(dt_raw, dt_pasef, by = "FRAME", all = TRUE,
    allow.cartesian = TRUE) %>%
    as_tibble()

tbl_raw_group <- tbl_raw_pasef %>%
    filter(SCANNUMBEGIN <= SCAN & SCAN <= SCANNUMEND) %>%
    distinct() %>%
    group_by(PRECURSOR, FRAME) %>%
    arrange(MZ) %>%
    mutate(INTENSITIES = paste(INTENSITY, collapse = ";")) %>%
    mutate(MZ = paste(MZ, collapse = ";")) %>%
    ungroup() %>%
    select(PRECURSOR, INTENSITIES, MZ, RETENTION_TIME,
        COLLISION_ENERGY, INV_ION_MOBILITY) %>%
    distinct()

dt_scan_precursor_map <- data.table(
    scan_precursor_map %>% mutate(PRECURSOR = as.integer(PRECURSOR)))
dt_raw_group <- data.table(tbl_raw_group)

tbl_raw_scan <- merge(dt_raw_group, dt_scan_precursor_map) %>% as_tibble()

tbl_raw_scan_group <- tbl_raw_scan %>%
    group_by(SCAN_NUMBER) %>%
    mutate(median_CE = median(COLLISION_ENERGY)) %>%
    mutate(combined_INTENSITIES = paste0(INTENSITIES, collapse = ";")) %>%
    mutate(combined_MZ = paste0(MZ, collapse = ";")) %>%
    mutate(median_RETENTION_TIME = median(RETENTION_TIME)) %>%
    mutate(median_INV_ION_MOBILITY = median(INV_ION_MOBILITY)) %>%
    ungroup() %>%
    select(-c(COLLISION_ENERGY, INTENSITIES, MZ, RETENTION_TIME,
        INV_ION_MOBILITY, FRAME, PRECURSOR)) %>%
    distinct()

dt_msms <- data.table(tbl_msms)
dt_raw_scan_group <- data.table(tbl_raw_scan_group)

tbl_msms_raw <- merge(dt_raw_scan_group, dt_msms, by = "SCAN_NUMBER",
    all = TRUE, allow.cartesian = TRUE) %>%
    as_tibble() %>%
    select(RAW_FILE, SCAN_NUMBER, combined_INTENSITIES,
        combined_MZ, MASS_ANALYZER, FRAGMENTATION,
        median_RETENTION_TIME, median_INV_ION_MOBILITY,
        median_CE, CHARGE)

file_name <- tbl_msms_raw %>% select(RAW_FILE) %>% distinct() %>% pull(RAW_FILE)
dt_msms_raw <- data.table(tbl_msms_raw)

output_path <- paste(txt_path, "/", file_name, ".csv", sep = "")
fwrite(dt_msms_raw, output_path)