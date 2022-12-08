# ssh cadams@10.152.135.57

library(tidyverse)
library(data.table)
library(timsr)

# -----------------------------------------------------------------------------
path_to_bruker_dll <- "/media/kusterlab/users_files/charlotte_adams/libtimsdata.so"
setup_bruker_so(path_to_bruker_dll)

accept_Bruker_EULA_and_on_Windows_or_Linux <- TRUE

if (accept_Bruker_EULA_and_on_Windows_or_Linux) {
    path_to_bruker_dll <- "/media/kusterlab/users_files/charlotte_adams/libtimsdata.so"
    setup_bruker_so(path_to_bruker_dll)
    all_columns <- c("frame", "scan", "tof", "intensity", "mz",
        "inv_ion_mobility", "retention_time")
} else {
    all_columns <- c("frame", "scan", "tof", "intensity", "retention_time")
}
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
raw_path <- args[1]
txt_path <- args[2]
output_path <- args[3]

# base_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/" # nolint
base_path <- "/Users/adams/Projects/300K/" # nolint
# raw_path <- paste(base_path, "PXD030334/190926_TIMSiDE_LCMS02_sample-3_90min_R2_Slot1-30_01_3622.d", sep = "") # nolint
# txt_path <- paste(base_path, "PXD030334/combined/txt/", sep = "")
# output_path <- paste(base_path, "PXD030334/precursor-mapped/S3_R2.csv", sep = "")

raw_path <- paste(base_path, "2022-library-run/raw-folders/HLAI_1_96_p2-A3_S2-A3_1_6928.d", sep = "") # nolint
txt_path <- paste(base_path, "2022-library-run/MaxQuant/TUM_HLA_3-txt/", sep = "") # nolint
all_columns <- c('frame','scan','tof','intensity')
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
frame_min
frame_max
df_tims <- TimsR(raw_path) # get data handle
df_raw <- df_tims[all_columns]
print(df_tims)
X = df_tims[c(1,5,67)]

tbl_raw <- df_raw %>%
    as_tibble() %>%
    distinct()
tbl_raw
tbl_raw_pasef <- merge(tbl_raw, tbl_pasef, by = "frame") %>% as_tibble()
tbl_raw_pasef
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

fwrite(tbl_msms_raw, output_path)