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
output_path <- args[2]
pasef_path <- args[3]

tbl_pasef <- read.csv(pasef_path) %>%
    as_tibble() %>%
    mutate(frame = Frame)

df_tims <- TimsR(raw_path) # get data handle
df_raw <- df_tims[tbl_pasef$frame, all_columns]

tbl_raw <- df_raw %>%
    as_tibble() %>%
    distinct()

tbl_raw_pasef <- merge(tbl_raw, tbl_pasef, by = "frame") %>% as_tibble()

tbl_raw_group <- tbl_raw_pasef %>%
    filter(ScanNumBegin <= scan & scan <= ScanNumEnd) %>%
    distinct() %>%
    group_by(Precursor) %>%
    arrange(mz) %>%
    mutate(intensities = paste(intensity, collapse = ";")) %>%
    mutate(mz = paste(mz, collapse = ";")) %>%
    ungroup() %>%
    select(Precursor, frame, intensities, mz, retention_time,
        CollisionEnergy, inv_ion_mobility) %>%
    distinct()

fwrite(tbl_raw_group, output_path)