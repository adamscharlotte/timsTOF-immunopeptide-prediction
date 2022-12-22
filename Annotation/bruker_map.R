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
frame_min <- args[3]
frame_max <- args[4]

df_tims <- TimsR(raw_path) # get data handle
df_raw <- df_tims[frame_min : frame_max, all_columns]

tbl_raw <- df_raw %>%
    as_tibble() %>%
    distinct()

fwrite(tbl_raw, output_path)