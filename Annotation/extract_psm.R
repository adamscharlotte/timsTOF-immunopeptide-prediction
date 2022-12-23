# ssh cadams@10.152.135.57

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
# pool <- "TUM_first_pool_5"

base_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/" # nolint
txt_path <- paste(base_path, "TUM-Searches/first_pool-unsp/", pool, "/combined/txt/", sep = "") # nolint
extract_path <- paste(base_path, "Annotation/extract-d/", pool, ".csv", sep = "") # nolint

tbl_msms <- fread(paste(txt_path, "msms.txt", sep = "")) %>% as_tibble()

tbl_precursors <- fread(paste(txt_path, "accumulatedMsmsScans.txt", sep = "")) %>% # nolint
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

tbl_msms_precursor <- merge(tbl_msms, scan_precursor_map) %>% as_tibble()

tbl_d_extract <- fread(extract_path) %>% as_tibble()

tbl_msms_raw <- tbl_msms_precursor %>%
    select(-c(Intensities)) %>%
    merge(tbl_d_extract, by = "Precursor") %>%
    as_tibble()

output_path <- paste(base_path, "Annotation/extract-psm/", pool, ".csv", sep = "") # nolint
fwrite(tbl_msms_raw, output_path)