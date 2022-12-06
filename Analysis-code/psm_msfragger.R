library(tidyverse)
library(data.table)

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

meta_path <- paste(base_path, "Metadata/full-meta-map.txt", sep = "")
tbl_meta <- fread(meta_path) %>%
    as_tibble() %>%
    mutate(folder_names = substring(folder_names, 1, nchar(folder_names)-2))

result_path <- "/Users/adams/Projects/300K/2022-library-run/MSFragger-results/Kuster_Lib_HLA1_plate2_MSFraggerResults/psm.tsv" # nolint
tbl_psm <- fread(result_path) %>%
    as_tibble() %>%
    spaceless()

colnames(tbl_psm)

# 331399
# 5477 psm are found in folders that should be empty



tbl_psm_folder <- tbl_psm %>%
    select(Spectrum, Peptide, Protein) %>%
    mutate(folder_names = str_sub(Spectrum, end = -15)) %>%
    mutate(Scan_Number = str_sub(Spectrum, start = -13, end = -9)) %>%
    select(-Spectrum) %>%
    distinct()

tbl_map <- tbl_meta %>%
    select(pool_name, folder_names)

tbl_merged <- merge(tbl_psm_folder, tbl_map, by = "folder_names") %>%
    as_tibble()

tbl_psm_folder %>%
    filter(!folder_names %in% tbl_meta$folder_names) %>%
    count(Protein)
