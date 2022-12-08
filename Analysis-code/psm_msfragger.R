library(tidyverse)
library(data.table)

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

pool <- "TUM_HLA_3"
folder_name <- "HLAI_1_96_p2-A3_S2-A3_1_6928"
base_path <- "/Users/adams/Projects/300K/2022-library-run/"
meta_path <- paste(base_path, "Metadata/full-meta-map.txt", sep = "")
pool_path <- paste(base_path, "Metadata/full-pool-sequence.txt", sep = "")
msfragger_path <-  paste(base_path, "msfragger-results/test/d/HLAI_1_96_p2-A3_S2-A3_1_6928.tsv", sep = "") # nolint
result_path <- "/Users/adams/Projects/300K/2022-library-run/MSFragger-results/Kuster_Lib_HLA1_plate2_MSFraggerResults/psm.tsv" # nolint

# Expected peptide sequences
tbl_folder_pool <- fread(meta_path) %>%
    as_tibble() %>%
    mutate(folder_names = substring(folder_names, 1,
        nchar(folder_names) - 2)) %>%
    select(folder_names, pool_name)

# path_trunc_full <- paste(base_path, "Annotation/full-length-map/", pool, ".csv", sep = "") # nolint

tbl_pool_all <- fread(pool_path) %>%
    as_tibble() %>%
    filter(pool_name == pool) %>%
    dplyr::rename(PROTEINS = pool_name) %>%
    rename_with(toupper) %>%
    select(SEQUENCE, PROTEINS) %>%
    distinct()

# Load results
# tbl_msfragger <- fread(msfragger_path) %>%
#     as_tibble() %>%
#     mutate(PSM = paste(peptide, "_", scannum, "_", charge))
# tbl_msfragger %>% count(protein)

tbl_psm <- fread(result_path) %>%
    as_tibble() %>%
    spaceless()
tbl_psm %>% count(Protein) %>% print(n = 250)

# 331399
# 5477 psm are found in folders that should be empty

tbl_psm_folder <- tbl_psm %>%
    select(Spectrum, Peptide, Protein) %>%
    mutate(folder_names = str_sub(Spectrum, end = -15)) %>%
    mutate(Scan_Number = str_sub(Spectrum, start = -13, end = -9)) %>%
    select(-Spectrum) %>%
    distinct()

tbl_psm_folder %>%
    count(folder_names) %>%
    print(n = 250)

tbl_psm_pool <- tbl_psm_folder %>%
    merge(tbl_folder_pool, all = TRUE) %>%
    as_tibble() %>%
    filter(!is.na(Peptide)) %>%
    distinct()

# 85178 psm
nrow(tbl_psm_pool %>% filter(!Protein == pool_name)) -
    nrow(tbl_psm_pool %>% filter(startsWith(Protein, "QC")))
nrow(tbl_psm_pool) - nrow(tbl_psm_pool %>% filter(Protein == pool_name)) -
    nrow(tbl_psm_pool %>% filter(startsWith(Protein, "QC")))

tbl_psm_pool %>%
    # filter(!Protein == pool_name) %>%
    select(Protein, pool_name) %>%
    distinct() %>%
    count(pool_name) %>%
    arrange(desc(n)) %>%
    print(n = 200)


tbl_psm_pool %>%
    filter(is.na(pool_name)) %>%
    select(Protein, pool_name, Scan_Number, folder_names) %>%
    distinct() %>%
    count(Protein)
    count(folder_names)

tbl_psm_pool %>%
    # filter(!pool_name == Protein) %>%
    filter(pool_name == "TUM_HLA_73") %>%
    distinct() %>%
    select(Protein, pool_name, Scan_Number) %>%
    count(Protein)
    select(Protein, pool_name) %>%
    distinct() %>%
    count(pool_name) %>%
    arrange(desc(n)) %>%
    print(n = 200)

tbl_map <- tbl_meta %>%
    select(pool_name, folder_names)

tbl_merged <- merge(tbl_psm_folder, tbl_map, by = "folder_names") %>%
    as_tibble()

tbl_psm_folder %>%
    filter(!folder_names %in% tbl_meta$folder_names) %>%
    count(Protein)
