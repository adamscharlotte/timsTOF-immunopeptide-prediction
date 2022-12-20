library(tidyverse)
library(data.table)

base_path <- "/Users/adams/Projects/300K/Tryptic/"
msms_path <- paste(base_path, "ProteomeTools_Proteotypic_result/msms.txt", sep = "") # nolint
pool_path <- paste(base_path, "Metadata/Pool_Sequence.csv", sep = "") # nolint

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

tbl_msms <- fread(msms_path) %>%
    as_tibble() %>%
    rename_with(toupper) %>%
    spaceless()

colnames(tbl_msms)

tbl_msms %>%
    count(GENE_NAMES) %>%
    arrange(desc(n))

tbl_msms %>%
    filter(GENE_NAMES == "") %>%
    select(SEQUENCE, GENE_NAMES) %>%
    add_count(SEQUENCE) %>%
    arrange(desc(n)) %>%
    distinct()

tbl_pool <- fread(pool_path) %>% as_tibble() %>% rename_with(toupper)
tbl_map <- tbl_pool %>%
    filter(QC_TYPE == "") %>%
    select(POOL_NAME, SEQUENCE) %>%
    # filter(!GENE_NAMES == "") %>%
    distinct()
tbl_merge <- merge(tbl_msms, tbl_map, by = "SEQUENCE") %>% as_tibble()

tbl_merge %>%
    select(RAW_FILE, POOL_NAME) %>%
    distinct() %>%
    count(POOL_NAME) %>%
    arrange(desc(n)) %>%
    count(n) %>%
    arrange(desc(nn)) %>%
    distinct()

tbl_map %>%
    select(SEQUENCE, POOL_NAME) %>%
    add_count(SEQUENCE) %>%
    arrange(SEQUENCE) %>%
    arrange(desc(n)) %>%
    count(n)

tbl_msms %>% select(RAW_FILE) %>% distinct()

tbl_msms %>%
    select(RAW_FILE) %>%
    distinct() %>%
    filter(str_detect(RAW_FILE, "Plate2"))

tbl_pool %>% select(POOL_NAME) %>% distinct()