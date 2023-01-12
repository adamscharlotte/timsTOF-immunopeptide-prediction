library(tidyverse)
library(data.table)

base_path <- "/Users/adams/Projects/300K/2022-library-run/"
filtered_path <- paste(base_path, "Annotation/total-scan-consensus/filtered/", sep = "")  # nolint
meta_path <- paste(base_path, "Metadata/full-meta-map.txt", sep = "")

meta_tbl <- fread(meta_path) %>%
    as_tibble() %>%
    mutate(RAW_FILE = str_remove(folder_names, ".d"))

# ------- get test list --------
meta_tbl %>%
    # filter(startsWith(pool_name, "TUM_first_pool")) %>%
    count(pool_name) %>%
    arrange(desc(n)) %>%
    print(n = 600)

tbl_map <- fread(meta_path) %>%
    as_tibble() %>%
    select(plate, pool_name) %>%
    distinct()

tbl_meta_map <- meta_tbl %>%
    select(pool_name) %>%
    distinct() %>%
    merge(tbl_map, all = TRUE) %>%
    as_tibble() %>%
    mutate(plate = replace_na(plate, "first_pool"))

tbl_meta_map %>%
    count(plate) %>%
    distinct()

set.seed(7)

result <- tbl_meta_map %>%
    add_count(plate) %>%
    mutate(sample_size = n / 10) %>%
    group_by(plate) %>%
    sample_n(sample_size)

# ---------- End of test split ----

setwd(filtered_path)
files <- dir(pattern = "*.csv")
filtered_tbl <- files %>%
    map(read_csv) %>%
    reduce(rbind)

mapped_tbl <- merge(filtered_tbl, meta_tbl,
    by = "RAW_FILE", all = TRUE) %>%
    as_tibble() %>%
    distinct()

counts_tbl <- mapped_tbl %>%
    group_by(plate) %>%
    filter(!is.na(plate)) %>%
    filter(!CHARGE == 4) %>%
    filter(!pool_name %in% result$pool_name) %>%
    count(CHARGE) %>%
    print(n = 20)

median(counts_tbl$n)

mapped_tbl %>%
    filter(plate == "20220621_AspNlysN") %>%
    filter(is.na(CHARGE))

colnames(mapped_tbl)

min(filtered_df$median_CE)
max(filtered_df$median_CE)
colnames(meta_tbl)

mapped_tbl %>%
    filter(!is.na(plate)) %>%
    filter(!CHARGE == 4) %>%
    filter(!pool_name %in% result$pool_name) %>%
    count(CHARGE) %>%
    print(n = 20)