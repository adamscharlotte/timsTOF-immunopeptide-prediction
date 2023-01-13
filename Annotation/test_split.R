library(tidyverse)
library(data.table)

base <- "/Users/adams/Projects/300K/2022-library-run/"
msms_path <- paste(base, "msms-txtx/", sep = "")
meta_path <- paste(base, "Metadata/full-pool-sequence.txt", sep = "")
map_path <- paste(base, "Metadata/full-meta-map.txt", sep = "")

tbl_meta <- fread(meta_path) %>% as_tibble()
tbl_meta %>%
    # filter(startsWith(pool_name, "TUM_first_pool")) %>%
    count(pool_name) %>%
    arrange(desc(n)) %>%
    print(n = 600)

tbl_map <- fread(map_path) %>%
    as_tibble() %>%
    select(plate, pool_name) %>%
    distinct()

tbl_meta_map <- tbl_meta %>%
    select(pool_name) %>%
    distinct() %>%
    merge(tbl_map, all = TRUE) %>%
    as_tibble() %>%
    mutate(plate = replace_na(plate, "first_pool"))

tbl_meta_map %>%
    count(plate) %>%
    distinct()

set.seed(7)

test_set <- tbl_meta_map %>%
    add_count(plate) %>%
    mutate(sample_size = n / 10) %>%
    group_by(plate) %>%
    sample_n(sample_size)

test_set %>% count(plate)

validation_set <- tbl_meta_map %>%
    add_count(plate) %>%
    mutate(sample_size = n / 10) %>%
    filter(!pool_name %in% test_set$pool_name) %>%
    group_by(plate) %>%
    sample_n(sample_size)

validation_set %>% count(plate)

validation_set %>%
    filter(pool_name %in% test_set$pool_name)

# -------------------- WRITE THE POOL NAMES TO A TXT FILE --------------------


# -------------------- WRITE THE POOL NAMES TO A TXT FILE --------------------
test_set_t <- test_set$pool_name

write.table(test_set_t,
    file = "/Users/adams/Code/timsTOF-immunopeptide-prediction/Names/test-set.txt", # nolint
    col.names = FALSE,
    row.names = FALSE,
    sep = "")

validation_set_t <- validation_set$pool_name

write.table(validation_set_t,
    file = "/Users/adams/Code/timsTOF-immunopeptide-prediction/Names/validation-set.txt", # nolint
    col.names = FALSE,
    row.names = FALSE,
    sep = "")