library(tidyverse)

base <- "/Users/adams/Projects/300K/2022-library-run/Metadata/"
meta_path <- paste(base, "full-pool-sequence.txt", sep = "")
map_path <- paste(base, "full-meta-map.txt", sep = "")

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

result <- tbl_meta_map %>%
    add_count(plate) %>%
    mutate(sample_size = n / 10) %>%
    group_by(plate) %>%
    sample_n(sample_size)

result

result %>% count(plate)

tbl_meta %>%
    filter(pool_name %in% result$pool_name) %>%
    filter(startsWith(pool_name, "TUM_first_pool")) %>%
    count(pool_name)
