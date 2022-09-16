library(tidyverse)
library(data.table)
library(stringr)

# pool <- "HLA-II"
# pool <- "HLA-I"
# pool <- "LysN"
pool <- "AspN"

pool_name_output_path <- paste("/Users/adams/Projects/300K/2022-library-run/metadata/poolnames/", pool, "_names.txt", sep="")
meta_output_path <- paste("/Users/adams/Projects/300K/2022-library-run/metadata/", pool, "-meta.txt", sep="")
meta_path <- paste("/Users/adams/Projects/300K/2022-library-run/metadata/", pool, ".txt", sep="")
full_meta_map_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/full-meta-map.txt"
full_pool_sequence_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/full-pool-sequence.txt"

tbl_meta <- fread(meta_path) %>% as_tibble()
tbl_pool_names <- tbl_meta %>% select(`Pool name`) %>% unique()
fwrite(tbl_pool_names, pool_name_output_path, col.names = FALSE)

tbl_meta_length <- tbl_meta %>% mutate(length = str_count(Sequence))
tbl_meta_minmax <- tbl_meta_length %>%
    mutate(pool_name = `Pool name`) %>%
    group_by(`Pool name`) %>%
    mutate(min_length = min(length), max_length = max(length)) %>%
    ungroup()
fwrite(tbl_meta_minmax, meta_output_path, col.names = TRUE, sep = "\t")

tbl_meta_minmax %>% select(min_length, max_length) %>% unique()

# ------------------------------ File mapper ----------------------------------

folder_names_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/Folder-names.txt"
meta_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/"
list_files <- dir(meta_path, pattern = "*.txt")
list_map <- list_files[grep(pattern = "location-map", list_files)]
list_meta <- list_files[grep(pattern = "meta", list_files)]
list_pool_sequence <- list_files[-grep(pattern = "meta|location-map|-names", list_files)]

files_dir <- paste(meta_path, list_map, sep = "")
tbl_map <- files_dir %>%
    map(read_tsv) %>%       # read in all the files individually, using
                            # the function read_tsv() from the readr package
    purrr::reduce(rbind)           # reduce with rbind into one dataframe

files_dir <- paste(meta_path, list_meta, sep = "")
tbl_meta <- files_dir %>%
    map(read_tsv) %>%       # read in all the files individually, using
                            # the function read_tsv() from the readr package
    purrr::reduce(rbind)           # reduce with rbind into one dataframe

tbl_map %>% select(pool_name) %>% unique()
tbl_meta %>% select(pool_name) %>% unique()

tbl_folder_names <- fread(folder_names_path) %>% as_tibble
tbl_folder_names %>% select(plate) %>% unique()

tbl_folder_map <- merge(tbl_map, tbl_folder_names) %>%
    as_tibble() %>%
    filter(str_detect(folder_names, paste(location, "_", sep = "")))

files_dir <- paste(meta_path, list_pool_sequence, sep = "")
tbl_pool_sequence <- files_dir %>%
    map(read_tsv) %>%       # read in all the files individually, using
                            # the function read_tsv() from the readr package
    purrr::reduce(rbind)           # reduce with rbind into one dataframe

fwrite(tbl_pool_sequence, full_pool_sequence_path, col.names = TRUE, sep = "\t")

# -----------------------------------------------------------------------------
# Some folders are missing
# I assume these were empty wells that were analyzed anyway. We could test
# to see whether we identify anything, but I think we won't.
tbl_folder_names %>%
    filter(!folder_names %in% tbl_folder_map$folder_names) %>%
    print(n = 60)

# -----------------------------------------------------------------------------

tbl_folder_meta_map <- merge(tbl_folder_map, tbl_meta, by = "pool_name") %>%
    as_tibble() %>%
    select(-c(length, Sequence)) %>%
    distinct()

fwrite(tbl_folder_meta_map, full_meta_map_path, col.names = TRUE, sep = "\t")