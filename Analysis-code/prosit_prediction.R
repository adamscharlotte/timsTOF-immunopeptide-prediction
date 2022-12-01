library(tidyverse)
library(data.table)

pool <- args[1]
precursor <- args[2]

pool <- "TUM_HLA2_31"

base_path <- "/Users/adams/Projects/300K/2022-library-run/"
meta_path <- paste(base_path, "Metadata/full-meta-map.txt", sep = "")
top_path <- paste(base_path, "Metadata/top-identifications.txt", sep = "")
sum_path <- paste(base_path, "Annotation/precursor-consensus/summed", sep = "")

tbl_meta <- fread(meta_path) %>%
    as_tibble() %>%
    mutate(folder_names = substring(folder_names, 1, nchar(folder_names) - 2))

tbl_top <- fread(top_path) %>%
    as_tibble() %>%
    unite(pool_precursor, V1:V2, sep = "_")

setwd(sum_path)
files <- dir(pattern = "*.csv")
scans_txt <- files %>%
    map(read_csv) %>%        # read in all the files individually, using
                            # the function read_csv() from the readr package
    reduce(rbind)           # reduce with rbind into one dataframe

tbl_sum <- scans_txt %>%
    merge(tbl_meta, by.x = "RAW_FILE", by.y = "folder_names") %>%
    as_tibble() %>%
    select(pool_name, PRECURSOR, everything()) %>%
    unite(pool_precursor, pool_name:PRECURSOR, sep = "_") %>%
    filter(pool_precursor %in% tbl_top$pool_precursor) %>%
    select(MODIFIED_SEQUENCE, COLLISION_ENERGY, PRECURSOR_CHARGE)


prosit_path <- paste(base_path, "Metadata/top-identifications-prosit.csv", sep = "")
write_csv(tbl_sum, prosit_path)

# 6A760C626CB0508BE6F08C49B34094FD
