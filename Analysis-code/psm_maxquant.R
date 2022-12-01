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

path <- "/Users/adams/Projects/300K/2022-library-run/accumulatedMsmsScans"
setwd(path)
files <- dir(pattern = "*.txt")
scans_txt <- files %>%
    map(read_tsv) %>%        # read in all the files individually, using
                            # the function read_csv() from the readr package
    reduce(rbind)           # reduce with rbind into one dataframe

tbl_scans_txt <- scans_txt %>%
    spaceless() %>%
    rename(folder_names = Raw_file) %>%
    merge(tbl_meta) %>%
    as_tibble

tbl_scans_txt %>%
    select(Score) %>%
    drop_na() %>%
    max()

tbl_scans_txt %>%
    select(Score, plate, pool_name, PASEF_precursor_IDs, Sequence) %>%
    group_by(plate) %>%
    top_n(10, Score) %>%
    ungroup() %>%
    arrange(plate) %>%
    select(-c(Score, plate, pool_name)) %>%
    print(n = 50)

