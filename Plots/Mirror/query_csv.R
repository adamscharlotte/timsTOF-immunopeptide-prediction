library(tidyverse)
library(data.table)

base_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/"

# ---------------------------------- FRAMES -----------------------------------
sa_path <- paste(base_path, "spectral-angle/pairwise/frames/all", sep = "")
psm_path <- paste(base_path, "full-truncated-qc/annotated-20ppm", sep = "")

setwd(sa_path)
files <- dir(pattern = "*.csv")
df_sa <- map_df(files, ~read.csv(.x) %>% mutate(pool_name = basename(.x)))

setwd(psm_path)
files_psm <- paste(substring(files, 1, nchar(files) - 6), ".csv", sep = "")
df_psm <- map_df(files_psm, ~read.csv(.x) %>% mutate(pool_name = basename(.x)))
tbl_psm <- df_psm %>%
    as_tibble() %>%
    mutate(pool_name = substring(pool_name, 1, nchar(pool_name) - 4)) %>%
    select(pool_name, PRECURSOR, PRECURSOR_CHARGE)

tbl_sa <- df_sa %>%
    mutate(pool_name = substring(pool_name, 1, nchar(pool_name) - 6)) %>%
    merge(tbl_psm, by = c("pool_name", "PRECURSOR")) %>%
    as_tibble() %>%
    distinct() %>%
    select(pool_name, A, B, PRECURSOR_CHARGE, PRECURSOR, SPECTRAL_ANGLE)

query_path <- "/Users/adams/Projects/300K/Results/Figures/Mirror/query/"
save_path <- paste(query_path, "frames.csv")
write_csv(tbl_sa, save_path)

# --------------------------------- PRECURSOR ---------------------------------
sa_path <- paste(base_path, "spectral-angle/pairwise/precursors/all", sep = "")
psm_path <- paste(base_path, "precursor-consensus/annotated-20ppm", sep = "")

setwd(sa_path)
df_sa <- map_df(files, ~read.csv(.x) %>% mutate(pool_name = basename(.x)))

setwd(psm_path)
files_psm <- paste(substring(files, 1, nchar(files) - 6), ".csv", sep = "")
df_psm <- map_df(files_psm, ~read.csv(.x) %>% mutate(pool_name = basename(.x)))
tbl_psm <- df_psm %>%
    as_tibble() %>%
    mutate(pool_name = substring(pool_name, 1, nchar(pool_name) - 4)) %>%
    select(pool_name, PRECURSOR, PRECURSOR_CHARGE)

tbl_sa %>%
    filter(pool_name == "TUM_HLA2_88") %>%
    filter(PRECURSOR_CHARGE == 1) %>%
    filter(SEQUENCE == "LDGFSVK")
    count(PRECURSOR_CHARGE)
tbl_sa <- df_sa %>%
    mutate(pool_name = substring(pool_name, 1, nchar(pool_name) - 6)) %>%
    mutate(PRECURSOR = A) %>%
    merge(tbl_psm, by = c("pool_name", "PRECURSOR")) %>%
    as_tibble() %>%
    distinct() %>%
    select(pool_name, A, B, PRECURSOR_CHARGE, PRECURSOR, SEQUENCE,
        SPECTRAL_ANGLE)

query_path <- "/Users/adams/Projects/300K/Results/Figures/Mirror/query/"
save_path <- paste(query_path, "precursors.csv")
write_csv(tbl_sa, save_path)
