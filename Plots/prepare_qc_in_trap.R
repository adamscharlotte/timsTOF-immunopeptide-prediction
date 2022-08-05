# ssh cadams@10.152.135.57

library(tidyverse)
library(data.table)
library(timsr)

args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
# pool <- "TUM_aspn_1"

base_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/" # nolint
meta_path <- paste(base_path, "Metadata/full-pool-sequence.txt", sep = "")
mapped_precursor_path <- paste(base_path, "Annotation/precursor-mapped/", pool, ".csv", sep = "")
# mapped_precursor_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/precursor-mapped/TUM_aspn_1.csv"
rt_output_path <- paste(base_path, "Annotation/QC-retention-time.csv", sep = "")
# ce_output_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/QC-collission-energy.csv"

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

tbl_mapped_precursor <- fread(mapped_precursor_path) %>% as_tibble() %>% spaceless()
tbl_qc_rt <- tbl_mapped_precursor %>%
    filter(stringr::str_starts(Proteins, "QC")) %>%
    mutate(pool_name = pool) %>%
    select(pool_name, Scan_number, Proteins, Sequence, Charge, Modifications, Score, retention_time, Precursor, frame) %>%
    distinct()

fwrite(tbl_qc_rt, rt_output_path, append = TRUE)
