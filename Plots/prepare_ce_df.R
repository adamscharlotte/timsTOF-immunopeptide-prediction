# ssh cadams@10.152.135.57

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
# pool <- "TUM_aspn_8"

base_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/Annotation/" # nolint
un_annotated_ce_path <- paste(base_path, "full-truncated-qc/un-annotated-ce/", pool, ".csv", sep = "") # nolint
mapped_precursor_path <- paste(base_path, "precursor-mapped/", pool, ".csv", sep = "") # nolint
# mapped_precursor_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/precursor-mapped/TUM_aspn_1.csv" # nolint
ce_output_path <- paste(base_path, "QC-collission-energies.csv", sep = "")
# ce_output_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/QC-collission-energy.csv" # nolint

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

tbl_mapped_precursor <- fread(mapped_precursor_path) %>%
    as_tibble() %>%
    spaceless() %>%
    select(Scan_number, Proteins) %>%
    filter(startsWith(Proteins, "QC_")) %>%
    distinct()

tbl_un_annotated_ce <- fread(un_annotated_ce_path) %>%
    mutate(Scan_number = SCAN_NUMBER) %>%
    merge(tbl_mapped_precursor) %>%
    as_tibble() %>%
    distinct()

fwrite(tbl_un_annotated_ce, ce_output_path, append = TRUE)
