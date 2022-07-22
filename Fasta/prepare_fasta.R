library(tidyverse)
library(data.table)

base_path <- "/Users/adams/Projects/300K/2022-library-run/"
output_path <- paste(base_path, "metadata/poolnames/All.txt", sep = "")
peptide_path <- paste(base_path, "metadata/full-pool-sequence.txt", sep = "")
tbl_peptides <- fread(peptide_path) %>% as_tibble

tbl_poolnames <- tbl_peptides %>% select(pool_name) %>% unique()

fwrite(tbl_poolnames, output_path, col.names = FALSE)
