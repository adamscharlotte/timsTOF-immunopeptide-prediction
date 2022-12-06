library(tidyverse)
library(data.table)

txt_path <- paste(base_path, "Metadata/identification-summary.txt", sep = "") # nolint

tbl_summary <- fread(txt_path) %>%
    as_tibble()

colnames(tbl_summary)

tbl_summary %>% filter(pool_wrong == 1056)

sum(tbl_summary$qc_wrong)
min(tbl_summary$qc_wrong)
max(tbl_summary$qc_wrong)

mean(tbl_summary$pool_trunc/tbl_summary$pool_peptides)

tbl_summary %>%
    mutate(wrong_percentage = pool_wrong/total_psm) %>%
    select(pool_name, pool_peptides, total_psm, pool_wrong, wrong_percentage) %>%
    arrange(desc(wrong_percentage))

tbl_summary %>%
    select(pool_name, pool_peptides, total_psm, pool_wrong, pool_cterm, pool_internal, qc_cterm) %>%
    arrange(desc(pool_cterm))
