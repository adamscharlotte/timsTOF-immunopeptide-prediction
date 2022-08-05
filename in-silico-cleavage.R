# BiocManager::install("cleaver") # nolint
library("cleaver")
library(tidyverse)
library(data.table)

jpt_qc_peptide <- "PEPTIDEKYVGDSYDSSAKHTAYSDFLSDKYGFSSEDIFTKVDTFLDGFSVKSILDYVSLVEKKPEPTIDE" # nolint
jpt_rt_peptide <- "PEPTIDEKYSAHEEHHYDKHEHISSDYAGKTFAHTESHISKISLGEHEGGGKLSSGYDGTSYKTASGVGGFSTKFGTGTYAGGEKVGASTGYSGLKSYASDFGSSAKTLIAYDDSSTKLYSYYSSTESKLYTGAGYDEVKHLTGLTFDTYKGFLDYESTGAKFLASSEGGFTKYALDSYSLSSKFVGTEYDGLAKHFALFSTDVTKHDTVFGSYLYKALFSSITDSEKYFGYTSDTFGKTFTGTTDSFFKVSGFSDISIYKTFGTETFDTFKFLFTGYDTSVKASDLLSGYYIKGIFGAFTDDYKTSIDSFIDSYKVYAETLSGFIKGFVIDDGLITKGASDFLSFAVKVSSIFFDTFDKFFLTGTSIFVKGDFTFFIDTFKIDVYILALLLKFLISLLEEYFKLFISALVDFFKSILAFLYLYFKSLFFIIDGFVKSLIFFLSTLLKKPEPTIDE" # nolint

jpt_qc_peptide_c <- cleave(jpt_qc_peptide, enzym = "trypsin") %>% unname() %>% unlist() # nolint
jpt_rt_peptide_c <- cleave(jpt_rt_peptide, enzym = "trypsin") %>% unname() %>% unlist() # nolint

list_all_qc <- c(jpt_qc_peptide_c, jpt_rt_peptide_c)

tbl_qc <- tibble(Sequence = list_all_qc) %>%
    mutate(length = str_length(Sequence)) %>%
    arrange(desc(length))

min(tbl_qc$length)
max(tbl_qc$length)

tbl_full_qc<- tbl_qc %>% filter(!str_detect(Sequence, "PEPTIDE"))

qc_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/qc-peptides.txt"
fwrite(tbl_full_qc, qc_path, col.names=T, sep="\t")
