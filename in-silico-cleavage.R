# BiocManager::install("cleaver") # nolint
library("cleaver")
library(tidyverse)
library(data.table)

jpt_qc_peptide <- "KYVGDSYDSSAKHTAYSDFLSDKYGFSSEDIFTKVDTFLDGFSVKSILDYVSLVEKK" # nolint
jpt_rt_peptide <- "KYSAHEEHHYDKHEHISSDYAGKTFAHTESHISKISLGEHEGGGKLSSGYDGTSYKTASGVGGFSTKFGTGTYAGGEKVGASTGYSGLKSYASDFGSSAKTLIAYDDSSTKLYSYYSSTESKLYTGAGYDEVKHLTGLTFDTYKGFLDYESTGAKFLASSEGGFTKYALDSYSLSSKFVGTEYDGLAKHFALFSTDVTKHDTVFGSYLYKALFSSITDSEKYFGYTSDTFGKTFTGTTDSFFKVSGFSDISIYKTFGTETFDTFKFLFTGYDTSVKASDLLSGYYIKGIFGAFTDDYKTSIDSFIDSYKVYAETLSGFIKGFVIDDGLITKGASDFLSFAVKVSSIFFDTFDKFFLTGTSIFVKGDFTFFIDTFKIDVYILALLLKFLISLLEEYFKLFISALVDFFKSILAFLYLYFKSLFFIIDGFVKSLIFFLSTLLKK" # nolint

jpt_qc_peptide_c <- cleave(jpt_qc_peptide, enzym = "trypsin") %>% unname() %>% unlist() # nolint
jpt_rt_peptide_c <- cleave(jpt_rt_peptide, enzym = "trypsin") %>% unname() %>% unlist() # nolint

tbl_qc_qc <- tibble(Sequence = jpt_qc_peptide_c) %>%
    mutate(length = str_length(Sequence)) %>%
    mutate(Proteins = "QC_JPT_QC_Peptide")

tbl_qc_rt <- tibble(Sequence = jpt_rt_peptide_c) %>%
    mutate(length = str_length(Sequence)) %>%
    mutate(Proteins = "QC_JPT_RT_Peptide")

tbl_full_qc <- rbind(tbl_qc_qc, tbl_qc_rt) %>%
    filter(length > 10) %>%
    mutate(trunc_sequence = Sequence) %>%
    select(Protein, Sequence, trunc_sequence)

qc_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/qc-peptides.txt" # nolint
fwrite(tbl_full_qc, qc_path, col.names = TRUE, sep = "\t")
