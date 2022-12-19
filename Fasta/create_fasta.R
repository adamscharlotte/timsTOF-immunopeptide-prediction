library(tidyverse)
library(data.table)
library(seqinr)

args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
# pool <- "TUM_first_pool_12"
# base_path <- "/Users/adams/Projects/300K/2022-library-run/"
base_path <- "/Users/adams/Projects/300K/Tryptic/"
output <- paste(base_path, "fasta/", pool, ".fasta", sep = "")
peptide_path <- paste(base_path, "metadata/full-pool-sequence.txt", sep = "")

tbl_peptides <- fread(peptide_path) %>% as_tibble
list_sequences <- tbl_peptides %>%
    filter(QC_type == "") %>%
    filter(pool_name == pool) %>%
    pull(Sequence)
mock_protein <- paste(list_sequences, collapse = "")
mock_protein_flank <- paste("PEPTIDEK", mock_protein, "KPEPTIDE", sep = "")
JPT_QC_Peptide <- "PEPTIDEKYVGDSYDSSAKHTAYSDFLSDKYGFSSEDIFTKVDTFLDGFSVKSILDYVSLVEKKPEPTIDE" # nolint
JPT_RT_Peptide <- "PEPTIDEKYSAHEEHHYDKHEHISSDYAGKTFAHTESHISKISLGEHEGGGKLSSGYDGTSYKTASGVGGFSTKFGTGTYAGGEKVGASTGYSGLKSYASDFGSSAKTLIAYDDSSTKLYSYYSSTESKLYTGAGYDEVKHLTGLTFDTYKGFLDYESTGAKFLASSEGGFTKYALDSYSLSSKFVGTEYDGLAKHFALFSTDVTKHDTVFGSYLYKALFSSITDSEKYFGYTSDTFGKTFTGTTDSFFKVSGFSDISIYKTFGTETFDTFKFLFTGYDTSVKASDLLSGYYIKGIFGAFTDDYKTSIDSFIDSYKVYAETLSGFIKGFVIDDGLITKGASDFLSFAVKVSSIFFDTFDKFFLTGTSIFVKGDFTFFIDTFKIDVYILALLLKFLISLLEEYFKLFISALVDFFKSILAFLYLYFKSLFFIIDGFVKSLIFFLSTLLKKPEPTIDE" # nolint

list_names <- list(pool, "QC_JPT_QC_Peptide", "QC_JPT_RT_Peptide")
list_mock_protein <- list(mock_protein_flank, JPT_QC_Peptide, JPT_RT_Peptide)
write.fasta(list_mock_protein, list_names, output, open = "w", nbchar = 60, as.string = FALSE) # nolint