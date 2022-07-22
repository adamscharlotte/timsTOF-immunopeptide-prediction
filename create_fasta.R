library(tidyverse)
library(data.table)
library(seqinr)

args <- commandArgs(trailingOnly = TRUE)
pool <- args[1]
type <- args[2]
output <- paste("/Users/adams/Projects/300K/2022-library-run/fasta/", pool, ".fasta", sep = "")
peptide_path <- paste("/Users/adams/Projects/300K/2022-library-run/metadata/", type, ".txt", sep = "")
# pool <- "CFP_aspn_1"

mass_tolerance <- as.numeric(args[3])

tbl_peptides <- fread(peptide_path) %>% as_tibble
tbl_RT_QC <- fread("/Users/adams/Projects/300K/2022-library-run/fasta/QC/RT_QC.fasta") %>% as_tibble

list_sequences <- tbl_peptides %>% filter(`Pool name` == pool) %>% pull(Sequence) 
mock_protein <- paste(list_sequences, collapse = '')
mock_protein_flank <- paste("PEPTIDEK", mock_protein, "KPEPTIDE", sep = "")
JPT_QC_Peptide <- "PEPTIDEKYVGDSYDSSAKHTAYSDFLSDKYGFSSEDIFTKVDTFLDGFSVKSILDYVSLVEKKPEPTIDE"
JPT_RT_Peptide <- "PEPTIDEKYSAHEEHHYDKHEHISSDYAGKTFAHTESHISKISLGEHEGGGKLSSGYDGTSYKTASGVGGFSTKFGTGTYAGGEKVGASTGYSGLKSYASDFGSSAKTLIAYDDSSTKLYSYYSSTESKLYTGAGYDEVKHLTGLTFDTYKGFLDYESTGAKFLASSEGGFTKYALDSYSLSSKFVGTEYDGLAKHFALFSTDVTKHDTVFGSYLYKALFSSITDSEKYFGYTSDTFGKTFTGTTDSFFKVSGFSDISIYKTFGTETFDTFKFLFTGYDTSVKASDLLSGYYIKGIFGAFTDDYKTSIDSFIDSYKVYAETLSGFIKGFVIDDGLITKGASDFLSFAVKVSSIFFDTFDKFFLTGTSIFVKGDFTFFIDTFKIDVYILALLLKFLISLLEEYFKLFISALVDFFKSILAFLYLYFKSLFFIIDGFVKSLIFFLSTLLKKPEPTIDE"
PRTC_RT_Peptide <- "PEPTIDEKSSAAPPPPPRGISNEGQNASIKHVLTSIGEKDIPVPKPKIGDYAGIKTASEFDSAIAQDKSAAGAFGPELSRELGQSGVDTYLQTKGLILVGGYGTRGILFVGSGVSGGEEGARSFANQPLEVVYSKLTILEELRNGFILDGFPRELASGLSFPVGFKLSSEAPALFQFDLKKPEPTIDE"

list_names <- list(pool, "QC: JPT_QC_Peptide", "QC: JPT_RT_Peptide", "QC: PRTC_RT_Peptide")
list_mock_protein <- list(mock_protein_flank, JPT_QC_Peptide, JPT_RT_Peptide, PRTC_RT_Peptide)
write.fasta(list_mock_protein, list_names, output, open = "w", nbchar = 60, as.string = FALSE)