library(tidyverse)
library(data.table)

# setwd("/Users/adams/Projects/300K/2022-library-run/MaxQuant/AspN/aspn_1_2-txt/")
# setwd("/Users/adams/Projects/300K/2022-library-run/MaxQuant/HLA-I/HLA_16-txt/")
setwd("/Users/adams/Projects/300K/2022-library-run/MaxQuant/LysN/LysN_5-txt/")
# setwd("/Users/adams/Projects/300K/2022-library-run/MaxQuant/LysN/LysN_5_3-txt/") -> fractions 
# setwd("/Users/adams/Projects/300K/2022-library-run/MaxQuant/LysN/LysN_5_no_min_max-txt/") -> setting the min and max peptide length is important
QC_path <- paste("/Users/adams/Projects/300K/2022-library-run/metadata/QC.txt", sep="")

list.files()

tbl_pasef <- fread("pasefMsmsScans.txt") %>% as_tibble
tbl_pasef %>% select(-c(`Raw file`, IsolationWidth))
tbl_pasef %>% group_by(Precursor) %>% mutate(difference=IsolationMz - lag(IsolationMz)) %>% ungroup() %>% select(difference) %>% unique
tbl_single_scan <- tbl_pasef %>% count(Precursor) %>% filter(n==1)

tbl_msms <- fread("msms.txt") %>% as_tibble
tbl_msms %>% select(Proteins) %>% unique()
tbl_msms %>% select(Sequence) %>% unique()
tbl_msms %>% filter(`Scan number`==5255)
colnames(tbl_msms)
tbl_msms %>% select(Sequence, `Missed cleavages`) %>% unique() %>% count(`Missed cleavages`) 


tbl_scans <- fread("accumulatedMsmsScans.txt") %>% as_tibble
tbl_scans %>% filter(`PASEF precursor IDs`==tbl_single_scan$Precursor)

meta_path <- paste("/Users/adams/Projects/300K/2022-library-run/metadata/LysN-meta.txt", sep="")
tbl_meta <- fread(meta_path) %>% as_tibble %>% select(pool_name, Sequence)
tbl_QC <- fread(QC_path) %>% as_tibble

tbl_pool <- tbl_meta %>% filter(pool_name=="TUM_lysn_5") %>% unique()
tbl_pool %>% mutate(trunc_sequence = substr(Sequence, 1, 7)) %>% select(trunc_sequence) %>% unique()

tbl_truncated <- tbl_msms %>% select(Sequence) %>% unique() %>% filter(!Sequence %in% tbl_pool$Sequence)
list_truncated <- tbl_truncated  %>% pull(Sequence) %>% paste(collapse="|")

tbl_pool %>% filter(stringr::str_starts(Sequence, tbl_truncated$Sequence))

tbl_found_truncated <- (tbl_pool, grepl(list_truncated, Sequence))
tbl_found_truncated <- (tbl_truncated, grepl(list_truncated, Sequence))

match <- str_match_all(tbl_pool$Sequence, list_truncated) 
head(match, n=1000)
tbl_truncated <- tbl_pool %>% select(Sequence) %>% unique() %>% filter(str_detect(Sequence, tbl_msms$Sequence))

tbl_found_truncated
tbl_not_found <- tbl_truncated %>%

tbl_pool %>% 

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(timsr)
path = "/Users/adams/Projects/300K/2022-library-run/MaxQuant/LysN/AspNLysN-A5_S1-A5_1_6546.d"
all_columns = c('frame','scan','tof','intensity','retention_time')
D = TimsR(path) # get data handle
print(D)
D[2347]

X = D[c(1,5,67), all_columns]
X = D[13:354, 355]


print(X)