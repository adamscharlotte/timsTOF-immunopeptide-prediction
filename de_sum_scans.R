library(tidyverse)
library(data.table)

base_path <- "/Users/adams/Projects/300K/2022-library-run/"
pool <- "TUM_HLA_16"    # "TUM_aspn_1_2", "TUM_HLA_16", "TUM_lysn_5"
plate <- "20220624_HLAI_1_96"    # 20220624_HLAI_1_96, 20220621_AspNlysN

setwd(paste(base_path, "MaxQuant/", plate, "/", pool, "-txt/", sep = ""))
list.files()

tbl_info <- tibble(pool_name = pool, pool_PSMs = nrow(tbl_msms),
    fasta_size = nrow(tbl_pool))

tbl_msms <- fread("msms.txt") %>% as_tibble
tbl_scans <- fread("accumulatedMsmsScans.txt") %>%
    as_tibble %>%
    filter(`Scan number` %in% tbl_msms$`Scan number`) %>%
    mutate(Precursor = `PASEF precursor IDs`) %>%
    # Split the precursor ID and create a new row for each
    mutate(Precursor = strsplit(as.character(Precursor), ";")) %>%
    unnest(Precursor) %>%
    mutate(pool_name = pool) %>%
    unique()

# The following resulted in more PSMs than in msms.txt
# there seems to be another filter that is applied
# filter(!Sequence == "" & is.na(Reverse)) %>%

meta_path <- paste(base_path, "metadata/full-pool-sequence.txt", sep = "")
tbl_meta <- fread(meta_path) %>% as_tibble
tbl_pool <- tbl_meta %>%
    filter(pool_name == pool) %>%
    mutate(trunc_sequence = str_sub(Sequence, start = -7)) %>%
    rename(Sequence_pool = Sequence) %>%
    unique()

tbl_filtered_scans <- merge(tbl_scans, tbl_pool, by = "pool_name") %>%
    as_tibble() %>%
    filter(stringr::str_ends(Sequence, trunc_sequence))

tbl_pasef <- fread("pasefMsmsScans.txt") %>% as_tibble
tbl_scans_pasef <- merge(tbl_filtered_scans, tbl_pasef, by = "Precursor") %>%
    as_tibble()

tbl_scans_pasef %>% select(Sequence) %>% unique()

# Number of full-length sequences
tbl_filtered_scans %>%
    filter(Sequence == Sequence_pool) %>%
    select(Sequence) %>%
    unique()

all_sequences <- tbl_filtered_scans %>%
    pull(Sequence) %>%
    unique()

tbl_scans %>%
    filter(!Sequence %in% all_sequences) %>%
    select(Sequence) %>%
    unique()

# ----------------------------------------------------------------------------
tbl_pool %>% filter(stringr::str_starts(Sequence, tbl_truncated$Sequence))

# The frame and scan beginning and end can be found in tbl_scans_pasef
# and can be mapped to get the raw scans.

qc_path <- paste(base_path, "metadata/QC.txt", sep = "")
tbl_qc <- fread(qc_path) %>% as_tibble

tbl_scans_pasef %>% select(Precursor) %>% unique()
tbl_scans_pasef %>% pull(`PASEF precursor IDs`)
tbl_pasef %>% filter(!Precursor %in% tbl_scans$`PASEF precursor IDs`)
tbl_msms_scans %>% filter(Sequence.x == Sequence.y)

tbl_msms_scans %>%
    select()
tbl_single_scan <- tbl_pasef %>% count(Precursor) %>% filter(n == 1)

tbl_msms %>% select(Proteins) %>% unique()
tbl_msms %>% select(Sequence) %>% unique()
tbl_msms %>% filter(`Scan number` == 5255)
colnames(tbl_msms)
tbl_msms %>%
    select(Sequence, `Missed cleavages`) %>%
    unique() %>%
    count(`Missed cleavages`)

tbl_scans %>% filter(`PASEF precursor IDs` == tbl_single_scan$Precursor)



tbl_truncated <- tbl_msms %>%
    select(Sequence) %>%
    unique() %>%
    filter(!Sequence %in% tbl_pool$Sequence)

list_truncated <- tbl_truncated  %>% pull(Sequence) %>% paste(collapse = "|")


# tbl_found_truncated <- (tbl_pool, grepl(list_truncated, Sequence))
# tbl_found_truncated <- (tbl_truncated, grepl(list_truncated, Sequence))

match <- str_match_all(tbl_pool$Sequence, list_truncated)
head(match, n = 1000)
tbl_truncated <- tbl_pool %>%
    select(Sequence) %>%
    unique() %>%
    filter(str_detect(Sequence, tbl_msms$Sequence))

tbl_found_truncated
tbl_not_found <- tbl_truncated %>%

tbl_pool

# -----------------------------------------------------------------------------
library(timsr)
base_path <- "/Users/adams/Projects/300K/2022-library-run/MaxQuant/LysN/"
path <- paste(base_path, "AspNLysN-A5_S1-A5_1_6546.d", sep = "")
all_columns <- c("frame", "scan", "tof", "intensity", "retention_time")
df <- TimsR(path) # get data handle
print(df)
df[2347]

df_x <- df[c(1,5,67), all_columns]
df_x <- df[13:354, 355]

print(df_x)