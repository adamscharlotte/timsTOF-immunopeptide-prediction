library(tidyverse)
library(data.table)

base_path <- "/Users/adams/Projects/300K/2022-library-run/"
pool <- "TUM_HLA2_140"
# "TUM_aspn_1_2", "TUM_HLA_16", TUM_lysn_5_no_enzyme, TUM_HLA2_122
plate <- "20220623_HLAII_p2"
# 20220624_HLAI_1_96, 20220621_AspNlysN, 20220623_HLAII_p2
info_path <- paste(base_path, "MaxQuant/info.txt", sep = "")

setwd(paste(base_path, "MaxQuant/", plate, "/", pool, "-txt/", sep = ""))
list.files()
getwd()

tbl_msms <- fread("msms.txt") %>% as_tibble
tbl_pool_msms <- tbl_msms %>% filter(Proteins == pool)
tbl_qc_msms <- tbl_msms %>% filter(stringr::str_starts(Proteins, "QC"))
tbl_reverse_msms <- tbl_msms %>%
    filter(!stringr::str_starts(Proteins, "QC")) %>%
    filter(!Proteins == pool) %>%
    select(Proteins, Sequence, Reverse)

tbl_pool_scans <- fread("accumulatedMsmsScans.txt") %>%
    as_tibble %>%
    filter(`Scan number` %in% tbl_pool_msms$`Scan number`) %>%
    mutate(Precursor = `PASEF precursor IDs`) %>%
    # Split the precursor ID and create a new row for each
    mutate(Precursor = strsplit(as.character(Precursor), ";")) %>%
    unnest(Precursor) %>%
    mutate(pool_name = pool) %>%
    unique()

tbl_unique_pool_peptides <- tbl_pool_scans %>% select(Sequence) %>% unique()

# This is also an option
# tbl_pool_msms %>% select(Sequence) %>% unique() # nolint

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

tbl_filtered_scans <- merge(tbl_pool_scans, tbl_pool, by = "pool_name") %>%
    as_tibble() %>%
    filter(stringr::str_ends(Sequence, trunc_sequence))

tbl_scans_contain <- merge(tbl_pool_scans, tbl_pool, by = "pool_name") %>%
    as_tibble() %>%
    filter(str_detect(Sequence_pool, Sequence))
colnames(tbl_filtered_scans)

tbl_scans_contain %>%
    filter(!Sequence %in% tbl_filtered_scans$Sequence) %>%
    select(Sequence, Sequence_pool, `Precursor retention time`) %>%
    unique() %>%
    print(n = 40)

tbl_filtered_scans %>%
    filter(Sequence_pool == "LTPEEEEILNKK") %>%
    select(Sequence, Sequence_pool, `Precursor retention time`,
    Modifications) %>%
    unique()

tbl_full_length <- tbl_filtered_scans %>%
    filter(Sequence == Sequence_pool) %>%
    select(Sequence, trunc_sequence) %>%
    unique()

tbl_truncated <- tbl_filtered_scans %>%
    filter(!Sequence %in% tbl_full_length$Sequence) %>%
    select(Sequence) %>%
    unique()

all_sequences <- tbl_filtered_scans %>%
    pull(Sequence) %>%
    unique()

tbl_wrong <- tbl_pool_scans %>%
    filter(!Sequence %in% all_sequences) %>%
    select(Sequence) %>%
    unique()

tbl_info <- tibble(pool_name = pool, PSMs = nrow(tbl_msms),
    pool_PSMs = nrow(tbl_pool_msms), qc_PSMs = nrow(tbl_qc_msms),
    reverse_PSMs = nrow(tbl_reverse_msms),
    pool_unique_peptides = nrow(tbl_unique_pool_peptides),
    full_length = nrow(tbl_full_length),
    truncated = nrow(tbl_truncated),
    wrong_identification = nrow(tbl_wrong),
    fasta_size = nrow(tbl_pool),
    percentage = full_length / fasta_size * 100)

fwrite(tbl_info, info_path, append = TRUE)
# ----------------------------------------------------------------------------

tbl_single_scan <- tbl_pasef %>% count(Precursor) %>% filter(n == 1)

# -----------------------------------------------------------------------------
library(timsr)

tbl_double_precursor <- tbl_pool_scans %>%
    select(Precursor) %>%
    count(Precursor) %>%
    filter(n > 1)

tbl_pool_scans %>%
    filter(Precursor %in% tbl_double_precursor$Precursor) %>%
    select(Sequence, Precursor, `Scan number`, `Precursor retention time`) %>%
    arrange(desc(Precursor))

tbl_pasef <- fread("pasefMsmsScans.txt") %>% as_tibble
tbl_scans_pasef <- merge(tbl_pool_scans, tbl_pasef, by = "Precursor") %>%
    as_tibble()

tbl_scans_pasef %>%
    filter(Precursor %in% tbl_double_precursor$Precursor) %>%
    select(Precursor, `Scan number`, Frame, ScanNumBegin, ScanNumEnd) %>%
    arrange(desc(Precursor)) %>%
    unique()

base_path <- "/Users/adams/Projects/300K/2022-library-run/MaxQuant/"
path <- paste(base_path, plate, "/HLAII_p2-D8_S1-D8_1_6777.d", sep = "")
all_columns <- c("frame", "scan", "tof", "intensity", "retention_time")
df <- TimsR(path) # get data handle
print(df)
tbl_raw <- df[c(2449, 2451, 2115, 2114, 2610, 3247),
    all_columns] %>% as_tibble()
tbl_raw %>% pull(retention_time) %>% unique()
tbl_raw %>% select(retention_time, frame) %>% unique()

tbl_raw %>% filter(frame == 3247) %>% pull(retention_time)

