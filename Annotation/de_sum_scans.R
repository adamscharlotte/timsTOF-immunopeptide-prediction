library(tidyverse)
library(data.table)
library(timsr)

base_path <- "/Users/adams/Projects/300K/2022-library-run/"
pool <- "TUM_HLA_71"
# "TUM_aspn_1_2", "TUM_HLA_16", TUM_lysn_5_no_enzyme, TUM_HLA2_122
# TUM_HLA_71, TUM_lysn_12
plate <- paste("20220624_HLAI_1_96", "/", sep = "")
# 20220624_HLAI_1_96, 20220621_AspNlysN, 20220623_HLAII_p2
raw_folder <- "HLAI_1_96_p2-F11_S2-F11_1_6996.d"
# HLAII_p2-D8_S1-D8_1_6777.d, AspNLysN-A12_S1-A12_1_6553
meta_path <- paste(base_path, "metadata/full-pool-sequence.txt", sep = "")
info_path <- paste(base_path, "MaxQuant/info.txt", sep = "")
raw_folder_path <- paste(base_path, "MaxQuant/", plate, raw_folder, sep = "")

setwd(paste(base_path, "MaxQuant/", plate, pool, "-txt/", sep = ""))
list.files()
getwd()

# Load PSMs
tbl_psms %>% select(Proteins) %>% distinct()
tbl_psms <- fread("msms.txt") %>% as_tibble
tbl_pool_psms <- tbl_psms %>% filter(Proteins == pool)
tbl_qc_psms <- tbl_psms %>% filter(stringr::str_starts(Proteins, "QC"))
tbl_reverse_psms <- tbl_psms %>%
    filter(!stringr::str_starts(Proteins, "QC")) %>%
    filter(!Proteins == pool) %>%
    select(Proteins, Sequence, Reverse)
tbl_unique_pool_peptides <- tbl_pool_psms %>% select(Sequence) %>% distinct()

# Load precursors
tbl_pool_precursors <- fread("accumulatedMsmsScans.txt") %>%
    as_tibble %>%
    filter(`Scan number` %in% tbl_pool_psms$`Scan number`) %>%
    mutate(Precursor = `PASEF precursor IDs`) %>%
    # Split the precursor ID and create a new row for each
    mutate(Precursor = strsplit(as.character(Precursor), ";")) %>%
    unnest(Precursor) %>%
    mutate(pool_name = pool) %>%
    distinct()

# Retreive all expected peptide sequences in the pool
tbl_meta <- fread(meta_path) %>% as_tibble
tbl_pool <- tbl_meta %>%
    filter(pool_name == pool) %>%
    mutate(trunc_sequence = str_sub(Sequence, start = -7)) %>%
    rename(Sequence_pool = Sequence) %>%
    distinct()

# Filter the identifications and retreive truncated and full-length peptides
tbl_filtered_precursors <- merge(tbl_pool_precursors,
    tbl_pool, by = "pool_name") %>%
    as_tibble() %>%
    filter(stringr::str_ends(Sequence, trunc_sequence))

# Load the precursor frames (filter identified scans)
tbl_frames <- fread("pasefMsmsScans.txt") %>%
    as_tibble %>%
    filter(Precursor %in% tbl_pool_precursors$Precursor) %>%
    mutate(frame = Frame) %>%
    select(Precursor, frame)

frame_min <- min(tbl_frames$frame)
frame_max <- max(tbl_frames$frame)

# Retreive the retention time for each scan (frame)
all_columns <- c("frame", "scan", "tof", "intensity", "retention_time")
some_columns <- c("frame", "retention_time")
df <- TimsR(raw_folder_path) # get data handle
tbl_raw <- df[frame_min : frame_max, some_columns] %>%
    as_tibble() %>%
    distinct()

tbl_pool_scans <- merge(tbl_pool_precursors, tbl_frames, by = "Precursor") %>%
    merge(tbl_raw, by = "frame") %>%
    as_tibble() %>%
    distinct()

colnames(tbl_pool_scans)

tbl_pool_scans %>%
    group_by(Precursor) %>%
    mutate(mean_precursor_rt = mean(retention_time)) %>%
    ungroup() %>%
    group_by(`Scan number`) %>%
    mutate(mean_rt = mean(mean_precursor_rt)) %>%
    ungroup() %>%
    mutate(mean_rt_min = mean_rt / 60) %>%
    mutate(rt_difference = mean_rt_min - `Precursor retention time`) %>%
    arrange(desc(rt_difference)) %>%
    select(rt_difference) %>%
    distinct() %>%
    print(n = 1000)

# ----------------------------------------------------------------------------

tbl_precursor_contain <- merge(tbl_pool_precursors,
    tbl_pool, by = "pool_name") %>%
    as_tibble() %>%
    filter(str_detect(Sequence_pool, Sequence)) %>%
    select(Sequence, Sequence_pool) %>%
    distinct()

tbl_full_length %>% select(Sequence) %>% distinct()

tbl_precursor_contain %>% select(Sequence) %>% distinct()
tbl_mapped_wrong <- tbl_precursor_contain %>%
    filter(!Sequence %in% all_sequences) %>%
    select(Sequence, Sequence_pool) %>%
    distinct()

tbl_full_length <- tbl_filtered_precursors %>%
    filter(Sequence == Sequence_pool) %>%
    select(Sequence, trunc_sequence) %>%
    distinct()

tbl_truncated <- tbl_filtered_precursors %>%
    filter(!Sequence %in% tbl_full_length$Sequence) %>%
    select(Sequence) %>%
    distinct()

all_sequences <- tbl_filtered_precursors %>%
    pull(Sequence) %>%
    unique()

tbl_wrong <- tbl_pool_precursors %>%
    filter(!Sequence %in% all_sequences) %>%
    select(Sequence) %>%
    unique()

tbl_wrong$Sequence %in% tbl_filtered_precursors$Sequence

# ----------------------------------------------------------------------------

tbl_info <- tibble(pool_name = pool, PSMs = nrow(tbl_psms),
    pool_PSMs = nrow(tbl_pool_psms), qc_PSMs = nrow(tbl_qc_psms),
    reverse_PSMs = nrow(tbl_reverse_psms),
    pool_unique_peptides = nrow(tbl_unique_pool_peptides),
    full_length = nrow(tbl_full_length),
    truncated = nrow(tbl_truncated),
    wrong_identification = nrow(tbl_wrong),
    fasta_size = nrow(tbl_pool),
    percentage = full_length / fasta_size * 100)

fwrite(tbl_info, info_path, append = TRUE)
# ----------------------------------------------------------------------------

tbl_mapped_wrong
tbl_pool_scans_mapped <- tbl_pool_scans %>%
    merge(tbl_precursor_contain, by = "Sequence") %>%
    as_tibble() %>%
    distinct()

colnames(tbl_pool_scans)

tbl_pool_scans_mapped %>%
    filter(Sequence_pool %in% tbl_mapped_wrong$Sequence_pool) %>%


tbl_precursor_contain

# ----------------------------------------------------------------------------

tbl_single_scan <- tbl_pasef %>% count(Precursor) %>% filter(n == 1)

# -----------------------------------------------------------------------------

tbl_double_precursor <- tbl_pool_precursors %>%
    select(Precursor) %>%
    count(Precursor) %>%
    filter(n > 1)

tbl_pool_precursors %>%
    filter(Precursor %in% tbl_double_precursor$Precursor) %>%
    select(Sequence, Precursor, `Scan number`, `Precursor retention time`) %>%
    arrange(desc(Precursor))


tbl_scans_pasef %>%
    filter(Precursor %in% tbl_double_precursor$Precursor) %>%
    select(Precursor, `Scan number`, Frame, ScanNumBegin, ScanNumEnd) %>%
    arrange(desc(Precursor)) %>%
    unique()
