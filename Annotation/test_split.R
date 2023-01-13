library(tidyverse)
library(data.table)
library(readr)
library(ggplot2)
# install.packages("hrbrthemes")
library(hrbrthemes)
library(viridis)

base <- "/Users/adams/Projects/300K/2022-library-run/"
meta_path <- paste(base, "Metadata/full-pool-sequence.txt", sep = "")
map_path <- paste(base, "Metadata/full-meta-map.txt", sep = "")
msms_path <- paste(base, "msms-txt/", sep = "")
# meta_qc_path <- paste(base, "Metadata/qc-peptides.txt", sep = "")
filter_path <- paste(base, "Annotation/total-scan-consensus/filtered/", sep = "") # nolint

tbl_meta <- fread(meta_path) %>% as_tibble()
tbl_meta %>%
    # filter(startsWith(pool_name, "TUM_first_pool")) %>%
    count(pool_name) %>%
    arrange(desc(n)) %>%
    print(n = 600)

tbl_map <- fread(map_path) %>%
    as_tibble() %>%
    select(plate, pool_name) %>%
    distinct()

tbl_meta_map <- tbl_meta %>%
    select(pool_name) %>%
    distinct() %>%
    merge(tbl_map, all = TRUE) %>%
    as_tibble() %>%
    mutate(plate = replace_na(plate, "first_pool"))

tbl_meta_map %>%
    count(plate) %>%
    distinct()

set.seed(7)

test_set <- tbl_meta_map %>%
    add_count(plate) %>%
    mutate(sample_size = n / 10) %>%
    group_by(plate) %>%
    sample_n(sample_size)

test_set %>% count(plate)

validation_set <- tbl_meta_map %>%
    add_count(plate) %>%
    mutate(sample_size = n / 10) %>%
    filter(!pool_name %in% test_set$pool_name) %>%
    group_by(plate) %>%
    sample_n(sample_size)

validation_set %>% count(plate)

validation_set %>%
    filter(pool_name %in% test_set$pool_name)

# ---------------------------- TEST THE SETS --------------------------------
spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

setwd(msms_path)
list_csv_files <- list.files(path = msms_path)
df_msms <- readr::read_tsv(list_csv_files, id = "file_name")

setwd(filter_path)
list_csv_files <- list.files(path = filter_path)
df_filtered <- readr::read_csv(list_csv_files, id = "file_name")

tbl_filtered <- df_filtered %>%
    mutate(pool_name = str_remove(file_name, ".csv"))

tbl_map <- df_msms %>%
    spaceless() %>%
    select(Proteins, Scan_number) %>%
    rename_with(toupper) %>%
    mutate(Proteins = PROTEINS)

tbl_filtered_mapped <- merge(tbl_filtered, tbl_map,
    by.x = c("SCAN_NUMBER", "pool_name"),
    by.y = c("SCAN_NUMBER", "PROTEINS"), all.x = TRUE) %>%
    as_tibble()

tbl_filtered_mapped %>% count(Proteins) %>% arrange(desc(n)) #%>% print(n=300)

tbl_test <- tbl_filtered_mapped %>%
    filter(pool_name %in% test_set$pool_name) %>%
    filter(!is.na(Proteins))

tbl_validation <- tbl_filtered_mapped %>%
    filter(pool_name %in% validation_set$pool_name) %>%
    filter(!is.na(Proteins))

tbl_train <- tbl_filtered_mapped %>%
    filter(!pool_name %in% validation_set$pool_name) %>%
    filter(!pool_name %in% tbl_test$pool_name)

nrow(tbl_test)
nrow(tbl_validation)
nrow(tbl_train)

tbl_test %>% filter(OBS_SEQUENCE %in% tbl_validation$OBS_SEQUENCE)

tbl_test_unique <- tbl_test %>%
    filter(!OBS_SEQUENCE %in% tbl_train$OBS_SEQUENCE) %>%
    filter(!OBS_SEQUENCE %in% tbl_validation$OBS_SEQUENCE)

tbl_validation_unique <- tbl_validation %>%
    filter(!OBS_SEQUENCE %in% tbl_train$OBS_SEQUENCE) %>%
    filter(!OBS_SEQUENCE %in% tbl_test$OBS_SEQUENCE)

tbl_test %>%
    filter(OBS_SEQUENCE %in% tbl_train$OBS_SEQUENCE) %>%
    filter(OBS_SEQUENCE == SEQUENCE) %>%
    count(Proteins) %>%
    arrange(desc(n))

tbl_train %>%
    filter(OBS_SEQUENCE %in% tbl_test$OBS_SEQUENCE) %>%
    filter(OBS_SEQUENCE == SEQUENCE) %>%
    count(Proteins) %>%
    arrange(desc(n))

tbl_validation %>%
    filter(OBS_SEQUENCE %in% tbl_train$OBS_SEQUENCE) %>%
    filter(OBS_SEQUENCE == SEQUENCE)

8319 + 8287 + 78522
8319/95128
8287/95128
78522/95128

7818 + 7765  + 78522
7818/94105
7765/94105
78522/94105

94105/97216
tbl_pool_peptides <- fread(meta_path) %>%
    as_tibble() %>%
    dplyr::rename(Proteins = pool_name) %>%
    distinct()

tbl_mapped <- merge(tbl_msms, tbl_pool_peptides, by = "Proteins") %>%
    as_tibble()

tbl_filtered <- tbl_mapped %>%
    filter(endsWith(Sequence, obs_sequence)) %>%
    filter(str_length(obs_sequence) <= str_length(Sequence)) %>%
    group_by(obs_sequence) %>%
    filter(str_length(Sequence) == min(str_length(Sequence))) %>%
    ungroup()

colnames(tbl_msms)
# -------------------- WRITE THE POOL NAMES TO A TXT FILE --------------------
test_set_t <- test_set$pool_name

write.table(test_set_t,
    file = "/Users/adams/Code/timsTOF-immunopeptide-prediction/Names/test-set.txt", # nolint
    col.names = FALSE,
    row.names = FALSE,
    sep = "")

validation_set_t <- validation_set$pool_name

write.table(validation_set_t,
    file = "/Users/adams/Code/timsTOF-immunopeptide-prediction/Names/validation-set.txt", # nolint
    col.names = FALSE,
    row.names = FALSE,
    sep = "")

# -----------------------------

c_overlap <- function(group1, group2) {
    tbl_1 <- tbl_meta %>% filter(pool_name == group1)
    tbl_2 <- tbl_meta %>%
        filter(pool_name == group2) %>%
        rename(pool_name2 = pool_name) %>%
        rename(Sequence2 = Sequence)
    tbl_overlap <- merge(tbl_1, tbl_2) %>%
        as_tibble() %>%
        filter(endsWith(Sequence, Sequence2))
    return(nrow(tbl_overlap))
}

tbl_meta_map %>% count(plate)
hla2_p2 <- tbl_meta_map %>%
    filter(plate == "20220623_HLAII_p2") %>%
    select(pool_name)

hla2_p2_compare <- hla2_p2 %>%
    rename(group1 = pool_name) %>%
    merge(hla2_p2) %>%
    as_tibble() %>%
    filter(!group1 == pool_name)

hla2_p2_overlap <- hla2_p2_compare %>%
    rowwise() %>%
    mutate(overlap = c_overlap(group1, pool_name))

# hla2_p2_overlap_B <- hla2_p2_compare %>%
#     rowwise() %>%
#     mutate(overlap = c_overlap(pool_name, group1))

# hla2_p2_overlap %>% arrange(desc(overlap))
# hla2_p2_overlap_B %>% arrange(desc(overlap))

ggplot(hla2_p2_overlap, aes(group1, pool_name, fill = overlap)) +
    geom_tile() +
    scale_fill_viridis(discrete = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

hla_p2_both_sides <- hla2_p2_overlap_B %>%
    rename(pool_name_t = group1) %>%
    rename(group1 = pool_name) %>%
    rename(pool_name = pool_name_t) %>%
    rbind(hla2_p2_overlap) %>%
    distinct()

hla2_p2 <- tbl_meta_map %>%
    filter(plate == "20220623_HLAII_p2") %>%
    select(pool_name)

hla2_p2_compare <- hla2_p2 %>%
    rename(group1 = pool_name) %>%
    merge(hla2_p2) %>%
    as_tibble() %>%
    filter(!group1 == pool_name)

hla2_p2_overlap <- hla2_p2_compare %>%
    rowwise() %>%
    mutate(overlap = c_overlap(group1, pool_name))

