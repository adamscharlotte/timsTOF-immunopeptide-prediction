library(tidyverse)
library(data.table)
library(readr)
# library(ggplot2)
# # install.packages("hrbrthemes")
# library(hrbrthemes)
# library(viridis)
# library(gtools)

base <- "/Users/adams/Projects/300K/2022-library-run/"
meta_path <- paste(base, "Metadata/full-pool-sequence.txt", sep = "")
map_path <- paste(base, "Metadata/full-meta-map.txt", sep = "")
msms_path <- paste(base, "msms-txt/", sep = "")
# meta_qc_path <- paste(base, "Metadata/qc-peptides.txt", sep = "")
filter_path <- paste(base, "Annotation/total-scan-consensus/annotated-40-ppm/", sep = "") # nolint

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
    filter(!is.na(Proteins)) %>%
    mutate(set = "test")

tbl_validation <- tbl_filtered_mapped %>%
    filter(pool_name %in% validation_set$pool_name) %>%
    filter(!is.na(Proteins)) %>%
    mutate(set = "validation")

tbl_train <- tbl_filtered_mapped %>%
    filter(!pool_name %in% validation_set$pool_name) %>%
    filter(!pool_name %in% tbl_test$pool_name) %>%
    mutate(set = "train")

nrow(tbl_test)
nrow(tbl_validation)
nrow(tbl_train)

24060 + 26297 + 240689

all_sets <- bind_rows(tbl_test, tbl_validation, tbl_train)

peptide_counts <- all_sets %>%
    # filter(PRECURSOR_CHARGE < 4) %>%
    group_by(set) %>%
    count(OBS_SEQUENCE, name = "count") %>%
    ungroup()

top_peptides <- peptide_counts %>%
    group_by(OBS_SEQUENCE) %>%
    top_n(1, count) %>%
    ungroup()

peptides_to_keep <- top_peptides %>%
    add_count(OBS_SEQUENCE, name = "pep_count") %>%
    filter(pep_count == 1)

tbl_test_filtered <- tbl_test %>%
    # filter(PRECURSOR_CHARGE < 4) %>%
    filter(OBS_SEQUENCE %in%
        peptides_to_keep[peptides_to_keep$set == "test", ]$OBS_SEQUENCE)
tbl_train_filtered <- tbl_train %>%
    # filter(PRECURSOR_CHARGE < 4) %>%
    filter(OBS_SEQUENCE %in%
        peptides_to_keep[peptides_to_keep$set == "train", ]$OBS_SEQUENCE)
tbl_validation_filtered <- tbl_validation %>%
    # filter(PRECURSOR_CHARGE < 4) %>%
    filter(OBS_SEQUENCE %in%
        peptides_to_keep[peptides_to_keep$set == "validation", ]$OBS_SEQUENCE)

nrow(tbl_test)
nrow(tbl_validation)
nrow(tbl_train)

nrow(tbl_test_filtered)
nrow(tbl_validation_filtered)
nrow(tbl_train_filtered)

22749 + 25014 + 236873
23062 + 25241 + 238466
23568 + 25803 + 239404
23717 + 25968 + 239595

# ------------------ WRITE THE SUMMED PEPTIDES FOR EACH SET ------------------
test_path <- paste(base, "Annotation/total-scan-consensus/split/tryptic/annotated-40-ppm-test.csv", sep = "") # nolint
train_path <- paste(base, "Annotation/total-scan-consensus/split/tryptic/annotated-40-ppm-train.csv", sep = "") # nolint
validation_path <- paste(base, "Annotation/total-scan-consensus/split/tryptic/annotated-40-ppm-validation.csv", sep = "") # nolint

tbl_test_filtered %>%
    distinct() %>%
    filter(startsWith(pool_name, "TUM_first_pool")) %>%
    filter(PRECURSOR_CHARGE %in% c(2,3)) %>%
    count(PRECURSOR_CHARGE)
    select(OBS_SEQUENCE) %>%
    distinct()
    filter(!SEQUENCE == OBS_SEQUENCE)

    %>%
    fwrite(., test_path)

tbl_train_filtered %>%
    filter(startsWith(pool_name, "TUM_first_pool")) %>%
    filter(PRECURSOR_CHARGE %in% c(2,3)) %>%
    select(OBS_SEQUENCE) %>%
    distinct()

    filter(!SEQUENCE == OBS_SEQUENCE)

    fwrite(., train_path)

tbl_validation_filtered %>%
    filter(startsWith(pool_name, "TUM_first_pool")) %>%
    filter(PRECURSOR_CHARGE %in% c(2,3)) %>%
    select(OBS_SEQUENCE) %>%
    distinct()
    filter(!SEQUENCE == OBS_SEQUENCE)



    fwrite(., validation_path)

test_path <- paste(base, "Annotation/total-scan-consensus/split/non-tryptic/annotated-40-ppm-test.csv", sep = "") # nolint
train_path <- paste(base, "Annotation/total-scan-consensus/split/non-tryptic/annotated-40-ppm-train.csv", sep = "") # nolint
validation_path <- paste(base, "Annotation/total-scan-consensus/split/non-tryptic/annotated-40-ppm-validation.csv", sep = "") # nolint

tbl_test_filtered %>%
    filter(!startsWith(pool_name, "TUM_first_pool")) %>%
    fwrite(., test_path)
tbl_train_filtered %>%
    filter(!startsWith(pool_name, "TUM_first_pool")) %>%
    fwrite(., train_path)
tbl_validation_filtered %>%
    filter(!startsWith(pool_name, "TUM_first_pool")) %>%
    fwrite(., validation_path)

tbl_train_filtered %>%
    filter(!startsWith(pool_name, "TUM_first_pool")) %>%
    filter(PRECURSOR_CHARGE == 4)

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

# -------------------------------- QC CE PLOT ---------------------------------

tbl_train_filtered %>%
    filter(is.na(Proteins)) %>%
    count(PRECURSOR_CHARGE)

# 121 distinct QC peptides
# 177 distinct QC peptides-charge combination

# 1  3383
# 2 10284
# 3   678
