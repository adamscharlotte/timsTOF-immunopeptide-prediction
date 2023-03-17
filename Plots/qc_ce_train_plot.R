library(tidyverse)
library(data.table)
library(readr)
library(plotly)

base <- "/Users/adams/Projects/300K/2022-library-run/"
train_non_tryptic_path <- paste(base, "Annotation/total-scan-consensus/split/non-tryptic/annotated-40-ppm-train.csv", sep = "") # nolint
train_tryptic_path <- paste(base, "Annotation/total-scan-consensus/split/tryptic/annotated-40-ppm-train.csv", sep = "") # nolint
test_non_tryptic_path <- paste(base, "Annotation/total-scan-consensus/split/non-tryptic/annotated-40-ppm-test.csv", sep = "") # nolint
test_tryptic_path <- paste(base, "Annotation/total-scan-consensus/split/tryptic/annotated-40-ppm-test.csv", sep = "") # nolint

tbl_train_non_tryptic <- fread(train_non_tryptic_path) %>% as_tibble()
tbl_train_tryptic <- fread(train_tryptic_path) %>% as_tibble()
tbl_test_non_tryptic <- fread(test_non_tryptic_path) %>% as_tibble()
tbl_test_tryptic <- fread(test_tryptic_path) %>% as_tibble()

tbl_test_tryptic %>%
    filter(!OBS_SEQUENCE == SEQUENCE)

    select(OBS_SEQUENCE) %>%
    distinct()





tbl_train_non_tryptic %>%
    filter(!OBS_SEQUENCE == SEQUENCE)

    select(OBS_SEQUENCE) %>%
    distinct()

tbl_train_tryptic %>%
    filter(!OBS_SEQUENCE == SEQUENCE)

    select(OBS_SEQUENCE) %>%
    distinct()

# -------------------------------- QC CE PLOT ---------------------------------

tbl_non_tryptic_qc <- tbl_train_non_tryptic %>% filter(Proteins == "")
tbl_tryptic_qc <- tbl_train_tryptic %>% filter(Proteins == "")

tbl_non_tryptic_qc %>%
    filter(PRECURSOR_CHARGE == 1) %>%
    add_count(OBS_SEQUENCE) %>%
    select(OBS_SEQUENCE, PRECURSOR_CHARGE, n) %>%
    arrange(desc(n)) %>%
    distinct()

# non-tryptic
# 99 distinct QC peptides
# 125 distinct QC peptides-charge combination

# 1  2160 -> "SYASDFGSSAK", "VGASTGYSGLK", "DTFLDGFSVK"
# 2  7001 -> "HDTVFGSYLYK", "GDFTFFIDTFK", "HFALFSTDVTK"
# 3   443 -> TFAHTESHISK, HEHISSDYAGK, YSAHEEHHYDK

# tryptic
# 93 distinct QC peptides
# 127 distinct QC peptides-charge combination

# 1  1223 -> c("TSIDSFIDSYK", "GIFGAFTDDYK", "GDFTFFIDTFK")
# 2  3283 -> "SILDYVSLVEK", "HTAYSDFLSDK", "YGFSSEDIFTK"
# c("SILDYVSLVEK", "HTAYSDFLSDK", "HTAYSDFLSDK",
#             "HDTVFGSYLYK", "GDFTFFIDTFK", "HFALFSTDVTK")
# 3   235 -> c("YSAHEEHHYDK", "HEHISSDYAGK", "TFAHTESHISK")

tbl_tryptic_qc_filtered <- tbl_tryptic_qc %>%
    filter(PRECURSOR_CHARGE == 1 &
        OBS_SEQUENCE == c("TSIDSFIDSYK", "GIFGAFTDDYK", "GDFTFFIDTFK",
            "SYASDFGSSAK", "VGASTGYSGLK", "DTFLDGFSVK")) %>%
    select(OBS_SEQUENCE, median_CE, pool_name) %>%
    mutate(Type = "Tryptic")

tbl_non_tryptic_qc_filtered <- tbl_non_tryptic_qc %>%
    filter(PRECURSOR_CHARGE == 1 &
        OBS_SEQUENCE == c("TSIDSFIDSYK", "GIFGAFTDDYK", "GDFTFFIDTFK",
            "SYASDFGSSAK", "VGASTGYSGLK", "DTFLDGFSVK")) %>%
    select(OBS_SEQUENCE, median_CE, pool_name) %>%
    mutate(Type = "Non-tryptic")

qc_filtered <- rbind(tbl_tryptic_qc_filtered, tbl_non_tryptic_qc_filtered)

fig <- plot_ly(qc_filtered, x = ~OBS_SEQUENCE, y = ~median_CE,
    color = ~Type, type = "box") %>%
    layout(boxmode = "group",
        title = 'Precursor charge = 1')

fig