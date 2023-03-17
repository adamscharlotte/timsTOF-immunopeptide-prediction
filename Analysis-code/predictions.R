library(tidyverse)
library(data.table)
library(ggplot2)

base_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/"  # nolint
file_1_path <- paste(base_path, "20230315104417_prediction.csv", sep = "")
file_2_path <- paste(base_path, "20230315111004_prediction.csv", sep = "")

tbl_file_1 <- fread(file_1_path) %>%
    as_tibble() %>%
    mutate(label = "HCD Prosit 2020") %>%
    mutate(group = substr(rawfile, 1, 7))

tbl_file_2 <- fread(file_2_path) %>%
    as_tibble() %>%
    mutate(label = "HCD Bruker Prosit 2023") %>%
    mutate(group = substr(rawfile, 1, 7))

mean(tbl_file_1$spectral_angle)
median(tbl_file_1$spectral_angle)

tbl_file_1 %>%
    filter(str_detect(rawfile, "first_pool")) %>%
    filter(spectral_angle >= 0.9)
319/22122
# 1%
192/7872
# 2.4%
127/14260 # 2 are missing!
# 0.9%

tbl_file_2 %>%
    filter(str_detect(rawfile, "TUM_first_pool")) %>%
    filter(spectral_angle >= 0.9)
8067/22122
# 36.4%
2074/7872
# 26.3%
6003/14260 # 2 are missing!
# 42.1%

median_group <- function(tbl_file, str) {
    tbl_tmp <- tbl_file %>% filter(startsWith(rawfile, str))
    median(tbl_tmp$spectral_angle)
}

median_group(tbl_file_1, "b'HLAI_")
median_group(tbl_file_2, "b'HLAI_")
median_group(tbl_file_1, "b'HLAII")
median_group(tbl_file_2, "b'HLAII")
median_group(tbl_file_1, "b'AspNLysN")
median_group(tbl_file_2, "b'AspNLysN")
median_group(tbl_file_1, "b'230719_f1-TUM_first_pool")
median_group(tbl_file_2, "b'230719_f1-TUM_first_pool")




mean(tbl_file_2$spectral_angle)
median(tbl_file_2$spectral_angle)

tbl_file_1 %>% select(rawfile) %>% distinct() %>% print(n = 50)