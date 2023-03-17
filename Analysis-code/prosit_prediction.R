library(tidyverse)
library(data.table)

base_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/spectral-angle-comparison" # nolint

csv_names <- list.files(path = base_path, pattern = "*.csv", full.names = TRUE)
csv_names2 <- data.frame(month = csv_names, id = as.character(1:length(csv_names)))

data <- csv_names %>%
    lapply(read_csv) %>% # read all the files at once
    bind_rows(.id = "id") %>% # bind all tables into one object, and give id for each
    left_join(csv_names2) # join month column created earlier

data_1 <- data %>%
    distinct() %>%
    filter(PRECURSOR_CHARGE == 1)

data_2 <- data %>%
    distinct() %>%
    filter(PRECURSOR_CHARGE == 2)

data_3 <- data %>%
    distinct() %>%
    filter(PRECURSOR_CHARGE == 3)

mean(data$SPECTRAL_ANGLE)
mean(data_1$SPECTRAL_ANGLE)
mean(data_2$SPECTRAL_ANGLE)
mean(data_3$SPECTRAL_ANGLE)
