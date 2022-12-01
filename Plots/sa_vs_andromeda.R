library(tidyverse)
library(data.table)
library(ggplot2)

ppm <- "20"
base_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/"
sum_path <- paste(base_path, "precursor-consensus/annotated-", ppm, "-ppm/", sep = "") # nolint
sa_path <- paste(base_path, "spectral-angle/pairwise/precursors/all-", ppm, "-ppm/", sep = "") # nolint
setwd(sa_path)
files <- dir(pattern = "*.csv")
df_sa <- map_df(files, ~read.csv(.x) %>%
    mutate(pool_charge = str_sub(basename(.x), start = 1, end = -5)))
setwd(sum_path)
sum_files <- str_sub(files, start = 1, end = -7) %>% unique()
sum_files <- paste(sum_files, ".csv", sep = "")
df_sum <- map_df(sum_files, ~read.csv(.x) %>%
    mutate(pool_name = str_sub(basename(.x), start = 1, end = -5)))
tbl_sum <- df_sum %>%
    as_tibble() %>%
    select(pool_name, PRECURSOR_CHARGE, SEQUENCE, SCORE, PRECURSOR) %>%
    unite("pool_charge", pool_name:PRECURSOR_CHARGE, sep = "_") %>%
    distinct()

tbl_sa <- df_sa %>%
    as_tibble() %>%
    merge(tbl_sum, by.x = c("pool_charge", "A", "SEQUENCE"),
        by.y = c("pool_charge", "PRECURSOR", "SEQUENCE")) %>%
    as_tibble() %>%
    rename(score_A = SCORE) %>%
    merge(tbl_sum, by.x = c("pool_charge", "B", "SEQUENCE"),
        by.y = c("pool_charge", "PRECURSOR", "SEQUENCE")) %>%
    as_tibble() %>%
    rename(score_B = SCORE) %>%
    # only look at SAs of precursors with the same score
    filter(score_A == score_B) %>%
    group_by(pool_charge, SEQUENCE) %>%
    mutate(mean_sa = mean(SPECTRAL_ANGLE)) %>%
    mutate(sd_sa = sd(SPECTRAL_ANGLE)) %>%
    replace(is.na(.), 0) %>%
    ungroup() %>%
    select(pool_charge, SEQUENCE, mean_sa, score_A) %>%
    distinct() %>%
    as_tibble() %>%
    distinct()

plot <- ggplot(tbl_sa, aes(x = score_A, y = mean_sa)) +
    geom_point(alpha = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste("Mean SA VS Andromeda Score - ", ppm, "ppm")) +
    xlab("Andromeda Score") +
    ylab("mean SA")

plot_path <- paste("/Users/adams/Projects/300K/Results/Figures/SA/sa-vs-score/",ppm, "-ppm.pdf", sep = "") # nolint
ggsave(plot_path, plot, width = 32, height = 14, units = "cm")
