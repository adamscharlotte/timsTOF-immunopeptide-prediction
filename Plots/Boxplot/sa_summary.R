
library(tidyverse)
library(data.table)
library(ggplot2)

base_path <- "/Users/adams/Projects/300K/2022-library-run/"
frame_path <- paste(base_path, "Annotation/spectral-angle/pairwise/frames/all/", sep = "") # nolint
setwd(frame_path)
files <- dir(pattern = "*.csv")
df_frame <- map_df(files, ~read.csv(.x) %>% mutate(pool_name = basename(.x)))
precursor_path <- paste(base_path, "Annotation/spectral-angle/pairwise/precursors/all-20-ppm/", sep = "") # nolint
setwd(precursor_path)
files <- dir(pattern = "*.csv")
df_precursor <- map_df(files, ~read.csv(.x) %>% mutate(pool_name = basename(.x))) # nolint

tbl_frame <- df_frame %>%
    as_tibble() %>%
    group_by(PRECURSOR) %>%
    mutate(mean_SA = mean(SPECTRAL_ANGLE)) %>%
    ungroup() %>%
    select(PRECURSOR, mean_SA) %>%
    distinct() %>%
    select(mean_SA) %>%
    # select(SPECTRAL_ANGLE) %>%
    mutate(level = "Frame")

tbl_precursor <- df_precursor %>%
    as_tibble() %>%
    group_by(SEQUENCE) %>%
    mutate(mean_SA = mean(SPECTRAL_ANGLE)) %>%
    ungroup() %>%
    select(SEQUENCE, mean_SA) %>%
    distinct() %>%
    select(mean_SA) %>%
    # select(SPECTRAL_ANGLE) %>%
    mutate(level = "Precursor")

tbl_plot <- rbind(tbl_frame, tbl_precursor)
# plot <- ggplot(tbl_plot, aes(x = level, y = SPECTRAL_ANGLE)) +
plot <- ggplot(tbl_plot, aes(x = level, y = mean_SA)) +
    geom_boxplot(alpha = 0.7) +
    # geom_jitter(shape = 16, width = 0.13) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste(pool,
        " - Pairwise comparison - ", charge)) +
    xlab("Precursor") +
    ylab("SA")
plot

plot_path <- "/Users/adams/Projects/300K/Results/Figures/SA/pairwise/mean-all-levels.pdf" # nolint
ggsave(plot_path, plot, width = 30, height = 14, units = "cm")