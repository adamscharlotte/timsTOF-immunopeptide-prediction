
library(tidyverse)
library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

pool_charge <- args[1]
# pool_charge <- "TUM_HLA_47_1"

pool <- str_sub(pool_charge, start = 1, end = -3)
charge <- str_sub(pool_charge, start = -1)

base_path <- "/Users/adams/Projects/300K/2022-library-run/"

# ---------------------------------- FRAME ------------------------------------
sa_path <- paste(base_path, "Annotation/spectral-angle/pairwise/frames/all/", pool, "_", charge, ".csv", sep = "") # nolint

tbl_cd_sa <- fread(sa_path) %>%
    as_tibble() %>%
    mutate(PRECURSOR = as.character(PRECURSOR))

plot <- ggplot(tbl_cd_sa, aes(x = PRECURSOR, y = SPECTRAL_ANGLE)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(shape = 16, width = 0.13) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste(pool,
        " - Pairwise comparison - ", charge)) +
    xlab("Precursor") +
    ylab("SA")
plot

plot_path <- paste("/Users/adams/Projects/300K/Results/Figures/SA/pairwise/frames/all/", pool, "_", charge, ".pdf", sep = "") # nolint
ggsave(plot_path, plot, width = 47, height = 14, units = "cm")

# -------------------------------- PRECURSOR ----------------------------------

sa_path <- paste(base_path, "Annotation/spectral-angle/pairwise/precursors/all/", pool, "_", charge, ".csv", sep = "") # nolint

tbl_cd_sa <- fread(sa_path) %>%
    as_tibble()

plot <- ggplot(tbl_cd_sa, aes(x = SEQUENCE, y = SPECTRAL_ANGLE)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(shape = 16, width = 0.13) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste(pool,
        " - Pairwise comparison - ", charge)) +
    xlab("Sequence") +
    ylab("SA")
plot

plot_path <- paste("/Users/adams/Projects/300K/Results/Figures/SA/pairwise/precursors/all/", pool, "_", charge, ".pdf", sep = "") # nolint
ggsave(plot_path, plot, width = 32, height = 14, units = "cm")
