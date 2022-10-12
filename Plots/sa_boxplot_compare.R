
library(tidyverse)
library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

pool_charge <- args[1]
# pool_charge <- "TUM_HLA2_88_1"

pool <- str_sub(pool_charge, start = 1, end = -3)
charge <- str_sub(pool_charge, start = -1)

base_path <- "/Users/adams/Projects/300K/2022-library-run/"

# -------------------------------- PRECURSOR ----------------------------------

sa_20_path <- paste(base_path, "Annotation/spectral-angle/pairwise/precursors/all-20-ppm/", pool, "_", charge, ".csv", sep = "") # nolint
sa_40_path <- paste(base_path, "Annotation/spectral-angle/pairwise/precursors/all-40-ppm/", pool, "_", charge, ".csv", sep = "") # nolint

tbl_sa_20 <- fread(sa_20_path) %>%
    as_tibble() %>%
    select(combination, SPECTRAL_ANGLE, SEQUENCE) %>%
    mutate(Accuracy = "20 ppm")

tbl_sa_40 <- fread(sa_40_path) %>%
    as_tibble() %>%
    select(combination, SPECTRAL_ANGLE, SEQUENCE) %>%
    mutate(Accuracy = "40 ppm")

tbl_plot <- rbind(tbl_sa_20, tbl_sa_40)

plot <- ggplot(tbl_plot, aes(x = SEQUENCE,
    y = SPECTRAL_ANGLE, fill = Accuracy)) +
    geom_boxplot(alpha = 0.7) +
    geom_point(position = position_jitterdodge(), shape = 16,
        aes(colour = Accuracy)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste(pool,
        " - Spectral angle compared with prosit prediction - ", charge)) +
    xlab("Sequence") +
    ylab("SA") +
    scale_fill_manual(values = c("#f33b16", "#a2bce0")) +
    scale_colour_manual(values = c("#f33b16", "#a2bce0"))

plot_path <- paste("/Users/adams/Projects/300K/Results/Figures/SA/pairwise/precursors/20-vs-40-ppm/", pool, "_", charge, ".pdf", sep = "") # nolint
ggsave(plot_path, plot, width = 32, height = 14, units = "cm")
