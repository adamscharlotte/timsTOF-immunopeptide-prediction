# install.packages("viridis")
# install.packages("reshape2")
library(reshape2)
library(tidyverse)
library(data.table)
library(ggplot2)
# library(viridis)
args <- commandArgs(trailingOnly = TRUE)

pool_charge <- args[1]
# pool_charge <- "TUM_HLA_47_1"

pool <- str_sub(pool_charge, start = 1, end = -3)
charge <- str_sub(pool_charge, start = -1)

base_path <- "/Users/adams/Projects/300K/2022-library-run/"

# ---------------------------------- FRAME ------------------------------------
sa_path <- paste(base_path, "Annotation/spectral-angle/frames/40-ppm/", pool, "_", charge, ".csv", sep = "") # nolint

tbl_cd_sa <- fread(sa_path) %>%
    as_tibble() %>%
    mutate(PRECURSOR = as.character(PRECURSOR))

tbl_cd_sa_melt <- melt(tbl_cd_sa, id.vars = 'FRAME',
    measure.vars = c('SPECTRAL_ANGLE', 'SPECTRAL_ANGLE_CAL'))
tbl_plot <- tbl_cd_sa %>%
    select(FRAME, RETENTION_TIME, ORIG_COLLISION_ENERGY, PRECURSOR) %>%
    merge(tbl_cd_sa_melt) %>%
    as_tibble() %>%
    mutate(Calibration = str_replace(variable, "SPECTRAL_ANGLE_CAL", "TRUE")) %>%
    mutate(Calibration = str_replace(Calibration, "SPECTRAL_ANGLE", "FALSE"))

plot <- ggplot(tbl_plot, aes(x = PRECURSOR,
    y = value, fill = Calibration)) +
    geom_boxplot(alpha = 0.7) +
    geom_point(position = position_jitterdodge(), shape = 16,
        aes(colour = Calibration)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste(pool,
        " - Spectral angle compared with prosit prediction - ", charge)) +
    xlab("Precursor") +
    ylab("SA") +
    scale_fill_manual(values = c("#f33b16", "#a2bce0")) +
    scale_colour_manual(values = c("#f33b16", "#a2bce0"))

plot_path <- paste("/Users/adams/Projects/300K/Results/Figures/SA/frames/40-ppm/", pool, "_", charge, ".pdf", sep = "") # nolint
ggsave(plot_path, plot)

# -------------------------------- PRECURSOR ----------------------------------

sa_path <- paste(base_path, "Annotation/spectral-angle/precursors/40-ppm/", pool, "_", charge, ".csv", sep = "") # nolint

tbl_cd_sa <- fread(sa_path) %>%
    as_tibble() %>%
    mutate(PRECURSOR = as.character(PRECURSOR))

tbl_cd_sa_melt <- melt(tbl_cd_sa, id.vars = 'PRECURSOR',
    measure.vars = c('SPECTRAL_ANGLE', 'SPECTRAL_ANGLE_CAL'))
tbl_plot <- tbl_cd_sa %>%
    select(PRECURSOR, RETENTION_TIME, ORIG_COLLISION_ENERGY, SEQUENCE) %>%
    merge(tbl_cd_sa_melt) %>%
    as_tibble() %>%
    mutate(Calibration = str_replace(variable, "SPECTRAL_ANGLE_CAL", "TRUE")) %>%
    mutate(Calibration = str_replace(Calibration, "SPECTRAL_ANGLE", "FALSE"))

plot <- ggplot(tbl_plot, aes(x = SEQUENCE,
    y = value, fill = Calibration)) +
    geom_boxplot(alpha = 0.7) +
    geom_point(position = position_jitterdodge(), shape = 16,
        aes(colour = Calibration)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste(pool,
        " - Spectral angle compared with prosit prediction - ", charge)) +
    xlab("Sequence") +
    ylab("SA") +
    scale_fill_manual(values = c("#f33b16", "#a2bce0")) +
    scale_colour_manual(values = c("#f33b16", "#a2bce0"))

plot_path <- paste("/Users/adams/Projects/300K/Results/Figures/SA/precursors/40-ppm/", pool, "_", charge, ".pdf", sep = "") # nolint
ggsave(plot_path, plot)
