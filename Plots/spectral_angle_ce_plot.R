install.packages("viridis")
library(tidyverse)
library(data.table)
library(ggplot2)
library(viridis)

# pool <- "TUM_HLA_17"
pool <- "TUM_HLA2_7"
charge <- 2
base_path <- "/Users/adams/Projects/300K/2022-library-run/"
ce_sa_path <- paste(base_path, "Annotation/full-truncated-qc/spectral-angle/", pool, "_", charge, ".csv", sep = "") # nolint

tbl_cd_sa <- fread(ce_sa_path) %>%
    as_tibble() %>%
    mutate(ce_bin = as.character(round(ORIG_COLLISION_ENERGY, digits = 0))) %>%
    group_by(COLLISION_ENERGY, ce_bin) %>%
    mutate(mean_sa = mean(SPECTRAL_ANGLE)) %>%
    ungroup()

plot <- ggplot(tbl_cd_sa, aes(x = COLLISION_ENERGY,
    y = mean_sa, group = ce_bin, color = ce_bin)) +
    geom_line() +
    scale_color_viridis(discrete = TRUE) +
    labs(title = paste("Collision energy calibration - charge =", charge))

plot_path <- paste(base_path, "Figures/Calibration/", pool, "_", charge, ".pdf", sep = "") # nolint
ggsave(plot_path, plot)
