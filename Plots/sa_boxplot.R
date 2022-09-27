# install.packages("viridis")
library(tidyverse)
library(data.table)
library(ggplot2)
# library(viridis)
args <- commandArgs(trailingOnly = TRUE)

pool_charge <- args[1]
# pool_charge <- "TUM_HLA2_7_1"

pool <- str_sub(pool_charge, start = 1, end = -3)
charge <- str_sub(pool_charge, start = -1)

base_path <- "/Users/adams/Projects/300K/2022-library-run/"
# sa_path <- paste(base_path, "Annotation/spectral-angle/frames/40-ppm/", pool, "_", charge, ".csv", sep = "") # nolint
sa_path <- paste(base_path, "Annotation/spectral-angle/precursors/40-ppm/", pool, "_", charge, ".csv", sep = "") # nolint

tbl_cd_sa <- fread(sa_path) %>%
    as_tibble() %>%
    mutate(PRECURSOR = as.character(PRECURSOR))

colnames(tbl_cd_sa)
tbl_cd_sa %>% select(SCORE, PREDICTED_INTENSITY, aligned_collision_energy, SEQUENCE) %>% distinct()

tbl_cd_sa %>% select(SCORE, SPECTRAL_ANGLE, ORIG_COLLISION_ENERGY) %>% arrange(SPECTRAL_ANGLE) %>% print(n = 35)

plot <- ggplot(tbl_cd_sa, aes(x = SEQUENCE,
    y = aligned_collision_energy)) +
    geom_boxplot() +
    geom_jitter(shape = 16, position = position_jitter(0.2), aes(colour = RETENTION_TIME)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste(pool,
        " - Spectral angle compared with prosit prediction - ", charge))

# plot_path <- paste("/Users/adams/Projects/300K/Results/Figures/SA/frames/", pool, "_", charge, ".pdf", sep = "") # nolint
plot_path <- paste("/Users/adams/Projects/300K/Results/Figures/SA/precursors/40-ppm/", pool, "_", charge, ".pdf", sep = "") # nolint
ggsave(plot_path, plot)
