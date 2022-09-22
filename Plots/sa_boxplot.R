# install.packages("viridis")
library(tidyverse)
library(data.table)
library(ggplot2)
# library(viridis)

pool <- "TUM_aspn_28"
# pool <- "TUM_HLA2_7"
charge <- 2
base_path <- "/Users/adams/Projects/300K/2022-library-run/"
# sa_path <- paste(base_path, "Annotation/spectral-angle/frames/", pool, "_", charge, ".csv", sep = "") # nolint
sa_path <- paste(base_path, "Annotation/spectral-angle/precursors/calibrated_", pool, "_", charge, ".csv", sep = "") # nolint

tbl_cd_sa <- fread(sa_path) %>%
    as_tibble() %>%
    mutate(PRECURSOR = as.character(PRECURSOR))

plot <- ggplot(tbl_cd_sa, aes(x = SEQUENCE,
    y = SPECTRAL_ANGLE)) +
    geom_boxplot() +
    geom_jitter(shape = 16, position = position_jitter(0.2)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste(pool,
        " - Spectral angle compared with prosit prediction - ", charge))

# plot_path <- paste("/Users/adams/Projects/300K/Results/Figures/SA/frames/", pool, "_", charge, ".pdf", sep = "") # nolint
plot_path <- paste("/Users/adams/Projects/300K/Results/Figures/SA/precursors/calibrated_", pool, "_", charge, ".pdf", sep = "") # nolint
ggsave(plot_path, plot)
