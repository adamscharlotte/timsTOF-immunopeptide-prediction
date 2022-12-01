library(tidyverse)
library(data.table)
library(ggplot2)

file_name <- "log_scan_40ppm_calibrated_frozence"
plot_title <- "Scan - Calibrated CE - Refinement"                      # Precursor Scratch Refinement # nolint

# file_name <- "un-summed/log_fixed_scratch"
# plot_title <- "Frames - Fixed CE - Scratch"                      #  Scratch Refinement # nolint
log_base <- "/Users/adams/Projects/300K/2022-library-run/Annotation/log/"
log_path <- paste(log_base, file_name, ".csv", sep = "")

# -------------- Mean SA between precursors of the same sequence --------------

base_path <- "/Users/adams/Projects/300K/2022-library-run/"
# sa_path <- paste(base_path, "Annotation/spectral-angle/pairwise/precursors/all-40-ppm/", sep = "") # nolint
sa_path <- paste(base_path, "Annotation/spectral-angle/pairwise/frames/all/", sep = "") # nolint
setwd(sa_path)
files <- dir(pattern = "*.csv")
df_sa <- map_df(files, ~read.csv(.x) %>% mutate(pool_name = basename(.x)))

tbl_sa <- df_sa %>%
    as_tibble() %>%
    # group_by(pool_name, SEQUENCE) %>%
    group_by(pool_name, PRECURSOR) %>%
    mutate(mean_sa = mean(SPECTRAL_ANGLE)) %>%
    replace(is.na(.), 0) %>%
    ungroup() %>%
    # select(pool_name, SEQUENCE, mean_sa) %>%
    select(pool_name, PRECURSOR, mean_sa) %>%
    distinct()

sa_mean <- mean(tbl_sa$mean_sa)

# ---------------------- Plot the losses and the mean SA ----------------------

tbl_log <- fread(log_path) %>% as_tibble()

plot_log <- ggplot(tbl_log, aes(x = epoch)) +
    geom_line(aes(y = loss, colour = "loss")) +
    geom_line(aes(y = val_loss, colour = "val_loss")) +
    labs(title = plot_title) +
    theme_minimal() +
    # geom_hline(yintercept = (1 - sa_mean), linetype = "dashed") +
    scale_color_manual(values = c("#a2bce0", "#f33b16")) +
    labs(color = "")
plot_log

plot_base <- "/Users/adams/Projects/300K/Results/Figures/Log/"
plot_path <- paste(plot_base, file_name, ".pdf", sep = "")
ggsave(plot_path, plot_log)