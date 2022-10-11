library(tidyverse)
library(data.table)
library(ggplot2)

file_name <- "log_precursor_20ppm_calibrated_scratch"
plot_title <- "Precursor - Calibrated CE - Scratch"                      #  Scratch Refinement # nolint
log_base <- "/Users/adams/Projects/300K/2022-library-run/Annotation/log/"
log_path <- paste(log_base, file_name, ".csv", sep = "")

tbl_log <- fread(log_path) %>% as_tibble()

plot_log <- ggplot(tbl_log, aes(x = epoch)) +
    geom_line(aes(y = loss, colour = "loss")) +
    geom_line(aes(y = val_loss, colour = "val_loss")) +
    labs(title = plot_title) +
    theme_minimal() +
    scale_color_manual(values = c("#a2bce0", "#f33b16")) +
    labs(color = "")
plot_log

plot_base <- "/Users/adams/Projects/300K/Results/Figures/Log/"
plot_path <- paste(plot_base, file_name, ".pdf", sep = "")
ggsave(plot_path, plot_log)