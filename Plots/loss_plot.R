library(tidyverse)
library(data.table)
library(ggplot2)

file_name <- "log_precursor_20ppm_calibrated_refinement"
plot_title <- "Precursor - Calibrated CE - Refinement"
log_base <- "/Users/adams/Projects/300K/2022-library-run/Annotation/log/"
log_path <- paste(log_base, file_name, ".csv", sep = "")

tbl_log <- fread(log_path) %>% as_tibble()

plot_log <- ggplot(tbl_log, aes(x = epoch, y = loss)) +
    geom_line() +
    labs(title = plot_title)

plot_base <- "/Users/adams/Projects/300K/Results/Figures/Log/"
plot_path <- paste(plot_base, file_name, ".pdf", sep = "")
ggsave(plot_path, plot_log)