library(tidyverse)
library(data.table)
library(ggplot2)

fixed_s_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/log/log_fixed_scratch.csv"
fixed_r_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/log/log_fixed_refinement.csv"
un_callibrated_r_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/log/log_un_callibrated_refinement.csv"

tbl_fixed_s <- fread(fixed_s_path) %>% as_tibble()
tbl_fixed_r <- fread(fixed_r_path) %>% as_tibble()
tbl_un_callibrated_r <- fread(un_callibrated_r_path) %>% as_tibble()

plot_s <- ggplot(tbl_fixed_s, aes(x=epoch, y=loss)) + 
    geom_line() +
    labs(title = "Training model from scratch - fixed CE")

plot_r <- ggplot(tbl_fixed_r, aes(x=epoch, y=loss)) + 
    geom_line() +
    labs(title = "Refining model - fixed CE")

plot_uc_r <- ggplot(tbl_un_callibrated_r, aes(x=epoch, y=loss)) + 
    geom_line() +
    labs(title = "Refining model - original CE")


ggsave("/Users/adams/Projects/300K/2022-library-run/Figures/log_fixed_scratch_plot.pdf", plot_s)
ggsave("/Users/adams/Projects/300K/2022-library-run/Figures/log_fixed_refinement_plot.pdf", plot_r)
ggsave("/Users/adams/Projects/300K/2022-library-run/Figures/log_un_callibrated_refinement_plot.pdf", plot_uc_r)