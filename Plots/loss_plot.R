library(tidyverse)
library(data.table)
library(ggplot2)

file_name <- "log_precursor_20ppm_calibrated_frozence_scratch"
plot_title <- "Precursor - Calibrated CE - Scratch"                      #  Scratch Refinement # nolint
log_base <- "/Users/adams/Projects/300K/2022-library-run/Annotation/log/"
log_path <- paste(log_base, file_name, ".csv", sep = "")

# -------------- Mean SA between precursors of the same sequence --------------

base_path <- "/Users/adams/Projects/300K/2022-library-run/"
sa_path <- paste(base_path, "Annotation/spectral-angle/pairwise/precursors/all-20-ppm/", sep = "") # nolint
setwd(sa_path)
files <- dir(pattern = "*.csv")
df_sa <- map_df(files, ~read.csv(.x) %>% mutate(pool_name = basename(.x)))

tbl_sa_ncomb <- df_sa %>%
    as_tibble() %>%
    # mutate(n_comb = count(pool_name, SEQUENCE)) %>%
    group_by(pool_name, SEQUENCE) %>%
    summarise(n_comb = n()) %>%
    ungroup()

tbl_sa <- df_sa %>%
    as_tibble() %>%
    # mutate(n_comb = count(pool_name, SEQUENCE)) %>%
    group_by(pool_name, SEQUENCE) %>%
    mutate(mean_sa = mean(SPECTRAL_ANGLE)) %>%
    mutate(sd_sa = sd(SPECTRAL_ANGLE)) %>%
    replace(is.na(.), 0) %>%
    ungroup() %>%
    merge(tbl_sa_ncomb, by = c("pool_name", "SEQUENCE")) %>%
    as_tibble() %>%
    select(pool_name, SEQUENCE, sd_sa, mean_sa, n_comb) %>%
    distinct()

n <- tbl_sa %>% nrow()

# sa_mean <- (sum(tbl_sa$mean_sa * tbl_sa$n_comb) / sum(tbl_sa$n_comb))
# sa_sd <- sqrt(sum(tbl_sa$sd_sa * tbl_sa$n_comb) / sum(tbl_sa$n_comb))
# sa_1_mean <- 1 - mean(tbl_sa$mean_sa)
# sa_sd <- sum(tbl_sa$sd_sa) / sqrt(n)

tbl_log <- fread(log_path) %>% as_tibble()

plot_log <- ggplot(tbl_log, aes(x = epoch)) +
    geom_line(aes(y = loss, colour = "loss")) +
    geom_line(aes(y = val_loss, colour = "val_loss")) +
    labs(title = plot_title) +
    theme_minimal() +
    geom_hline(yintercept = sa_1_mean, linetype = "dashed") +
    geom_hline(yintercept = (sa_1_mean + sa_sd), color = "lightgrey") +
    geom_hline(yintercept = (sa_1_mean - sa_sd), color = "lightgrey") +
    scale_color_manual(values = c("#a2bce0", "#f33b16")) +
    labs(color = "")
plot_log

plot_base <- "/Users/adams/Projects/300K/Results/Figures/Log/"
plot_path <- paste(plot_base, file_name, ".pdf", sep = "")
ggsave(plot_path, plot_log)