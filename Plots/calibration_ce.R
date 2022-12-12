library(tidyverse)
library(data.table)
library(ggplot2)
library(cowplot)
# install.packages("ggpubr")
library(ggpubr)
# install.packages("patchwork")
# library(patchwork)

base_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/"
sa_path <- paste(base_path, "scan-consensus/spectral-angle/", sep = "")
setwd(sa_path)

tbl_charge_1 <-
    list.files(pattern = "*_1.csv") %>%
    map(read_csv) %>%       # read in all the files individually, using
    reduce(rbind)

tbl_charge_2 <-
    list.files(pattern = "*_2.csv") %>%
    map(read_csv) %>%
    reduce(rbind)

tbl_charge_3 <-
    list.files(pattern = "*_3.csv") %>%
    map(read_csv) %>%
    reduce(rbind)

tbl_charge_4 <-
    list.files(pattern = "*_4.csv") %>%
    map(read_csv) %>%
    reduce(rbind)

tbl_top_sa_1 <- tbl_charge_1 %>%
    select(SCAN_NUMBER, SPECTRAL_ANGLE, MASS,
        ORIG_COLLISION_ENERGY, COLLISION_ENERGY) %>%
    group_by(SCAN_NUMBER) %>%
    top_n(1, SPECTRAL_ANGLE) %>%
    ungroup() %>%
    mutate(CE_diff = ORIG_COLLISION_ENERGY - COLLISION_ENERGY) %>%
    mutate(MySpecificBins =
        cut(SPECTRAL_ANGLE, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))

tbl_top_sa_1 %>% arrange(SPECTRAL_ANGLE)

# -----
tbl_top_sa_2 <- tbl_charge_2 %>%
    select(SCAN_NUMBER, SPECTRAL_ANGLE, MASS,
        ORIG_COLLISION_ENERGY, COLLISION_ENERGY) %>%
    group_by(SCAN_NUMBER) %>%
    top_n(1, SPECTRAL_ANGLE) %>%
    ungroup() %>%
    mutate(ORIG_COLLISION_ENERGY = as.double(ORIG_COLLISION_ENERGY)) %>%
    mutate(CE_diff = ORIG_COLLISION_ENERGY - COLLISION_ENERGY) %>%
    mutate(MySpecificBins =
        cut(SPECTRAL_ANGLE, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))

tbl_top_sa_2 %>% arrange(SPECTRAL_ANGLE)

# -----
tbl_top_sa_3 <- tbl_charge_3 %>%
    select(SCAN_NUMBER, SPECTRAL_ANGLE, MASS,
        ORIG_COLLISION_ENERGY, COLLISION_ENERGY) %>%
    group_by(SCAN_NUMBER) %>%
    top_n(1, SPECTRAL_ANGLE) %>%
    ungroup() %>%
    filter(COLLISION_ENERGY > 4) %>%
    mutate(SPECTRAL_ANGLE = as.double(SPECTRAL_ANGLE)) %>%
    mutate(CE_diff = ORIG_COLLISION_ENERGY - COLLISION_ENERGY) %>%
    mutate(MySpecificBins =
        cut(SPECTRAL_ANGLE, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))

tbl_top_sa_3 %>% arrange(SPECTRAL_ANGLE)
tbl_top_sa_3 %>%
    mutate(ORIG_COLLISION_ENERGY = as.double(SPECTRAL_ANGLE)) %>%
    filter(is.na(SPECTRAL_ANGLE))
tbl_top_sa_3 %>% arrange(SCAN_NUMBER)

# -----
tbl_top_sa_4 <- tbl_charge_4 %>%
    select(SCAN_NUMBER, SPECTRAL_ANGLE, MASS,
        ORIG_COLLISION_ENERGY, COLLISION_ENERGY) %>%
    group_by(SCAN_NUMBER) %>%
    top_n(1, SPECTRAL_ANGLE) %>%
    ungroup() %>%
    mutate(CE_diff = ORIG_COLLISION_ENERGY - COLLISION_ENERGY) %>%
    mutate(MySpecificBins =
        cut(SPECTRAL_ANGLE, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))

tbl_top_sa_4 %>% arrange(SPECTRAL_ANGLE)

tbl_top_sa_1$MySpecificBins <- factor(tbl_top_sa_1$MySpecificBins,
    levels = rev(levels(tbl_top_sa_1$MySpecificBins)))

tbl_top_sa_2$MySpecificBins <- factor(tbl_top_sa_2$MySpecificBins,
    levels = rev(levels(tbl_top_sa_2$MySpecificBins)))

tbl_top_sa_3$MySpecificBins <- factor(tbl_top_sa_3$MySpecificBins,
    levels = rev(levels(tbl_top_sa_3$MySpecificBins)))

tbl_top_sa_4$MySpecificBins <- factor(tbl_top_sa_4$MySpecificBins,
    levels = rev(levels(tbl_top_sa_4$MySpecificBins)))

plot_1 <- ggplot(tbl_top_sa_1,
    aes(x = MASS, y = CE_diff, fill = MySpecificBins)) +
    geom_point(aes(colour = MySpecificBins), alpha = 0.7) +
    theme_minimal() +
    xlab("Peptide mass") +
    ylab("original CE - top CE") +
    scale_colour_manual(values =
        c("#FDE624", "#5DC763", "#20908C", "#3B538B", "#430254"),
        name = "Spectral Angle")

plot_2 <- ggplot(tbl_top_sa_2,
    aes(x = MASS, y = CE_diff)) +
    geom_point(aes(colour = MySpecificBins), alpha = 0.7) +
    theme_minimal() +
    xlab("Peptide mass") +
    ylab("original CE - top CE") +
    scale_colour_manual(values =
        c("#FDE624", "#5DC763", "#20908C", "#3B538B", "#430254")
        name = "Spectral Angle")

plot_3 <- ggplot(tbl_top_sa_3,
    aes(x = MASS, y = CE_diff)) +
    geom_point(aes(colour = MySpecificBins), alpha = 0.7) +
    theme_minimal() +
    # theme(legend.position = "bottom") +
    xlab("Peptide mass") +
    ylab("original CE - top CE") +
    scale_colour_manual(values =
        c("#FDE624", "#5DC763", "#20908C", "#3B538B", "#430254"),
        name = "Spectral Angle")

plot_4 <- ggplot(tbl_top_sa_4,
    aes(x = MASS, y = CE_diff)) +
    geom_point(aes(colour = MySpecificBins), alpha = 0.7) +
    theme_minimal() +
    xlab("Peptide mass") +
    ylab("original CE - top CE") +
    scale_colour_manual(values =
        c("#FDE624", "#5DC763", "#20908C", "#3B538B", "#430254")
        name = "Spectral Angle")

leg <- get_legend(plot_3)
# Convert to a ggplot and print
legend_plot <- as_ggplot(leg)

full_plot <- ggarrange(plot_1, plot_2, plot_3, plot_4,
    ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom",
    labels = c("A", "B", "C", "D"))

plot_path <- "/Users/adams/Projects/300K/Results/Figures/CE/calibration_difference.png" # nolint
ggsave(plot_path, full_plot, width = 32, height = 22, units = "cm")

plot_path <- "/Users/adams/Projects/300K/Results/Figures/CE/calibration_difference_legend.png" # nolint
ggsave(plot_path, legend_plot, width = 32, height = 14, units = "cm")


# --------------
