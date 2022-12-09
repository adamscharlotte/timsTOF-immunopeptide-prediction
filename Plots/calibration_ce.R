library(tidyverse)
library(data.table)
library(ggplot2)
# install.packages("ggpubr")
library(ggpubr)
# install.packages("patchwork")
# library(patchwork)

base_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/"
sa_path <- paste(base_path, "scan-consensus/spectral-angle/", sep = "")
setwd(sa_path)

tbl_charge_1 <-
    list.files(pattern = "*_1.csv") %>%
    map_df(~read_csv(.))

tbl_charge_2 <-
    list.files(pattern = "*_2.csv") %>%
    map_df(~read_csv(.))

tbl_charge_3 <-
    list.files(pattern = "*_3.csv") %>%
    map_df(~read_csv(.))

tbl_charge_4 <-
    list.files(pattern = "*_4.csv") %>%
    map_df(~read_csv(.))

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

ggplot(tbl_top_sa_1, aes(ORIG_COLLISION_ENERGY, COLLISION_ENERGY)) +
    geom_point()

ggplot(tbl_top_sa_1, aes(MASS, CE_diff)) +
    geom_point()

# -----
tbl_top_sa_2 <- tbl_charge_2 %>%
    select(SCAN_NUMBER, SPECTRAL_ANGLE, MASS,
        ORIG_COLLISION_ENERGY, COLLISION_ENERGY) %>%
    group_by(SCAN_NUMBER) %>%
    top_n(1, SPECTRAL_ANGLE) %>%
    ungroup() %>%
    mutate(CE_diff = ORIG_COLLISION_ENERGY - COLLISION_ENERGY) %>%
    mutate(MySpecificBins =
        cut(SPECTRAL_ANGLE, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))

tbl_top_sa_2 %>% arrange(SPECTRAL_ANGLE)

ggplot(tbl_top_sa_2, aes(ORIG_COLLISION_ENERGY, COLLISION_ENERGY)) +
    geom_point()

ggplot(tbl_top_sa_2, aes(MASS, CE_diff)) +
    geom_point()

# -----
tbl_top_sa_3 <- tbl_charge_3 %>%
    select(SCAN_NUMBER, SPECTRAL_ANGLE, MASS,
        ORIG_COLLISION_ENERGY, COLLISION_ENERGY) %>%
    group_by(SCAN_NUMBER) %>%
    top_n(1, SPECTRAL_ANGLE) %>%
    ungroup() %>%
    mutate(CE_diff = ORIG_COLLISION_ENERGY - COLLISION_ENERGY) %>%
    mutate(MySpecificBins =
        cut(SPECTRAL_ANGLE, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))

tbl_top_sa_3 %>% arrange(SPECTRAL_ANGLE)

ggplot(tbl_top_sa_3, aes(ORIG_COLLISION_ENERGY, COLLISION_ENERGY)) +
    geom_point()



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

ggplot(tbl_top_sa_4, aes(ORIG_COLLISION_ENERGY, COLLISION_ENERGY)) +
    geom_point()

plot_1 <- ggplot(tbl_top_sa_1,
    aes(x = MASS, y = CE_diff)) +
    geom_point(aes(colour = MySpecificBins), alpha = 0.8) +
    scale_colour_manual(values =
        c("#430254", "#3B538B", "#20908C", "#5DC763", "#FDE624"))

plot_2 <- ggplot(tbl_top_sa_2,
    aes(x = MASS, y = CE_diff)) +
    geom_point(aes(colour = MySpecificBins), alpha = 0.8) +
    scale_colour_manual(values =
        c("#430254", "#3B538B", "#20908C", "#5DC763", "#FDE624"))

plot_3 <- ggplot(tbl_top_sa_3,
    aes(x = MASS, y = CE_diff)) +
    geom_point(aes(colour = MySpecificBins), alpha = 0.8) +
    scale_colour_manual(values =
        c("#430254", "#3B538B", "#20908C", "#5DC763", "#FDE624"))

plot_4 <- ggplot(tbl_top_sa_4, aes(x = MASS, y = CE_diff)) +
    geom_point(aes(colour = MySpecificBins), alpha = 0.8) +
    scale_colour_manual(values =
        c("#430254", "#3B538B", "#20908C", "#5DC763", "#FDE624"))

ggarrange(plot_1, plot_4, plot_2, plot_3, ncol = 2, nrow = 2,
    common.legend = TRUE, legend = "bottom",
    labels = c("A", "B", "C", "D"))


# plot_3 + plot_4 + plot_layout(guides = "collect")