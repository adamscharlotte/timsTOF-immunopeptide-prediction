# install.packages("ggtext")
# install.packages("scales")
library(tidyverse)
library(data.table)
library(ggtext)
library(scales)


hla_16_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/full-truncated-qc/un-annotated/TUM_HLA_16.csv"
tbl_hla_16 <- fread(hla_16_path) %>% as_tibble()
colnames(tbl_hla_16)

tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    select(COLLISION_ENERGY, MODIFIED_SEQUENCE, PRECURSOR, FRAME) %>%
    pull(COLLISION_ENERGY)


tbl_1565_mz <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1565") %>%
    # count(FRAME) %>% distinct()
    select(MZ) %>%
    mutate(MZ = strsplit(as.character(MZ), ";")) %>%
    unnest(MZ)

tbl_1565_i <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1565") %>%
    select(INTENSITIES) %>%
    mutate(INTENSITIES = strsplit(as.character(INTENSITIES), ";")) %>%
    unnest(INTENSITIES)

tbl_1565 <- cbind(tbl_1565_mz, tbl_1565_i) %>%
    as_tibble() %>%
    mutate(MZ = as.integer(MZ)) %>%
    mutate(INTENSITIES = as.double(INTENSITIES))

plot_1565 <- ggplot(tbl_1565) +
    geom_col(aes(x = MZ, y = INTENSITIES)) +
    # theme_bw() +
    # theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
    theme_minimal() +
    labs(x = "*m/Z*", y = "Intensity") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = ggtext::element_markdown()) +
    scale_x_continuous(limits = c(-0, 1110)) +
    geom_hline(yintercept = 0, size = 0.1) +
    geom_vline(xintercept = 0, size = 0.1)



tbl_1565_1803_mz <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1565") %>%
    filter(FRAME == "1803") %>%
    select(MZ) %>%
    mutate(MZ = strsplit(as.character(MZ), ";")) %>%
    unnest(MZ)

tbl_1565_1803_i <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1565") %>%
    filter(FRAME == "1803") %>%
    select(INTENSITIES) %>%
    mutate(INTENSITIES = strsplit(as.character(INTENSITIES), ";")) %>%
    unnest(INTENSITIES)

tbl_1565_1803 <- cbind(tbl_1565_1803_mz, tbl_1565_1803_i) %>%
    as_tibble() %>%
    mutate(MZ = as.integer(MZ)) %>%
    mutate(INTENSITIES = as.double(INTENSITIES))

plot_1565_1803 <- ggplot(tbl_1565_1803) +
    geom_col(aes(x = MZ, y = INTENSITIES)) +
    # theme_bw() +
    # theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
    theme_minimal() +
    labs(x = "*m/Z*", y = "Intensity") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = ggtext::element_markdown()) +
    scale_x_continuous(limits = c(-0, 1110)) +
    geom_hline(yintercept = 0, size = 0.1) +
    geom_vline(xintercept = 0, size = 0.1)

# -----------------------------------------------------------------------------

tbl_1565_1804_mz <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1565") %>%
    filter(FRAME == "1804") %>%
    select(MZ) %>%
    mutate(MZ = strsplit(as.character(MZ), ";")) %>%
    unnest(MZ)

tbl_1565_1804_i <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1565") %>%
    filter(FRAME == "1804") %>%
    select(INTENSITIES) %>%
    mutate(INTENSITIES = strsplit(as.character(INTENSITIES), ";")) %>%
    unnest(INTENSITIES)

tbl_1565_1804 <- cbind(tbl_1565_1804_mz, tbl_1565_1804_i) %>%
    as_tibble() %>%
    mutate(MZ = as.integer(MZ)) %>%
    mutate(INTENSITIES = as.double(INTENSITIES))

plot_1565_1804 <- ggplot(tbl_1565_1804) +
    geom_col(aes(x = MZ, y = INTENSITIES)) +
    # theme_bw() +
    # theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
    theme_minimal() +
    labs(x = "*m/Z*", y = "Intensity") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = ggtext::element_markdown()) +
    scale_x_continuous(limits = c(-0, 1110)) +
    geom_hline(yintercept = 0, size = 0.1) +
    geom_vline(xintercept = 0, size = 0.1)

# -----------------------------------------------------------------------------

tbl_1565_1805_mz <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1565") %>%
    filter(FRAME == "1805") %>%
    select(MZ) %>%
    mutate(MZ = strsplit(as.character(MZ), ";")) %>%
    unnest(MZ)

tbl_1565_1805_i <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1565") %>%
    filter(FRAME == "1805") %>%
    select(INTENSITIES) %>%
    mutate(INTENSITIES = strsplit(as.character(INTENSITIES), ";")) %>%
    unnest(INTENSITIES)

tbl_1565_1805 <- cbind(tbl_1565_1805_mz, tbl_1565_1805_i) %>%
    as_tibble() %>%
    mutate(MZ = as.integer(MZ)) %>%
    mutate(INTENSITIES = as.double(INTENSITIES))

plot_1565_1805 <- ggplot(tbl_1565_1805) +
    geom_col(aes(x = MZ, y = INTENSITIES)) +
    # theme_bw() +
    # theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
    theme_minimal() +
    labs(x = "*m/Z*", y = "Intensity") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = ggtext::element_markdown()) +
    scale_x_continuous(limits = c(-0, 1110)) +
    geom_hline(yintercept = 0, size = 0.1) +
    geom_vline(xintercept = 0, size = 0.1)

# -----------------------------------------------------------------------------

tbl_1602_1817_mz <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1602") %>%
    filter(FRAME == "1817") %>%
    select(MZ) %>%
    mutate(MZ = strsplit(as.character(MZ), ";")) %>%
    unnest(MZ)

tbl_1602_1817_i <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1602") %>%
    filter(FRAME == "1817") %>%
    select(INTENSITIES) %>%
    mutate(INTENSITIES = strsplit(as.character(INTENSITIES), ";")) %>%
    unnest(INTENSITIES)

tbl_1602_1817 <- cbind(tbl_1602_1817_mz, tbl_1602_1817_i) %>%
    as_tibble() %>%
    mutate(MZ = as.integer(MZ)) %>%
    mutate(INTENSITIES = as.double(INTENSITIES))

plot_1602_1817 <- ggplot(tbl_1602_1817) +
    geom_col(aes(x = MZ, y = INTENSITIES)) +
    # theme_bw() +
    # theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
    theme_minimal() +
    labs(x = "*m/Z*", y = "Intensity") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = ggtext::element_markdown()) +
    scale_x_continuous(limits = c(-0, 1110)) +
    geom_hline(yintercept = 0, size = 0.1) +
    geom_vline(xintercept = 0, size = 0.1)


path_plot <- "/Users/adams/Projects/300K/Results/Figures/Spectra/HLA_16_1565_1803.png" # nolint
ggsave(path_plot, plot_1565_1803, width = 10.5, height = 6, units = "cm")
path_plot <- "/Users/adams/Projects/300K/Results/Figures/Spectra/HLA_16_1565_1804.png" # nolint
ggsave(path_plot, plot_1565_1804, width = 10.5, height = 6, units = "cm")
path_plot <- "/Users/adams/Projects/300K/Results/Figures/Spectra/HLA_16_1565_1805.png" # nolint
ggsave(path_plot, plot_1565_1805, width = 10.5, height = 6, units = "cm")
path_plot <- "/Users/adams/Projects/300K/Results/Figures/Spectra/HLA_16_1602_1817.png" # nolint
ggsave(path_plot, plot_1602_1817, width = 10.5, height = 6, units = "cm")
path_plot <- "/Users/adams/Projects/300K/Results/Figures/Spectra/HLA_16_1565.png" # nolint
ggsave(path_plot, plot_1565, width = 10.5, height = 6, units = "cm")




# -----------------------------------------------------------------------------
spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

msms_path <- "/Users/adams/Projects/300K/2022-library-run/MaxQuant/20220624_HLAI_1_96/TUM_HLA_16-txt/msms.txt"
tbl_msms <- fread(msms_path) %>% as_tibble() %>% spaceless()

tbl_msms %>% filter(Scan_number == 2502) %>% select(Matches)
colnames(tbl_msms)

# -----------------------------------------------------------------------------











##############################################################################

tbl_1565_1803_mz <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1565") %>%
    filter(FRAME == "1803") %>%
    select(MZ) %>%
    mutate(MZ = strsplit(as.character(MZ), ";")) %>%
    unnest(MZ)

tbl_1565_1803_i <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1565") %>%
    filter(FRAME == "1803") %>%
    select(INTENSITIES) %>%
    mutate(INTENSITIES = strsplit(as.character(INTENSITIES), ";")) %>%
    unnest(INTENSITIES) %>%
    mutate(INTENSITIES = rescale(as.integer(INTENSITIES), to = c(0, 100)))

tbl_1565_1803 <- cbind(tbl_1565_1803_mz, tbl_1565_1803_i) %>%
    as_tibble() %>%
    mutate(MZ = as.integer(MZ))

max(tbl_1565_1803$INTENSITIES)

plot_1565_1803 <- ggplot(tbl_1565_1803) +
    geom_col(aes(x = MZ, y = INTENSITIES)) +
    # theme_bw() +
    # theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
    theme_minimal() +
    labs(x = "*m/Z*", y = "Intensity") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = ggtext::element_markdown()) +
    scale_x_continuous(limits = c(-0, 1550)) +
    geom_hline(yintercept = 0, size = 0.1) +
    geom_vline(xintercept = 0, size = 0.1)

# -----------------------------------------------------------------------------

tbl_1565_1804_mz <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1565") %>%
    filter(FRAME == "1804") %>%
    select(MZ) %>%
    mutate(MZ = strsplit(as.character(MZ), ";")) %>%
    unnest(MZ)

tbl_1565_1804_i <- tbl_hla_16 %>%
    filter(MODIFIED_SEQUENCE == "_KQGPTGIHV_") %>%
    filter(PRECURSOR == "1565") %>%
    filter(FRAME == "1804") %>%
    select(INTENSITIES) %>%
    mutate(INTENSITIES = strsplit(as.character(INTENSITIES), ";")) %>%
    unnest(INTENSITIES) %>%
    mutate(INTENSITIES = rescale(as.integer(INTENSITIES), to = c(0, 100)))

tbl_1565_1804 <- cbind(tbl_1565_1804_mz, tbl_1565_1804_i) %>%
    as_tibble() %>%
    mutate(MZ = as.integer(MZ))

plot_1565_1804 <- ggplot(tbl_1565_1804) +
    geom_col(aes(x = MZ, y = INTENSITIES)) +
    # theme_bw() +
    # theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
    theme_minimal() +
    labs(x = "*m/Z*", y = "Intensity") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = ggtext::element_markdown()) +
    scale_x_continuous(limits = c(-0, 1550)) +
    geom_hline(yintercept = 0, size = 0.1) +
    geom_vline(xintercept = 0, size = 0.1)