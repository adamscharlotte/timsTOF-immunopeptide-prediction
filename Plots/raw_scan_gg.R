# install.packages("ggtext")
# install.packages("scales")
library(tidyverse)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

pool <- args[1]
precursor <- args[2]

# pool <- "TUM_HLA2_31"
# precursor <- "12063"

base_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/"
sum_path <- paste(base_path, "precursor-consensus/summed/", pool, ".csv", sep = "") # nolint
frame_path <- paste(base_path, "full-truncated-qc/un-annotated/", pool, ".csv", sep = "") # nolint

# -------------------------- Load the summed spectrum -------------------------
tbl_sum <- fread(sum_path) %>% as_tibble()

tbl_sum_mz <- tbl_sum %>%
    filter(PRECURSOR == precursor) %>%
    select(MZ) %>%
    mutate(MZ = str_sub(MZ, 2, -2)) %>%
    mutate(MZ = strsplit(as.character(MZ), ", ")) %>%
    unnest(MZ)

tbl_sum_i <- tbl_sum %>%
    filter(PRECURSOR == precursor) %>%
    select(INTENSITIES) %>%
    mutate(INTENSITIES = str_sub(INTENSITIES, 2, -2)) %>%
    mutate(INTENSITIES = strsplit(as.character(INTENSITIES), ", ")) %>%
    unnest(INTENSITIES) %>%
    mutate(INTENSITIES = as.double(INTENSITIES))

tbl_sum_normi <- tbl_sum_i %>%
    mutate(INTENSITIES = (INTENSITIES) / (max(tbl_sum_i$INTENSITIES)))

tbl_sum_mzi <- cbind(tbl_sum_mz, tbl_sum_normi) %>%
    as_tibble() %>%
    mutate(MZ = as.integer(MZ)) %>%
    mutate(INTENSITIES = as.double(INTENSITIES))

# ------------------------ Load spectrum for each frame -----------------------

tbl_frame <- fread(frame_path) %>% as_tibble()

list_frame <- tbl_frame %>%
    filter(PRECURSOR == precursor) %>%
    pull(FRAME)

for (i in list_frame) {
    frame <- i
    tbl_frame_mz <- tbl_frame %>%
        filter(PRECURSOR == precursor) %>%
        filter(FRAME == frame) %>%
        select(MZ) %>%
        mutate(MZ = strsplit(as.character(MZ), ";")) %>%
        unnest(MZ)
    tbl_frame_i <- tbl_frame %>%
        filter(PRECURSOR == precursor) %>%
        filter(FRAME == frame) %>%
        select(INTENSITIES) %>%
        mutate(INTENSITIES = strsplit(as.character(INTENSITIES), ";")) %>%
        unnest(INTENSITIES) %>%
        mutate(INTENSITIES = as.double(INTENSITIES))
    tbl_frame_normi <- tbl_frame_i %>%
        mutate(INTENSITIES = (INTENSITIES) / (max(tbl_frame_i$INTENSITIES)))
    tbl_frame_mzi <- cbind(tbl_frame_mz, tbl_frame_normi) %>%
        as_tibble() %>%
        mutate(MZ = as.integer(MZ)) %>%
        mutate(INTENSITIES = as.double(INTENSITIES))
    tbl_sum_frame_mzi <- tbl_sum_mzi %>%
        rename(sum_INTENSITIES = INTENSITIES) %>%
        left_join(tbl_frame_mzi) %>%
        mutate(INTENSITIES = replace_na(INTENSITIES, 0))
    plot_frame_sum <- ggplot(tbl_sum_frame_mzi) +
        geom_col(aes(x = MZ, y = INTENSITIES, width = 1)) +
        geom_col(aes(x = MZ, y = -sum_INTENSITIES, width = 1)) +
        # theme_bw() +
        # theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme_minimal() +
        labs(x = "*m/Z*", y = "Intensity") +
        theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = ggtext::element_markdown()) +
        scale_y_continuous(limits = c(-60, 60)) +
        # scale_x_continuous(limits = c(-0, 1110)) +
        geom_hline(yintercept = 0, size = 0.1) +
        geom_vline(xintercept = 0, size = 0.1)
    result_path <- "/Users/adams/Projects/300K/Results/Figures/Mirror/"
    plot_path <- paste(result_path, pool, "_", precursor, "_", frame, ".pdf", sep = "")
    ggsave(plot_path, plot_frame_sum, width = 22, height = 14, units = "cm")
}
