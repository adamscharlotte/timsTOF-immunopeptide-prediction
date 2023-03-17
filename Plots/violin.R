library(tidyverse)
library(data.table)
library(ggplot2)

# -----------------------------------------------------------------------------

geom_split <- ggproto("GeomSplitViolin", GeomViolin,
    draw_group = function(self, data, ..., draw_quantiles = NULL) {
        data <- transform(data, xminv = x - violinwidth * (x - xmin),
            xmaxv = x + violinwidth * (xmax - x))
        grp <- data[1, "group"]
        new_df <- plyr::arrange(transform(data,
            x = if (grp %% 2 == 1) xminv else xmaxv),
            if (grp %% 2 == 1) y else -y)
        new_df <- rbind(new_df[1, ], new_df, new_df[nrow(new_df), ],
            new_df[1, ])
    new_df[c(1, nrow(new_df) - 1, nrow(new_df)), "x"] <- round(new_df[1, "x"])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
        quantiles <- ggplot2:::create_quantile_segment_frame(data,
            draw_quantiles)
        aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data),
            c("x", "y")), drop = FALSE]
        aesthetics$alpha <- rep(1, nrow(quantiles))
        both <- cbind(quantiles, aesthetics)
        quantile_grob <- GeomPath$draw_panel(both, ...)
        ggplot2:::ggname("geom_split_violin",
            grid::grobTree(GeomPolygon$draw_panel(new_df, ...), quantile_grob))
    } else {
        ggplot2:::ggname("geom_split_violin",
            GeomPolygon$draw_panel(new_df, ...))
    }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
    position = "identity", ..., draw_quantiles = NULL, trim = TRUE,
    scale = "area", na_rm = FALSE, show_legend = NA, inherit_aes = TRUE) {
        layer(data = data, mapping = mapping, stat = stat, geom = geom_split,
        position = position, show.legend = show_legend,
        inherit.aes = inherit_aes, params = list(trim = trim,
            scale = scale, draw_quantiles = draw_quantiles, na.rm = na_rm, ...))
}

# -----------------------------------------------------------------------------

base_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/"  # nolint
file_1_path <- paste(base_path, "20230315104417_prediction.csv", sep = "")
file_2_path <- paste(base_path, "20230315111004_prediction.csv", sep = "")
file_3_path <- paste(base_path, "20230316105504_prediction.csv", sep = "")

test_non_path <- paste(base, "Annotation/total-scan-consensus/split/non-tryptic/annotated-40-ppm-test.csv", sep = "") # nolint
train_non_path <- paste(base, "Annotation/total-scan-consensus/split/non-tryptic/annotated-40-ppm-train.csv", sep = "") # nolint
test_tryp_path <- paste(base, "Annotation/total-scan-consensus/split/tryptic/annotated-40-ppm-test.csv", sep = "") # nolint
train_tryp_path <- paste(base, "Annotation/total-scan-consensus/split/tryptic/annotated-40-ppm-train.csv", sep = "") # nolint

test_non <- fread(test_non_path) %>% as_tibble()
train_non <- fread(train_non_path) %>% as_tibble()
test_tryp <- fread(test_tryp_path) %>% as_tibble()
train_tryp <- fread(train_tryp_path) %>% as_tibble()

tbl_test_filtered <- rbind(test_non, test_tryp)
tbl_train_filtered <- rbind(train_non, train_tryp)

# 2020 VS 2023 Model
tbl_file_1 <- fread(file_1_path) %>%
    as_tibble() %>%
    mutate(label = "HCD Prosit 2020") %>%
    mutate(group = substr(rawfile, 1, 7))

tbl_file_2 <- fread(file_2_path) %>%
    as_tibble() %>%
    mutate(label = "HCD Bruker Prosit 2023") %>%
    mutate(group = substr(rawfile, 1, 7))

tbl_prediction <- rbind(tbl_file_1, tbl_file_2)
tbl_test_filtered <- tbl_test_filtered %>% mutate(scan_number = SCAN_NUMBER)
df_prediction_map <- left_join(tbl_prediction, tbl_test_filtered, by = "scan_number")

tbl_prediction_map <- df_prediction_map %>%
    as_tibble() %>%
    filter(str_detect(rawfile, RAW_FILE)) %>%
    distinct() %>%
    mutate(type = substr(pool_name, 5, 8)) %>%
    mutate(new_type = case_when(
        type == "HLA_" ~ "MHC-I",
        type == "HLA2" ~ "MHC-II",
        type == "lysn" ~ "LysN",
        type == "aspn" ~ "AspN",
        type == "firs" ~ "Tryptic",
        TRUE ~ NA_character_
    )) %>%
    mutate(charge = as.character(PRECURSOR_CHARGE)) %>%
    mutate(n_term = substr(OBS_SEQUENCE, 1, 1)) %>%
    mutate(c_term = str_sub(OBS_SEQUENCE, start = -1)) %>%
    mutate(length = as.character(nchar(OBS_SEQUENCE)))

# Test VS Train
tbl_file_2 <- fread(file_2_path) %>%
    as_tibble() %>%
    mutate(label = "Test set") %>%
    mutate(group = substr(rawfile, 1, 7))
tbl_test_filtered <- tbl_test_filtered %>% mutate(scan_number = SCAN_NUMBER)
tbl_test <- left_join(tbl_file_2, tbl_test_filtered, by = "scan_number")

tbl_file_3 <- fread(file_3_path) %>%
    as_tibble() %>%
    mutate(label = "Training set") %>%
    mutate(group = substr(rawfile, 1, 7))
tbl_train_filtered <- tbl_train_filtered %>% mutate(scan_number = SCAN_NUMBER)
tbl_train <- left_join(tbl_file_3, tbl_train_filtered, by = "scan_number")

tbl_test_train <- rbind(tbl_test, tbl_train)

tbl_test_train_map <- tbl_test_train %>%
    as_tibble() %>%
    filter(str_detect(rawfile, RAW_FILE)) %>%
    distinct() %>%
    mutate(type = substr(pool_name, 5, 8)) %>%
    mutate(new_type = case_when(
        type == "HLA_" ~ "MHC-I",
        type == "HLA2" ~ "MHC-II",
        type == "lysn" ~ "LysN",
        type == "aspn" ~ "AspN",
        type == "firs" ~ "Tryptic",
        TRUE ~ NA_character_
    )) %>%
    mutate(charge = as.character(PRECURSOR_CHARGE)) %>%
    mutate(n_term = substr(OBS_SEQUENCE, 1, 1)) %>%
    mutate(c_term = str_sub(OBS_SEQUENCE, start = -1)) %>%
    mutate(length = as.character(nchar(OBS_SEQUENCE)))


group_violin_plot <- ggplot(tbl_test_train_map,
    aes(length, spectral_angle, fill = label)) +
    geom_split_violin(colour = "white", size = 0) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.2),
        legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
    scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
        0.9, 1.0)) +
    scale_x_discrete(limits = c("7", "8", "9", "10", "11", "12", "13", "14",
        "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26",
        "27", "28", "29", "30", "31")) +
    # scale_x_discrete(limits = c("MHC-I", "MHC-II", "LysN",
    #     "AspN", "Tryptic")) +
    # scale_x_discrete(limits = c("1", "2", "3")) +
    scale_fill_manual(values = c("#7A8DB3", "#022873")) +
    # scale_fill_manual(labels = c("HCD TOF Prosit 2023",
    #         expression(paste("HCD Prosit 2020 (Wilhelm ",
    #             italic("et al."), ")"))),
                # values = c("#0e1c36", "#cdeac0")) +
    ggtitle("2023 HCD TOF Prosit Accuracy") +
    labs(y = "Spectral angle",
        # x = "")
        x = "Peptide length")
        # x = "N-terminus")
        # x = "Precursor charge")

group_violin_plot

# plot_path <- "/Users/adams/Projects/300K/Results/Figures/model-performance/violin/test-vs-train/group-violin-plot.png" # nolint
plot_path <- "/Users/adams/Projects/300K/Results/Figures/model-performance/violin/predicted-vs-synthesized/lentgh-violin-plot.png" # nolint
ggsave(plot_path, group_violin_plot, width = 18, height = 8, units = "cm")
