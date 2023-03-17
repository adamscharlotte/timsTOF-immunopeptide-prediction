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

tbl_file_1 <- fread(file_1_path) %>%
    as_tibble() %>%
    mutate(label = "HCD Prosit 2020") %>%
    mutate(group = substr(rawfile, 1, 7))

tbl_file_2 <- fread(file_2_path) %>%
    as_tibble() %>%
    mutate(label = "HCD Bruker Prosit 2023") %>%
    mutate(group = substr(rawfile, 1, 7))

tbl_prediction <- rbind(tbl_file_1, tbl_file_2)

# base_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/total-scan-consensus/prediction/"  # nolint
# file_1_path <- paste(base_path, "20230315104417_prediction.csv", sep = "")
# file_2_path <- paste(base_path, "20230315111004_prediction.csv", sep = "")

tbl_test_filtered <- tbl_test_filtered %>% mutate(scan_number = SCAN_NUMBER)
df_prediction_map <- left_join(tbl_prediction, tbl_test_filtered, by = "scan_number")

tbl_prediction_map <- df_prediction_map %>%
    as_tibble() %>%
    filter(str_detect(rawfile, RAW_FILE)) %>%
    distinct() %>%
    mutate(type = substr(pool_name, 5, 8)) %>%
    mutate(threshold = type) %>%
    mutate(threshold = case_when(
        threshold == "HLA_" ~ 0.5,
        threshold == "HLA2" ~ 0.6,
        threshold == "lysn" ~ 0.8,
        threshold == "aspn" ~ 0.4,
        threshold == "firs" ~ 0.5,
        TRUE ~ as.numeric(as.character(threshold))
    )) %>%
    mutate(new_type = case_when(
        type == "HLA_" ~ "Class I",
        type == "HLA2" ~ "Class II",
        type == "lysn" ~ "LysN",
        type == "aspn" ~ "AspN",
        type == "firs" ~ "Tryptic",
        TRUE ~ NA_character_
    ))

tbl_prediction_map %>% select(threshold)
colnames(tbl_prediction_map)
tbl_prediction_map %>%
    count(type)

level_order <- c("HLA_", "HLA2", "lysn", "aspn", "firs")
        (labels=c("HLA_" = "Class I", "HLA2" = "Class II", "lysn" = "LysN",
            "aspn" = "AspN", "firs" = "Tryptic"))
        (labels = c("Class I", "Class II", "LysN", "AspN", "Tryptic"))

group_violin_plot <- ggplot(tbl_prediction_map,
    aes(new_type, spectral_angle, fill = label)) +
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
    scale_x_discrete(limits = c("Class I", "Class II", "LysN",
        "AspN", "Tryptic")) +
    scale_fill_manual(labels = c("HCD Bruker Prosit 2023",
            expression(paste("HCD Prosit 2020 (Wilhelm ",
                italic("et al."), ")"))),
                values = c("#0e1c36", "#cdeac0")) +
    ggtitle("2023 HCD Bruker Prosit Accuracy") +
    labs(y = "Spectral angle\nPrediction vs synthesized peptide spectrum",
        x = "")

plot_path <- "/Users/adams/Projects/300K/Results/Figures/model-performance/violin/group-violin-plot.png" # nolint
ggsave(plot_path, group_violin_plot)
threshold <- c(0.5,0.6,0.8,0.7,0.6)

group_violin_plot +
    stat_summary(aes(y = 0.25, group = new_type), fun = median, geom = "errorbar",
            width = 0.5, color = "red")


    stat_summary(aes(y = threshold, group = new_type), fun = median, geom = "hline", color = "red", alpha = 0.5, yintercept = 1)

    stat_summary(fun = median, geom = "errorbar", width = 0.5, color = "black", size = 1, alpha = 0.7, 
                 position = position_dodge(width=0.75)) +
    stat_summary(aes(y = threshold), fun = median, geom = "errorbar", width = 0.5, col = "black", alpha = 0.7,
                 position = position_dodge(width=0.75), yintercept = "middle")

    stat_summary(fun.data = "median_hilow", geom = "errorbar", width = 0.5, col = "black", alpha = 0.7)

    stat_summary(fun = median, geom = "hline", aes(yintercept = 0.5))
    stat_summary(fun.data = "median_hilow", geom = "errorbar", width = 0.5, col = "red", alpha = 0.7, fun.args = list(y = c(0.3, 0.5, 0.2, 0.4, 0.6)))


    geom_hline(aes(yintercept = med_threshold), data = tbl_prediction_map, color = "red", size = 1.5)


    facet_wrap(~new_type, scales = "free_x", ncol = 5) +
    geom_hline(aes(yintercept = threshold), data = tbl_prediction_map %>% distinct(new_type, threshold), color = "red", size = 1.5)

    geom_hline(aes(yintercept = threshold), data = tbl_prediction_map %>% distinct(group, threshold), color = "red", size = 1.5)

    geom_segment(aes(x = new_type, xend = new_type, y = thresholds, yend = thresholds), color = "red", size = 1.5)


        stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar",
               width = 0.25,
               position = position_dodge(width = .25),
               colour = "#f33b16")







# create sample data
set.seed(123)
df <- data.frame(
  new_type = rep(c("a I", "a II", "a", "b", "c"), each = 20),
  spectral_angle = runif(100, 0, 1)
)

# calculate median threshold per new_type group
med_threshold <- aggregate(spectral_angle ~ new_type, data = df, FUN = median)

# plot split violin with median threshold line for each new_type group
ggplot(df, aes(new_type, spectral_angle, fill = new_type)) +
  geom_split_violin(colour = "white", size = 0) +
  geom_hline(aes(yintercept = spectral_angle), data = med_threshold, color = "red", size = 1) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.2),
        legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_discrete(limits = c("a I", "a II", "a", "b", "c")) +
  labs(fill = "new_type") +
  ggtitle("Split Violin Plot with Median Threshold Line for Each Group")

ggplot(tbl_prediction_map,
       aes(new_type, spectral_angle, fill = label)) +
  geom_split_violin(colour = "white", size = 0) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.2),
        legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  scale_x_discrete(limits = c("a I", "a II", "a", "b", "c")) +
  stat_summary(aes(y = 0.25, group = new_type), fun = median, geom = "errorbar", 
               width = 0, color = "black") +
  stat_summary(aes(y = 0.45, group = new_type), fun = median, geom = "errorbar", 
               width = 0, color = "black") +
  stat_summary(aes(y = 0.15, group = new_type), fun = median, geom = "errorbar", 
               width = 0, color = "black") +
  stat_summary(aes(y = 0.35, group = new_type), fun = median, geom = "errorbar", 
               width = 0, color = "black") +
  stat_summary(aes(y = 0.55, group = new_type), fun = median, geom = "errorbar", 
               width = 0, color = "black")