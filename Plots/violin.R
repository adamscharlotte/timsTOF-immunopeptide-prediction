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

tbl_plot <- rbind(tbl_file_1, tbl_file_2)

ggplot(tbl_plot,
    aes(group, spectral_angle, fill = label)) +
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
    scale_fill_manual(labels = c("HCD Bruker Prosit 2023",
            expression(paste("HCD Prosit 2020 (Wilhelm ",
                italic("et al."), ")"))),
                values = c("#0e1c36", "#cdeac0")) +
    ggtitle("2023 HCD Bruker Prosit Accuracy") +
    labs(y = "Spectral angle\nPrediction vs synthesized peptide spectrum",
        x = "")

plot_path <- "/Users/adams/Projects/300K/Results/Figures/model-performance/violin/group-violin-plot.png" # nolint
ggsave(plot_path, group_violin_plot)