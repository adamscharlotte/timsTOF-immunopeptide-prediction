library(tidyverse)
library(data.table)
library(ggplot2)
# library(plotly)

rt_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/QC-retention-time.csv"
meta_qc_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/qc-peptides.txt"
full_meta_map_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/full-meta-map.txt"

tbl_qc_peptides <- fread(meta_qc_path) %>%
    as_tibble() %>%
    rename(full_sequence = Sequence)
tbl_meta <- fread(full_meta_map_path) %>%
    as_tibble() %>%
    select(pool_name, folder_names)

tbl_rt <- fread(rt_path) %>% as_tibble

tbl_rt_filtered <- tbl_rt %>%
    filter(Sequence %in% tbl_qc_peptides$full_sequence) %>%
    # unite(Sequence_Charge, Sequence:Charge, remove = FALSE) %>%
    # mutate(Charge = as.character(Charge))
    # filter(Charge == 2) %>%
    as_tibble() %>%
    distinct()

tbl_rt %>% count(Proteins) %>% distinct()
tbl_rt_filtered %>% count(Proteins) %>% distinct()

tbl_qc_map <- merge(tbl_rt, tbl_qc_peptides) %>%
    as_tibble() %>%
    filter(str_detect(full_sequence, Sequence))

tbl_qc_map %>%
    filter(!Sequence %in% tbl_qc_peptides$full_sequence) %>%
    filter(Proteins == "QC_JPT_RT_Peptide") %>%
    count(Sequence) %>%
    arrange(desc(n))

tbl_qc_map %>%
    filter(Sequence == "GDFTFFID") %>%
    select(Sequence, full_sequence)

tbl_dot_plot <- tbl_qc_map %>%
    filter(full_sequence == "GDFTFFIDTFK") %>%
    select(pool_name, Sequence, everything()) %>%
    unite("pool_seq", pool_name:Sequence, remove = FALSE) %>%
    group_by(pool_seq) %>%
    top_n(1, Score) %>%
    ungroup()

tbl_dot_plot_ann <- merge(tbl_dot_plot, tbl_meta) %>%
    as_tibble %>%
    arrange(folder_names)

plot_scatter <- ggplot(tbl_dot_plot_ann, aes(x = folder_names, y = retention_time, color = Sequence)) +
    geom_point() +
    scale_color_manual(values = c(
        "#0E1C36", "#878C9B", "#79B473", "#014D02", "#A2BCE0", "#EA2B37"
    )) +
    theme_minimal() +
    ylab("Retention time") +
    theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank())

path_plot <- "/Users/adams/Projects/300K/Results/Figures/RT scatter.png" # nolint
ggsave(path_plot, plot_scatter, width = 27.5, height = 8, units = "cm")

# -----------------------------------------------------------------------------

# tbl_ce_filtered %>% count(Proteins) %>% distinct()
# tbl_ce_filtered %>% select(Charge) %>% unique()
# tbl_ce_filtered %>% select(Modifications) %>% unique()

# fig <- plot_ly(tbl_ce_filtered, x = ~Sequence, y = ~collision_energy, color = ~Charge, type = "box")
# fig <- fig %>% layout(boxmode = "group")
# fig


fig <- plot_ly(tbl_ce_filtered, x = ~Sequence, y = ~collision_energy, color = ~plate, type = "box")
fig <- fig %>% layout(boxmode = "group")

fig
