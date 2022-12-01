library(tidyverse)
library(data.table)
# library(ggplot2)
# # install.packages("plotly")
# library(plotly)

ce_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/QC-collission-energies.csv"
meta_qc_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/qc-peptides.txt"
full_meta_map_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/full-meta-map.txt"

tbl_qc_peptides <- fread(meta_qc_path) %>% as_tibble()
tbl_meta <- fread(full_meta_map_path) %>%
    as_tibble() %>%
    select(folder_names, plate) %>%
    mutate(RAW_FILE = str_sub(folder_names, 1, nchar(folder_names) - 2))

tbl_ce <- fread(ce_path) %>% as_tibble

tbl_ce %>%
    select(MODIFIED_SEQUENCE, RAW_FILE) %>%
    distinct() %>%
    count(MODIFIED_SEQUENCE) %>%
    arrange(desc(n))

tbl_ce %>%
    select(RAW_FILE) %>%
    distinct()

qc_present_in_most_files <- "_SYASDFGSSAK_"

filtered_ce <- tbl_ce %>%
    # select(RAW_FILE) %>%
    mutate(time = str_sub(RAW_FILE, start = -4)) %>%
    filter(MODIFIED_SEQUENCE == "_HDTVFGSYLYK_") %>%
    filter(CHARGE == 2) %>%
    select(RAW_FILE, time, MODIFIED_SEQUENCE, SCORE,
    COLLISION_ENERGY) %>%
    distinct() %>%
    arrange(time) %>%
    merge(tbl_meta) %>%
    as_tibble() %>%
    distinct()

mean(filtered_ce$SCORE)

fig <- plot_ly(filtered_ce, x = ~time, y = ~COLLISION_ENERGY, color = ~plate, type = "box")
fig <- fig %>% layout(boxmode = "group")

fig

tbl_mean_ce <- tbl_ce %>%
    # select(RAW_FILE) %>%
    mutate(time = str_sub(RAW_FILE, start = -4)) %>%
    filter(MODIFIED_SEQUENCE == "_HDTVFGSYLYK_") %>%
    filter(CHARGE == 2) %>%
    select(RAW_FILE, time, MODIFIED_SEQUENCE, SCORE,
    COLLISION_ENERGY) %>%
    arrange(time) %>%
    merge(tbl_meta) %>%
    as_tibble() %>%
    distinct() %>%
    group_by(RAW_FILE) %>%
    mutate(mean_ce = mean(COLLISION_ENERGY)) %>%
    ungroup() %>%
    select(mean_ce, time, plate) %>%
    distinct()

mean(filtered_ce$SCORE)

fig <- plot_ly(tbl_mean_ce, x = ~time, y = ~mean_ce, color = ~plate,
    type = "scatter")
fig <- fig %>% layout(boxmode = "group")

fig

# -----------------------------------------------------------------------------


sampled_pools <- tbl_ce_filtered %>%
    select(pool_name, plate) %>%
    group_by(plate) %>%
    sample_n(2) %>%
    ungroup()

tbl_ce_sampled <- tbl_ce_filtered %>%
    filter(pool_name %in% sampled_pools$pool_name) %>%
    group_by(pool_name, Sequence) %>%
    mutate(mean_ce = mean(collision_energy)) %>%
    ungroup() %>%
    select(Sequence, pool_name, mean_ce) %>%
    distinct() %>%
    group_by(Sequence) %>%
    mutate(pool_count = n_distinct(pool_name)) %>%
    ungroup() %>%
    filter(pool_count == 10)

ggplot(data = tbl_ce_sampled, aes(x=Sequence, y=mean_ce, group = pool_name, colour = as.factor(pool_name))) +
    geom_line() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


fig <- plot_ly(tbl_ce_sampled, x = ~Sequence, y = ~collision_energy, color = ~pool_name, type = "lines")

fig <- fig %>% layout(boxmode = "group")
fig <- plot_ly(data, x = ~x, y = ~random_y, type = 'scatter', mode = 'lines')

fig <- plot_ly(tbl_ce_filtered, x = ~Sequence)
fig <- fig %>% add_trace(y = ~trace_0, name = 'trace 0',mode = 'lines')
fig <- fig %>% add_trace(y = ~trace_1, name = 'trace 1', mode = 'lines+markers')
fig <- fig %>% add_trace(y = ~trace_2, name = 'trace 2', mode = 'markers')

fig


# -----------------------------------------------------------------------------

p <- ggplot(tbl_ce_filtered, aes(x=Sequence_Charge, y=collision_energy)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

tbl_ce %>% filter(Sequence == "GDFTFFIDTFK") %>% select(Score) %>% max()

ggsave("/Users/adams/Projects/300K/2022-library-run/Figures/QC-CE-boxplot.pdf", p)