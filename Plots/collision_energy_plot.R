library(tidyverse)
library(data.table)
library(ggplot2)
# install.packages("plotly")
library(plotly)

ce_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/QC-collission-energy.csv"
meta_qc_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/qc-peptides.txt"
full_meta_map_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/full-meta-map.txt"

tbl_qc_peptides <- fread(meta_qc_path) %>% as_tibble()
tbl_meta <- fread(full_meta_map_path) %>%
    as_tibble() %>%
    select(pool_name, plate)

tbl_ce <- fread(ce_path) %>% as_tibble
tbl_ce_filtered <- tbl_ce %>%
    filter(Sequence %in% tbl_qc_peptides$Sequence) %>%
    # unite(Sequence_Charge, Sequence:Charge, remove = FALSE) %>%
    # mutate(Charge = as.character(Charge))
    # filter(Charge == 2) %>%
    merge(tbl_meta) %>%
    as_tibble() %>%
    distinct()

# tbl_ce_filtered %>% count(Proteins) %>% distinct()
# tbl_ce_filtered %>% select(Charge) %>% unique()
# tbl_ce_filtered %>% select(Modifications) %>% unique()

# fig <- plot_ly(tbl_ce_filtered, x = ~Sequence, y = ~collision_energy, color = ~Charge, type = "box")
# fig <- fig %>% layout(boxmode = "group")
# fig


fig <- plot_ly(tbl_ce_filtered, x = ~Sequence, y = ~collision_energy, color = ~plate, type = "box")
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