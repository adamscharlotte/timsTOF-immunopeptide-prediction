library(tidyverse)
library(data.table)
library(ggplot2)
# library(plotly)

rt_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/QC-retention-time.csv"
meta_qc_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/qc-peptides.txt"
# full_meta_map_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/full-meta-map.txt"

tbl_qc_peptides <- fread(meta_qc_path) %>% as_tibble()
# tbl_meta <- fread(full_meta_map_path) %>%
#     as_tibble() %>%
#     select(pool_name, plate)

tbl_rt <- fread(rt_path) %>% as_tibble
tbl_rt_filtered <- tbl_rt %>%
    filter(Sequence %in% tbl_qc_peptides$Sequence) %>%
    # unite(Sequence_Charge, Sequence:Charge, remove = FALSE) %>%
    # mutate(Charge = as.character(Charge))
    # filter(Charge == 2) %>%
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

mapped_precursor_path <- "/Users/adams/Projects/300K/2022-library-run/Annotation/precursor-mapped/TUM_aspn_1.csv"

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

tbl_mapped_precursor <- fread(mapped_precursor_path) %>% as_tibble() %>% spaceless()


# -----------------------------------------------------------------------------
# Info.txt
txt_output_path <- paste(base_path, "Metadata/info.txt", sep = "")

tbl_qc_psms <- tbl_mapped_precursor %>% filter(stringr::str_starts(Proteins, "QC"))
meta_qc_path <- paste(base_path, "Metadata/qc-peptides.txt", sep = "")

tbl_meta <- fread(meta_path) %>% as_tibble
tbl_pool <- tbl_meta %>%
    filter(pool_name == pool) %>%
    mutate(trunc_sequence = str_sub(Sequence, start = -7)) %>%
    rename(Sequence_pool = Sequence) %>%
    distinct()
tbl_qc_peptides <- fread(meta_qc_path) %>% as_tibble()

tbl_filtered_qc_frames <- tbl_mapped_precursor %>%
    filter(stringr::str_starts(Proteins, "QC")) %>%
    filter(Sequence %in% tbl_qc_peptides$Sequence) %>%
    select(Scan_number, frame, Precursor, Sequence) %>%
    distinct()

tbl_pool_psms <- tbl_mapped_precursor %>%
    mutate(pool_name = pool) %>%
    filter(Proteins == pool) %>%
    merge(tbl_pool, by = "pool_name") %>%
    as_tibble()

tbl_filtered_pool_frames <- tbl_pool_psms %>%
    filter(stringr::str_ends(Sequence, trunc_sequence)) %>%
    select(Scan_number, frame, Precursor, Sequence) %>%
    distinct()

tbl_psms <- tbl_mapped_precursor %>% select(Scan_number, Sequence) %>% distinct()
tbl_pool_psms <- tbl_pool_psms %>% select(Scan_number, Sequence) %>% distinct()
tbl_qc_psms <- tbl_qc_psms %>% select(Scan_number, Sequence) %>% distinct()
tbl_filtered_pool_psms <- tbl_filtered_pool_frames %>% select(Scan_number, Sequence) %>% distinct()
tbl_filtered_qc_psms <- tbl_filtered_qc_frames %>% select(Scan_number, Sequence) %>% distinct()
tbl_unique_pool_peptides <- tbl_filtered_pool_frames %>% select(Sequence) %>% distinct()
tbl_full_length <- tbl_filtered_pool_frames %>% filter(Sequence %in% tbl_pool$Sequence_pool) %>%
    select(Sequence) %>% distinct
tbl_filtered_qc_precursors <- tbl_filtered_qc_frames %>%
    select(Precursor, frame, Sequence) %>%
    distinct() %>%
    add_count(Precursor, frame) %>%
    filter(n==1)

tbl_filtered_pool_precursors <- tbl_filtered_pool_frames %>%
    select(Precursor, frame, Sequence) %>%
    distinct() %>%
    add_count(Precursor, frame) %>%
    filter(n==1)

tbl_filtered_precursors <- rbind(tbl_filtered_qc_precursors, tbl_filtered_pool_precursors)

# -----------------------------------------------------------------------------

tbl_annotation <- tbl_mapped_precursor %>%
    # add_count(Precursor, frame) %>%
    # filter(n == 1) %>%    # Remove precursors selected twice for scan summing
    # select(-c(Precursor, frame, n)) %>%
    # select(-c(Precursor, frame)) %>%
    # filter(across(c(Precursor, frame) ~.x %in% tbl_filtered_precursors[, c("Precursor", "frame")])) %>%
    filter(Sequence %in% tbl_filtered_precursors$Sequence &
        Precursor %in% tbl_filtered_precursors$Precursor) %>%
    # filter(Precursor %in% tbl_filtered_qc_psms$Precursor |
    #     Precursor %in% tbl_filtered_precursors$Precursor) %>%
    rename_with(toupper) %>%
    select(RAW_FILE, SCAN_NUMBER, MODIFIED_SEQUENCE, CHARGE, FRAGMENTATION,
    MASS_ANALYZER, SCAN_EVENT_NUMBER, MASS, SCORE, REVERSE, RETENTION_TIME,
    MZ, INTENSITIES, COLLISION_ENERGY, PRECURSOR, FRAME)


output_path <- paste(base_path, "Annotation/full-truncated-qc/un-annotated/", pool, ".csv", sep = "")
fwrite(tbl_annotation, output_path)

# -----------------------------------------------------------------------------

tbl_info <- tibble(pool_name = pool, PSMs = nrow(tbl_psms),
    max_length = max(tbl_pool_psms$Length),
    pool_PSMs = nrow(tbl_pool_psms), qc_PSMs = nrow(tbl_qc_psms),
    filtered_pool_psms = nrow(tbl_filtered_pool_psms), 
    filtered_qc_psms = nrow(tbl_filtered_qc_psms), 
    pool_unique_peptides = nrow(tbl_unique_pool_peptides),
    full_length = nrow(tbl_full_length),
    # truncated = nrow(tbl_truncated),
    precursor = nrow(tbl_filtered_precursors),
    frame = nrow(tbl_filtered_pool_frames),
    # wrong_identification = nrow(tbl_wrong),
    fasta_size = nrow(tbl_pool),
    percentage = full_length / fasta_size * 100,
    min_score = min(tbl_mapped_precursor$Score),
    max_score = max(tbl_mapped_precursor$Charge),
    total_annotation = nrow(tbl_annotation)
)

fwrite(tbl_info, txt_output_path, append = TRUE)

# -----------------------------------------------------------------------------

