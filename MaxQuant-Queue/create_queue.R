library(tidyverse)
library(data.table)

base <- "Z:\\\\internal_projects\\\\active\\\\ProteomeTools\\\\ProteomeTools\\\\" # nolint
path <- "External_data\\\\Bruker\\\\UA-TimsTOF-300K\\\\"
local_path <- "/Users/adams/Projects/300K/Tryptic/"
# local_path <- "/Users/adams/Projects/300K/2022-library-run/"
full_meta_map_path <- paste(local_path, "Metadata/full-meta-map.txt", sep = "")
queue_path <- paste(local_path, "MaxQuant/Queue/queue.csv", sep = "")
min_max_length <- 11
tbl_full_meta_map <- fread(full_meta_map_path) %>% as_tibble

tbl_queue <- tbl_full_meta_map %>%
    mutate(pool = pool_name, MQ_version = "2.1.2.0",
    pre_payload = "generateMQpar", post_payload = "cleanPostSearch",
    threads = 4,
    fasta_file = paste(base, path, "TUM-fasta\\\\",
        pool_name, ".fasta", sep = ""),
    raw_folder = paste(base, "Analysis\\\\", pool_name,
        "_01_01\\\\TIMS-30min-R1-TIMS_semitryptic\\\\230719_f1-",
        pool_name, "_01_01-TIMS-30min-R1.d", sep = ""),
    base_mqpar = "mqpar-tryptic.xml",
    output_folder = paste(base, path, "TUM-Searches\\\\",
        pool_name, sep = "")) %>%
    mutate(max_length_qc = ifelse(max_length < min_max_length,
        min_max_length, max_length)) %>%
    select(pool, MQ_version, pre_payload, post_payload, threads, plate,
        fasta_file, raw_folder, max_length_qc, base_mqpar, output_folder)
fwrite(tbl_queue, queue_path)