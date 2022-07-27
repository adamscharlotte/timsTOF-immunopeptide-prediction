library(tidyverse)
library(data.table)

base <- "Z:\\\\internal_projects\\\\active\\\\ProteomeTools\\\\"
path <- "ProteomeTools\\\\External_data\\\\Bruker\\\\UA-TimsTOF-300K\\\\"
local_path <- "/Users/adams/Projects/300K/2022-library-run/"
full_meta_map_path <- paste(local_path, "metadata/full-meta-map.txt", sep = "")
queue_path <- paste(local_path, "MaxQuant/Queue/queue.csv", sep = "")
min_max_length <- 17
tbl_full_meta_map <- fread(full_meta_map_path) %>% as_tibble

tbl_queue <- tbl_full_meta_map %>%
    mutate(pool = pool_name, MQ_version = "2.1.2.0",
    pre_payload = "generateMQpar", post_payload = "cleanPostSearch",
    threads = 1,
    fasta_file = paste(base, path, "fasta\\\\", pool_name, ".fasta", sep = ""),
    raw_folder = paste(base, path, "raw-library-run\\\\",
    plate, "\\\\", folder_names, sep = ""),
    base_mqpar = "mqpar-unspecific.xml",
    output_folder = paste(base, path, "Searches\\\\", plate, "\\\\",
    pool_name, sep = "")) %>%
    mutate(max_length_qc = ifelse(max_length < min_max_length,
    min_max_length, max_length)) %>%
    select(pool, MQ_version, pre_payload, post_payload, threads, plate,
    fasta_file, raw_folder, max_length_qc, base_mqpar, output_folder)
fwrite(tbl_queue, queue_path)

min(tbl_queue$max_length_qc)