library(tidyverse)
library(data.table)

base <- "Z:\\\\internal_projects\\\\active\\\\ProteomeTools\\\\ProteomeTools\\\\" # nolint
# path <- "External_data\\\\Bruker\\\\UA-TimsTOF-300K\\\\"
path <- "External_data\\\\Bruker\\\\PXD038782-comparison\\\\"
local_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/" # nolint
# d_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/d-folder" # nolint
d_path <- "/media/kusterlab/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/PXD038782-comparison/raw" # nolint
# queue_path <- paste(local_path, "Queue/Annotation/queue-raw-ii.csv", sep = "")
queue_path <- paste(local_path, "Queue/Annotation/queue-raw.csv", sep = "")
# queue_path <- paste(local_path, "Queue/Annotation/queue-hla-i-d.csv", sep = "")
# local_path <- "/Users/adams/Projects/300K/Tryptic/"
# local_path <- "/Users/adams/Projects/300K/2022-library-run/"
# full_meta_map_path <- paste(local_path, "Metadata/full-meta-map.txt", sep = "")
# queue_path <- paste(local_path, "MaxQuant/Queue/queue-unsp.csv", sep = "")
# min_max_length <- 11
# tbl_full_meta_map <- fread(full_meta_map_path) %>% as_tibble

# hla_i_d <- list.files(d_path, pattern = "W6-32_17%_DDA_Rep")
list_raw <- list.files(d_path, pattern = ".raw")

# tbl_queue <- tbl_full_meta_map %>%
# tbl_queue <- tibble(folder_name = hla_i_d) %>%
tbl_queue <- tibble(folder_name = list_raw) %>%
    filter(str_detect(folder_name, "W6-32")) %>%
    # filter(str_detect(folder_name, "Tue39")) %>%
    # mutate(file_name = substr(folder_name, 1, nchar(folder_name) - 2),
    mutate(file_name = substr(folder_name, 1, nchar(folder_name) - 4),
        MQ_version = "2.1.2.0", pre_payload = "generateMQpar",
        post_payload = "cleanPostSearch", threads = 4,
        # plate = "hla-i",
        plate = "raw-hla-i",
        # plate = "raw-hla-ii",
        max_length_qc = 16,
        # max_length_qc = 30,
        fasta_file = paste(base, path,
            "fasta\\\\uniprot_sprot-230323.fasta", sep = ""),
        # raw_folder = paste(base, path, "d-folder\\\\", folder_name, sep = ""),
        raw_folder = paste(base, path, "raw\\\\", folder_name, sep = ""),
        base_mqpar = "mqpar-raw.xml",
        output_folder = paste(base, path, "Searches\\\\", file_name, sep = "")
    ) %>%
    # mutate(pool = pool_name, MQ_version = "2.1.2.0",
    #     pre_payload = "generateMQpar", post_payload = "cleanPostSearch",
    #     threads = 4,
    #     fasta_file = paste(base, path, "TUM-fasta\\\\",
    #     pool_name, ".fasta", sep = ""),
        # raw_folder = paste(base, "Analysis\\\\", pool_name,
        #     "_01_01\\\\TIMS-30min-R1-TIMS_semitryptic\\\\230719_f1-",
        #     pool_name, "_01_01-TIMS-30min-R1.d", sep = ""),
        # base_mqpar = "mqpar-unspecific.xml",
        # output_folder = paste(base, path, "TUM-Searches\\\\first_pool\\\\",
        #     pool_name, sep = "")) %>%
    # mutate(max_length_qc = ifelse(max_length < min_max_length,
    #     min_max_length, max_length)) %>%
    select(file_name, MQ_version, pre_payload, post_payload, threads, plate,
        fasta_file, raw_folder, max_length_qc, base_mqpar, output_folder)
    # select(pool, MQ_version, pre_payload, post_payload, threads, plate,
    #     fasta_file, raw_folder, max_length_qc, base_mqpar, output_folder)
fwrite(tbl_queue, queue_path)