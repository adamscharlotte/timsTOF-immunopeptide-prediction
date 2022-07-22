library(tidyverse)
library(data.table)

example_queue_path <- "/Users/adams/Library/CloudStorage/OneDrive-UniversiteitAntwerpen/Code/300K/queue/queue_tbd.csv"
full_meta_map_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/full-meta-map.txt"
lysn_queue_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/queue/LysN-queue.csv"
lysn_annotation_path <- "/Users/adams/Projects/300K/2022-library-run/metadata/queue/annotation/LysN-RawFilesAnnotation.csv"

tbl_full_meta_map <- fread(full_meta_map_path) %>% as_tibble
tbl_example_queue <- fread(example_queue_path) %>% as_tibble

tbl_lysn_queue <- tbl_full_meta_map %>% filter(str_detect(pool_name, "lysn")) %>% 
	mutate(pool = pool_name, MQ_version = "2.1.2.0", pre_payload = "generateMQpar", 
	post_payload = "cleanPostSearch", threads = 1, base = plate, 
	fasta_file = paste("/cygdrive/z/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/fasta/", pool_name, ".fasta", sep=""), 
	raw_folder = paste("/cygdrive/z/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/raw-library-run/", plate, "/", folder_names, sep="")) %>%
	select(pool, MQ_version, pre_payload, post_payload, threads, base, fasta_file, raw_folder, max_length)
fwrite(tbl_lysn_queue, lysn_queue_path)

tbl_lysn_annotation <- tbl_full_meta_map %>% filter(str_detect(pool_name, "lysn")) %>% 
	mutate(File = str_remove(folder_names, ".d"), Output_folder = paste("/cygdrive/z/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/raw-library-run/Searches", plate, "/", pool_name, sep=""), 
	Experiment = "", MS_Plate = plate, Fraction = "", Parameter_group = "", PTM = "", Reference_channel = "", 
	Compound_order = "", Project = "", Comments = "") %>%
	mutate(raw_folder = paste("/cygdrive/z/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/raw-library-run/", plate, "/", File, sep=""),
	fasta_file = paste("/cygdrive/z/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/fasta/", pool_name, ".fasta", sep="")) %>%
	separate(pool_name, into=c("tum", "enzyme", "File_order"), sep="_") %>%
	select(raw_folder,fasta_file,max_length,Output_folder,MS_Plate)
	# select(File,Output_folder,Experiment,MS_Plate,Fraction,Parameter_group,PTM,Reference_channel,File_order,Compound_order,Project,Comments)
	
fwrite(tbl_lysn_annotation, lysn_annotation_path)

1. `File`: Raw file name without full path (the path is specified in `queue.csv`) and without `.raw` extension
2. `Output_folder`: Should correspond to `_name` field of the corresponding row in `queue.csv`. This (partly) specifies the MQ processing folder.
3. `Experiment`: `Experiment` field in `mqpar.xml`. Fractionated samples should have the same experiment name
4. `MS_Plate`: Should be the same for the RAW files searched together in one search, e.g. can typically be set to a fixed string for an entire project or, to be more safe, taken as the same as `Output_folder`. This column can be used to reuse the same RAW file in multiple searches, e.g. use the same DMSO for different `Output_folder`s.
5. `Fraction`: `Fraction` field in `mqpar.xml`
6. `Parameter_group`: `Parameter group` field in `mqpar.xml`
7. `PTM`: `PTM` field in `mqpar.xml`. If set to `True` this run will not be used for protein quantification.
8. `Reference_channel`: `Reference channel` field in `mqpar.xml`. Used to specify the bridge channel for isobaric match between runs. Leave empty if you don't want to use the isobaric match between runs normalization.
9. `File_order`: Not used
10. `Compound_order`: Not used
11. `Project`: Should correspond to `_project` field of the corresponding row in `queue.csv`. This (partly) specifies the MQ processing folder and acts as a filter for which Raw files to add to the `mqpar.xml`.
12. `Comments`: Comments



tbl_full_meta_map %>% select(min_length) %>% count(min_length) %>% arrange(min_length) %>% unique


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# tbl_lysn_queue <- tbl_full_meta_map %>% filter(str_detect(pool_name, "lysn")) %>% 
# 	mutate(pool = pool_name, MQ_version = "2.1.2.0", pre_payload = "cpAndGenerateMQpar", 
# 	post_payload = "cleanPostSearch", threads = 1, project = "UA-TimsTOF-300K", 
# 	fasta_file = paste("/cygdrive/z/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/fasta/", pool_name, ".fasta", sep=""), 
# 	raw_folder = paste("/cygdrive/z/internal_projects/active/ProteomeTools/ProteomeTools/External_data/Bruker/UA-TimsTOF-300K/raw-library-run/", plate, "/", folder_names, sep="")) %>%
# 	select(pool, MQ_version, pre_payload, post_payload, threads, project, fasta_file, raw_folder, max_length)
# fwrite(tbl_lysn_queue, lysn_queue_path)
