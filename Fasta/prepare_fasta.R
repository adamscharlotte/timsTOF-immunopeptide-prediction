library(tidyverse)
library(data.table)

tbl_AspN <- fread("/Users/adams/Projects/300K/20220613-library-test/metadata/AspN.txt") %>% as_tibble
tbl_AspN_poolname <- tbl_AspN %>% select(`Pool name`) %>% unique()
fwrite(tbl_AspN_poolname, "/Users/adams/Projects/300K/20220613-library-test/metadata/poolnames/AspN_names.txt", col.names=FALSE)
