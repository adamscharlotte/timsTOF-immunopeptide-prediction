library(tidyverse)
library(data.table)

base_path <- "/Users/adams/Projects/300K/2022-library-run/"
meta_path <- paste(base_path, "Metadata/full-meta-map.txt", sep = "")
sa_path <- paste(base_path, "Annotation/spectral-angle/frames", sep = "")

setwd(sa_path)
files <- dir(pattern = "*.csv")
df_sa <- files %>%
    map(read_csv) %>%        # read in all the files individually, using
                            # the function read_csv() from the readr package
    reduce(rbind)           # reduce with rbind into one dataframe

tbl_meta <- fread(meta_path) %>%
    as_tibble() %>%
    select(pool_name, folder_names) %>%
    mutate(folder_names = substring(folder_names, 1, nchar(folder_names) - 2))

tbl_sa <- merge(df_sa, tbl_meta, by.x = "RAW_FILE", by.y = "folder_names") %>%
    as_tibble()
tbl_sa %>% group_by(SEQUENCE) %>% do(data.frame(t(combn(.$PRECURSOR, 2))))

fruits <- tibble(
  type   = c("apple", "orange", "apple", "orange", "orange", "orange"),
  year   = c(2010, 2010, 2012, 2010, 2011, 2012),
  size  =  factor(
    c("XS", "S",  "M", "S", "S", "M"),
    levels = c("XS", "S", "M", "L")
  ),
  weights = rnorm(6, as.numeric(size) + 2)
)

install.packages("tictoc")
library(tictoc)
fruits %>% group_by(type) %>% expand(year, full_seq(year, 1)) %>% ungroup() #%>% colnames()

tic("group fast")
tbl_sa %>%
    group_by(SEQUENCE) %>%
    summarise(PRECURSOR = combn(PRECURSOR, 2, simplify = FALSE)) %>%
    unnest_wider(PRECURSOR, names_sep = "_")
toc()

tic("group faster")
tbl_sa %>% group_by(SEQUENCE) %>% do(data.frame(t(combn(.$PRECURSOR, 2))))
toc()

tic("group even faster")
setDT(tbl_sa)[, {i1 <-  combn(PRECURSOR, 2)
                   list(i1[1,], i1[2,]) }, by =  SEQUENCE]
toc()

tic("group even even faster")
setDT(tbl_sa)[, transpose(combn(PRECURSOR, 2, FUN = list)), by = SEQUENCE]
toc()

tic("group fastest")
lst <- by(tbl_sa$PRECURSOR, tbl_sa$SEQUENCE, FUN = combn, m= 2)
data.frame(SEQUENCE = rep(unique(as.character(tbl_sa$SEQUENCE)),
                    sapply(lst, ncol)), t(do.call(cbind, lst)))
toc()


sample <- data.frame(
  group=c("a","a","a","a","b","b","b"),
  number=c(1,2,3,2,4,5,3)
) %>% as_tibble()

tic("group fast")
sample %>%
    group_by(group) %>%
    summarise(number = combn(number, 2, simplify = FALSE)) %>%
    unnest_wider(number, names_sep = "_")
toc()

tic("group faster")
sample %>% group_by(group) %>% do(data.frame(t(combn(.$number, 2))))
toc()

tic("group even faster")
setDT(sample)[, {i1 <-  combn(number, 2)
                   list(i1[1,], i1[2,]) }, by =  group]
toc()

tic("group even even faster")
setDT(sample)[, transpose(combn(number, 2, FUN = list)), by = group]
toc()

tic("group fastest")
lst <- by(sample$number, sample$group, FUN = combn, m= 2)
data.frame(group = rep(unique(as.character(sample$group)), 
                    sapply(lst, ncol)), t(do.call(cbind, lst)))
toc()



    select(pool_name, PRECURSOR, )

colnames(tbl_sa)
query_path <- "/Users/adams/Projects/300K/Results/Figures/Mirror/"
