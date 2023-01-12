library(tidyverse)
library(data.table)

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

path <- "/Users/adams/Projects/300K/2022-library-run/accumulatedMsmsScansTryptic" # nolint
setwd(path)
files <- dir(pattern = "*.txt")
scans_txt <- files %>%
    map(read_tsv) %>%       # read in all the files individually, using
                            # the function read_csv() from the readr package
    reduce(rbind)			# reduce with rbind into one dataframe

tbl_barplot <- scans_txt %>%
    spaceless() %>%
    select(Raw_file, Identified, Scan_number) %>%
    distinct() %>%
    add_count(Raw_file) %>%
    rename(total = n) %>%
    group_by(Raw_file) %>%
    add_count(Identified) %>%
    ungroup() %>%
    select(Raw_file, Identified, n, total) %>%
    mutate(Identified = replace_na(Identified, "nul")) %>%
    distinct()

tbl_barplot$Identified <- factor(tbl_barplot$Identified, levels = c("+", "nul"))

plot_bar <- ggplot(data = tbl_barplot,
    aes(x = Raw_file, y = n, fill = Identified)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(name = "", labels = c("Identified", "Not identified"),
    values = c("#EA2B37", "#0E1C36")) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # legend.position = "top",
        # legend.position = c(0.3, 0.8)
        legend.position = "none"
        ) +
    xlab("") +
    ylab("#Spectra")

path_plot <- "/Users/adams/Projects/300K/Results/Figures/Identification rate unspecific enzyme run proteotypic peptides.png" # nolint
ggsave(path_plot, plot_bar, width = 21.5, height = 6, units = "cm")

tbl_identified <- tbl_barplot %>% filter(Identified == "+")

sum(tbl_identified$total)
sum(tbl_identified$n)

# tryptic-unsp
381237 / 3155208
# tryptic tryptic
344354 / 2890141
# non-tryptic
379709 / 2175627