library(tidyverse)
library(data.table)
library(ComplexHeatmap)

tbl_fasta <- fread("/Users/adams/Downloads/Kuster/TUM_HLA_138.fasta", header=F) %>% as_tibble

tbl_target_sep_search <- fread("/Users/adams/Downloads/Kuster/crux-output-sep-search/comet.target.txt") %>% as_tibble
tbl_decoy_sep_search <- fread("/Users/adams/Downloads/Kuster/crux-output-sep-search/comet.decoy.txt") %>% as_tibble
tbl_target_conc_search <- fread("/Users/adams/Downloads/Kuster/crux-output-conc-search/comet.target.txt") %>% as_tibble
# tbl_target_sep_no_enz <- fread("/Users/adams/Downloads/Kuster/crux-output-sep-no-enz/comet.target.txt") %>% as_tibble
# tbl_decoy_sep_no_enz <- fread("/Users/adams/Downloads/Kuster/crux-output-sep-no-enz/comet.decoy.txt") %>% as_tibble

tbl_target_sep_search %>% count(`sp rank`)
tbl_decoy_sep_search %>% count(`sp rank`)
tbl_target_conc_search %>% count(`sp rank`)

# -------------------------------------------------------------------------------------------------------------------------
tbl_sep_confidence <- fread("/Users/adams/Downloads/Kuster/crux-output-sep-search/crux-output/assign-confidence.target.txt") %>% as_tibble
tbl_conc_confidence <- fread("/Users/adams/Downloads/Kuster/crux-output-conc-search/crux-output/assign-confidence.target.txt") %>% as_tibble

tbl_sep_confidence_filtered <- tbl_sep_confidence %>% filter(`mix-max q-value` < 0.0001) %>%
	separate(`protein id`, into=c("sp", "sequence_id", "experiment")) %>%
	filter(sequence_id==sequence)

tbl_conc_confidence_filtered <- tbl_conc_confidence %>% filter(`tdc q-value` < 0.001) %>%
	separate(`protein id`, into=c("sp", "sequence_id", "experiment")) %>%
	filter(sequence_id==sequence)


tbl_conc_confidence_filtered %>% count(charge)
	select(sequence, charge) %>% unique() %>% count(charge)

tbl_conc_confidence_filtered %>% filter(charge == 4) %>% select(sequence)

#	Load results -----------------------------------------------------------------------------------------------------------------
list_all_sequences <- tbl_fasta %>% filter(!str_detect(V1, ">sp")) %>% unique() %>% pull(V1)

mix_max_sequences <- tbl_sep_confidence_filtered %>% pull(sequence) %>% unique()
tdc_sequences <- tbl_conc_confidence_filtered %>% pull(sequence) %>% unique()

#	Create a binary matrix -----------------------------------------------------------------------------------------------------------
list_sequences <- list(
	"all possible sequences" = list_all_sequences,
	"mix-max 0.01%" = mix_max_sequences,
	"TDC 0.1%" = tdc_sequences
)

matrix_sequences <- list_to_matrix(list_sequences)

#	Create a combinarion matrix ------------------------------------------------------------------------------------------------------
matrix_sequences_combination <- make_comb_mat(matrix_sequences, top_n_sets = 3) # mode = "intersect")
set_size(matrix_sequences_combination)
set_name(matrix_sequences_combination)
matrix_sequences_combination <- matrix_sequences_combination[comb_degree(matrix_sequences_combination) > 0]

pdf(file = "/Users/adams/Downloads/Kuster/Figures/Upset/identified-sequences-mix-max-vs-tdc.pdf", height = 3.4, width = 4)
ss = set_size(matrix_sequences_combination)
cs = comb_size(matrix_sequences_combination)
ht = UpSet(matrix_sequences_combination,
	set_order = order(ss),
	comb_order = order(comb_degree(matrix_sequences_combination), -cs),
	top_annotation = HeatmapAnnotation(
		"PPI intersection size" = anno_barplot(cs,
			ylim = c(0, max(cs)*1.1),
			border = FALSE,
			gp = gpar(fill = "black"),
			height = unit(4, "cm")
		),
		annotation_name_gp= gpar(fontsize = 9),
		annotation_name_side = "left",
		annotation_name_rot = 90),
	left_annotation = rowAnnotation(
		"Set size" = anno_barplot(-ss,
			baseline = 0,
			axis_param = list(
				at = c(0, -1000, -2000),
				labels = c(0, 1000, 2000),
				labels_rot = 0),
			border = FALSE,
			gp = gpar(fill = "black"),
			width = unit(3, "cm")),
		set_name = anno_text(set_name(matrix_sequences_combination),
			location = 0.5,
			just = "center",
			width = max_text_width(set_name(matrix_sequences_combination)) + unit(-11, "mm"),
			gp= gpar(fontsize = 8)
			),
			annotation_name_gp= gpar(fontsize = 9)
			),
	right_annotation = NULL,
	show_row_names = FALSE)

ht = draw(ht)
od = column_order(ht)
decorate_annotation("PPI intersection size", {
	grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
		default.units = "native", just = c("center", "bottom"),
		gp = gpar(fontsize = 6, col = "#404040"))
})
dev.off()


# --------------------------------------------------------------------------------------------------------------------------------
conc_counts <- tbl_conc_confidence_filtered %>%
	select(sequence, scan) %>% unique() %>%
	add_count(sequence) %>% select(-scan) %>% unique() %>%
	count(n) %>% mutate(n=as.character(n))

f_01 <- ggplot(conc_counts, aes(n, nn))

plot_number_of_identifications <- f_01 + geom_col() +
	scale_x_discrete(name="# identifications",
	limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
	"12", "13", "14", "15", "16", "17", "18", "19", "20", "22", "24", "25", "26", "29", "30", "32", "52")) +
	geom_text(aes(label=nn), vjust=-0.25) +
	labs(y="# peptides") +
	theme_minimal() +
	theme(panel.grid.major=element_blank(), 
	panel.grid.minor=element_blank())

ggsave("/Users/adams/Downloads/Kuster/Figures/Upset/number-of-identifications-0.1-fdr.png", plot_number_of_identifications, width = 12, height = 9, units = "cm")

charge_counts <- tbl_conc_confidence_filtered %>%
	select(charge, sequence) %>% unique() %>%
	add_count(charge) %>% select(-sequence) %>% 
	mutate(charge=as.character(charge)) %>%
	# mutate(n=as.character(n)) %>%
	unique()

f_02 <- ggplot(charge_counts, aes(charge, n))

plot_number_of_charges <- f_02 + geom_col() +
	# scale_x_discrete(name="charge") +
	geom_text(aes(label=n), vjust=-0.25, color="white", fontface = "bold") +
	labs(y="# peptides", x ="charge") +
	theme_minimal() +
	theme(panel.grid.major=element_blank(), 
	panel.grid.minor=element_blank(),
	axis.title.x=element_text(color="white"),
	axis.title.y=element_text(color="white"))

ggsave("/Users/adams/Downloads/Kuster/Figures/Upset/number-of-charges-0.1-fdr-peptides.png", plot_number_of_charges, width = 5, height = 9, units = "cm")



tbl_conc_confidence_filtered %>% count(charge)
	select(sequence, charge) %>% unique() %>% count(charge)

tbl_conc_confidence_filtered %>% filter(sequence == "RLREHLVRF") %>% select(sequence, charge) %>% unique()

sep_counts <- tbl_conc_confidence_filtered %>%
	select(sequence, scan) %>% unique() %>%
	add_count(sequence) %>% select(-scan) %>% unique() %>%
	count(n) %>% mutate(n=as.character(n))

f_01 <- ggplot(sep_counts, aes(n, nn))

plot_number_of_identifications <- f_01 + geom_col() +
	scale_x_discrete(name="# identifications",
	limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
	"12", "13", "14", "15", "16", "17", "18", "19", "20", "22", "24", "25", "26", "29", "30", "32", "52")) +
	geom_text(aes(label=nn), vjust=-0.25) +
	labs(y="# peptides") +
	theme_minimal()

ggsave("/Users/adams/Downloads/Kuster/Figures/Upset/number-of-edentifications-0.1-fdr.png", plot_number_of_identifications, width = 12, height = 9, units = "cm")

# --------------------------------------------------------------------------------------------------------------------------------


tbl_peprec <- tbl_conc_confidence_filtered %>%
	mutate(spec = "controllerType=0 controllerNumber=1 scan=") %>%
	select(spec, scan, sequence, charge, `sp score`) %>% unique() %>% 
	unite(spec_id, spec:scan, sep="") %>%
	unite(seq_charge, sequence:charge, sep="_", remove=F) %>%
	group_by(seq_charge) %>%
	slice_max(order_by = `sp score`, n = 1) %>%
	ungroup() %>%
	select(-c(seq_charge, `sp score`))


fwrite(tbl_peprec, file="/Users/adams/Downloads/Kuster/peprec.tsv", sep='\t')
"controllerType=0 controllerNumber=1 scan="