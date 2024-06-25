###### script info: script for plotting frequency histogram of sequences length data ######
source("scripts/plotting_scripts/1_utils/01.1_libs.R")
source("scripts/plotting_scripts/1_utils/01.2_own_funct_softw.R")
source("scripts/plotting_scripts/1_utils/01.3_graphics.R")

args = commandArgs(trailingOnly=TRUE)
load(args[1])
output_dir <- args[2]
sample_order <- EvoTraceR_object$sample_order



###### input data (from EvoTraceR package) ######
df_to_plot_perf_match <- 
  EvoTraceR_object$statistics$all_asv_statistics %>% 
  mutate(sample=fct_relevel(sample, sample_order)) # adjust columns names and order


###### data outputs ######
## output dir: for graphs analysis
graphs_analysis_dir <- output_dir
if (!dir.exists(graphs_analysis_dir)) 
{dir.create(graphs_analysis_dir, recursive = TRUE)}


### samples order ###
sample_columns <- setdiff(colnames(EvoTraceR_object$asv_prefilter), c("seq_names", "seq"))


### adjust data ###
data_for_hist <- 
  df_to_plot_perf_match %>% 
  dplyr::filter(stringr::str_detect(asv_names, "ASV")) # remove non-marked barcode (e.g., BC10v0, depending on the experiment)


###### plot: parameters ######
## automatic scaling of y axis
## maximum value for y
perc_in_sample_max <- 
  data_for_hist %>% 
  dplyr::select(seq_n, sample, perc_in_sample) %>% 
  group_by(seq_n, sample) %>% 
  mutate(perc_in_sample_max = sum(perc_in_sample)) %>% 
  pull(perc_in_sample_max) %>% 
  max()
## y axis ticks: <10% -> 2% | > 10% <= 25% -> 5% | > 25% <= 50% -> 10% | > 50% -> 20%
y_axis_ticks <- ifelse(perc_in_sample_max <= 10, 2, 
                       ifelse(perc_in_sample_max > 10 & perc_in_sample_max <= 25, 5, 
                              ifelse(perc_in_sample_max > 25 & perc_in_sample_max <= 50, 10, 
                                     ifelse(perc_in_sample_max > 50, 20))))


###### plot: plotting histogram data of sequences length, hist_freq_seq_length.pdf ######
hist_seq_count <- 
	ggplot(data=data_for_hist, aes(y=perc_in_sample, x=seq_n)) + 
	geom_bar(stat="identity", position = "stack", fill="#B484A9", width=1) +
	scale_x_continuous(labels=scales::comma,
	                   limits=c(0, nchar(EvoTraceR_object$reference$ref_seq)*2), 
	                   breaks=c(1, seq(floor(nchar(EvoTraceR_object$reference$ref_seq)*2/10), nchar(EvoTraceR_object$reference$ref_seq)*2,  floor(nchar(EvoTraceR_object$reference$ref_seq)*2/10))), 
						         expand = c(0.01, 0.01)) +
  scale_y_continuous(labels=function(x) paste0(x, "%"),
                     limits=c(0, plyr::round_any(perc_in_sample_max, y_axis_ticks, f=ceiling)), # automatic
                     breaks=seq(from=0, to=plyr::round_any(perc_in_sample_max, y_axis_ticks, f=ceiling), by=y_axis_ticks),
                     expand = c(0.01, 0.01)) +
	xlab("ASV Length") +
	ylab("ASV Frequency") +
	geom_vline(xintercept=nchar(EvoTraceR_object$reference$ref_seq), linetype="dashed", linewidth=0.25, col="#84B48F") + # expected size
	lemon::coord_capped_cart(left="both", bottom="left") + # axis with lemon
	lemon::facet_rep_grid(rows = vars(sample), repeat.tick.labels = TRUE) + 
	# add theme
	barplot_nowaklab_theme() +
	theme(aspect.ratio = 1/4,
		plot.margin = unit(c(2, 0, 0, 0), "mm"),
		legend.position = "None")
## save pdf
ggsave(plot=hist_seq_count,
       filename=file.path(graphs_analysis_dir, "/hist_freq_seq_length.pdf"),
       width=20, height=5*length(sample_columns), units = "cm")


###### csv: data_for_hist.csv ######
write.csv(data_for_hist,  
			file.path(graphs_analysis_dir, "/hist_freq_seq_length.csv"),
			row.names = FALSE, quote = FALSE)
  