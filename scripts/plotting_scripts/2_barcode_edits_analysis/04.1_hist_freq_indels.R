###### script info: plotting frequency histogram of deletions and insertions on the barcode scale ######
## Plots show peaks centered around edit sites (green dashed line). Edits will be deletions (blue) and insertions (red).
## Created Files: hist_freq_indels.pdf 
source("~/snakemake_evotracer_machina/scripts/plotting_scripts/1_utils/01.1_libs.R")
source("~/snakemake_evotracer_machina/scripts/1_utils/01.2_own_funct_softw.R")
source("~/snakemake_evotracer_machina/scripts/1_utils/01.3_graphics.R")

args = commandArgs(trailingOnly=TRUE)
load(args[1])
output_dir <- args[2]
sample_order <- EvoTraceR_object$sample_order

###### input data (from EvoTraceR package) ######
del_sub_ins_df <- 
  EvoTraceR_object$alignment$mutations_df %>% 
  mutate(sample=fct_relevel(sample, sample_order)) # adjust columns names and order


###### data outputs ######
# output dir: for graphs analysis
graphs_analysis_dir <- paste0(output_dir, "/graphs_analysis")
if (!dir.exists(graphs_analysis_dir)) 
{dir.create(graphs_analysis_dir, recursive = TRUE)}

###### adjusts input data: "del_sub_ins_df" ######
## summarize stat for "del" and "sub" -> position is not stacked but added as one; i.e. pos 10 & freq: 12%, 10%, will be pos: 10 freq: 22%
del_sub_df_data_to_plot_sum_perc <-
  del_sub_ins_df %>% 
  ungroup() %>%
  group_by(sample, alt, position_bc260) %>% #
  dplyr::summarise(sum_perc = sum(perc_in_sample), .groups = "drop") %>%
  dplyr::filter(alt != "w" & alt != "i") # don't plot "wt" and "ins"
## "ins" only -> position is not stacked but added as one; i.e. pos 10 & freq: 12%, 10%, will be pos: 10 freq: 22%
ins_df_data_to_plot_sum_perc <-
  del_sub_ins_df %>%
  dplyr::select(asv_names, sample, position_bc260, alt, perc_in_sample) %>%
  dplyr::filter(alt == "i") %>% # get only "ins"
  unique() %>%
  group_by(sample, alt, position_bc260) %>%
  dplyr::summarise(sum_perc = sum(perc_in_sample), .groups = "drop")
## bind "ins" after recalculation with "del_sub"
del_sub_ins_df_data_to_plot_sum_perc <- rbind(del_sub_df_data_to_plot_sum_perc, ins_df_data_to_plot_sum_perc)
# add full names for labeling of plots
del_sub_ins_df_data_to_plot_sum_perc <-
  del_sub_ins_df_data_to_plot_sum_perc %>%
  mutate(alt_long_names = ifelse(alt == "i", "Insertions", 
                                 ifelse(alt == "d", "Deletions", 
                                        ifelse(alt == "s", "Substitutions", "No Edits"))))


###### auxiliary data for plotting histogram graph  ######
# position of PAM in guides
pam_pos <- EvoTraceR_object$reference$ref_cut_sites
# length of barcode
bc_len <- nchar(EvoTraceR_object$reference$ref_seq)
# annotating rectangles for target sites (gray color) = 26 bp -> 20x bp (target site) + 3x bp (PAM) + 3x bp (spacer)
del_sub_ins_df_data_to_plot_sum_perc <- 
  del_sub_ins_df_data_to_plot_sum_perc %>% 
  dplyr::filter(alt_long_names != "No Edits")
# add rectangles for target sites (gray/white) 
alt_count_annot_rect <-
  ggplot(data = del_sub_ins_df_data_to_plot_sum_perc, aes(x=position_bc260, y=sum_perc, fill=alt_long_names, group=sample))
annot_rect <- 
  EvoTraceR_object$reference$ref_border_sites
for (i in seq(1, length(annot_rect), by = 2)) {
  alt_count_annot_rect = alt_count_annot_rect + annotate("rect", xmin=annot_rect[i], xmax=annot_rect[i+1], ymin=-Inf, max=Inf, fill="black", alpha=0.1) 
}


###### plot: parameters ######
## automatic scaling of y axis
## maximum value for y
sum_perc_max <- max(del_sub_ins_df_data_to_plot_sum_perc$sum_perc) 
## y axis ticks: <10% -> 2% | > 10% <= 25% -> 5% | > 25% <= 50% -> 10% | > 50% -> 20%
y_axis_ticks <- ifelse(sum_perc_max <= 10, 2, 
                       ifelse(sum_perc_max > 10 & sum_perc_max <= 25, 5, 
                              ifelse(sum_perc_max > 25 & sum_perc_max <= 50, 10, 
                                     ifelse(sum_perc_max > 50, 20))))

###### plot: plotting histogram data of sequences length, hist_freq_seq_length.pdf ######
alt_count_bc <-
  alt_count_annot_rect + # add earlier prepared rectangles to graph
  # geom bar
  geom_vline(xintercept = pam_pos, linetype = "dashed", linewidth=0.3, col = "#59C74C") + # Cas9 Cleavage
  geom_bar(stat = "identity", width = 1) + #, size=1
  scale_fill_manual(values=c("Insertions" = "#FF0033", "Deletions" = "#3366FF"), breaks=c("Deletions", "Insertions")) +
  scale_x_continuous(labels=scales::comma, 
                     breaks=c(1, seq(ceiling(bc_len/10), bc_len, ceiling(bc_len/10))), 
                     limits=c(-4, bc_len + 5), 
                     expand = c(0.001, 0.001)) +
  scale_y_continuous(labels=function(x) paste0(x, "%"),
                     limits=c(0, plyr::round_any(max(del_sub_ins_df_data_to_plot_sum_perc$sum_perc), y_axis_ticks, f = ceiling)), # automatic
                     breaks=seq(from=0, to=plyr::round_any(max(del_sub_ins_df_data_to_plot_sum_perc$sum_perc), y_axis_ticks, f = ceiling), by=y_axis_ticks),
                     expand = c(0, 0)) +
  lemon::coord_capped_cart(left="both", bottom="both") +
  lemon::facet_rep_grid(rows = vars(sample), cols=vars(alt_long_names), scales="free_y", repeat.tick.labels = TRUE) + # repeat.tick.labels = TRUE -> the same y scales on all
  labs(y = "Marking Frequency", x = "Position of Nucleotides", fill = "Type of Marking:") +
# add theme
  theme(plot.margin = unit(c(1, 1, 1, 1), "mm"),
        axis.ticks = element_blank(), # disable ticks lines
        axis.line.y = element_line(colour="black", size=0.3), # axis y line only
        axis.line.x = element_line(colour="black", size=0.3), # axis x line only
        panel.border = element_blank(), # disable panel border
        panel.grid.major = element_blank(), # disable lines in grid on X-axis
        panel.grid.minor = element_blank(), # disable lines in grid on X-axis
        axis.text.y = element_text(size=8, angle=0, hjust=1, vjust=0.5),
        axis.text.x = element_text(size=8, angle=0, hjust=0.5, vjust=0.5),
        axis.ticks.x = element_line(colour="black", linewidth=0.3),
        axis.ticks.y = element_line(colour="black", linewidth=0.3),
        strip.background=element_blank(), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12),
        legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, -5, -5, -5),
        panel.spacing.y = unit(7.5, "mm"), 
        panel.background = element_rect(fill="white")) 
## save pdf
ggsave(plot=alt_count_bc, 
       filename=file.path(graphs_analysis_dir, "hist_freq_indels.pdf"),
       width=20, height=5*length(sample_order), units = "cm") 


##### csv data: save "hist_freq_indels.csv" ######
write.csv(del_sub_ins_df_data_to_plot_sum_perc,  
          file.path(graphs_analysis_dir, "hist_freq_indels.csv"),
          row.names = FALSE, quote = FALSE)


### tuniec ###
