###### script for plotting counts edits/marks size, type and assigned affected sites ######
source("~/snakemake_evotracer_machina/scripts/plotting_scripts/1_utils/01.1_libs.R")
source("~/snakemake_evotracer_machina/scripts/plotting_scripts/1_utils/01.2_own_funct_softw.R")
source("~/snakemake_evotracer_machina/scripts/plotting_scripts/1_utils/01.3_graphics.R")

args = commandArgs(trailingOnly=TRUE)
load(args[1])
output_dir <- args[2]
sample_order <- EvoTraceR_object$sample_order

###### download data ######
histo_sites_start <- 
  drop_na(EvoTraceR_object$plot_summary$df_to_plot_final) %>% # drops unedited barcode, usually CP00
  mutate(sample=fct_relevel(sample, sample_order)) # adjust columns names and order


###### data outputs ######
# output dir: for graphs analysis
graphs_analysis_dir <- paste0(output_dir, "/graphs_analysis")
if (!dir.exists(graphs_analysis_dir)) 
{dir.create(graphs_analysis_dir, recursive = TRUE)}


###### prepare data for plotting ######
# adjust data
histo_sites <- 
  histo_sites_start %>%
  mutate(n_sites_num = as.numeric(n_sites)) %>%
  #dplyr::select(asv_names, sample, count, mut_id, mutation_type, start, end, n_nucleotides, n_sites_num) %>%
  dplyr::select(sample, count, mutation_type, n_nucleotides, n_sites_num) %>%
  group_by(sample, n_nucleotides, mutation_type) %>%
  mutate(count=sum(count)) %>%
  dplyr::select(sample, n_nucleotides, mutation_type, n_sites_num, count) %>%
  unique() %>%
  mutate(count_log10=log(count, 10)) %>%
  mutate(count_log10=if_else(mutation_type == "d", -1*count_log10, count_log10)) %>%
  mutate(mutation_type_n_sites_num = if_else(mutation_type == "d" & n_sites_num == 1, "1xTS affected, deletion",
                                             if_else(mutation_type == "d" & n_sites_num == 2, "2xTS affected, deletion",
                                                     if_else(mutation_type == "d" & n_sites_num == 3, "3xTS affected, deletion",
                                                             if_else(mutation_type == "d" & n_sites_num >= 4, ">= 4xTS affected, deletion",
                                                                     if_else(mutation_type == "i" & n_sites_num == 1, "1xTS affected, insertion",
                                                                             if_else(mutation_type == "i" & n_sites_num == 2, "2xTS affected, insertion",
                                                                                     if_else(mutation_type == "i" & n_sites_num == 3, "3xTS affected, insertion",
                                                                                             if_else(mutation_type == "i" & n_sites_num >= 4, ">= 4xTS affected, insertion", "other")))))))))


# set order of affected sites
histo_sites$mutation_type_n_sites_num <- factor(histo_sites$mutation_type_n_sites_num, levels = indels_site_names)
# altewrnativel: drop levels
# histo_sites$mutation_type_n_sites_num <- droplevels(histo_sites$mutation_type_n_sites_num)


###### plotting aux data ######
# automatic scaling of y axis
histo_sites_max <- 
  histo_sites %>%
  ungroup() %>%
  dplyr::select(n_nucleotides, sample, mutation_type, count_log10) %>% 
  group_by(n_nucleotides, sample, mutation_type) %>% 
  mutate(histo_sites_count = sum(count_log10)) %>% 
  pull(histo_sites_count)
# set max and min for y-axis
max_axis_limit_y <- max(histo_sites_max)
min_axis_limit_y <- min(histo_sites_max) 

###### plot ######
bargraph_histo_sites <-
  ggplot(histo_sites, aes(x=n_nucleotides, y=count_log10, group=mutation_type_n_sites_num, fill=mutation_type_n_sites_num)) +
  geom_vline(xintercept = c(1, seq(20, 240, 20)), linetype = "dashed", color = "gray", size = 0.25) +
  geom_bar(stat = "unique", width = 0.8) +
  scale_fill_manual(values=indels_site_color[levels(histo_sites$mutation_type_n_sites_num)], drop = FALSE) +
  scale_x_continuous(labels=comma, 
                     limits=c(0, 241),
                     breaks=c(1, seq(20, 240, 20)), 
                     expand = c(0, 0)) +
  scale_y_continuous(limits=c(plyr::round_any(min_axis_limit_y, -2, f=ceiling), plyr::round_any(max_axis_limit_y, 2, f=ceiling)), # automatic
                     breaks=seq(from=plyr::round_any(min_axis_limit_y, -2, f=ceiling), to=plyr::round_any(max_axis_limit_y, 2, f=ceiling), by=2),
                     expand = c(0.01, 0.01),
                     labels = abs(seq(from=plyr::round_any(min_axis_limit_y, -2, f=ceiling), to=plyr::round_any(max_axis_limit_y, 2, f=ceiling), by=2))) +
  
  ### WiP second x-axis ###
  # x-axis description
  # scale_x_continuous(breaks=c(1, seq(26, 260, 26)), expand = c(0,0), limits=c(0, 261),
  #                    labels=c("1 \n 0%", "26 \n 10%", "52 \n 20%", "78 \n 30%", "104 \n 40%", "130 \n 50%",
  #                             "156 \n 60%", "182 \n 70%", "208 \n 80%", "234 \n 90%", "260 \n 100%")) +
  # # second x-axis
  # scale_x_continuous(labels=comma, breaks=c(1, seq(26, 260, 26)), expand = c(0,0), limits=c(1, 260),
  # sec.axis=sec_axis(~.,
  # breaks=c(1, seq(26, 260, 26)),
  # labels=c("0%", "10%","20%","30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"))) +
  ### WiP second x-axis ###
  
  labs(x = "Size of Indels (bp)", y=expression("log"[10]*"[frequency]"), fill = "Number of Affected Sites:") + #title = paste0("Number of affected Sites", ": ", group_id_select), 
  geom_hline(yintercept=0, linetype="solid", size=0.25, col="grey75") +
  lemon::coord_capped_cart(left="both", bottom="both") +
  lemon::facet_rep_grid(rows = vars(sample), repeat.tick.labels = TRUE) +
  guides(fill=guide_legend(ncol=4, byrow=TRUE))
  
  # add theme
  bargraph_histo_sites <- 
    bargraph_histo_sites + 
    theme(plot.margin = unit(c(2, 2, 2, 2), "mm"),
          axis.ticks = element_blank(), # disable ticks lines
          axis.line.y = element_line(colour="black", size=0.3), # axis y line only
          axis.line.x = element_line(colour="black", size=0.3), # axis x line only
          ### WiP ###
          #axis.line.x.top = element_blank(), # axis x line only
          #axis.ticks.x.top = element_blank(),
          #axis.text.x.top = element_text(colour="#6BAED6", size=8, angle=0, hjust=0.5, vjust=0.5),
          ### WiP ###
          panel.border = element_blank(), # disable panel border
          panel.grid.major = element_blank(), # disable lines in grid on X-axis
          panel.grid.minor = element_blank(), # disable lines in grid on X-axis
          axis.text.y = element_text(size=10, angle=0, hjust=1, vjust=0.5),
          axis.text.x = element_text(size=10, angle=0, hjust=0.5, vjust=0.5),
          axis.ticks.x = element_line(colour="black", size=0.3),
          axis.ticks.y = element_line(colour="black", size=0.3),
          strip.background=element_blank(), 
          strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12),
          legend.position="bottom", legend.box = "horizontal", legend.text=element_text(size=12), legend.title=element_text(size=12),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-5,-5,-5,-5),
          panel.spacing.y = unit(7.5, "mm"), 
          panel.background = element_rect(fill="white")) 
#bargraph_histo_sites <- bargraph_histo_sites + guides(fill = guide_legend(nrow = 1, byrow = TRUE))
# save pdf  
ggsave(plot=bargraph_histo_sites, 
       filename=file.path(graphs_analysis_dir, "hist_freq_site_affected.pdf"), 
       width=30, height=6*length(sample_order), units = "cm") 

# save csv
write.csv(histo_sites, #pg$data[[1]],  
          file.path(graphs_analysis_dir, "hist_freq_site_affected.csv"),
          row.names = FALSE, quote = FALSE)

### tuniec ###


