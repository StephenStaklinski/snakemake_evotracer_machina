##########################################################################################################################
# script to plot tree_msa_bubble_all_clonal_populations - graphs with tree, msa, bubble for all clonal populations (CPs) #
##########################################################################################################################
source("scripts/plotting_scripts/1_utils/01.1_libs.R")
source("scripts/plotting_scripts/1_utils/01.2_own_funct_softw.R")
source("scripts/plotting_scripts/1_utils/01.3_graphics.R")

args = commandArgs(trailingOnly=TRUE)

evo_object <- args[1]
output_dir <- args[2]

load(evo_object)
sample_order <- EvoTraceR_object$sample_order

###### plotting/calculations parameters ######
# space for additional scale on x for a space vertical bard fo CP ## MMUS1469=3
offset_for_cp_bar <- 4.75
## seed
set.seed(1980) 
## CP non-marked
cp_nmbc <- "CP00"

## CP non-marked automatic: create auotmatically CP with 0s
#df_to_plot_final <-
#cp_nmbc <- "CP0"

###### input ######
# ## create cp00 ASVs
# cp00_asv <- 
#   EvoTraceR_object$plot_summary$df_to_plot_final %>% 
#   filter(group == cp_nmbc) %>% 
#   dplyr::select(asv_names) %>% 
#   unique() %>% 
#   deframe()
# # remove non-marked barcode
# cp00_asv <- cp00_asv[! cp00_asv %in% "BC10v0"]  


## download asv stat data
df_to_plot_final <- 
  EvoTraceR_object$plot_summary$df_to_plot_final %>% 
  mutate(sample=fct_relevel(sample, sample_order)) # adjust columns names and order

#df_to_plot_final <- 
#  EvoTraceR_object$plot_summary$df_to_plot_final

## load trees for all CPs
tree_df <- 
  EvoTraceR_object$phylogeny$tree

## load tree phylo
tree_phylo <- 
  EvoTraceR_object$phylogeny$tree_phylo


###### output ######
# create output subdir: for phylogeny analysis (alt. output_dir -> EvoTraceR_object$output_directory)
if (!dir.exists(output_dir))
{dir.create(output_dir, recursive = TRUE)}

###### (1) adjust data for bubble plot (quantiles based) ######
## adjust based on quanties
df_to_plot_final <-
  df_to_plot_final %>% 
  ## create quantiles
  mutate(perc_asv_qnt = ifelse(perc_asv >= 0 & perc_asv < 25, 0.25,
                                      ifelse(perc_asv >= 25 & perc_asv < 50, 0.50, 
                                             ifelse(perc_asv >= 50 & perc_asv < 75, 0.75, 1.00)))) %>% 
  mutate(perc_asv_min_rank = dplyr::min_rank(-count)) %>% 
  mutate(perc_asv_ntile = ntile(-count, 10))


###### (2) adjust data for the phylogenetic tree plot ######
### data for non-marked barcode
barcode_tip <- 
  tree_df %>% 
  dplyr::filter(label == EvoTraceR_object$reference$ref_name)

### adjust position of clades
# Cassiopeia puts first the sequences that are not assigned to any cluster: it puts them at the bottom of the tree
# So, if the barcode is not at the bottom I should put it there, swapping it with another sequence
if (nrow(barcode_tip) > 0) {
  first_tip <- 
    tree_df %>% 
    filter(y == 1)
  if (barcode_tip$y != 1) {
    current_barcode_y = barcode_tip$y
    tree_df <- 
      tree_df %>% 
      dplyr::mutate(label = ifelse(y == 1, barcode_tip$label, label)) %>%
      dplyr::mutate(label = ifelse(y == current_barcode_y, first_tip$label, label))      
  }
}


### order position of clades in CPs so the top ones are the most edited in each clonal population (CP)
## create order of CPs for sapply
# order strings containing embedded numbers so that the numbers are numerically sorted rather than sorted by character value
cps <- as.character(gtools::mixedsort(unique(tree_df$group), decreasing = TRUE)) 
cp0 <- cps[length(cps)]
cps <- c(cp0, head(cps, -1))
## get tips/leaves for the given CP
rotated_trees <- sapply(cps, function(cp_of_interest) {
  cp_tips <- 
    tree_df %>% 
    dplyr::filter((group == cp_of_interest) & isTip == TRUE) %>% 
    pull(label)
  ## single CP tree to plot, remove tips in a phylogenetic tree based on tips names (all tips droped apart of CP of interest) 
  cp_tree <- 
    ape::drop.tip(phy=tree_phylo, tip=setdiff(tree_df %>% dplyr::filter(isTip) %>% pull(label), cp_tips))
  ## parents order per CP using x (number of edits per ASV)  
  cp_parent_ord <- 
    tree_df %>% 
    dplyr::filter(group == cp_of_interest) %>% 
    dplyr::filter(!is.na(label)) %>% 
    arrange(desc(x)) %>%
    dplyr::select(parent) %>%
    deframe() %>% 
    unique() %>%
    as.character()
  ## tip/leaves order per CP using x (number of edits per ASV)  
  cp_parent_x_ord <- 
    tree_df %>% 
    dplyr::filter(group == cp_of_interest) %>% 
    dplyr::filter(!is.na(label)) %>%
    mutate(parent=factor(parent)) %>% 
    mutate(parent = fct_relevel(parent, cp_parent_ord)) %>% 
    arrange(match(parent, cp_parent_ord)) %>%
    group_by(parent) %>% 
    arrange(desc(x), .by_group = TRUE) %>%
    pull(label)
  ## rotating internal branches giving a constraint on the order of the tips  
  cp_tree_ord <- 
    ape::rotateConstr(cp_tree, constraint = rev(cp_parent_x_ord))
  ## create tree data frame based on new tip order
  cp_tree_df_ord <- fortify(cp_tree_ord)
  ## plot tree
  ggtree_mp <-
    ggtree::ggtree(cp_tree_ord, ladderize = FALSE) #+ #right = TRUE
    #ggtree::geom_tiplab(geom = "text", align=TRUE, linesize=0.5, linetype="dotted", size=10)
  return(rev(cp_parent_x_ord))
})

### order the full tree | TEMP: it work until here
tree_phylo_ord <- 
  ape::rotateConstr(EvoTraceR_object$phylogeny$tree_phylo, unlist(rotated_trees))


###### (2A) plot phylogenetic tree ######
## phylogenetic tree colors (vertical bars) for clonal populations (CPs)
# parameter space: it controls additional scale on x for a space vertical bard fo CP ## MMUS1469=3
offset_for_cp_bar <- offset_for_cp_bar

## plot phylogenetic tree
ggtree_mp <- 
  ggtree::ggtree(tree_phylo_ord, ladderize = FALSE) + 
  ggtree::geom_tiplab(geom = "text", align=TRUE, linesize=0.5, linetype="dotted", size=7) +
  xlab("Phylogenetic Tree \n Cassiopeia Greedy")


## assign colors
if (("group" %in% colnames(tree_df))) {
  colors <- c(set3_12xcols, sample(rainbow(n = length(unique(tree_df$group))-length(set3_12xcols)))) # 12x set 3
  names(colors) <- unique(tree_df$group)
  colors[[cp_nmbc]] = "black"
  for (c in unique(tree_df$group)) {
    cluster_nodes <- 
      tree_df %>% 
      filter(group == c)
    if (c != cp_nmbc) { # & isTip %>% arrange(y) %>% pull(node)
      min_parent <- min(setdiff(cluster_nodes$parent, c(min(tree_df$parent))))
      ggtree_mp <- 
        ggtree_mp +
        ##
        geom_cladelab(node = min_parent, label = "", # width and position of CP bar
                      textcolor = NA, 
                      barcolor = colors[[c]], extend = 0.4,
                      offset = offset_for_cp_bar, angle = 90, barsize = 12.5, align = TRUE) + # offset: distance between tree and bar, default is 0, barsize: width of bar
        ##
        geom_cladelab(node = min_parent, label=c, fontsize = 8, fontface = "bold", hjust = "center", vjust = "center",
                      textcolor = "gray25", offset.text = -0.125, # controls text positioning the CP bar, smaller to the left (e.g. negative)
                      barcolor = NA, extend = NA,
                      offset = offset_for_cp_bar, angle = 90, barsize = NA, align = TRUE)
    } else {
      ggtree_mp <- 
        ggtree_mp +
        geom_strip(cluster_nodes[1], cluster_nodes[length(cluster_nodes)], barsize=15, color=colors[[c]],
                  #label= c,
                  #offset.text=0.1,
                  offset = 0.85,
                  #angle = 90,
                  align = T)
    }
  }
} 

## phylogenetic tree plotting (final tree)
# plot
ggtree_mp <-
  ggtree_mp +
  scale_x_continuous(expand = c(0, 0), limits=c(0, max(tree_df$x)+offset_for_cp_bar+1), breaks=c(0:max(tree_df$x))) + # limits=c(0, max(tree_df$x)+1)
  lemon::coord_capped_cart(bottom="both") +
# add theme  
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        axis.text.y = element_blank(), # disable y axis text
        #axis.title.x = element_text(size=8, angle=0),
        axis.ticks.x = element_line(colour="black", linewidth=0.5),
        axis.line.x = element_line(colour="black", linewidth=0.5))


###### (3) adjust data for the plot msa ######
## parameters
cleaned_deletions <- FALSE
subset_asvs <- unique(df_to_plot_final$asv_names)
## choose what mutations were included
if (cleaned_deletions == "del") {
  to_plot_df = EvoTraceR_object$alignment$asv_barcode_alignment %>% 
    mutate(alt = ifelse(alt == "d", alt, "w"))
} else if (cleaned_deletions == "del_ins") {
  to_plot_df = EvoTraceR_object$alignment$asv_barcode_alignment %>% 
    mutate(alt = ifelse(alt %in% c("d","i"), alt, "w"))
} else {
  to_plot_df = EvoTraceR_object$alignment$asv_barcode_alignment
}
## put insertions after wt for visualization
to_plot_df <- 
  to_plot_df %>% 
  dplyr::arrange(asv_names, position_bc260, desc(alt))
## add full names for labeling of plots
to_plot_df <-
  to_plot_df %>%
  mutate(alt_long_names = ifelse(alt == "i", "Insertion", 
                                 ifelse(alt == "d", "Deletion", 
                                        ifelse(alt == "s", "Substitutions", "No Edits"))))
## position of PAM in guides
pam_pos <- EvoTraceR_object$reference$ref_cut_sites
bc_len <- nchar(EvoTraceR_object$reference$ref_seq)
# adjust for height of tiles -> "sub" smaller
to_plot_df$tile_height <- ifelse(to_plot_df$alt == "s", 0.3, 0.75)


###### frames around msa apart ORG ######
msa_frame <- data.frame(xmin= 0.1, xmax = bc_len+0.9, 
                        ymin = seq(from = 1.6, to = nlevels(as.factor(to_plot_df$asv_names)), by=1), 
                        ymax = seq(from = 2.4, to = nlevels(as.factor(to_plot_df$asv_names))+1, by=1))
#to_plot_df$alt = factor(x = to_plot_df$alt, levels = c('sub', 'del', 'wt', 'ins'))
to_plot_df$alt_long_names = factor(x = to_plot_df$alt_long_names)#, levels = c('Deletion', 'w', 'i'))
to_plot_df$asv_names = factor(x = to_plot_df$asv_names, levels = subset_asvs)

## calculate insertion size
to_plot_df <- 
  to_plot_df %>% 
  group_by(asv_names, position_bc260) %>% 
  mutate(ins_size = length(position_bc260)-1) %>%  # count size of insertions | -1 to account for the next "wt"
  group_by(asv_names, alt) %>% 
  mutate(del_size = length(position_bc260)-1) # count size of deletions | -1 to account for the next "wt"


###### (3A) plot final msa ###### 
## plot
msa_cna_bc <- 
  ggplot(data=to_plot_df, aes(x=position_bc260, y=asv_names)) +
  geom_vline(xintercept=pam_pos, linetype="solid", size=0.7, col="#59C74C") + # Cas9 Cleavage "#59C74C"
  geom_tile(aes(fill=alt_long_names, width=0.75, height=tile_height), colour = NA) +
  scale_fill_manual(values=c("Deletion"="#3366FF", "Insertion"="#FF0033", "No Edits"="#f2f2f2"), breaks=c("Insertion", "No Edits", "Deletion")) + 
  scale_x_continuous(labels=scales::comma, breaks=c(1, seq(26, 260, 26)), expand = c(0.014, 0.014)) +
  
  ## filled diamond: shape=18; data specific for "ins", size = pos_diff based on size of insertion but check as diamonds size different
  geom_point(data=to_plot_df %>% filter(alt == "i"), colour = "red", aes(size=ins_size), shape=18) + 
  scale_size_continuous(range = c(2.5, 7.5), breaks = c(1, 5, 10, 25, 50)) + # range size of shape
  
  ## rectangles around marked barcodes
  geom_rect(data=msa_frame, mapping=aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"), colour = "grey75", fill=NA, inherit.aes = FALSE, size=0.6) + # around ASVs
  ## rectangle around non-marked barcode
  geom_rect(xmin=0.1, xmax=260+0.9, ymin=0.6, ymax=1.4, colour = "#65A7F3", fill=NA, size=0.6) + ## around non-marked BC10
  xlab("Barcode (BC10) Nucleotides \n (1-260)") +
  lemon::coord_capped_cart(bottom="both") +
  ## guides/legend
  guides(fill = guide_legend("Type of\nEditing"), order = 1) + # , order = 1
  guides(size = guide_legend("Insertion\nSize"), order = 2) + #, order = 2
## add theme
  barplot_nowaklab_theme() +
  scale_y_discrete(c(0, length(subset_asvs)), expand = c(0, 0.6)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        axis.line.y = element_blank(), # disable y axis lines
        axis.ticks.y = element_blank(), # disable y axis ticks lines
        axis.title.y = element_blank(), # disable y axis lines
        axis.text.y = element_blank()) # disable y axis text


###### (4) plot mutations width ######
## bar graph
bar_ins_del_sub_width <- 
  ggplot(data=dplyr::select(df_to_plot_final, asv_names, width_total_d, width_total_i) %>% #, width_total_s) %>% 
           unique() %>% 
           tidyr::gather(key="alt", "width_total_d", "width_total_i",  value="width_sum") %>% #"width_total_s",
           filter(!stringr::str_detect(asv_names, "NMBC"))) +
  geom_bar(aes(x=width_sum, y=asv_names, fill=alt), position="stack", stat="identity", width=0.8, size=0.2) +
  scale_fill_manual(values=c("width_total_d" = "#3366FF", "width_total_i" = "#FF0033"), #"width_total_s" = "#329932" 
                    breaks=c("width_total_d", "width_total_i")) + #, "width_total_s"
  geom_vline(xintercept=0, linetype="solid", size=0.5, col="black") + # no indels
  scale_x_continuous(labels=scales::comma, expand = c(0.01, 0.01), breaks=c(0, 130, 260))+#, limits=c(0, 286)) +
  geom_vline(xintercept=130, linetype="dotted", size=0.5, col="gray50") +
  geom_vline(xintercept=260, linetype="dotted", size=0.5, col="gray50") +
  labs(x = "Cummulative \n Widths of Markings", fill = "Width of Marking") +
# add Theme  
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        axis.line.y = element_blank(), # disable y axis lines
        axis.ticks.y = element_blank(), # disable y axis ticks lines
        axis.title.y = element_blank(), # disable y axis lines
        axis.text.y = element_blank()) # disable y axis text


###### (5) plot length of ASVs ######
## bar graph
bar_seq_n <- 
  ggplot(data=dplyr::select(df_to_plot_final, asv_names, seq_n) %>% 
           unique() %>%
           mutate(seq_n_col=ifelse(seq_n == 260, "260", ifelse(seq_n < 260, "< 260", "> 260"))) %>%
           filter(stringr::str_detect(asv_names, "ASV"))) +
  geom_bar(aes(x=seq_n, y=asv_names, fill=seq_n_col), position = "dodge", stat = "identity", width=0.8, size=0.2) +
  geom_text(aes(x=1, y=asv_names, label = seq_n), col="white", size=4, hjust="left") + #, position=position_dodge(width=0.9)
  scale_fill_manual(values=c("260" = "#84B48F", "< 260" = "#377EB8", "> 260" = "#F39B7FFF"), breaks=c("260", "< 260", "> 260")) +
  geom_vline(xintercept=0, linetype="dotted", size=0.5, col="#377EB8") + # 0 
  geom_vline(xintercept=260, linetype="dotted", size=0.5, col="#84B48F") + # 260 bp = original length
  geom_vline(xintercept=520, linetype="dotted", size=0.5, col="#F39B7FFF") + # 0 
  scale_x_continuous(labels=scales::comma, expand = c(0.01, 0.01), limits=c(0, 572), breaks=c(0, 130, 260, 390, 520)) + # 572 = 1.1 * 520 -> nice separation between graphs
  labs(x = "ASV\nLength", fill = "ASV Length") +
  lemon::coord_capped_cart(bottom="both") # axis with lemon
# add theme  
bar_seq_n <- 
  bar_seq_n + 
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        axis.line.y = element_blank(), # disable y axis lines
        axis.ticks.y = element_blank(), # disable y axis ticks lines
        axis.title.y = element_blank(), # disable y axis lines
        axis.text.y = element_blank()) # disable y axis text


###### (6) plot_similarity ######
bar_pid <- 
  ggplot(data=dplyr::select(df_to_plot_final, asv_names, pid) %>% unique() %>%
           filter(stringr::str_detect(asv_names, "ASV")) %>%
           dplyr::mutate_if(is.numeric, round, 0), 
         aes(x=pid, y=asv_names)) +
  geom_segment(aes(y=asv_names, yend=asv_names, x=0, xend=100), size=2, color="grey75") +
  geom_vline(xintercept=c(0, 25, 50, 75, 100), linetype="solid", size=0.5, col="white") + # 100 similar
  geom_point(size=9, fill="#ce5a57", col="white", shape=21, stroke=0.5) +
  geom_text(aes(label=pid), size=4, col="white") + 
  scale_x_continuous(labels=scales::comma, expand = c(0.01, 0.01), limits=c(0, 110), breaks=c(0, 25, 50, 75, 100)) +
  xlab("Sequence \n Identity (%)") +
  lemon::coord_capped_cart(bottom="both") # axis with lemon
# add Theme  
bar_pid <- 
  bar_pid + 
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        axis.line.y = element_blank(), # disable y axis lines
        axis.ticks.y = element_blank(), # disable y axis ticks lines
        axis.title.y = element_blank(), # disable y axis lines
        axis.text.y = element_blank())


###### (7) plot organs quantity bubble - drawing ######
## data adjustments for plot bubble size/color
#
subset_asvs <- unique(df_to_plot_final$asv_names)
#
df_to_plot_final$asv_names <- factor(df_to_plot_final$asv_names, levels = subset_asvs)
#
scale_bubble <- 
  ceiling(max(df_to_plot_final %>% dplyr::pull(perc_in_sample), na.rm = T))
#
scale_bubble_nonmbc <- ceiling(max(df_to_plot_final %>% 
                                     filter(stringr::str_detect(asv_names, "ASV")) %>% 
                                     dplyr::pull(perc_in_sample), na.rm = T))

###### (8) plot bubble size/color #####
bubble <- 
  ggplot(data=df_to_plot_final, aes(x=sample, y=asv_names, scale="globalminmax")) +
  geom_point(aes(size=perc_in_sample, fill=perc_fold_to_max), shape=21, stroke=0.5, col="black") +
  colorspace::scale_fill_continuous_sequential(palette = "Plasma", limits=c(0, 100), rev = F, na.value = "grey") +
  scale_size_area(limits = c(0, scale_bubble), trans = "log1p", breaks = c(1, scale_bubble_nonmbc/2, scale_bubble_nonmbc, scale_bubble)) + # manual
  labs(x = "Analyzed \n Samples", size = "Frequency in\nSample", fill = "Frequency\nNormalized to Max") +
  scale_y_discrete(c(0, length(subset_asvs))) +
  lemon::coord_capped_cart(bottom="both") +
# add theme
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        #legend.title.align = 0.5,
        axis.line.y = element_blank(), # disable y axis lines
        axis.ticks.y = element_blank(), # disable y axis ticks lines
        axis.title.y = element_blank(), # disable y axis lines
        axis.text.y = element_blank(), # disable y axis text
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(colour="grey75", linewidth=0.5, "dotted"), # y grid line 
        panel.grid.major.x = element_line(colour="grey75", linewidth=0.5, "dotted")) # disable lines in grid on X-axi


##### (9) bubble size/color with quartiles ######
bubble_qnt <- 
  ggplot(data=df_to_plot_final, aes(x=sample, y=asv_names)) + # scale="globalminmax"
  geom_point(aes(size=perc_asv_qnt, color=sample), shape=19) +# stroke=0.5
  scale_size_area(limits = c(0.25, 1.00), breaks = c(0.25, 0.50, 0.75, 1.00), max_size = 9.5) + # manual
  #scale_color_manual(values = sample_color[sample_order]) +
  labs(x = "Analyzed \n Samples", size = "Frequency in Sample (Organ)\nAdjusted Based on Quartiles", color = "Sample (Organ)\nPrimary/Metastases") +
  scale_y_discrete(c(0, length(subset_asvs))) +
  lemon::coord_capped_cart(bottom="both") +
  
  guides(color = guide_legend(override.aes = list(size = 5))) + # size in color circels
  
  # add theme
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        #legend.title.align = 0.5,
        axis.line.y = element_blank(), # disable y axis lines
        axis.ticks.y = element_blank(), # disable y axis ticks lines
        axis.title.y = element_blank(), # disable y axis lines
        axis.text.y = element_blank(), # disable y axis text
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(colour="grey75", linewidth=0.6, linetype="solid")) # y grid line 
        #panel.grid.major.x = element_line(colour="grey75", size=0.5, linetype="dotted")) # disable lines in grid on X-axi


###### (10) plot rtiles organs present/absent ######
# plot
rtile <- 
  ggplot(data=df_to_plot_final, aes(x=sample, y=asv_names, fill=sample)) +
  geom_rtile(color = "white", lwd = 2.0, radius = unit(1.5, "mm")) +
  #scale_fill_manual(values = sample_color[sample_order]) +
  #scale_x_discrete(limits = rev(topology_sort), labels = topology_labels_line, expand=c(0, 0)) +
  xlab("Organs") +
  labs(fill = "Organs") +
# add theme
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        #legend.title.align = 0.5,
        axis.line.y = element_blank(), # disable y axis lines
        axis.ticks.y = element_blank(), # disable y axis ticks lines
        axis.title.y = element_blank(), # disable y axis lines
        axis.text.y = element_blank(), # disable y axis text
        panel.grid.major.y = element_line(colour="grey75", linewidth=0.5, linetype="dotted"), # y grid line 
        #legend.position="none", 
        panel.background = element_rect(fill="white"))
        

###### (i) combine all graphs together using aplot package ######
# add bubble
# msa <- aplot::insert_right(msa_cna_bc, bubble, width = 0.12)
# # add alterations width
# msa_cna_bc.bar_ins_del_sub_width <- aplot::insert_right(msa, bar_ins_del_sub_width, width = 0.25) 
# # add ggtree_mp
# msa_cna_bc.bar_ins_del_sub_width.ggtree_mp <- aplot::insert_left(msa_cna_bc.bar_ins_del_sub_width, ggtree_mp, width = 1)#0.75)
# # add bar_seq_n
# msa_cna_bc.bar_ins_del_sub_width.ggtree_mp.bar_seq_n <- aplot::insert_right(msa_cna_bc.bar_ins_del_sub_width.ggtree_mp, bar_seq_n, width = 0.2)
# # add bar_pid
# msa_cna_bc.bar_ins_del_sub_width.ggtree_mp.bar_seq_n.bar_pid <- aplot::insert_right(msa_cna_bc.bar_ins_del_sub_width.ggtree_mp.bar_seq_n, bar_pid, width = 0.2)
# # add alterations width
# msa_cna_bc.bar_ins_del_sub_width <- aplot::insert_right(msa, bar_ins_del_sub_width, width = 0.25) 
# # add ggtree_mp
# msa_cna_bc.bar_ins_del_sub_width.ggtree_mp <- aplot::insert_left(msa_cna_bc.bar_ins_del_sub_width, ggtree_mp, width = 1)#0.75)
# # temp
# ggsave(filename=file.path("~/Desktop/cp_tree_all.pdf"), 
#        plot=msa_cna_bc.bar_ins_del_sub_width.ggtree_mp.bar_seq_n.bar_pid, 
#        width=60, height=dim(tree_df)[1]*0.7, units = "cm", limitsize = FALSE)


###### (ii) combine all graphs together using aplot package ######

#msa_cna_bc_bubble_qnt <- aplot::insert_right(msa_cna_bc, bubble_qnt)

msa_cna_bc_bubble_qnt_ggtree_mp <- aplot::insert_left(msa_cna_bc, ggtree_mp, width = 0.4) # width = 0.95 before 03/26/23

## add legend at the bottom of the plot
msa_cna_bc_bubble_qnt_ggtree_mp <- print(msa_cna_bc_bubble_qnt_ggtree_mp) & theme(legend.position = "top", legend.box = "horizontal")


## for publication and presentations
ggsave(filename=output_dir + "/cp_tree_msa_cna_bc_bubble_qnt_ggtree_mp.pdf",plot=msa_cna_bc_bubble_qnt_ggtree_mp,width=50, height=dim(tree_df)[1]*0.7, units = "cm", limitsize = FALSE)} else {

### tuniec ###

