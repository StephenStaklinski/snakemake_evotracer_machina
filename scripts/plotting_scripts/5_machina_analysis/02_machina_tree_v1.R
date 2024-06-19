########################################################
###### script drawing tree graph based on MACHINA ######
########################################################

source("~/snakemake_evotracer_machina/scripts/plotting_scripts/1_utils/01.1_libs.R")
source("~/snakemake_evotracer_machina/scripts/plotting_scripts/1_utils/01.2_own_funct_softw.R")
source("~/snakemake_evotracer_machina/scripts/plotting_scripts/1_utils/01.3_graphics.R")

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]
output_dir <- args[2]
sample_order <- c("PRL","FML","FMR","FMR1","FMR2","HMR","LGL","LGR","LNDL","LNSL","LNSR","LVL","LVM","LVR","RBL","TBL","UNL","UNR")
##### Data Comments ######
# COMMENT 1: Warning Messages Explanation: Removed ## rows containing missing values (`geom_point()`). -> inferred states transtions have NA so no number
# COMMENT 2: Warning Messages Explanation: Removed ## rows containing missing values (`geom_point()`). -> when inferred states plotted rest is set to NA

###### input: machina_tree_graph ######
## automatic data frame download
machina_tree_df <- read_delim(infile)
## output dir: for machina analysis
machina_tree_dir <- paste0(output_dir, "/", "machina_tree_plots")
if (!dir.exists(machina_tree_dir)) 
{dir.create(machina_tree_dir, recursive = TRUE)}

###### input parameter for clonal populations ######
## clonal populations to follow
#cp <- c("CP01")
# choose all cp
cp <- unique(machina_tree_df$CP)

###### adjust input ######
# loop to plot all trees
for (cp in cp) {
# get colors corresponding to numbers
machina_colors <-
  machina_tree_df %>%
  dplyr::filter(CP %in% cp, group == "color") %>% # get CPs and tree component 
  dplyr::select(CP, value1, value2) %>%
  pull(value1, value2) %>% 
  enframe() #%>% 
  #rename("name" = "sample")
  
# retrieve labels for single clonal population (cp)
machina_labels <-
  machina_tree_df %>%
  dplyr::filter(CP %in% cp, group == "label") %>% # get CPs and tree component
  dplyr::filter(value2 != "ASVXXX") %>% # remove dummy labels
  dplyr::select(CP, value1, value2) %>%
  dplyr::rename(from = value1, to = value2) %>%
  dplyr::rename(name = to, sample = from) %>% 
  mutate_if(is.character, as.factor) %>% 
  merge(x=., y=machina_colors, by.x="sample", by.y="name") %>% 
  dplyr::select(CP, value, name)
colnames(machina_labels) <- c("CP", "sample", "name")

# create data retrieve tree data for single clonal population (cp)
machina_tree <-
  machina_tree_df %>%
  dplyr::filter(CP %in% cp, group == "tree") %>% # get CPS and tree component
  dplyr::filter(value2 != "ASVXXX") %>% # remove dummy labels
  dplyr::select(CP, value1, value2, asv_freq_count, transition_probability) %>%
  #dplyr::select(CP, value1, value2) %>%
  dplyr::rename(from = value1, to = value2) %>%
  left_join(., machina_labels, by = c("from" = "name")) %>% 
  dplyr::select(!CP.y) %>% 
  dplyr::rename(CP = CP.x) 

###### create graph data from data frame ######
machina_graph_df <- 
  as_tbl_graph(machina_tree)
  
###### add data to edges and nodes ######
# add data about samples to nodes
machina_graph_df <- left_join(machina_graph_df, machina_labels, by = "name") 
# add ASV count freq data
machina_graph_df <- left_join(machina_graph_df, machina_tree[,3:4], by = c("name" = "to"))

# remove after "_"
machina_graph_df <-
  machina_graph_df %>% 
  ## nodes
  activate(nodes) %>%
  mutate(name = gsub("_.*","", name)) %>% 
  mutate(asv_freq_count_inferred = ifelse(is.na(asv_freq_count) == TRUE, "Inferred Ancestral State", "Observed ASV State")) %>% 
  mutate_if(is.character, as.factor) %>% 
  mutate(name=as.factor(name)) %>% 
  ## edges
  activate(edges) %>%
  mutate(CP=as.factor(CP)) %>% 
  mutate(sample=as.factor(sample))

###### visualization of tree/graph ###### 
### adjustments for strength scale ("Reds")
# max scale
max_transition_probability <- machina_graph_df %>% activate(edges) %>% data.frame() %>% dplyr::select(transition_probability) %>% deframe() %>% max(., na.rm = TRUE) %>% plyr::round_any(., 0.2, ceiling)
# seq of cale
seq_transition_probability <- seq(0.2, from=0, to=max_transition_probability)

# tree ggraph
machina_tree_graph <-
  ggraph(machina_graph_df, layout="tree") + 
  
  ## edges (1): transitions weight assigned to edges
  geom_edge_diagonal(aes(color=transition_probability), flipped=FALSE,
                     edge_width=unit(0.25, "mm"), strength=0.8, arrow = arrow(length = unit(0.75, "mm"), ends = "last", type = "closed", angle = 45), start_cap = circle(1.0, "mm"), end_cap = circle(2.00, "mm")) +
  # map color to of strength to edges
  scale_edge_colour_distiller(type = "seq", palette = "Reds", direction = 1,
                              #breaks=c(0, 0.25, 0.50, 0.75, 1.00), limits=c(0.00, 0.50), # the same scales for all CPs
                              breaks=seq_transition_probability, limits=c(0.00, max_transition_probability), # specific scale for CP
                              na.value = "grey50", name = "Seeding Strength (Edges):   ",
                              guide = guide_edge_colorbar(direction = "horizontal",
                                                          title.position = "left",
                                                          title.hjust = 0, title.vjust = 0.5,
                                                          ticks = TRUE, label = TRUE, ticks.colour = "white",
                                                          label.hjust = 0.5, label.vjust = 0.5,
                                                          label.position = "bottom",
                                                          label.theme = element_text(colour = "black", size=7),
                                                          frame.colour = "transparent", barwidth = 7, barheight = 0.5)) +
  
  # ## edges (2A): organs assigned to edges, width of edge is the same for all
  # geom_edge_diagonal(aes(color=sample), edge_width=unit(0.25, "mm"), strength=1, arrow = arrow(length = unit(0.75, "mm"), ends = "last", type = "closed", angle = 45), start_cap = circle(1.0, "mm"), end_cap = circle(2.00, "mm")) +
  # scale_edge_color_manual(values=sample_color[sample_order], name="Seeding Direction From:") +
  
  # ## edges (2B): organs assigned to edges, width of edge is dependent on weight of connection
  # geom_edge_diagonal(aes(width=transition_probability, color=sample), arrow = arrow(length = unit(0.75, "mm"), ends = "last", type = "closed", angle = 45), start_cap = circle(1.0, "mm"), end_cap = circle(2.00, "mm")) +
  # scale_edge_width(range = c(0.1, 0.8), name = "Seeding Strenght (Edges):") +
  # scale_edge_color_manual(values=sample_color[sample_order], name = "Seeding Direction From (Edges):") +
  
  ## organs colors and ASVs sizes are assigned to nodes
  # nodes
  geom_node_point(aes(color = sample, size = asv_freq_count), shape = 20) + # alternatively: asv_freq_count_radius_log10 vs asv_freq_count
  # color of nodes 
  scale_color_manual(values=sample_color[sample_order], breaks= sample_order, name = "Organ/ASV (Nodes):") +
  # size of nodes
  scale_size_area(max_size = 6, breaks = c(100, 1000, 10000), labels= c(expression(10^2), expression(10^3), expression(10^4))) +
  ## shapes of nodes
  geom_node_point(aes(color = sample, shape = asv_freq_count_inferred), size=2) +
  scale_shape_manual(values=c("Inferred Ancestral State" = 18), limits=c("Inferred Ancestral State"), name = "Character State (Nodes):") +

  ## labs
  labs(x = NULL, y = NULL, title = paste0("Migration History of ", cp)) +
  #scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) + # scale_x_reverse
  scale_x_reverse(expand = expansion(mult = c(0.01, 0.01))) + # scale_x_reverse
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.03))) +
  
  ## guides  
  guides(col = guide_legend("Organ/ASV (Nodes):", order = 1), 
         size = guide_legend("ASV Freq. (Nodes):", order=2),
         col = guide_legend("Seeding Direction From (Edges):", order = 4),
         width = guide_legend("Seeding Strength (Edges):", order = 5),
         shape = guide_legend("Character State (Nodes):", order = 3)) +
  
  ## add theme
  barplot_nowaklab_theme() +
  theme(# axis
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        # panels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # legends
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin=margin(0, 0, 0, 0, unit = "mm"), # legend.margin=margin(0, 0, 0, 0, unit = "mm"),
        legend.title=element_text(size=10, margin = margin(t = 0, r = -1, b = 0, l = 2, unit = "mm")),
        legend.text=element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))) # legend.text=element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = -2, unit = "mm")))
# save
ggsave(paste0(machina_tree_dir, "/", cp, "_", "machina_tree_graph.pdf"), machina_tree_graph, width = 45, height = 10, units = "cm")
}

### tuniec ###
