################################################################
###### script drawing migration patterns based on MACHINA ######
################################################################

source("scripts/plotting_scripts/1_utils/01.1_libs.R")
source("scripts/plotting_scripts/1_utils/01.2_own_funct_softw.R")
source("scripts/plotting_scripts/1_utils/01.3_graphics.R")

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]
output_dir <- args[2]

###### input: machina_migration ######
machina_migration_df <- read_delim(infile)

## output dir: for machina analysis
machina_migration_dir <- paste0(output_dir, "/", "machina_migration_plots")
if (!dir.exists(machina_migration_dir)) 
{dir.create(machina_migration_dir, recursive = TRUE)}

###### input: parameter for clonal populations ######
## clonal populations to follow
# cp <- "all"
# cp <- "CP13"
# choose all cp
cp <- unique(machina_migration_df$CP)

###### adjust input ######
## analyze single clones
# get colors corresponding to numbers
machina_migration_cp_df <-
  machina_migration_df %>%
  dplyr::filter(CP %in% cp) %>% # get CPS and tree component 
  dplyr::select(CP, model, migrations, TreeMetRate, contains(":")) %>% 
  pivot_longer(cols=5:dim(.)[2], names_to = "connection", values_to = "total_weight") %>% # 4 -> after: 1) CP, 2) model, 3) migrations; dim(.)[2] = dim(machina_mx)[2]
  tidyr::separate(col=connection, into = c("source_of_seeding", "recipient_of_seeding"), sep=":") %>% 
  dplyr::filter(!CP == "all") %>% # temp: remove all as there is a bug in the original data frame
  group_by(CP) %>% 
  mutate(total_weight_perc = round(total_weight/sum(total_weight), 2)) %>%  # total_weight_perc - % to total
  group_by(CP, source_of_seeding) %>% 
  mutate(row_weight_perc = round(total_weight/sum(total_weight), 2)) %>%  # total_weight_perc - % to total
  mutate(source_of_seeding = as.factor(source_of_seeding)) %>% 
  mutate(recipient_of_seeding = as.factor(recipient_of_seeding)) %>% 
  mutate(CP = as.factor(CP)) %>% 
  mutate(model_long = ifelse(model == "pS", "Single-Source Seeding (pS)", 
                             ifelse(model == "pPS", "Parallel Single-Source Seeding (pPS)",
                                    ifelse(model == "mPS", "Multi-Source Seeding (mPS)", 
                                           ifelse(model == "pR", "Re-Seeding (pR)", NA)))))  
machina_migration_cp_df$row_weight_perc[is.nan(machina_migration_cp_df$row_weight_perc)] <- 0

###### summarize all ######
machina_machina_migration_all_df <-
  machina_migration_cp_df %>% 
  group_by(source_of_seeding, recipient_of_seeding) %>% 
  summarise(total_weight = sum(total_weight)) %>%  # total_weight_perc - % to total
  ungroup() %>% 
  mutate(total_weight_perc = round(total_weight/sum(total_weight), 2)) %>%  # total_weight_perc - % to total
  ungroup() %>% 
  group_by(source_of_seeding) %>% 
  mutate(row_weight_perc = round(total_weight/sum(total_weight), 2))
# count all migrations
migration_all <-
  machina_migration_cp_df %>%
  ungroup %>%
  dplyr::select(CP, migrations) %>%
  group_by(CP) %>%
  unique() %>% 
  ungroup %>%
  summarise(migrations=sum(migrations)) %>%
  deframe()
# add columns with data for all 
machina_machina_migration_all_df <-
  machina_machina_migration_all_df %>% 
  add_column(CP = "all", .before = "source_of_seeding") %>% 
  add_column(migrations = migration_all, .after = "CP")


###### merge data for final analysis ######
machina_migration_final_df <- rbind(machina_migration_cp_df, machina_machina_migration_all_df)


###### visualization: drawing weighted transition matrices ######
# loop to plot all matrices
for (cp in machina_migration_final_df$CP) {
  
### adjustments for strength scale ("Reds")
# max scale
max_migration_weight_perc <-
  machina_migration_final_df %>%
  filter(CP == cp) %>%
  ungroup() %>% 
  dplyr::select(total_weight_perc) %>% 
  deframe() %>% 
  max(., na.rm = TRUE) %>% 
  plyr::round_any(., 0.2, ceiling)
# seq of cale
seq_migration_weight_perc <- seq(0.1, from=0, to=max_migration_weight_perc)  
  

tissue_order <- unique(machina_migration_final_df[machina_migration_final_df$CP == cp,]$recipient_of_seeding)
tissue_order <- tissue_order[tissue_order  != "PRL"]
tissue_order <- as.character(droplevels(tissue_order))
tissue_order <- c("PRL",tissue_order [order(names(setNames(tissue_order , tissue_order )))])

## plot transtiom matrixes 
tree_trans_mx_weight_plot <-
  dplyr::filter(machina_migration_final_df, CP == cp) %>% # choose row to be plotted 
  #mutate(recipient_of_seeding = factor(recipient_of_seeding, levels=tissue_order)) %>%
  mutate(recipient_of_seeding = factor(recipient_of_seeding, levels=tissue_order),
         source_of_seeding = factor(source_of_seeding, levels=tissue_order)) %>%
  ggplot(aes(x = recipient_of_seeding, y = source_of_seeding, fill= total_weight_perc), height=0.15, width=0.15) + 
  geom_rtile(color = "white", size = 1, radius = unit(15, "mm"), linejoin = "mitre", show.legend = NA, inherit.aes = TRUE) + 
  
  #scale_fill_continuous_sequential(palette = "Reds", limits=c(0.00, 1.00), rev = T, na.value = "grey") +
  scale_fill_continuous_sequential(palette = "Reds", breaks=seq_migration_weight_perc, limits=c(0.00, max_migration_weight_perc), rev = T, na.value = "grey") +
  
  geom_text(aes(label = round(total_weight_perc, 2)), color = "black", size = 4) +
  #scale_x_discrete(limits = sample_order, expand=c(0, 0)) +
  #scale_y_discrete(limits = sample_order, expand=c(0, 0)) +
  labs(y = "Source of Seeding", x = "Recipient of Seeding", title = cp) +
  # add theme
  barplot_nowaklab_theme() +
  theme(aspect.ratio=1,
        plot.margin = unit(c(1, 1, 1, 1), "mm"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(margin=margin(-3,0,0,0), angle=0, hjust=0.5, vjust=1, size=12),
        axis.text.y = element_text(margin=margin(0,-3,0,0), angle=90, hjust=0.5, vjust=1, size=12),
        plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5, lineheight = 0.9),
        plot.subtitle = element_text(color = "black", size = 10, face = "plain", hjust = 0.5, lineheight = 0.9),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        legend.title = element_blank(), # change legend title font size
        legend.text = element_text(size=10)) #c hange legend text font size
# try make bigger legend that is the size of geom_tiles()
tree_trans_mx_weight_plot <- tree_trans_mx_weight_plot + guides(fill = guide_colourbar(barwidth = 1.0, barheight = 8.5)) # 10.238
# save pdf
ggsave(filename = file.path(paste0(machina_migration_dir, "/", "trans_mx_", cp, ".pdf")), plot = tree_trans_mx_weight_plot, height = 7.5, width = 10, units = "cm")}

### tuniec ####

