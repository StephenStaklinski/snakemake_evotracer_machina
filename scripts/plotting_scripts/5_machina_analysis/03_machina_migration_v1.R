source("scripts/plotting_scripts/1_utils/01.1_libs.R")
source("scripts/plotting_scripts/1_utils/01.2_own_funct_softw.R")
source("scripts/plotting_scripts/1_utils/01.3_graphics.R")

args = commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]
primary_tissue <- args[3]

# Read input data
tissue_matrix <- read.csv(infile, row.names = 1)

# Convert matrix to long format
machina_migration_df <- as.data.frame(as.table(as.matrix(tissue_matrix)))
colnames(machina_migration_df) <- c("source_of_seeding", "recipient_of_seeding", "total_weight")

# Calculate percentages
machina_migration_df <- machina_migration_df %>%
  group_by(source_of_seeding) %>%
  mutate(total_weight_perc = round(total_weight / sum(total_weight), 2)) %>%
  ungroup()

machina_migration_df <- machina_migration_df %>%
  mutate(total_weight = ifelse(is.na(total_weight), 0, total_weight),
          total_weight_perc = ifelse(is.na(total_weight_perc), 0, total_weight_perc))

# Visualization: drawing weighted transition matrices
max_migration_weight_perc <- machina_migration_df %>%
  ungroup() %>%
  dplyr::select(total_weight_perc) %>%
  deframe() %>%
  max(., na.rm = TRUE) %>%
  plyr::round_any(., 0.2, ceiling)

seq_migration_weight_perc <- seq(0.1, from=0, to=max_migration_weight_perc)

tissue_order <- unique(machina_migration_df$recipient_of_seeding)
tissue_order <- tissue_order[tissue_order != primary_tissue]
tissue_order <- as.character(droplevels(tissue_order))
tissue_order <- c(primary_tissue, tissue_order[order(names(setNames(tissue_order, tissue_order)))])

# Plot transition matrices
tree_trans_mx_weight_plot <- machina_migration_df %>%
  mutate(recipient_of_seeding = factor(recipient_of_seeding, levels=tissue_order),
         source_of_seeding = factor(source_of_seeding, levels=tissue_order)) %>%
  ggplot(aes(x = recipient_of_seeding, y = source_of_seeding, fill= total_weight_perc), height=0.15, width=0.15) + 
  geom_rtile(color = "white", size = 1, radius = unit(15, "mm"), linejoin = "mitre", show.legend = NA, inherit.aes = TRUE) + 
  scale_fill_continuous_sequential(palette = "Reds", breaks=seq_migration_weight_perc, limits=c(0.00, max_migration_weight_perc), rev = T, na.value = "grey") +
  geom_text(aes(label = round(total_weight_perc, 2)), color = "black", size = 4) +
  labs(y = "Source of Seeding", x = "Recipient of Seeding") +
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
        legend.title = element_blank(),
        legend.text = element_text(size=10))

tree_trans_mx_weight_plot <- tree_trans_mx_weight_plot + guides(fill = guide_colourbar(barwidth = 1.0, barheight = 8.5))

# Save plot as PDF
ggsave(filename = file.path(outfile), plot = tree_trans_mx_weight_plot, height = 12.5, width = 15, units = "cm")
