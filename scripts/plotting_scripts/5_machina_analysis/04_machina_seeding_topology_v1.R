source("scripts/plotting_scripts/1_utils/01.1_libs.R")
source("scripts/plotting_scripts/1_utils/01.2_own_funct_softw.R")
source("scripts/plotting_scripts/1_utils/01.3_graphics.R")

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]
output_dir <- args[2]

machina_seeding_topology_df <- read_delim(infile)

machina_seeding_topology_dir <- paste0(output_dir, "/", "machina_seeding_topology")
if (!dir.exists(machina_seeding_topology_dir)) {
  dir.create(machina_seeding_topology_dir, recursive = TRUE)
}

topology_color <- c("absent" = "white",
                    "primary_confined" = "#8FCEF2",
                    "primary_mono_seeding" = "#4486B6",
                    "primary_parallel_seeding" = "#444C5C",
                    "primary_re_seeding" = "#E1B16A",
                    "metastatic_confined" = "#D9A3A1",
                    "metastatic_mono_seeding" = "#CE5A57",
                    "metastatic_parallel_seeding" = "#80120D",
                    "metastatic_re_seeding" = "#F49564")

topology_order <- names(topology_color)

topology_labels <- c("absent" = "Absent",
                     "primary_confined" = "Primary\nConfined",
                     "primary_mono_seeding" = "Primary\nMono-Seeding",
                     "primary_parallel_seeding" = "Primary\nParallel\nSeeding",
                     "primary_re_seeding" = "Primary\nRe-Seeding",
                     "metastatic_confined" = "Metastatic\nConfined",
                     "metastatic_mono_seeding" = "Metastatic\nMono-Seeding",
                     "metastatic_parallel_seeding" = "Metastatic\nParallel\nSeeding",
                     "metastatic_re_seeding" = "Metastatic\nRe-Seeding")

cp <- unique(machina_seeding_topology_df$CP)

# Process data for each clonal population (CP)
machina_per_cp_df <- machina_seeding_topology_df %>%
  pivot_longer(cols=2:dim(.)[2], names_to = "topology", values_to = "freq") %>%
  mutate(CP = as.factor(CP)) %>% 
  mutate(topology = as.factor(topology)) %>% 
  group_by(CP) %>% 
  mutate(total_freq = sum(freq)) %>% 
  mutate(perc = round(freq/total_freq*100, 1)) %>% 
  mutate(perc_absent = 100-perc) %>%
  ungroup() %>% 
  dplyr::filter(CP %in% cp)

# Process data for all mice
machina_per_mouse_df <- machina_per_cp_df %>% 
  dplyr::select(!CP) %>% 
  group_by(topology) %>%
  summarise(freq = sum(freq)) %>%
  mutate(total_freq = sum(freq)) %>% 
  mutate(perc = round(freq/total_freq*100, 1)) %>% 
  ungroup() %>% 
  add_column(.before="topology", CP="all") %>%
  mutate(perc_absent = 100-perc) %>%
  mutate(CP = as.factor(CP)) %>% 
  mutate(topology = fct_relevel(topology, topology_order)) %>% 
  arrange(factor(topology, levels = topology_order))

# Combine data for visualization
machina_seeding_topology_combined <- rbind(machina_per_cp_df, machina_per_mouse_df) %>% 
  mutate(bin_topology = ifelse(freq >= 1, 1, 0)) %>%
  mutate(bin_name_topology = ifelse(freq >= 1, paste0(topology), "absent"))

# Summarize rows present/absent
machina_seeding_per_topology_bin <- machina_seeding_topology_combined %>% 
  dplyr::filter(!CP == "all") %>%
  group_by(topology) %>% 
  summarise(freq = sum(bin_topology)) %>% 
  mutate(total_freq = length(unique(machina_per_cp_df$CP))) %>% 
  mutate(perc = round(freq/total_freq*100, 1)) %>% 
  mutate(perc_absent = 100-perc) %>%
  arrange(-perc) %>% 
  mutate(topology=factor(topology, levels=topology))

topology_sort <- as.vector(machina_seeding_per_topology_bin$topology)

# Summarize columns present/absent
machina_seeding_topology_per_cp_bin <- machina_seeding_topology_combined %>% 
  dplyr::filter(!CP == "all") %>%
  group_by(CP) %>% 
  summarise(freq = sum(bin_topology)) %>% 
  arrange(-freq) %>% 
  mutate(CP=factor(CP, levels=CP))

cp_ord <- as.vector(machina_seeding_topology_per_cp_bin$CP)

# Reorder and relevel based on the number of topologies per CPs
machina_seeding_topology_combined <- machina_seeding_topology_combined %>%
  dplyr::mutate(CP = fct_relevel(CP, cp_ord))

# Prepare data for pie charts
pie_per_all_cp_df <- rbind(machina_per_mouse_df, machina_per_cp_df) %>% 
  dplyr::select(CP, topology, perc, perc_absent) %>%
  pivot_longer(cols= c(perc, perc_absent), names_to = "perc_pie") %>%
  left_join(x=., y=enframe(topology_color, name = "topology", value = "color"), by="topology") %>% 
  left_join(x=., y=enframe(topology_labels, name = "topology", value = "labels"), "topology" ) %>% 
  mutate(color = ifelse(perc_pie == "perc_absent", "#DDDDDD", color)) %>% 
  mutate_if(is.character, as.factor) %>% 
  mutate(topology = fct_relevel(topology, topology_order[-1])) %>% 
  mutate(labels = fct_relevel(labels, topology_labels[-1])) %>%
  arrange(CP, topology) %>% 
  mutate(color = as.character(color))

# Draw pie charts
for (cp in unique(pie_per_all_cp_df$CP)) {
  pie_per_cp_df <- dplyr::filter(pie_per_all_cp_df, CP == cp)
  
  pie_per_cp_chart <- ggplot(data=pie_per_cp_df, aes(x = 2, y = value)) +
    geom_bar(fill=pie_per_cp_df$color, width = 0.4, stat = "identity", size=0.4, col="white") +
    geom_text(data = dplyr::filter(pie_per_cp_df, perc_pie == "perc"),
              aes(x=0.25, y=0, label = paste0(dplyr::filter(pie_per_cp_df, perc_pie == "perc")$value, " %")),
              col=dplyr::filter(pie_per_cp_df, perc_pie == "perc")$color, size=6.5, fontface="bold") +
    coord_polar(theta="y", start = 0, direction = -1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), limit=c(0.25, 2.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0, 100)) +
    facet_wrap(~labels, ncol=length(unique(pie_per_cp_df$topology))) +
    labs(x = cp) +
    guides(fill = "none") +
    barplot_nowaklab_theme() +
    theme(panel.spacing.x = unit(-12.5, "mm" ),
          panel.spacing.y = unit(-15, "mm" ),
          plot.margin = unit(c(-10, -5, -15, -5), "mm"),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(vjust = -5),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(vjust = -0.05, face="bold", size=12),
          strip.background.x = element_blank(),
          legend.position="none")
  
  ggsave(plot = pie_per_cp_chart,
         filename = file.path(paste0(machina_seeding_topology_dir, "/", "seed_topology_pie_per_cp", "_", cp, ".pdf")),
         height = 5.5, width = 25, units = "cm")
}

# Prepare data for tiles plot
tiles_topology_perc_per_cp <- ggplot(data=machina_seeding_topology_combined %>% dplyr::filter(!CP == "all"), aes(x=CP, y=topology, fill=bin_name_topology)) +
  geom_rtile(color = "white", lwd = 0.5, radius = unit(7.5, "mm")) +
  scale_fill_manual(values=topology_color[c("absent", topology_sort)]) +
  scale_y_discrete(limits = rev(topology_sort), labels = topology_labels_line, expand=c(0, 0)) +
  labs(x="Clonal Population (CP)") +
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=12),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(margin=margin(-3, 0, 0, 0), size=8, angle=30, hjust=1, vjust=1),
        axis.text.y = element_text(size=12, angle=0, hjust=1.0, vjust=0.5),
        axis.title.y = element_blank(),
        legend.position="none", 
        panel.background = element_rect(fill="white"))

# Prepare data for rows bargraph
machina_per_topology_bin_stat_bargraph_row <- dplyr::filter(machina_seeding_per_topology_bin) %>% 
  ggplot(aes(x = topology, y = perc, fill = topology)) +
  geom_chicklet(radius = grid::unit(1, "mm"), width = 0.8, color = "white", size=0.5) +
  scale_fill_manual(values=topology_color[topology_order]) +
  scale_x_discrete(limits = topology_sort, expand=c(0, 0)) +
  scale_y_continuous(labels=function(x) paste0(x, "%"), breaks=c(0, 25, 50, 75, 100), limits = c(-2.5, 125), expand = c(0, 0)) +
  geom_text(aes(x=topology, y=perc, label=paste0(round(perc, digits=0), "%"), color=topology), size=4, vjust = "middle", hjust=-0.25, fontface = "plain") +
  scale_color_manual(values=topology_color[topology_order]) +
  labs(y = "Topologies Presence\nin All CPs (%)") +
  barplot_nowaklab_theme() +
  lemon::coord_capped_flip(bottom="both") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x= element_text(margin=margin(0, 0, 0, 0), size=10, angle=0, hjust=0.5, vjust=0.5),
        legend.position = "none", 
        panel.background = element_rect(fill=NA))

# Prepare data for columns lollipop
machina_per_cp_bin_stat_lollipop_col <- ggplot(machina_seeding_topology_per_cp_bin, aes(y=CP, x=freq)) +
  geom_segment(aes(y=CP, yend=CP, x=0, xend=freq), color="grey", linewidth = 0.25) +
  geom_point(color="#59C74C", size=3) +
  scale_x_continuous(breaks=seq(0, max(machina_seeding_topology_per_cp_bin$freq)), limits = c(0, max(machina_seeding_topology_per_cp_bin$freq)+1), expand = c(0, 0)) +
  labs(x = "Topologies\nPresence\nin Single CP") +
  lemon::coord_capped_flip(left="both") +
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.y= element_text(margin=margin(0, 0, 0, 0), size=10, angle=0, hjust=0.5, vjust=0.5),
        legend.position="none", 
        panel.background = element_rect(fill="white"))

# Combine plots
machina_per_topology_bin_stat_tiles_bargraph_row <- aplot::insert_right(tiles_topology_perc_per_cp, machina_per_topology_bin_stat_bargraph_row, 0.2)
machina_per_topology_bin_stat_tiles_bargraph_row_col <- aplot::insert_top(machina_per_topology_bin_stat_tiles_bargraph_row, machina_per_cp_bin_stat_lollipop_col, 0.6)
