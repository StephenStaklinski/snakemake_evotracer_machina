source("scripts/plotting_scripts/1_utils/01.1_libs.R")
source("scripts/plotting_scripts/1_utils/01.2_own_funct_softw.R")
source("scripts/plotting_scripts/1_utils/01.3_graphics.R")

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]
outfile <- args[2]

machina_seeding_topology_df <- read_delim(infile)

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

# Normalize the values to proportions
machina_seeding_topology_df <- machina_seeding_topology_df %>%
  pivot_longer(cols = -CP, names_to = "topology", values_to = "freq") %>%
  mutate(freq = freq / sum(freq) * 100)

# Process data for visualization
machina_seeding_topology_combined <- machina_seeding_topology_df %>%
  mutate(bin_topology = ifelse(freq >= 1, 1, 0)) %>%
  mutate(bin_name_topology = ifelse(freq >= 1, paste0(topology), "absent"))

# Prepare data for pie chart
pie_df <- machina_seeding_topology_combined %>%
  pivot_longer(cols = c(freq), names_to = "perc_pie") %>%
  left_join(x = ., y = enframe(topology_color, name = "topology", value = "color"), by = "topology") %>%
  left_join(x = ., y = enframe(topology_labels, name = "topology", value = "labels"), "topology") %>%
  mutate(color = ifelse(perc_pie == "perc_absent", "#DDDDDD", color)) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(topology = fct_relevel(topology, topology_order[-1])) %>%
  mutate(labels = fct_relevel(labels, topology_labels[-1])) %>%
  arrange(CP, topology) %>%
  mutate(color = as.character(color))

# Draw pie chart
pie_chart <- ggplot(data = pie_df, aes(x = 2, y = value)) +
  geom_bar(aes(x = 2, y = 100), fill = "#F0F0F0", width = 0.4, stat = "identity", size = 0.4, col = "white") +
  geom_bar(fill = pie_df$color, width = 0.4, stat = "identity", size = 0.4, col = "white") +
  geom_text(data = dplyr::filter(pie_df, perc_pie == "freq"),
            aes(x = 0.25, y = 0, label = paste0(round(dplyr::filter(pie_df, perc_pie == "freq")$value, 1), " %")),
            col = dplyr::filter(pie_df, perc_pie == "freq")$color, size = 6.5, fontface = "bold") +
  coord_polar(theta = "y", start = 0, direction = -1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limit = c(0.25, 2.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit = c(0, 100)) +
  facet_wrap(~labels, ncol = length(unique(pie_df$topology))) +
  labs(x = unique(pie_df$CP)) +
  guides(fill = "none") +
  barplot_nowaklab_theme() +
  theme(panel.spacing.x = unit(-12.5, "mm"),
        panel.spacing.y = unit(-15, "mm"),
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
        strip.text.x = element_text(vjust = -0.05, face = "bold", size = 12),
        strip.background.x = element_blank(),
        legend.position = "none")

ggsave(plot = pie_chart,
       filename = outfile,
       height = 5.5, width = 25, units = "cm")
