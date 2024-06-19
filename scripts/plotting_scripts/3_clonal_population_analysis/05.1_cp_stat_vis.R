###### Description: Script to Analyze Clonal Populations (CPs) Defined by Cassiopeia Based on the Truncal Mutation ###### 
source("~/snakemake_evotracer_machina/scripts/plotting_scripts/1_utils/01.1_libs.R")
source("~/snakemake_evotracer_machina/scripts/plotting_scripts/1_utils/01.2_own_funct_softw.R")
source("~/snakemake_evotracer_machina/scripts/plotting_scripts/1_utils/01.3_graphics.R")


args = commandArgs(trailingOnly=TRUE)
load(args[1])
output_dir <- args[2]
sample_order <- EvoTraceR_object$sample_order

###### Download Data (from EvoTraceR_object) ######
## take df_to_plot_final 
df_to_plot_final <- 
  EvoTraceR_object$plot_summary$df_to_plot_final %>%
  dplyr::mutate(sample = fct_relevel(sample, sample_order)) %>% 
  dplyr::rename(cp = group) ## rename to cp (clonal population)


###### Data Outputs ######
# output dir: for graphs analysis
graphs_analysis_dir <- paste0(output_dir, "/graphs_analysis")
if (!dir.exists(graphs_analysis_dir)) 
{dir.create(graphs_analysis_dir, recursive = TRUE)}


###### Set Parameters ###### 
## name of Clonal Population that is: "cp_single_asv" = "CP_NMBC" + "CP_SASV", usually "CP00"
# "CP_NMBC" - "Non-Marked BarCodes": clonal population; possibly not related; not efficient editing in vivo
# "CP_SASV" - "Single ASV": clonal population with only one ASV; possibly extinct ASVs and/or never expanded and edited
# "CP_MASV" - "Multiple ASV": clonal population with at least two ASV or more; used for multiple analysis  

## find "CP00" (or higher number of 0, "CP00" = "CP_NMBC" + "CP_SASV")
# % are being calculated based on all three CPs: "CP_MASV", "CP_NMBC", "CP_SASV"
cp_single_asv <- data.frame(cp=df_to_plot_final$cp, cp_digit = gsub("[^0-9.-]", "", df_to_plot_final$cp)) %>% 
  filter(cp_digit == min(cp_digit)) %>% 
  unique() %>% 
  dplyr::select(cp) %>% 
  deframe() %>%
  as.vector()


###### Calculations ######
## adjust input data df_long (used for dispersal calculations)
## data used for dispersal calculations
# asv_names - amplicon sequence variants names
# sample - sample where ASV has been found
# cp - clonal population, defined based on founder truncal mutation
# count - the counts of an ASV in a given sample normalized by the total number of reads in that sample then multiplied by a scale factor
## data frame

df_long <- 
  df_to_plot_final %>% 
  dplyr::select(c("cp", "sample", "asv_names", "count")) %>%
  unique() %>% 
  mutate(cp = ifelse(asv_names == "BC10v0", "CP_NMBC", 
                     ifelse(cp == cp_single_asv, "CP_SASV", as.character(cp)))) %>% ## assign general group
  mutate(cp_general = ifelse(cp == "CP_NMBC" | cp == "CP_SASV", cp, "CP_MASV")) %>% ## assign general group
  dplyr::select(c("cp", "cp_general", "sample", "asv_names", "count")) %>%
  mutate(cp = as.factor(cp)) %>% 
  mutate(cp_general = as.factor(cp_general))

       
###### (0) calculate "perc summary" per "sample and cp" (clonal populations) ######
df_long_per_cp_general <- 
  df_long %>%
  dplyr::select(cp_general, count) %>%
  group_by(cp_general) %>%
  summarize(sum_per_cp_general = sum(count)) %>% 
  ungroup() %>% 
  mutate(perc_per_cp_general = 100*round(sum_per_cp_general/sum(sum_per_cp_general), 4))
  

###### (1) calculate perc summary per "cp & sample" (clonal populations) ######
## used for tiles graphing; distribution of organs in cp (count_perc_per_sample_cp)
## sample & cp is grouped
df_long_perc_per_sample_cp <- 
  df_long %>% 
  dplyr::select(cp, sample, count) %>%
  # summarize counts of ASVs in specific CP per tissue
  group_by(cp, sample) %>% 
  summarize(count_per_sample_cp = sum(count)) %>%
  mutate(count_per_sample_cp_log10 = log(count_per_sample_cp, 10)) %>% 
  # summarize log 10 total counts of ASVs in specific CP
  group_by(cp) %>%
  mutate(sum_per_sample_cp_log10 = sum(count_per_sample_cp_log10)) %>%
  ungroup() %>% 
  # summarize total counts of ASVs in specific CP
  group_by(cp) %>%
  mutate(sum_per_sample_cp = sum(count_per_sample_cp)) %>%
  ungroup() %>% 
  # percentage of each ASVs
  mutate(count_perc_per_sample_cp = round((count_per_sample_cp/sum_per_sample_cp*100), 3)) %>%
  #
  mutate(sum_per_sample_cp_max_perc = round((count_per_sample_cp/sum(count_per_sample_cp)*100), 3)) %>% 
  #
  group_by(cp) %>% 
  mutate(total_per_sample_cp_max_perc=sum(sum_per_sample_cp_max_perc)) %>%
  #
  group_by(cp) %>% 
  mutate(label = ifelse(total_per_sample_cp_max_perc <= 1, 1, total_per_sample_cp_max_perc)) %>% 
  mutate(label = ifelse(duplicated(total_per_sample_cp_max_perc), NA, total_per_sample_cp_max_perc)) %>%  
  mutate(label = round(label, 1)) %>% 
  mutate(label = ifelse(label < 1, 1, label)) %>% 
  mutate(label = as.character(label)) %>%
  mutate(label = as.character(paste0(label, "%"))) %>%
  mutate(label = ifelse(label == "NA%", "", label)) %>% 
  mutate(label = ifelse(label == "1%", "", label)) %>% 
  arrange(cp) %>% 
  droplevels()
## skip clone "CP_NMBC" and clone "CP_SASV"
df_long_perc_per_sample_cp <-
  df_long_perc_per_sample_cp %>% 
  dplyr::filter(!cp == "CP_NMBC") %>% 
  dplyr::filter(!cp == "CP_SASV")

###### (2) calculate summary per "cp & asv_names" ######
## used for 
df_long_summ_per_asv_cp <- 
  df_long %>%
  dplyr::select(asv_names, sample, cp, count) %>%
  group_by(cp) %>% # asv_names, 
  mutate(count_sum_per_cp = sum(count)) %>%
  group_by(asv_names) %>%
  mutate(count_sum_per_asv = sum(count)) %>%
  mutate(count_perc = round(count_sum_per_asv/count_sum_per_cp, 2)) %>% 
  ungroup() %>% 
  arrange(cp, asv_names)

###### (3) calculate alpha diversity per "cp": 1) shannons_index_percluster, 2) richness_asv_percluster, 3) pielous_evenness_percluster ######
adiv_per_cp <-
  df_long_summ_per_asv_cp %>%
  dplyr::select(cp, asv_names, count_sum_per_asv) %>%
  unique() %>%
  group_by(cp) %>%
  mutate(richness_per_cp = length(count_sum_per_asv), # S
         shannons_per_cp = vegan::diversity(count_sum_per_asv, index = "shannon", MARGIN = 2), # H′abundance≡H′=−∑Si=1piln
         pielous_per_cp = shannons_per_cp/log(richness_per_cp), # J′=H′/ln(S).
         count_sum_per_asv=count_sum_per_asv) %>% 
  dplyr::select(cp, asv_names, richness_per_cp, shannons_per_cp, pielous_per_cp, count_sum_per_asv) %>% 
  unique() %>% 
  arrange(cp) 

###### (4) calculate alpha diversity per "sample and cp": : 1) shannons_index_percluster, 2) richness_asv_percluster, 3) pielous_evenness_percluster ====== 
df_long_adiv_per_sample_cp <- 
  df_long %>%
  group_by(sample, cp) %>%
  mutate(richness_per_sample_cp = length(cp),
         shannons_per_sample_cp = vegan::diversity(count, index = "shannon", MARGIN = 2),
         pielous_per_sample_cp = shannons_per_sample_cp/log(richness_per_sample_cp)) %>% 
  dplyr::select(sample, cp, richness_per_sample_cp, shannons_per_sample_cp, pielous_per_sample_cp) %>% 
  unique() %>% 
  arrange(cp)


###### (5) calculate dispersal based on df_long ######
# Rational: The level of tissue dispersal is a consequence of metastatic dissemination and thus can inform on the intensity of past metastatic events.
# Goal: understanding how dispersed the counts in each clone are across the different samples.
# Data: to compute this for cluster c, we need to compare the total counts in each sample for clone c versus the total counts in each sample in any other clone.
# Problem: It saturates at intermediate metastatic regimes

## calculate dispersal
# prepare vector for loops
v_list = c()
nasv_list = c()

# loop
for (cluster in unique(df_long$cp)) {
  # counts per clonal population (=cluster)
  count_cluster <- 
    df_long %>% 
    filter(cp == cluster)  
  nasv <- length(count_cluster %>% 
                   pull(asv_names) %>% 
                   unique)
  # summaries c counts per cp and asv 
  count_cluster <- 
    count_cluster %>% 
    dplyr::select(-cp) %>%
    group_by(sample) %>% 
    summarise(c = sum(count)) 
  # count not in cluster (summ all) 
  count_notcluster <- 
    df_long %>% 
    filter(cp != cluster) %>% 
    dplyr::select(-cp) %>% 
    group_by(sample) %>% 
    summarise(notc = sum(count)) 
  # data frame with c notc for each sample  
  count_final <- 
    dplyr::full_join(count_cluster, count_notcluster, by = 'sample') %>%
    replace(is.na(.), 0) %>% 
    column_to_rownames('sample')
  total_counts <- colSums(count_final)
  scale_factor <- total_counts['c'] / total_counts['notc']
  count_final$notc <- count_final$notc * scale_factor
  colSums(count_final)
  
  # chi-square test assesses whether the distribution of counts in cluster C is different than the one in the background (Not C).  
  chisq <- chisq.test(count_final)
  
  cramers_v = function(phi, N, k, r) {
    phi2 = phi/N
    
    phi2corr = max(c(0, phi2 - ((k-1) * (r-1)) / (N-1)))
    rcorr = r - ((r-1)**2)/(N-1)
    kcorr = k - ((k-1)**2)/(N-1)
    
    v = sqrt( phi2corr / min( (kcorr - 1), (rcorr - 1) ))
    return(v)
  }
  
  # from the chi-square test we obtain the chi-square statistic, which is then used to compute the tissue dispersal score.
  # tissue dispersal score: (inverted Cramer’s V == v) ranges from 0 to 1 and indicates how dispersed the counts in a cluster are across tissues.
  # when it is 0 it indicates that there is no deviation from the background, when it’s 1 it means that the distribution across samples in cluster C completely deviates from the background.  
  v = cramers_v(phi = chisq$statistic, N = sum(count_final), k = nrow(count_final), r = ncol(count_final))
  v_list = c(v_list, 1-v)
  nasv_list = c(nasv_list, nasv)
}

# combine final dispersal data
dispersal_per_cp <- 
  data.frame(cp = unique(df_long$cp), v_stat = v_list, stringsAsFactors = F) %>%
  tibble()
# final join
df_long_dispersal_adiv_per_cp <- 
  left_join(dispersal_per_cp, adiv_per_cp, by="cp", multiple = "all") %>%
  dplyr::filter(!cp == "CP_NMBC") %>% 
  dplyr::filter(!cp == "CP_SASV") %>% 
  dplyr::select(cp, v_stat, richness_per_cp, shannons_per_cp, pielous_per_cp) %>% 
  unique() %>% 
  arrange(cp)


###### Combine All Data ######
# df_long_perc_per_sample_cp
# df_long_dispersal_adiv_per_cp


###### Visualizations Per Clonal Population (CP) ######


###### (1) Richness_per_cp: number of ASVs per CP, represented as a bar plot ######
### plot: parameters
## automatic scaling of y axis
## maximum value for y
y_max <- 
  df_long_dispersal_adiv_per_cp %>% 
  pull(richness_per_cp) %>% 
  max()
## y axis ticks:
y_axis_ticks <- ifelse(y_max <= 10, 2, 
                       ifelse(y_max > 10 & y_max <= 25, 5, 
                              ifelse(y_max > 25 & y_max <= 50, 10,
                                     ifelse(y_max > 50 & y_max <= 100, 20,
                                            ifelse(y_max > 100 & y_max <= 250, 50,
                                                   ifelse(y_max > 250 & y_max <= 500, 100,
                                                          ifelse(y_max > 500, 200)))))))
## maximum value for y
y_max <- plyr::round_any(df_long_dispersal_adiv_per_cp %>% 
                           pull(richness_per_cp) %>% 
                           max(), y_axis_ticks, f = ceiling)
## limit is what axis go but not seen 
y_min_limit <- -0.25 * y_max
y_max_limit <-  1.15 * y_max

## plot
bars_richness_per_cp <-
  ggplot(df_long_dispersal_adiv_per_cp, aes(x=cp, y=richness_per_cp)) +
  geom_bar(stat = "unique", width=0.9, fill="#ce5a57") +
  geom_text(aes(label=richness_per_cp), position=position_dodge(width=0.9), size=5, fontface="bold", vjust=-0.25) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(expand = c(0, 0), limits= c(y_min_limit, y_max_limit), breaks = seq(0, y_max, y_axis_ticks)) +
  labs(x="Clonal Populations (CPs)", y="ASVs Richness") +
  lemon::coord_capped_cart(left="both") +
  # geom_hline(yintercept=median(df_long_dispersal_adiv_per_cp$richness_per_cp), linetype="dashed", color = "red") +
  geom_hline(yintercept=0, colour="black", size=0.3) +
# add theme
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "mm"),
        axis.line.y = element_line(colour="black", size=0.3), # axis y line only
        axis.line.x.top = element_line(colour="black", size=0.3), # axis x line only
        axis.line.x.bottom = element_line(colour="black", size=0.3), # axis x line only
        panel.border = element_blank(), # disable panel border
        panel.grid.major.x = element_line(size = 0.25, linetype = "dotted", colour = "#999999"),
        panel.grid.minor = element_blank(), # disable lines in grid on X-axis
        axis.title.y = element_text(size=16),
        axis.title.x= element_text(size=16),
        axis.text.y = element_text(size=16, angle=0, hjust=1, vjust=0.5),
        axis.text.x = element_text(size=12, angle=45, hjust=0, vjust=0),
        axis.ticks.x = element_line(colour="black", size=0.3),
        axis.ticks.y = element_line(colour="black", size=0.3),
        strip.background=element_blank(), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))


###### (2) ASVs Total Counts for CPs - bar ######
## ASVs counts in all CPs
# this graph shows the counts distribution of sites and the dominance of the major clones
# this graphs show wells smaller CPs with less counts as it is is in the log10 scale
# set minimum/maximum for graph
min_ticks_y <- 10^0 
max_ticks_y <- 10 * max(df_long_perc_per_sample_cp$count_per_sample_cp)

# add +1 for dealing with 0
df_long_perc_per_sample_cp$count_per_sample_cp <- df_long_perc_per_sample_cp$count_per_sample_cp+.1 
## plot
bars_count_per_sample_cp <- 
  ggplot(data = df_long_perc_per_sample_cp, aes(x = cp, y = count_per_sample_cp, fill = sample)) +
  ## stacked
  # geom_bar(stat = "identity", position = "stack", width=0.9, color = "white", size=0.1) + # "white" to have lines around 
  ## next to each other
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"), width=0.9, color = "white", size=0.1) + # "white" to have lines around 
  # geom_text(aes(y = sum_per_sample_cp, label=label), vjust=-0.35, fontface="bold") +
  scale_fill_manual(values = sample_color[sample_order]) +
  #scale_x_discrete(expand = c(0, 0)) +
  scale_y_log10(expand = c(0, 0),
                limits=c(min_ticks_y, max_ticks_y),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(x=" ", y="ASVs Total Count", fill = "Samples (Organs)") +
  lemon::coord_capped_cart(left="both") +
  # add theme
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "mm"),
        axis.line.y = element_line(colour="black", size=0.3), # axis y line only
        axis.line.x = element_line(colour="black", size=0.3), # axis x line only
        panel.border = element_blank(), # disable panel border
        panel.grid.major.x = element_line(size = 0.25, linetype = "dotted", colour = "#999999"),
        panel.grid.minor = element_blank(), # disable lines in grid on X-axis
        axis.title.y = element_text(size=16),
        # axis.title.x = element_blank(),
        axis.text.y = element_text(size=16, angle=0, hjust=1, vjust=0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour="black", size=0.3),
        strip.background=element_blank(), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))


###### (3) ASVs Freq. % for CPs ######
## ASVs counts as % of all ASVs in all CPs
## this graph shows the % distribution of sites and the dominance of the major clones
## this graphs doesn't show well smaller CPs with less % dominance

### automatic scaling of y axis
## maximum value for y
y_max <- 
  max(df_long_perc_per_sample_cp$total_per_sample_cp_max_perc)

## y axis ticks: max is 100%
y_axis_ticks <- ifelse(y_max <= 5, 1, 
                       ifelse(y_max > 5 & y_max <= 10, 2,  
                              ifelse(y_max > 10 & y_max <= 25, 5, 
                                     ifelse(y_max > 25 & y_max <= 50, 10,
                                            ifelse(y_max > 50 & y_max <= 100, 20)))))
## maximum value for y
y_max <- plyr::round_any(y_max, y_axis_ticks, f = ceiling)

## limit is what axis go but not seen 
# y_min_limit <- -0.25 * y_max
# y_max_limit <-  1.15 * y_max

## plot
bars_sum_per_sample_cp_max_perc <- 
  ggplot(data=df_long_perc_per_sample_cp, aes(x = cp, y = sum_per_sample_cp_max_perc, fill = sample)) +
  geom_bar(stat = "identity", width=0.9, color = "white", size=0.1) + # "white" to have lines around 
  geom_text(aes(y=total_per_sample_cp_max_perc, label=label), vjust=-0.35, size=4, fontface="bold") +
  scale_fill_manual(values = sample_color[sample_order]) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, y_max), breaks=seq(from=0, to=y_max, by=y_axis_ticks)) +
  labs(x=" ", y="ASVs Frequency (%)", fill = "Samples (Organs)") +
  lemon::coord_capped_cart(left="both") +
  # add theme
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "mm"),
        axis.line.y = element_line(colour="black", size=0.3), # axis y line only
        axis.line.x = element_line(colour="black", size=0.3), # axis x line only
        panel.border = element_blank(), # disable panel border
        panel.grid.major.x = element_line(size = 0.25, linetype = "dotted", colour = "#999999"),
        panel.grid.minor = element_blank(), # disable lines in grid on X-axis
        axis.title.y = element_text(size=16),
        #axis.title.x = element_blank(),
        axis.text.y = element_text(size=16, angle=0, hjust=1, vjust=0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour="black", size=0.3),
        strip.background=element_blank(), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))


###### (4) Shannons_per_cp - bar plot ###### 

### automatic scaling of y axis
## maximum value for y
y_max <- 
  df_long_dispersal_adiv_per_cp$shannons_per_cp %>%  
  max()

## y axis ticks: max is 100%
y_axis_ticks <- ifelse(y_max <= 5, 1, 
                       ifelse(y_max > 5 & y_max <= 10, 2,  
                              ifelse(y_max > 10 & y_max <= 25, 5, 
                                     ifelse(y_max > 25 & y_max <= 50, 10,
                                            ifelse(y_max > 50 & y_max <= 100, 20)))))
## maximum value for y
y_max <- plyr::round_any(y_max, y_axis_ticks, f = ceiling)

## limit is what axis go but not seen 
# y_min_limit <- -0.25 * y_max
# y_max_limit <-  1.15 * y_max

## plot
bars_shannons_per_cp <-
  ggplot(df_long_dispersal_adiv_per_cp, aes(x=cp, y=shannons_per_cp)) +
  geom_bar(stat = "unique", width=0.9, fill="#78a5a3") +
  # scale_y_continuous(expand = c(0, 0), limits= c(0, 4), breaks = c(0, 1, 2, 3, 4)) +
  scale_y_continuous(expand = c(0, 0), limits= c(0, y_max), breaks = seq(0, y_max, y_axis_ticks)) +
  labs(x=" ", y="Shannon's Entropy") +
  lemon::coord_capped_cart(left="both") +
  geom_hline(yintercept=median(df_long_dispersal_adiv_per_cp$shannons_per_cp), linetype="dashed", color = "red") +
  # add theme
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "mm"),
        axis.line.y = element_line(colour="black", size=0.3), # axis y line only
        axis.line.x = element_line(colour="black", size=0.3), # axis x line only
        panel.border = element_blank(), # disable panel border
        panel.grid.major.x = element_line(size = 0.25, linetype = "dotted", colour = "#999999"),
        panel.grid.minor = element_blank(), # disable lines in grid on X-axis
        axis.title.y = element_text(size=16),
        # axis.title.x = element_blank(),
        axis.text.y = element_text(size=16, angle=0, hjust=1, vjust=0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour="black", size=0.3),
        strip.background=element_blank(), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))


###### (5) Pielous_per_cp - bar plot ###### 
bars_pielous_per_cp <-
  ggplot(df_long_dispersal_adiv_per_cp, aes(x=cp, y=pielous_per_cp)) +
  geom_bar(stat = "unique", width=0.9, fill="#e1b16a") +
  scale_y_continuous(expand = c(0, 0), limits= c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
  labs(x=" ", y="Pielou's Eveness") +
  lemon::coord_capped_cart(left="both") +
  geom_hline(yintercept=median(df_long_dispersal_adiv_per_cp$pielous_per_cp), linetype="dashed", color = "red") +
# add theme
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "mm"),
        axis.line.y = element_line(colour="black", size=0.3), # axis y line only
        axis.line.x = element_line(colour="black", size=0.3), # axis x line only
        panel.border = element_blank(), # disable panel border
        panel.grid.major.x = element_line(size = 0.25, linetype = "dotted", colour = "#999999"),
        panel.grid.minor = element_blank(), # disable lines in grid on X-axis
        axis.title.y = element_text(size=16),
        # axis.title.x = element_blank(),
        axis.text.y = element_text(size=16, angle=0, hjust=1, vjust=0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour="black", size=0.3),
        strip.background=element_blank(), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))


###### (6) dispersal score - lollipop ###### 
bars_dipsersal_score <-
  ggplot(df_long_dispersal_adiv_per_cp, aes(x=cp, y=v_stat)) +
  #geom_bar(stat = "unique", width=0.9, fill="#e1b16a") +
  geom_segment(aes(x=cp, xend=cp, y=0, yend=v_stat), color="#E3A73F") +
  geom_point(color="#C45872", size=4, shape = 19, fill="white") +  
  scale_y_continuous(expand = c(0, 0), limits= c(0, 1.025), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
  #scale_x_discrete(expand = c(0, 0)) +
  labs(x="Clonal Populations (CPs)", y="Dispersal Score") +
  lemon::coord_capped_cart(left="both") +
  geom_hline(yintercept=median(df_long_dispersal_adiv_per_cp$v_stat), linetype="dashed", color = "red") +
# add theme
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "mm"),
        axis.ticks = element_blank(), # disable ticks lines
        axis.line.y = element_line(colour="black", size=0.3), # axis y line only
        axis.line.x = element_line(colour="black", size=0.3), # axis x line only
        panel.border = element_blank(), # disable panel border
        panel.grid.major.x = element_line(size = 0.25, linetype = "dotted", colour = "#999999"),
        panel.grid.minor = element_blank(), # disable lines in grid on X-axis
        axis.title.y = element_text(size=16),
        axis.title.x= element_text(size=16),
        axis.text.y = element_text(size=16, angle=0, hjust=1, vjust=0.5),
        axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1),
        axis.ticks.x = element_line(colour="black", size=0.3),
        axis.ticks.y = element_line(colour="black", size=0.3),
        strip.background=element_blank(), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))  

  
###### (7) tiles (not used currently) - perc of samples distribution in samples ###### 
tiles_perc_per_sample_cp <-
  ggplot(data=df_long_perc_per_sample_cp, aes(x=cp, y=sample)) +
  geom_tile(aes(fill=count_perc_per_sample_cp), colour = "black") +
  #geom_text(aes(label=count_sample_group_perc), size=1.5, col="white") +
  labs(y="Tissues", fill = "Freq. of Tissues \n per CP (%)") +
  scale_fill_continuous_sequential(palette = "Plasma", limits=c(0, 100), rev = TRUE) +
  scale_y_discrete(expand = c(0, 0), limits = rev(sample_order)) + # change order to have PRL up 
  lemon::coord_capped_cart(left="both") +
  coord_equal() +
  # add theme
  barplot_nowaklab_theme() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "mm"),
        axis.line.y = element_line(colour="black", size=0.3), # axis y line only
        axis.line.x = element_blank(), # axis x line only
        panel.border = element_blank(), # disable panel border
        panel.grid.major = element_blank(), # disable lines in grid on X-axis
        panel.grid.minor = element_blank(), # disable lines in grid on X-axis
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=16, angle=0, hjust=1, vjust=0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour="black", size=0.3),
        strip.background=element_blank(), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))

###### (8) combine graphs:  (bars) + tiles (tiles) ###### 
## aplot parameters
height = 0.95

## combine plots
# bars_richness_per_cp & bars_count_per_sample_cp
plot_1 <- insert_bottom(bars_richness_per_cp, bars_count_per_sample_cp, height = height) 
# bars_sum_per_sample_cp_max_perc
plot_2 <- insert_bottom(plot_1, bars_sum_per_sample_cp_max_perc, height = height) 
# bars_shannons_per_cp
plot_3 <- insert_bottom(plot_2, bars_shannons_per_cp, height = height) 
# bars_pielous_per_cp
plot_4 <- insert_bottom(plot_3, bars_pielous_per_cp, height = height)
# bars_dipsersal_score
plot_5 <- insert_bottom(plot_4, bars_dipsersal_score, height = height)

# # bars_dipsersal_score
# plot_3 <- insert_bottom(plot_2, bars_dipsersal_score, height = height) 
# # tiles_perc_per_sample_cp
# plot_4 <- insert_bottom(plot_3, tiles_perc_per_sample_cp, height = height/1.5) 

# final plot
plot_6 <- print(plot_5) & theme(legend.position = "bottom", legend.box = "horizontal")

## save
ggsave(filename=paste0(graphs_analysis_dir, "/stat_cps_dispersal_bargraph_hm.pdf"), 
       plot=plot_6, 
       width=length(unique(df_long_dispersal_adiv_per_cp$cp))/1, 
       height=42, 
       units = "cm", limitsize = FALSE)


###### (5) WiP: dot plot ###### 
# dispersal_dotgraph <-
#   ggplot(data= df_long_dispersal_adiv_per_cp, aes(x = log(count_sum_per_asv, 10), y=v_stat, size=shannons_per_cp, fill=richness_per_cp)) + 
#   geom_point(shape = 21, colour = "black", stroke = 0.5) +
#   scale_fill_continuous_sequential(palette = "Burg", rev = TRUE) +
#   scale_y_continuous(expand = c(0.01, 0), limits= c(0, 1.1), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) + # v-stat
#   #scale_x_continuous(expand = c(0.01, 0), limits= c(0, 165), breaks = c(0, 40, 80, 120, 160)) + # 
#   #scale_x_continuous(expand = c(0.01, 0), limits= c(0, 1.5), breaks = c(0, 0.5, 1, 1.5)) + # shannon
#   #scale_x_continuous(expand = c(0.01, 0), limits= c(0, 1.1), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) + # piellous
#   #labs(x="Number of ASVs per Clone", y="Dispersal Score") +
#   lemon::coord_capped_cart(left="both", bottom="both") +
#   coord_equal() 
# # add theme
# dispersal_dotgraph <-
#   dispersal_dotgraph +
#   barplot_nowaklab_theme() +
#   theme(aspect.ratio=1/1)
# # save -> dont plot
# ggsave(filename=paste0(graphs_analysis_dir, "/dispersal_dotgraph.pdf"), plot=dispersal_dotgraph, width=20, height=20, units = "cm", limitsize = FALSE)


# ###### Visualizations Per Sample ######
# 
# ## (2)  clones with shannons ======
# chart
#sample_adiv <-
#  ggplot(data=df_long_adiv_per_sample_cp, aes(x = sample, y = shannons_per_sample_cp, fill=sample)) +
#  geom_boxplot(position=position_dodge(0.9), width=0.75, size=0.5, outlier.shape = NA, colour="gray25", outlier.colour = #"gray", outlier.size=1.5) +
#  scale_fill_manual(values = sample_color[sample_order]) +
#  xlab("Tissues") +
#  ylab("Shannon H") +
#  coord_capped_cart(left="both") +
#  barplot_nowaklab_theme() +
#  theme(legend.position="bottom",
#        axis.line.x = element_blank(), # disable y axis lines
#        axis.ticks.x = element_blank(), # disable y axis ticks
#        axis.text.x = element_blank(),
#        strip.text.x = element_text( size = 12, color = "#cc0000", face = "bold.italic")) +
#  theme(axis.text.x = element_text(angle=0, vjust=1, hjust=0.5))

# save
#ggsave(filename=paste0(graphs_analysis_dir, "/stat_sample_cp_shanon.pdf"), plot=sample_adiv, width=length(sample_order)*3, height=12, units = "cm")




### tuniec ###

