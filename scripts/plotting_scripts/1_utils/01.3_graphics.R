###### Color Choices #######
### set colors for number of affected sites ###
indels_unique_char_ts_color <- c("single_affected_ts_ins_per_asv" = "#ce5a57",
                                 "single_affected_ts_del_per_asv" = "#8fcef2", "multi_affected_ts_del_per_asv" = "#4486B6", 
                                 "total_affected_ts_per_asv" = "#e1b16a", "total_unique_char_per_asv" = "#78a5a3")
  

indels_site_color <- c("1xTS affected, deletion" = "#C6DBEF", "2xTS affected, deletion" = "#6BAED6", "3xTS affected, deletion" = "#2171B5", ">= 4xTS affected, deletion" = "#08306B",
                       "1xTS affected, insertion" = "#F29E9C", "2xTS affected, insertion" = "#E47D72", "3xTS affected, insertion" = "#EC3323", ">= 4xTS affected, insertion" = "#912D20")

# indels_site_color <- c("1xTS affected, deletion" = "#BECCFB", "2xTS affected, deletion" = "#7D9AF8", "3xTS affected, deletion" = "#3D69F6", ">= 4xTS affected, deletion" = "#2846A4",
#                        "1xTS affected, insertion" = "#E4ACAB", "2xTS affected, insertion" = "#CF5E5A", "3xTS affected, insertion" = "#BD271A", ">= 4xTS affected, insertion" = "#7D160D")



# assign names
indels_site_names <- names(indels_site_color)

# Colors: set final colors for organs
sample_color <- c(## Prostate
                  "PRL" =  "#65A7F3", # Prostate Left
                  "PRL1" = "#65A7F3", # Prostate Left 1 (Clone 1)
                  "PRL2" = "#A6C9FE", # Prostate Left 2 (Clone 2) 
                  "PRL3" = "#E3EDFF", # Prostate Left 3 (Clone 3) 
                  "PRL4" = "#81B8FF", # Prostate Left 4 (Clone 4) 
                  "PRL5" = "#C6DBFF", # Prostate Left 5 (Clone 5)  
                  "PRR" =  "#65A7F3", # Prostate Right
                  ## Seminal Vesicle
                  "SVL" =  "#ADD8E6", # Seminal Vesicle Left
                  "SVR" =  "#ADD8E6", # Seminal Vesicle Right
                  ## Bladder
                  "BDR" =  "mediumpurple4", # Bladder 
                  "BDR1b" ="mediumpurple4", # Bladder (Clone 1b)
                  "BDR2" = "mediumorchid3", 
                  "BDR3" = "mediumpurple3", 
                  "BDR4" = "mediumorchid4",
                  ## Liver
                  "LVL" =  "#BA361E", # Liver Left Lobe 
                  "LVL1" = "#BA361E", # Liver Left Lobe 1
                  "LVL2" = "#FFBBE1", # Liver Left Lobe 2 
                  "LVR" =  "#DD356E", # Liver Right Lobe
                  "LVM" =  "#FC7FB6", # Liver Median Lobe, 
                  "LVC" =  "#FFBBE1", # Liver Caudate Lobe   
                  ## Lungs
                  "LGL" = "#95194f", # Lung Left
                  "LGR" = "#4D0D29", # Lung Right
                  "LGT" = "#71133c", # Lung Trachea 
                  ## Adrenal Cortex
                  "ADL" = "#4387B7", # Adrenal Left
                  "ADR" = "#89AECC", # Adrenal Right
                  ## Kidneys
                  "KDL" = "#E7943B", # Kidney Left
                  "KDR" = "#E7943B", # Kidney Right
                  ## Testis
                  "TSL" = "#EE82EE", # Testis Left
                  "TSR" = "#EE82EE", # Testis Right
                  ## Brain
                  "BRN" = "#F1D351", # Brain
                  ## Blood
                  "BLD" = "#DB8C9B", # Blood
                  ## Lymphatic System 
                  "SPN" = "#1A3F13", # Spleen
                  "THM" = "#e0ebe1", # Thymus
                  ## Lymph Nodes
                  "LNSL" = "#4D684B", # Lymph Node Superficial Cervical Left
                  "LNSR" = "#669F6C", # Lymph Node Superficial Cervical Right
                  "LNDL" = "#BEBE00", # Lymph Node Deep Cervical Left
                  "LNDR" = "#7E7E26", # Lymph Node Deep Cervical Right
                  "LNML" = "#6CD24B", # Lymph Node Mediastinal Left 
                  "LNMR" = "#D4EEC6", # Lymph Node Mediastinal Right
                  "LNIL" = "#8DD3C7", # Lymph Node Lumbar Left
                  "LNIR" = "#BBE0D9", # Lymph Node Lumbar Right
                  "LNLL" = "#3F3E0E", # Lymph Node Inguinal Left
                  "LNLR" = "#7F7D62", # Lymph Node Inguinal Right
                  
                  "LNRR"= "#3d5f40", "LNRL" = "#3d5f40", # Lymph Node, Lymph Node
                  "LNAL"= "#8DD3C7", "LNAR" = "#8DD3C7", # Lymph Node, Lymph Node
                  "LNBL"= "#BEBE00", "LNBR" = "#BEBE00", # Lymph Node, Lymph Node
                  "LNCA"= "#e0ebe1", # Lymph Node Caudal
                  ## Bones
                  # Upper Limbs
                  "SCL" = "#DEB887", # Scapula Left 
                  "SCR" = "#DEB887", # Scapula Right
                  "HML" = "#8B4513", "HML1" = "#8B4513", "HML2" = "#8B4513", # Humerus Left
                  "HMR" = "#F4A460", "HMR1" = "#F4A460", "HMR2" = "#F4A460", # Humerus Right
                  
                  "UNL" = "#875151", # Ulna Left
                  "UNL1"= "#875151", # Ulna Left 
                  "UNL2"= "#e4d2d2", # Ulna Left 

                  "UNR" = "#BC8F8F", # Ulna Right
                  "UNR1"= "#BC8F8F", # Ulna Right
                  "UNR2"= "#eee3e3", # Ulna Right
                  # Lower Limbs
                  "TBL" = "#DEB887", # Tibia Left 
                  "TBR" = "#DEB887", # Tibia Right
                  "FML" = "#8B4513", # Femur Left
                  "FMR" = "#8B4513", "FMR1" = "#8B4513", "FMR2" = "#8B4513", # Femur Right
                  "PVL" = "#F4A460", # Pelvis Left
                  "PVR" = "#F4A460", # Pelvis Right
                  # Ribs
                  "RBL" = "#E7B985", # Rib Left 
                  "RBR" = "#E7B985", # Rib Right
                  ## Spine
                  "SPI" = "#54585E", # spine temp
                  "SPC" = "#54585E", "SPC1" = "#54585E", "SPC2" = "#54585E", # Spine Cervical
                  #"SPT" =
                  #"SPL" =
                  #"SPS" =
                  ## Colon 
                  "CLT" = "#093CF5", "CLT1" = "#093CF5", "CLT2" = lighten("#093CF5", 0.2, space = "HCL"), # Colon Transverse (CLT) 
                  # Small Intestine
                  "SMI" = "#031E8D", "SMI1" = "#031E8D", "SMI2" = lighten("#031E8D", 0.2, space = "HCL"), # Small Intenstine (CLT)
                  ## Days
                  "D03" = "#78a5a3", "D07" = "#ce5a57", "D14" = "#e1b16a", "D28" = "#4486B6", # Days 
                  ## Genotypes
                  "WT" = "#78a5a3", # wild type genotype
                  "2P" = "#e1b16a", # Pten/Trp53-loss genotype
                  "2PR" = "#ce5a57", # Pten/Trp53/Rb1-loss genotype
                  "2PR" = "#54585E", # Pten/Trp53/Smad4-loss genotype 
                  "2PL" = "#4486B6") # # Pten/Trp53/Lrrc15-loss genotype 




# Color Option 1 
# create multiple colors
colony_25xcols <- c("#b2df8a", "#1F78C8", "#ff0000", "#ff7f00", "#36648B", "#FFD700", "#FB6496", "#a6cee3", "#33a02c", "#CAB2D6", 
                    "#FDBF6F", "#999999", "#EEE685", "#C8308C", "#FF83FA", "#C814FA", "#0000FF", "#6A33C2", "#00E2E5", "#00FF00", 
                    "#778B00", "#BEBE00", "#8B3B00", "#A52A3C", "#D9D9D9FF")
# set order of colors
colony_25xcols[-which(colony_25xcols %in% "#D9D9D9FF")]
# add greseqtab_df for cut-off
colony_col_max <- c("#D9D9D9FF", colony_25xcols)
#grid::grid.raster(colony_25xcols, interpolate = FALSE)

# color Option for multiple colors
#set3_col <- brewer.pal(n = 12, name = "Set3")
set3_12xcols <- c("#80B1D3", "#FDB462", "#8DD3C7", "#FCCDE5", "#FB8072", "#FFFFB3", "#B3DE69", "#BEBADA", "#FFED6F", "#BC80BD", "#CCEBC5", "#D9D9D9")



###### Themes for ggplot2 ######
# nowaklab theme for ggplot2
barplot_nowaklab_theme <- function(axis.title.font = "Helvetica", axis.title.col = "black",
                                   axis.text.font = "Helvetica", axis.text.col = "black",  
                                   legend.text.font = "Helvetica")
{
  # General settings
  theme(
    
    # plot
    plot.margin = unit(c(1, 2, 1, 2), "mm"),    
    plot.background = element_blank(), # alt. plot.background = element_rect(colour = NA, fill = "transparent"),
    plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5, lineheight = 0.9),
    plot.subtitle = element_text(color = "black", size = 10, face = "plain", hjust = 0.5, lineheight = 0.9),
    
    # axis lines
    #axis.line = element_line(colour = "black"),
    axis.line = element_line(), # for 
    # Tick axis x and y axes
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(1.5, "mm"),
    
    # X axis text
    # X text straight
    axis.text.x  = element_text(size = 10, color = axis.text.col, angle = 0, vjust = 0, hjust= 0.5, family = axis.text.font),
    axis.title.x = element_text(size = 12, color = axis.title.col, face = "bold", margin = unit(c(3, 3, 3, 3), "mm"), family = axis.title.font),
    # X text angled
    #axis.text.x  = element_text(color="black", angle=55, vjust=1, hjust=1, size=14),
    #axis.title.x = element_text(size=14, color="black", face = "plain"),
    
    # Y axis text  
    axis.text.y = element_text(size = 10, color = axis.text.col, angle = 0, vjust = 0.5, hjust= 1, family = axis.text.font),
    axis.title.y = element_text(size = 12, color = axis.title.col , face = "bold", margin = unit(c(3, 3, 3, 3), "mm"), family = axis.title.font),    
    
    # Legend   
    legend.title = element_text(size = 12, face = "bold"),
    legend.title.align = c(0),
    legend.text = element_text(size = 12, color="black", face = "plain", hjust= 1, margin = margin(l=0.1, r=0.2, unit="cm"), family = legend.text.font),
    legend.text.align = c(0),
    legend.position = "right", 
    legend.box = "vertical",
    legend.background = element_rect(colour = "transparent", fill="transparent"),
    legend.key = element_rect(colour = "transparent", fill="transparent"),
    #legend.key.height = unit(0.5, "cm"),
    #legend.key.width = unit(0.5, "cm"),
    #legend.key.size = unit(0.5, "cm"),
    #legend.spacing = unit(0.1, "cm"),
    #legend.spacing.x = unit(0.1, "cm"), # space between key and text in the legend
    #legend.spacing.y = unit(0.2,"cm"), # space between different legends
    #legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    #legend.box.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    
    # Guides
    
    # Strip
    #strip.background=element_blank(),
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text.x = element_text(size = 12, color="black", face = "bold"),
    strip.switch.pad.grid=unit(0, "cm"),
    strip.switch.pad.wrap=unit(0, "cm"),
    
    # panels: elements control the appearance of the plotting panels
    #panel.background=element_rect(colour = NA, fill="transparent", size = NA, linetype = NA), # background of panel
    panel.background=element_blank(),  # background of panel
    ## panel grid
    #panel.grid = element_blank(),
    #panel.grid.major = element_line(size = 0.5, linetype = "dashed", colour = "white"), 
    #panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "white"),    
    ## X and Y Major | alternative: panel.grid.major.y = element_line(size = 0.5, linetype = "dotted", colour = "#999999"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    ## X and Y Minor | alternative: panel.grid.major.y = element_line(size = 0.25, linetype = "dotted", colour = "#999999"),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
    #panel.spacing.y = unit(-0.5, "lines"),
    #panel.spacing.x = unit(-0.5, "lines")
  )
  
}

### tuniec ###

