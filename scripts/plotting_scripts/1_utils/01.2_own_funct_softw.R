###### Software Dependency ###### 

###### Own Functions To Download ######  
###### drawing rectangles with circularized angles (stadium like shapes) ######
### Inspired from: https://stackoverflow.com/questions/64355877/round-corners-in-ggplots-geom-tile-possible ###
## rounded rgtiles() -> geom_tiles
# GeomRtile
`%||%` <- function(a, b) {
  if(is.null(a)) b else a
}
## fuction
GeomRtile <- ggproto("GeomRtile", 
                     statebins:::GeomRrect, # 1) only change compared to ggplot2:::GeomTile
                     
                     extra_params = c("na.rm"),
                     setup_data = function(data, params) {
                       data$width <- data$width %||% params$width %||% resolution(data$x, FALSE)
                       data$height <- data$height %||% params$height %||% resolution(data$y, FALSE)
                       
                       transform(data,
                                 xmin = x - width / 2,  xmax = x + width / 2,  width = NULL,
                                 ymin = y - height / 2, ymax = y + height / 2, height = NULL
                       )
                     },
                     default_aes = aes(
                       fill = "grey20", colour = NA, size = 0.1, linetype = 1,
                       alpha = NA, width = NA, height = NA
                     ),
                     required_aes = c("x", "y"),
                     
                     # These aes columns are created by setup_data(). They need to be listed here so
                     # that GeomRect$handle_na() properly removes any bars that fall outside the defined
                     # limits, not just those for which x and y are outside the limits
                     non_missing_aes = c("xmin", "xmax", "ymin", "ymax"),
                     draw_key = draw_key_polygon
)

geom_rtile <- function(mapping = NULL, data = NULL,
                       stat = "identity", position = "identity",
                       radius = grid::unit(6, "pt"), # 2) add radius argument
                       ...,
                       linejoin = "mitre",
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRtile, # 3) use ggproto object here
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = rlang::list2(
      linejoin = linejoin,
      na.rm = na.rm,
      ...
    )
  )
}



###### Identify Mutations ######
identifyMutations <- function(seq, ref) {
  
  indicesToKeep <- as.vector(as.matrix(ref)!="-")
  seqSubset <- seq[indicesToKeep]
  refSubset <- ref[indicesToKeep]
  # vectors for translating the coordinates
  barcodpToBarcode <- rep(NA, length(seq))
  barcodpToBarcode[which(indicesToKeep)] <- 1:length(which(indicesToKeep))  
  barcodpToBarcode <- nafill(barcodpToBarcode, type = "locf")
  
  barcodeToBarcodep <- rep(NA, length(indicesToKeep))  
  barcodeToBarcodep <-  which(!is.na(barcodpToBarcode))
  
  # barcode + scale
  arle.inf.df <- data.frame(unclass(rle(as.numeric(seq == ref)))) 
  # find start and end coordinates of each segment
  arle.inf.df$end.barcodep <- cumsum(arle.inf.df$lengths)
  if (nrow(arle.inf.df)==1) {
    arle.inf.df$start.barcodep <- 1
  } else {
    arle.inf.df$start.barcodep <- c(0, arle.inf.df$end.barcodep[1:(nrow(arle.inf.df)-1)])+1
  }
  # mismatches mean mutations - subset to mutations
  arle.inf.df <- subset(arle.inf.df, values==0)
  
  # if any mutations found, add their sequence and type
  if (nrow(arle.inf.df)>0) {
    # annotate each mutations with alt and ref alleles, type and sequence ID
    for (mi in 1:nrow(arle.inf.df)) {
      arle.inf.df$ref[mi]<-paste(ref[ arle.inf.df$start.barcodep[mi]:arle.inf.df$end.barcodep[mi]], collapse='')
      arle.inf.df$alt[mi]<-paste(seq[ arle.inf.df$start.barcodep[mi]:arle.inf.df$end.barcodep[mi]], collapse='')
    }
    arle.inf.df$type <- "sub"
    arle.inf.df[arle.inf.df$lengths>1, "type"] <- "complex"
    # if alt sequence is made entirely of hyphens, it is a deletions
    arle.inf.df[gsub("-", "", arle.inf.df$alt)=="", "type"] <- "del"
    arle.inf.df[gsub("-", "", arle.inf.df$ref)=="", "type"] <- "ins"
    # add the barcode scale
    arle.inf.df <- subset(arle.inf.df, type=="ins")
    arle.inf.df$start.barcode <- barcodpToBarcode[arle.inf.df$start.barcodep]
    arle.inf.df$end.barcode <- arle.inf.df$start.barcode + 1
    # from the barcodep scale, only select insertions
    
  }    
  
  
  ### barcode scale with a scale ------------------------------------------------------ 
  arle.df <- data.frame(unclass(rle(as.numeric(as.vector(seqSubset) == as.vector(refSubset)))))      
  # find start and end coordinates of each segment
  arle.df$end.barcode <- cumsum(arle.df$lengths)
  if (nrow(arle.df)==1) {
    arle.df$start.barcode <- 1
  } else {
    arle.df$start.barcode <- c(0, arle.df$end.barcode[1:(nrow(arle.df)-1)])+1
  }   
  # mismatches mean mutations - subset to mutations
  arle.df <- subset(arle.df, values==0)
  
  # if any mutations found, add their sequence and type
  if (nrow(arle.df)>0) {
    # annotate each mutations with alt and ref alleles, type and sequence ID
    for (mi in 1:nrow(arle.df)) {
      arle.df$ref[mi]<-paste(refSubset[ arle.df$start.barcode[mi]:arle.df$end.barcode[mi]] ,collapse='')
      arle.df$alt[mi]<-paste(seqSubset[ arle.df$start.barcode[mi]:arle.df$end.barcode[mi]], collapse='')
    }
    arle.df$type <- "sub"
    arle.df[arle.df$lengths>1, "type"] <- "complex"
    # if alt sequence is made entirely of hyphens, it is a deletions
    arle.df[gsub("-", "", arle.df$alt)=="", "type"] <- "del"
    # there should not be any insertions here: they are dealt 
  }    
  # add the barcodep scale
  arle.df$start.barcodep <- barcodeToBarcodep[arle.df$start.barcode]
  arle.df$end.barcodep <- barcodeToBarcodep[arle.df$end.barcode]
  
  if ((nrow(arle.inf.df)>0) &&(nrow(arle.df)>0)) {
    all.muts <- rbind(arle.inf.df[,colnames(arle.df)], arle.df)  
  } else if ((nrow(arle.inf.df)>0)) {
    all.muts <- arle.inf.df
  }  else {
    all.muts <- arle.df
  }
  
  return(all.muts)
  
}

## Correlations Analysis   ------------------------------------------------------ 
# Get lower triangle of the correlation matrix
get_lower_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[upper.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[lower.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

# hc.order correlation matrix
hc_cormat_order <- function(cormat, hc.method = "complete") {
  dd <- stats::as.dist((1 - cormat) / 2)
  hc <- stats::hclust(dd, method = hc.method)
  cormat <- cormat[hc$order, hc$order]
}


# Phylogenetic functions ------------------------------------------------------ 
#as.multiPhylo.phylo    
as.multiPhylo.phylo<-function(x,...){
  obj<-list(x)
  class(obj)<-"multiPhylo"
  obj
}

as.multiPhylo<-function(x,...){
  if (identical(class(x),"multiPhylo")) return(x)
  UseMethod("as.multiPhylo")
}


# 1k axis formatting   ------------------------------------------------------ 
ks <- function (x) { number_format(accuracy = 1,
                                   scale = 1/1000,
                                   suffix = "k",
                                   big.mark = ",")(x) }


# binding columns function   ------------------------------------------------------ 
rbind.all.columns <- function(x, y) {
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  x[, c(as.character(y.diff))] <- NA
  y[, c(as.character(x.diff))] <- NA
  return(rbind(x, y))
}

# add labels angles for for pie plot    ------------------------------------------------------ 
compute_angle = function(perc_in_sample){
  angle = 1
  if(perc_in_sample < 0.25) # 1st q [90,0]
    angle = 90 - (perc_in_sample/0.25) * 90
  else if(perc_in_sample < 0.5) # 2nd q [0, -90]
    angle = (perc_in_sample - 0.25) / 0.25 * - 90
  else if(perc_in_sample < 0.75) # 3rd q [90, 0]
    angle = 90 - ((perc_in_sample - 0.5) / 0.25 * 90)
  else if(perc_in_sample < 1.00) # last q [0, -90]
    angle = ((perc_in_sample - 0.75)/0.25) * - 90
  
  # Or even more compact, but less readable
  #if(perc_in_sample < 0.5) # 1st half [90, -90]
  # angle = (180 - (perc_in_sample/0.5) * 180) - 90
  #else # 2nd half [90, -90]
  # angle = (90 - ((perc_in_sample - 0.5)/0.5) * 180)
  #return(angle)
}


### tuniec ###


