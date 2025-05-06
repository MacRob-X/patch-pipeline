# Create visuals showing patch colours on colour space ----
# Robert MacDonald
# 13th May 2024

# Load libraries ---- 
library(dplyr)
library(ggplot2)
library(extrafont)

# clear environment
rm(list=ls())

# Choose parameters ----
## EDITABLE CODE ##
## Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neoaves"
## Choose colour space mapping
## Note that if usmldbl or usmldblr, will have to manually modify date in filename below
space <- "lab"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
## UMAP or PCA space?
## FALSE - use the PCA space
## TRUE - load a UMAP space from a file
## "perform" - load a PCA space from a file, then perform UMAP on it
## Note that if umap == TRUE, user will have to manually set path to umap file
load_umap <- FALSE
umap_filepath <- here::here(
  "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
  paste(clade, sex_match, "patches", space, "pca", "canonUMAP.rds", sep = ".")
)
# umap_filepath <- here::here(
#   "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
#   "Neoaves.patches.pca.jndxyzlumr.UMAPs.iterations.240326.rds"
# )
# If load_umap == FALSE, can set perform_umap == TRUE
# to perform UMAP on the loaded PCA and then plot the UMAP axes
perform_umap <- FALSE
## Choose PCs to plot (if plotting PCA)
x_axis <- "PC2"
y_axis <- "PC3"
## Choose aspect ratio 
## -"square" fixes to size 2500x2500px 
## - "wrap" sets the axis with the larger range to 2500px and scales the other accordingly)
## "square" is better if you want to combine multiple plots into one figure, "wrap" 
## is better for a single plot as it wraps the size of the plot to the range sizes
## Note that the aspect ratio of the plot itself is always fixed to 1, this parameter is 
## solely to control the aspect ratio of the output png
asp_ratio <- "wrap"
## Choose font "default" or name font (e.g. "Century Gothic", "AvenirNext LT Pro Regular)
## Use extrafont::fonts() to see available fonts
font_par <- "AvenirNext LT Pro Regular"

# Load data ----

if(load_umap == TRUE) {
    plot_space <- readr::read_rds(umap_filepath) %>% 
      magrittr::extract2("layout") %>% 
      as.data.frame() %>% 
      rename(
        UMAP1 = V1, UMAP2 = V2
      )
    x_axis <- "UMAP1"
    y_axis <- "UMAP2"
} else if (load_umap == FALSE){
  # set path to file
  pca_filepath <- here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA", 
    paste(clade, sex_match, "patches.231030.PCAcolspaces", "rds", sep = ".")
#    paste0("Neoaves.patches.231030.PCAcolspaces.", space, ".240829.rds")
  )
  # load data from file
  pca_all <- readr::read_rds(pca_filepath)[[space]]
  # extract pc scores to plot
  plot_space <- pca_all %>% 
    magrittr::extract2("x") %>% 
    as.data.frame()
  # get proportions of variance for each axis
  pca_variance <- pca_all %>% 
    summary() %>% 
    magrittr::extract2("importance") %>% 
    as.data.frame()
  
  if(perform_umap == TRUE){
    plot_space <- plot_space %>% 
      umap::umap() %>% 
      magrittr::extract2("layout") %>% 
      as.data.frame() %>% 
      rename(
        UMAP1 = V1, UMAP2 = V2
      )
    x_axis <- "UMAP1"
    y_axis <- "UMAP2"
  }
}


## Plot space ----

if(asp_ratio == "wrap"){
  # determine axis with greatest range (to determine size of png)
  range_x_axis <- diff(range(plot_space[, x_axis]))
  range_y_axis <- diff(range(plot_space[, y_axis]))
  if(range_x_axis > range_y_axis){
    png_width <- 2500
    png_height <- png_width * (range_y_axis / range_x_axis)
  } else if(range_y_axis > range_x_axis){
    png_height <- 2500
    png_width <- png_height * (range_x_axis / range_y_axis)
  }
} else if (asp_ratio == "square"){
  png_width <- 2500
  png_height <- 2500
}


# get x and y ranges (note that I actually think this is obsolete when asp = T)
# x_5 <- (max(plot_space[, x_axis]) - min(plot_space[, x_axis])) / 20
# y_5 <- (max(plot_space[, y_axis]) - min(plot_space[, y_axis])) / 20
# xrange <- c((min(plot_space[, x_axis]) - x_5), (max(plot_space[, x_axis]) + x_5))
# yrange <- c((min(plot_space[, y_axis]) - y_5), (max(plot_space[, y_axis]) + y_5))

# set axis labels
if(load_umap == TRUE | perform_umap == TRUE){
  xlabel <- "UMAP1"
  ylabel <- "UMAP2"
} else {
  xlabel <- paste0(x_axis, " (", round(pca_variance[, x_axis][2], 2) * 100, "% of variance)")
  ylabel <- paste0(y_axis, " (", round(pca_variance[, y_axis][2], 2) * 100, "% of variance)")
}



# Initialise saving
# Edit the filepath here if you want to save to a different location
png(
  here::here(
    "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", space, 
    paste(clade, sex_match, "patches", space, x_axis, y_axis, "colour_grids_CG.png", sep = "_")
  ),
  width = png_width, height = png_height,
  units = "px",
  pointsize = 40,
  family = font_par
  
)

plot(plot_space[, y_axis] ~ plot_space[, x_axis], 
     asp = T, 
     type = "n", 
     xlab = "", ylab = "", las = 1,
     xaxt = "n", yaxt = "n"
#     xlim = xrange,  ylim = yrange
     )
title(xlab = xlabel, ylab = ylabel, line = 1)

box(lwd = 2)
rw <- diff(range(plot_space[,1]))/12

for(i in 1:nrow(plot_space)){
  fname <- rownames(plot_space)[i]
  fpng <- png::readPNG(paste0(
    "C:/Users/bop23rxm/Documents/colour_grids_repositioned/", 
    paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
  rasterImage(fpng, 
              xleft = plot_space[i, x_axis] - (rw/15), 
              ybottom = plot_space[i, y_axis]-(rw/9), 
              xright = plot_space[i, x_axis]+(rw/15), 
              ytop = plot_space[i, y_axis]+(rw/9))
  cat(paste0("\r", i, " of ", nrow(plot_space), " processed"))
  
  if(isTRUE(all.equal(i/500, as.integer(i/500)))){
    # Force an update to the graphics device after each 500 images
    # This helps prevent memory issues
    if(interactive()) {
      Sys.sleep(0.01)  # Small delay to allow graphics to update
    }
    
    # Free up memory
    gc()
  }
  
}

dev.off()

