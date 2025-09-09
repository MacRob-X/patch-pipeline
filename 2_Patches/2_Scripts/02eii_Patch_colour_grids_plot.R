# Create visuals showing patch colours on colour space ----
# Robert MacDonald
# 13th May 2024

# Load libraries ---- 
library(dplyr)
library(ggplot2)
library(extrafont)

# clear environment
rm(list=ls())

# load plotting functions ----
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "plotting.R"
  )
)

# Choose parameters ----
## EDITABLE CODE ##
## Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neognaths"
## Choose colour space mapping
## Note that if usmldbl or usmldblr, will have to manually modify date in filename below
space <- "lab"
# Restrict to only species for which we have male and female data?
mf_restrict <- TRUE
if(mf_restrict == TRUE){
  spec_sex <- "matchedsex"
} else{
  spec_sex <- "allspecimens"
}
# Plot single double axis plot (FALSE) or four double axis plot (TRUE)?
plot_multiple <- FALSE
## UMAP or PCA space?
## FALSE - use the PCA space
## TRUE - load a UMAP space from a file
## "perform" - load a PCA space from a file, then perform UMAP on it
## Note that if umap == TRUE, user will have to manually set path to umap file
load_umap <- TRUE
umap_filepath <- here::here(
  "2_Patches", "3_OutputData", clade, "2_PCA_ColourPattern_spaces", "2_UMAP",
  paste(clade, spec_sex, "patches", "nn.25.mindist.0.1", space, "UMAP.rds", sep = ".")
)
# umap_filepath <- here::here(
#   "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
#   "Neoaves.patches.pca.jndxyzlumr.UMAPs.iterations.240326.rds"
# )
# If load_umap == FALSE, can set perform_umap == TRUE
# to perform UMAP on the loaded PCA and then plot the UMAP axes
perform_umap <- FALSE
## Choose PCs to plot (if plotting PCA)
x_axis <- "PC1"
y_axis <- "PC2"
## Choose aspect ratio 
## -"square" fixes to size 2500x2500px 
## - "wrap" sets the axis with the larger range to 2500px and scales the other accordingly)
## "square" is better if you want to combine multiple plots into one figure, "wrap" 
## is better for a single plot as it wraps the size of the plot to the range sizes
## Note that the aspect ratio of the plot itself is always fixed to 1, this parameter is 
## solely to control the aspect ratio of the output png
asp_ratio <- "wrap"
## Choose font "default" or name font (e.g. "Century Gothic")
## Use extrafont::fonts() to see available fonts
font_par <- "Century Gothic"
# set path to colour grids folder
grid_path <- "C:/Users/bop23rxm/Documents/colour_grids_repositioned"

# Load data ----
if(plot_multiple == FALSE){
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
      "2_Patches", "3_OutputData", clade, "2_PCA_ColourPattern_spaces", "1_Raw_PCA", 
      paste(clade, "patches.250716.PCAcolspaces", "rds", sep = ".")
      #    paste0("Neoaves.patches.250716.PCAcolspaces.", space, ".240829.rds")
    )
    # load data from file
    pca_all <- readr::read_rds(pca_filepath) %>% 
      magrittr::extract2(space)
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
} else if(plot_multiple == TRUE){
  
  # load PCA space
  
  # set path to file
  pca_filepath <- here::here(
    "2_Patches", "3_OutputData", clade, "2_PCA_ColourPattern_spaces", "1_Raw_PCA", 
    paste(clade, spec_sex, "patches.250716.PCAcolspaces", "rds", sep = ".")
    #    paste0("Neoaves.patches.250716.PCAcolspaces.", space, ".240829.rds")
  )
  # load data from file
  pca_all <- readr::read_rds(pca_filepath) %>% 
    magrittr::extract2(space)
  # extract pc scores to plot
  pca_space <- pca_all %>% 
    magrittr::extract2("x") %>% 
    as.data.frame()
  # get proportions of variance for each axis
  pca_variance <- pca_all %>% 
    summary() %>% 
    magrittr::extract2("importance") %>% 
    as.data.frame()
  
  if(load_umap == TRUE) {
    
    umap_space <- readr::read_rds(umap_filepath) %>% 
      magrittr::extract2("layout") %>% 
      as.data.frame() %>% 
      rename(
        UMAP1 = V1, UMAP2 = V2
      )

  } else if (load_umap == FALSE){
    
    if(perform_umap == TRUE){
      umap_space <- pca_space %>% 
        umap::umap() %>% 
        magrittr::extract2("layout") %>% 
        as.data.frame() %>% 
        rename(
          UMAP1 = V1, UMAP2 = V2
        )

    }
  }
  
}



# Plot space

# for single plot
if(plot_multiple == FALSE){
  
  # if plotting UMAP
  if(load_umap == TRUE){
    # set save path for plot
    save_path <- here::here(
      "2_Patches", "4_OutputPlots", clade, "1_Colourspace_visualisation", space,
      paste(clade, spec_sex, "all", "patches_UMAPnn25mindist0.1.png", sep = "_")
    )
    # plot
    plot_patch_grids(
      data_matrix = plot_space,
      x_axis = x_axis, y_axis = y_axis,
      colour_grid_path = grid_path, 
      asp_ratio = "wrap",
      save_as = "png",
      save_path = save_path
    )
    
  } else if(perform_umap == TRUE){
    # set save path for plot
    save_path <- here::here(
      "2_Patches", "4_OutputPlots", clade, "1_Colourspace_visualisation", space,
      paste(clade, spec_sex, "all", "patches_defaultUMAP.png", sep = "_")
    )
    # plot
    plot_patch_grids(
      data_matrix = umap_space,
      colour_grid_path = grid_path, 
      asp_ratio = "wrap",
      save_as = "png",
      save_path = save_path
    )
  } else {
    # plot PC axes
    
    plot_patch_grids(
      data_matrix = pca_space,
      x_axis = x_axis, y_axis = y_axis,
      colour_grid_path = grid_path, 
      asp_ratio = "wrap",
      save_as = "png",
      save_path = save_path
    )
    
  }
  
}

# For multiple (four) plots
if(plot_multiple == TRUE){
  folder_path <- here::here(
    "2_Patches", "4_OutputPlots", clade, "1_Colourspace_visualisation", space
  )
  # if plotting PC1-6 and UMAP
  if(load_umap == TRUE){
    
    filename <- paste(clade, spec_sex, "all", "patches_PC1-PC6_UMAP.png", sep = "_")
    plot_four_cg(
      pca_all, "PC1", "PC2",
      pca_all, "PC3", "PC4",
      pca_all, "PC5", "PC6",
      umap_space, "UMAP1", "UMAP2",
      cg_path = grid_path,
      save_type = "png",
      write_folder = folder_path,
      #  thin_number = 500,
      file_name = filename
    )
  } else if(load_umap == FALSE){
    # if plotting PC1-8
    
    filename <- paste(clade, spec_sex, "all", "patches_PC1-PC8.png", sep = "_")
    plot_four_cg(
      pca_all, "PC1", "PC2",
      pca_all, "PC3", "PC4",
      pca_all, "PC5", "PC6",
      pca_all, "PC7", "PC8",
      cg_path = grid_path,
      save_type = "png",
      write_folder = folder_path,
      #  thin_number = 500,
      file_name = filename
    )
    
  }
  
}


# Heatmap plotting ----

# clear environment
rm(list=ls())

# load plotting functions ----
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "plotting.R"
  )
)

# load libraries
library(terra)
library(tidyterra)

## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neognaths"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
# N.B. will need to change the date in the pca_all filename if using usmldbl or usmldblr
# (from 240603 to 240806)
space <- "lab"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# select sex of interest ("all", "male_female", "male_only", "female_only", "unknown_only")
sex_interest <- "male_female"

# Choose PC axes for which to plot heatmaps
axes <- paste0("PC", 1:4)
pixel_type <- "lab"
# choose whether to plot heatmaps with constant unified colour scale for all PC axes ("constant"), individual 
# colour scales for each axis ("axis") or individual colour scales for each channel of each axis ("axis_channel")
col_scale_type <- "constant"

# Load data ----

# Load PCA data
pca_filename <- paste(clade, sex_match, "patches.250716.PCAcolspaces", "rds", sep = ".")
pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", clade, "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]]
pca_space <- pca_all


# Generate non-stacked (standard) heatmaps
p_standard <- patch_loading_heatmap(pca_all, axes, pixel_type, col_scale_type = col_scale_type)

# generate stacked heatmaps
p_stacked <- patch_loading_heatmap(pca_all, axes, pixel_type, col_scale_type = col_scale_type, stack_channels = TRUE)




### ARCHIVE ###

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
    paste(clade, spec_sex, "patches", space, x_axis, y_axis, "colour_grids.png", sep = "_")
  ),
  width = png_width, height = png_height,
  units = "px",
  pointsize = 50,
  family = font_par
  
)

plot(plot_space[, y_axis] ~ plot_space[, x_axis], 
     asp = T, 
     type = "n", 
     xlab = xlabel, ylab = ylabel, las = 1,
#     xlim = xrange,  ylim = yrange
     )

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
}

dev.off()

