## Calculate diversity metrics averages for ecoregions, biomes
## Robert MacDonald
## 17th October 2024


# clear environment
rm(list=ls())

## Load libraries ----
library(ggplot2)
library(dplyr)
#library(terra)
library(tidyterra)
library(dispRity)

# Load shared functions ----
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "mapping.R"
  )
)

source(
  here::here(
    "3_SharedScripts", "dispRity_metric_functions.R"
  )
)

## EDITABLE CODE ## ----
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neognaths"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
# N.B. will need to change the date in the pca_all filename if using usmldbl or usmldblr
# (from 240603 to 240806)
space <- "lab"
# select sex (""allspecimens" or "matchedsex")
# "matchedsex" will subset to only species for which we have at least one
# male and female specimen
sex_match <- "matchedsex"
# select sex of interest ("all", "male_female", "male_only", "female_only", "unknown_only")
sex_interest <- "male_female"
# select metric ("centr-dist", "nn-k", "nn-count")
# note that nn-k is EXTREMELY slow to run - needs parallelisation (but will probably still
# be too slow to run)
metric <- "centr-dist"
# select averaging type ("mean", "median", "mode")
avg_par <- "mean"
# select whether to calculate local or global diversity loss (i.e., mean distance to local or global centroid)
div_loss_type <- "local"
# select whether to exclude grid cells with species richness below a certain threshold (e.g. 5)
# set as 0 if no threshold wanted
sr_threshold <- 5
# select whther to use liberal, conservative, or nominate IUCN data
# "liberal" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the highest threat level of the multiple species
# "conservative" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the lowest threat level of the multiple species
# "nominate" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the threat level of the BL species that corresponds to the nominate subspecies
iucn_type <- "nominate"
# select whether to use liberal, conservative, or specified PAM
pam_type <- "conservative"
# clip PAM to land only? ("_clipped" or "")
pam_seas <- "_clipped"
# select PAM grid cell resolution
pam_res <- "50km"
# enter PAM files location
pams_filepath <- "X:/cooney_lab/Shared/Rob-MacDonald/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/PAMs"
#pams_filepath <- "X:/cooney_lab/Shared/Rob-MacDonald/SpatialData/BirdLife/BirdLife_Shapefiles_v9/PAMs/100km/Behrmann_cea/"
# use ecoregions or biomes?
regions <- "biomes"
# parallelise working (parallelisation doesn't currently work so set to FALSE)
parallelise <- FALSE

## Plotting parameters
# Select whether to use binned (based on quantiles) or continuous colour scale
col_scale_type <- "binned"
# If binned, choose the number of quantiles to use
nquants <- 10
# Also plot species richness map?
plot_sr <- TRUE
## Use viridis_c, viridis turbo palette or custom colour palette?
palette_choice <- "viridis"
## END EDITABLE CODE ##

# set dispRity metric
metric_get <- set_metric(metric)

## Load data ----

# load patch data (PCA of whichever colourspace - generated in 02_Patch_Analyse_features.R)
pca_filename <- paste(clade, sex_match, "patches.250716.PCAcolspaces", "rds", sep = ".")
pca_dat <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", clade, "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]][["x"]]

# load PAM
pam_filename <- paste0("PAM_birds_Behrman", pam_res, "_Pres12_Orig12_Seas12_", pam_type, pam_seas, ".rds")
# or
#pam_filename <- paste0(pam_res, "_Behrmann_cea_PAM_combined.rds")
pam_raw <- readRDS(
  paste(pams_filepath, pam_filename, sep = "/")
)

# load region data (either Dinerstein et al 2017 ecoregion data or biome data derived from same)
region_shapes <- load_regions(regions)



## Analysis ----

# Set spatial parameters
# extract null raster from PAM file
null_rast <- extract_null_rast(pam_raw)

# extract matrix values from PAM file
pam <- extract_pam_vals(pam_raw)

# remove raw PAM file (for RAM reasons)
rm(pam_raw)
# change "Genus species" style to "Genus_species" style in PAM colnames
colnames(pam) <- gsub(" ", "_", colnames(pam))

# convert ecoregions to Behrman CEA CRS
region_shapes <- sf::st_transform(region_shapes, crs = terra::crs(null_rast, proj = TRUE))


# get list of species included in diversity/PCA data
spp_list <- get_unique_spp(pca_dat)

# subset pam to only species in diversity data
pam <- subset_pam(pam, spp_list)

# subset diversity pca_dat to only species present in PAM
pca_dat <- subset_data_by_spp(pca_dat, colnames(pam))

# clear garbage
gc()

# generate individual list of sex-specific species values of metric of interest (keep or discard sexes
# depending on value of sex_match)
# first set sex list to iterate over
sexes <- set_sex_list(sex_interest)

# now apply a function to extract named vector lists of values for each sex individually
sexed_pca_list <- lapply(sexes, extract_sex_vals, div_data = pca_dat, assemblage_level = TRUE)
names(sexed_pca_list) <- sexes

# if calculating diversity loss vs global species pool, calculate distance to global centroid for individual
# sexed species
if(div_loss_type == "global"){
  
  # first calculate the metric
  sexed_global_metric <- lapply(sexed_pca_list, metric_vals, metric = metric)
  
  # now get the average metric for each region
  # get number of regions
  n_regions <- length(attr(region_shapes, "row.names"))
  # apply function across sexes and across rows of ecoregions
  avg_vals <- lapply(sexed_global_metric, function(sexed_metric){
    lapply(
      1:n_regions, calc_regions_global_metric,
      region_shapes = region_shapes, region_type = regions, metric = metric, metric_vals = sexed_metric, div_loss_type = div_loss_type, null_terrast = null_rast, avg_par = avg_par
    )
    }
  )
  
} else if(div_loss_type == "local"){
  # if calculating diversity loss vs local species pool, calculate distance to local centroid for individual
  # regions  
  
  
  # get number of regions
  n_regions <- length(attr(region_shapes, "row.names"))
  
  # need to apply to each sexed pca dataset
  avg_vals <- lapply(
    
    sexed_pca_list, function(sexed_pca_data){
      
      # non-parallel version
      if(parallelise == FALSE){
        
        # get the average metric for each region - feed sexed pca data into the function to calculate 
        # regional average metric (applied over each region)
        result <- lapply(
          1:n_regions, calc_regions_global_metric,
          region_shapes = region_shapes, region_type = regions, metric = metric, metric_vals = NULL, pca_data = sexed_pca_data, div_loss_type = div_loss_type, null_terrast = null_rast, avg_par = avg_par
        )
        
      } else if(parallelise == TRUE){
        
        #### NOT USE - doesn't currently work
        
        # set up cluster
        no_cores <- parallel::detectCores() - 10
        cl <- parallel::makeCluster(no_cores)
        
        # export calculation function
        parallel::clusterExport(cl, c("calc_regions_global_metric", "region_shapes", "regions", "metric", "sexed_pca_data", "div_loss_type", "null_rast", "avg_par", "n_regions"), envir = environment())
        

        # parallel apply over the regions
        result <- parallel::parLapply(cl, 1:n_regions, function(region_number){
          calc_regions_global_metric(
            region_index = region_number, region_shapes = region_shapes, region_type = regions, 
            metric = metric, metric_vals = NULL, pca_data = sexed_pca_data, 
            div_loss_type = div_loss_type, null_terrast = null_rast, avg_par = avg_par
          )
        })
        
        # stop cluster
        parallel::stopCluster(cl)
        
      }
      
      return(result)
      
    }
  )
    
  }


# convert output to 3-column matrix
results_df <- data.frame(
  M_values = sapply(avg_vals$M, function(x) x[1]),
  F_values = sapply(avg_vals$F, function(x) x[1]),
  species_richness = sapply(avg_vals$M, function(x) x[2])  # or use avg_vals$F
)

# if using a quantile colour scale, create a new column of classes based on metric value quantiles shared between sexes
if(col_scale_type == "binned"){
  
  # combine male and female values together in a list
  all_vals <- c(results_df$M_values, results_df$F_values)
  breaks <- quantile(all_vals, na.rm = T, probs = seq(0, 1, by = 1/nquants))
  # make names of quants a rounded version of the quantiles (for plotting)
  names(breaks) <- round(breaks, digits = 1)
  breaks_sr <- quantile(results_df[, "species_richness"], na.rm = TRUE, probs = seq(0, 1, by = 1/nquants))
  names(breaks_sr) <- round(breaks_sr, digits = 1)
  
  # create function to provide class names for single data column
  class_namer <- function(data_col, colour_breaks, n_quantiles){
    
    class_col <- rep(NA, times = length(data_col))
    
    for(quant in 1:n_quantiles){
      class_col[data_col <= colour_breaks[quant + 1] & is.na(class_col)] <- as.character(quant)
    }
    
    return(class_col)
    
  }
  
  # apply function to get quantile for each sex
  for(sex in sexes){
    
    results_df[paste("col_class", sex, sep = "_")] <- class_namer(results_df[, paste(sex, "values", sep = "_")], breaks, nquants)
    
  }
  # and for species richness
  results_df["col_class_sr"] <- class_namer(results_df[, "species_richness"], breaks_sr, nquants)
}


# ## For INDIVIDUAL sex quantile scales
# # if using a quantile colour scale, create a new column of classes based on metric value quantiles
# # for each sex
# if(col_scale_type == "binned"){
#   breaks <- lapply(sexes, function(sex){
#     
#     colname <- paste(sex, "values", sep = "_")
#     
#     bks <- quantile(results_df[, colname], na.rm = TRUE, probs = seq(0, 1, by = 1/nquants))
#     
#     # make names of quants a rounded version of the quantiles (for plotting)
#     names(bks) <- round(bks, digits = 1)
#     
#     return(bks)
#     
#   })
#   breaks[["sr"]] <- quantile(results_df[, "species_richness"], na.rm = TRUE, probs = seq(0, 1, by = 1/nquants))
#   names(breaks[["sr"]]) <- round(breaks[["sr"]], digits = 1)
#   names(breaks) <- c(sexes, "sr")
#   
#   # create function to provide class names for single data column
#   class_namer <- function(data_col, colour_breaks, n_quantiles){
#     
#     class_col <- rep(NA, times = length(data_col))
#     
#     for(quant in 1:n_quantiles){
#       class_col[data_col <= colour_breaks[quant + 1] & is.na(class_col)] <- as.character(quant)
#     }
#     
#     return(class_col)
#     
#   }
#   
#   for(sex in sexes){
#     
#     results_df[paste("col_class", sex, sep = "_")] <- class_namer(results_df[, paste(sex, "values", sep = "_")], breaks[[sex]], nquants)
#     
#   }
#   ### BELOW LINE DOESN'T WORK
#   results_df["col_class_sr"] <- class_namer(results_df[, "species_richness"], breaks[["sr"]], nquants)
# }


# bind to ecoregion data
region_shapes <- cbind(region_shapes, results_df)


# Plot ----

# clear irrelevant objects
rm(list=setdiff(ls(), c("clade", "space", "sex_match", "regions", "avg_par", "div_loss_type", "metric", "sr_threshold", "pam_res", "pam_type", "pam_seas", "nquants", "col_scale_type", "plot_sr", "palette_choice")))

# load results bound to ecoregion data
results_shapes_filename <- paste(clade, "patches", sex_match, regions, avg_par, div_loss_type, metric, "sr_thresh", sr_threshold, pam_res, "Behrman", pam_type, pam_seas, "shp", sep = ".")
output_folder <- here::here("2_Patches", "3_OutputData", clade, "6_Spatial_mapping", "4_Regional_diversity", pam_res)
region_shapes <- sf::st_read(paste(output_folder, results_shapes_filename, sep = "/")) %>% 
  rename(
    BIOME_NUM = BIOME_NU,
    BIOME_NAME = BIOME_NA,
    M_values = M_valus,
    F_values = F_valus,
    species_richness = spcs_rc,
    col_class_M = cl_cl_M,
    col_class_F = cl_cl_F,
    col_class_sr = cl_cls_
  )

# plot ecoregions coloured by species richness and by average metric
plots <- list()

# for continuous colour scale
if(col_scale_type == "continuous"){
  
  if(plot_sr == TRUE){
    # sr plot
    psr <- ggplot() + 
      geom_sf(data = region_shapes, aes(fill = species_richness, colour = species_richness)) +
      scale_fill_viridis_c() + 
      scale_colour_viridis_c() + 
      coord_sf(expand = FALSE) + 
      labs(fill = "species_richness", colour = "species_richness") + 
      theme_minimal() + 
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    plots[["species_richness"]] <- psr
  }
  
  # male plot
  pm <- ggplot() + 
    geom_sf(data = region_shapes, aes(fill = M_values, colour = M_values)) +
    scale_fill_viridis_c() + 
    scale_colour_viridis_c() + 
    coord_sf(expand = FALSE) + 
    labs(fill = paste0(avg_par, "_metric_M"), colour = paste0(avg_par, "_metric_M"))+ 
    theme_minimal() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  plots[["male_metric"]] <- pm
  
  # female plot
  pf <- ggplot() + 
    geom_sf(data = region_shapes, aes(fill = F_values, colour = F_values)) +
    scale_fill_viridis_c() + 
    scale_colour_viridis_c() + 
    coord_sf(expand = FALSE) + 
    labs(fill = paste0(avg_par, "_metric_F"), colour = paste0(avg_par, "_metric_F"))+ 
    theme_minimal() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  plots[["female_metric"]] <- pf
  
} else if(col_scale_type == "binned"){
  
  # create colour palette
  col_pal <- viridisLite::viridis(nquants)
  names(col_pal) <- 1:nquants
  na_col <- "grey50"
  
  # plots
  if(plot_sr == TRUE){
    # sr plot
    psr <- ggplot() + 
      geom_sf(data = region_shapes, aes(fill = col_class_sr, colour = col_class_sr)) +
      scale_fill_manual(values = col_pal, na.value = na_col) + 
      scale_colour_manual(values = col_pal, na.value = na_col) + 
      coord_sf(expand = FALSE) + 
      labs(fill = "species_richness", colour = "species_richness")+ 
      theme_minimal() + 
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    plots[["species_richness"]] <- psr
  }
  
  # male plot
  pm <- ggplot() + 
    geom_sf(data = region_shapes, aes(fill = col_class_M, colour = col_class_M)) +
    scale_fill_manual(values = col_pal, na.value = na_col) + 
    scale_colour_manual(values = col_pal, na.value = na_col) + 
    coord_sf(expand = FALSE) + 
    labs(fill = paste0(avg_par, "_metric_M"), colour = paste0(avg_par, "_metric_M"))+ 
    theme_minimal() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  plots[["male_metric"]] <- pm
  
  # female plot
  pf <- ggplot() + 
    geom_sf(data = region_shapes, aes(fill = col_class_F, colour = col_class_F)) +
    scale_fill_manual(values = col_pal, na.value = na_col) + 
    scale_colour_manual(values = col_pal, na.value = na_col) + 
    coord_sf(expand = FALSE) + 
    labs(fill = paste0(avg_par, "_metric_F"), colour = paste0(avg_par, "_metric_F"))+ 
    theme_minimal() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  plots[["female_metric"]] <- pf
 
  
}

# save
png_filename <- paste(clade, "patches", sex_match, regions, col_scale_type, avg_par, div_loss_type, metric, "sr_thresh", sr_threshold, pam_res, "Behrman", pam_type, pam_seas, "png", sep = ".")

space_plot_path <- here::here(
  "2_Patches", "4_OutputPlots", clade, "3_Spatial_mapping", "4_Bioregion_level_mapping", pam_res, space, pam_type
)
if(!dir.exists(space_plot_path)){
  dir.create(space_plot_path, recursive = TRUE)
}

if(plot_sr == TRUE){
  n_col <- 3
  png_width <- 600
} else {
  n_col <- 2
  png_width <- 400
}
png(
  paste(space_plot_path, png_filename, sep = "/"), 
  width = 500, height = 150, units = "mm", pointsize = 24, res = 100
)

gridExtra::grid.arrange(grobs = plots, ncol = n_col)
# or
# ggpubr::ggarrange(plotlist = plots, common.legend = T, ncol = n_col)

dev.off()



# save ecoregion data with diversity values attached
# first as shapes file
results_shapes_filename <- paste(clade, "patches", sex_match, regions, avg_par, div_loss_type, metric, "sr_thresh", sr_threshold, pam_res, "Behrman", pam_type, pam_seas, "shp", sep = ".")
output_folder <- here::here("2_Patches", "3_OutputData", clade, "6_Spatial_mapping", "4_Regional_diversity", pam_res)
if(!dir.exists(output_folder)){
  dir.create(output_folder, recursive = TRUE)
}
sf::st_write(region_shapes, paste(output_folder, results_shapes_filename, sep = "/"))


# inspect most diverse ecoregions
region_shapes %>% 
  arrange(desc(F_values)) %>% 
  head(n = 10)

# inspect ecoregions in top quantiles
region_shapes %>% 
  filter(
    col_class_F == as.character(nquants)
  )


# dev ----

# first pivot longer to get all values in one column
r_s_long <- region_shapes %>% 
  pivot_longer(c(M_values, F_values, -geometry), names_to = "sex") %>% 
  select(
    -c(col_class_M, col_class_F, col_class_sr)
  )



# add quantile column
value_breaks <- quantile(r_s_long$value, probs = seq(0, 1, 1/nquants), na.rm = TRUE)
break_labels <- rep(NA, times = nquants)
for(label_index in 1:nquants){
  
  break_labels[label_index] <- paste0(round(value_breaks[label_index], 1), "-", round(value_breaks[label_index + 1], 1))
  
}
r_s_long$val_quant <- r_s_long$value %>% 
  cut(
    breaks = value_breaks,
    include.lowest = TRUE,
    labels = break_labels
  )

# create a short version for dev purposes
r_s_long <- r_s_long[1:20, ]

# plot
ggplot() + 
  geom_sf(data = r_s_long, aes(fill = val_quant, colour = val_quant)) + 
  scale_fill_viridis_d(name = metric, guide = guide_legend()) + 
  scale_colour_viridis_d(guide = "none") +
  facet_wrap(~ sex, ncol = 1) + 
  theme_minimal() + 
  theme(legend.position="bottom",
        legend.spacing.x = unit(0, 'cm'))
  guides(fill = guide_legend(label.position = "bottom"))
  
  

