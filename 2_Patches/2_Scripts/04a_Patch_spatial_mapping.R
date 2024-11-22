# New version of PAM spatial mapping
# Takes advantage of separate generation of diversity metrics
# 4th November 2024
# Robert MacDonald

# Clear environment
rm(list=ls())

# Load libraries
library(dplyr)
library(raster)
library(ggplot2)
library(tidyterra)

## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Passeriformes"
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
# select metric ("centr-dist", "nn-k", "nn-count")
metric <- "centr-dist"
# select type of averaging to use ("mean" or "median")
avg_par <- "mean"
# select whether to exclude grid cells with species richness below a certain threshold (e.g. 5)
# set as 0 if no threshold wanted
sr_threshold <- 5
# select whether to exclude species with metric value below a certain percentile threshold (default is 75)
sift_div_data <- FALSE
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

## Plotting parameters
# select whether to filter data to only the top quartile of the metric of interest
# this makes it a bit easier to see trends when plott
sift_rast_data <- FALSE
# Select whether to use binned (based on quantiles) or continuous colour scale
col_scale_type <- "binned"
# If binned, choose the number of quantiles to use
nquants <- 10
# Also plot species richness map?
plot_sr <- TRUE
## Use viridis_c, viridis turbo palette or custom colour palette?
palette_choice <- "viridis"

## END EDITABLE CODE ##

# Functions ----

# filter diversity data by percentile (i.e. keep only species in the top nth percentile of the 
# chosen diversity metric)
filter_div_data <- function(div_dat, div_metric, separate_sexes = TRUE, pcile = 75){
  
  if(separate_sexes == TRUE){
    
    # get species and sex as extra columns
    div_dat <- cbind(div_dat, get_spp_sex(div_dat))
    
    # get percentile values for each sex
    sexes <- unique(div_dat[, "sex"])
    sexed_filtered_dat <- list()
    for(sex in sexes){
      
      # get percentile values for each sex
      pcile_val <- quantile(div_dat[div_dat[, "sex"] == sex, div_metric], prob = pcile/100, na.rm = TRUE)
      
      # set metric value to NA for the species which have metric values lower than the percentile threshold
      div_dat[div_dat[, "sex"] == sex, div_metric] <- ifelse(div_dat[div_dat[, "sex"] == sex, div_metric] < pcile_val, NA, div_dat[div_dat[, "sex"] == sex, div_metric])
      
      # remove all species for which the metric value is now set as NA
      div_dat <- div_dat[!is.na(div_dat[, metric]), ]
      
    }
    
    # remove species/sex columns
    div_dat <- div_dat[, colnames(div_dat) != "species" & colnames(div_dat) != "sex"]
    
    
  } else if(separate_sexes == FALSE){
    
    # get percentile values
    pcile_val <- quantile(div_dat[, div_metric], prob = pcile/100, na.rm = TRUE)
    
    # filter data to only values of metric which are above the percentile threshold
    div_dat <- div_dat[div_dat[, metric] >= pcile_val, ]
    
  }
 
  return(div_dat)
  
}

# extract null raster from PAM file
extract_null_rast <- function(pam_raw){
  
  return(terra::rast(pam_raw[[2]]))
  
}

# extract PAM values
extract_pam_vals <- function(pam_raw){
  
  return(pam_raw[[1]])
  
}

# get list of unique species (of any sex) in data
get_unique_spp <- function(div_data){
  
  # extract species
  species_list <- unique(sapply(strsplit(rownames(div_data), split = "-"), "[", 1))
  
  return(species_list)
  
}

# subset PAM by a list of species (must be a vector)
subset_pam <- function(pam, spp_list){
  
  spp_keep <- colnames(pam)[colnames(pam) %in% spp_list]
  
  return(pam[, spp_keep])
  
}

# subset dataframe/matrix (with species/sex as rownames) by a list of species (vector)
subset_data_by_spp <- function(df, spp_list){
  
  # extract list of species (non-unique)
  df_spp <- sapply(strsplit(rownames(df), split = "-"), "[", 1)
  
  # subset
  return(df[df_spp %in% spp_list, ])
}

# set list of sexes to iterate over
set_sex_list <- function(sex_interest){
  
  if(sex_interest == "all"){
    sexes <- c("M", "F", "U")
  } else if(sex_interest == "male_female"){
    sexes <- c("M", "F")
  } else if(sex_interest == "male"){
    sexes <- "M"
  } else if(sex_interest == "female"){
    sexes <- "F"
  } else if(sex_interest == "unknown"){
    sexes <- "U"
  }
  
  return(sexes)
  
}

# extract named vectors of diversity metric for an individual sex
extract_sex_vals <- function(div_data, metric, sex){
  
  # generate sex suffix
  sex_suffix <- paste0("-", sex)
  
  # filter to individual sex
  metric_sex <- div_data[grepl(sex_suffix, rownames(div_data)), metric, drop = FALSE]
  rownames(metric_sex) <- sapply(strsplit(rownames(metric_sex), split = "-"), "[", 1)
  
  return(metric_sex)
  
}

# calculate species richness for each grid cell
calc_species_richness <- function(pam_matrix){
  
  return(rowSums(pam_matrix, na.rm = TRUE))
  
}

# create species richness mask that can be applied to PAM
make_sr_mask <- function(sr_vector, sr_threshold){
  
  # create mask
  sr_mask <- ifelse(sr_vector >= sr_threshold, TRUE, NA)
  
  return(sr_mask)
  
}

# create species richness raster
make_sr_raster <- function(sr_vals, null_rast){
  
  # set values == 0 to NA (so they don't get plotted)
  sr_vals[sr_vals == 0] <- NA
  
  # create copy of null rast
  sr_raster <- null_rast
  
  # set raster values to species richness values
  terra::values(sr_raster) <- sr_vals
  
  return(sr_raster)
  
  
}

# create raster of diversity metrics
make_diversity_raster <- function(sexed_metric_vals, sex, pam, null_rast, avg_type = "mean"){
  
  # get diversity metric values for the sex of interest
  div_vals <- as.matrix(sexed_metric_vals[[sex]])
  
  # subset PAM to only the species in the sexed list of metric values
  # (this step is only necessary if filtering the diversity data to the upper percentile of values
  # separately for the different sexes - otherwise it should have no effect)
  spp <- get_unique_spp(div_vals)
  pam <- subset_pam(pam, spp)
  
  # make wrapper function for averaging (allows selection of average type e.g. mean, median)
  avg <- function(vals, avg_type, na.rm = TRUE){
    return(get(avg_type)(vals, na.rm = na.rm))
  }
  
  # calculate average value of metric across all species in each grid cell
  cell_avg <- apply(pam, 1, function(one_row) {
    species_present <- which(one_row == 1)
    if (length(species_present) > 0) {
      return(avg(div_vals[names(species_present), , drop = F], avg_type, na.rm = TRUE))
    } else {
      return(NA)  # If no species are present in the cell
    }
  })
  
  # create raster to populate with values
  div_raster <- null_rast
  
  # assign values to raster
  terra::values(div_raster) <- cell_avg
  
  return(div_raster)
  
}

# combine individual rasters stored in a list into a single multilayer raster
combine_raster_list <- function(raster_list){
  
  # get number of rasters in list
  n_rasters <- length(raster_list)
  
  # extract elements and combine to single raster
  # (unfortunately I think the only way to do this is recursively)
  multi_rast <- raster_list[[1]]
  for(n in 2:n_rasters){
    # combine multilayer raster with new layer
    multi_rast <- c(multi_rast, raster_list[[n]])
  }
  
  # set layer names as element names from list
  if(!is.null(names(raster_list))){
    names(multi_rast) <- names(raster_list)
  }
  
  return(multi_rast)
  
}

apply_sr_threshold <- function(diversity_raster, sr_mask, exclude_sr_lyr = TRUE){
  
  # exclude species richness values from masking, if requested
  if(exclude_sr_lyr == TRUE){
    
    # set aside species richness values
    sr_values <- terra::values(diversity_raster[["species_richness"]])
    
    # apply sr threshold mask to whole raster
    div_rast_masked <- diversity_raster
    terra::values(div_rast_masked) <- terra::values(diversity_raster) * sr_mask
    
    # reset sr values to their original
    terra::values(div_rast_masked[["species_richness"]]) <- sr_values
    
  } else {
    
    # apply sr threshold mask to whole raster
    div_rast_masked <- diversity_raster
    terra::values(div_rast_masked) <- terra::values(diversity_raster) * sr_mask
    
  }
  
  return(div_rast_masked)
  
  
}

# get list of species (with sex) included in diversity data
get_spp_sex <- function(div_data){
  
  # extract species and sex
  spp_sex <- t(sapply(strsplit(rownames(div_data), split = "-"), "[", 1:2))
  
  # set column names
  colnames(spp_sex) <- c("species", "sex")
  
  return(spp_sex)
  
}

# filter diversity raster by percentile
filter_div_raster <- function(div_rast, pcile = 75){
  
  # get percentile values
  pcile_vals <- terra::global(div_rast, fun = quantile, na.rm = TRUE, prob = pcile / 100)[[1]]
  names(pcile_vals) <- terra::names(div_rast)
  
  # set aside species richness values
  sr_values <- terra::values(div_rast[["species_richness"]])
  
  # apply percentile cutoff to each layer of raster
  for (layer in terra::names(div_rast)){
    terra::values(div_rast[[layer]])[terra::values(div_rast[[layer]]) < pcile_vals[layer]] <- NA
  }
  
  # reset sr values to their original
  terra::values(div_rast[["species_richness"]]) <- sr_values
  
  return(div_rast)
  
}

get_world <- function(new_crs, spatvec){
  
  # get world map from rnaturalearth
  world <- rnaturalearth::ne_countries(returnclass = "sf")
  
  # transform to correct crs
  world <- sf::st_transform(world, crs = new_crs)
  
  # transform to spatvector, if required (need to do this if using base R plot)
  if(spatvec == TRUE){
    
    # transform to spatvector
    world <- terra::vect(world)
    
  }
  
  return(world)
}

make_colour_breaks <- function(raster, nquants, plot_sr = FALSE){
  
  # calculate species richness quantiles, if required
  if(plot_sr == TRUE){
    sr_quants <- quantile(
      terra::values(raster[["species_richness"]]),
      na.rm = TRUE,
      probs = seq(0, 1, by = 1/nquants)
    )
    names(sr_quants) <- round(sr_quants, digits = 1)
  }
  
  # remove the species_richness layer
  layer_names <- terra::names(raster)
  sr_index <- which(layer_names == "species_richness")
  new_layer_names <- layer_names[-sr_index]
  raster <- raster[[new_layer_names]]
  
  # get quantiles from raster diversity values
  div_quants <- quantile(
    terra::values(raster), 
    na.rm = TRUE, 
    probs = seq(0, 1, by = 1/nquants)
    )
  
  # make names of quants a rounded version of the quantiles (for plotting)
  names(div_quants) <- round(div_quants, digits = 1)
  
  if(plot_sr == TRUE){
    quants <- list(div_breaks = div_quants, sr_breaks = sr_quants)
  } else {
    quants <- div_quants
  }
  
  return(quants)
  
}

# match extents of world and diversity rasters (necessary to plot in base R)
match_extents <- function(div_rast, world_spatvec){
  
  # match extent of diversity raster to that of world spatvector
  terra::ext(div_rast) <- terra::ext(world_spatvec)
  
  return(div_rast)
  
}

# plot individual layers of diversity raster on the same plot
plot_div_raster <- function(div_rast, 
                            world_spatvec,
                            div_metric,
                            scale_type = c("binned", "continuous"), 
                            div_breaks = NULL, sr_breaks = NULL, 
                            pal_choice = c("viridis", "turbo", "custom"), 
                            plot_sr = FALSE,
                            save_png = FALSE,
                            png_filename = NULL){
  
  # get names of raster layers
  layer_names <- terra::names(div_rast)
  
  # set up colour scale and palette
  if(scale_type == "binned"){
    
    nquants <- length(div_breaks) - 1
    
    # set colour palette for diversity metric
    if(pal_choice == "viridis"){
      divpal <- viridisLite::viridis(nquants)
    } else if(pal_choice == "turbo"){
      divpal <- viridisLite::turbo(nquants)
    } else if(pal_choice == "custom"){
      # create own colour palette (based on Cooney et al (2022) Fig. 2)
      cols <- c("#3e9eb5ff", "#eacc2cff", "#f82202ff")
      divpal <- colorRampPalette(cols)(nquants)
    }
    
    # reverse colour palette if using nn-count (as low numbers indicate more unusual colours for this metric)
    if(metric == "nn-count"){
      divpal <- rev(divpal)
    }
    
    # set colour palette for species richness, if plotting
    if(plot_sr == TRUE){
      
      nquants_sr <- length(sr_breaks) - 1
      
      if(palette_choice == "viridis"){
        srpal <- viridisLite::viridis(nquants_sr)
      } else if(palette_choice == "turbo"){
        srpal <- viridisLite::turbo(nquants_sr)
      } else if(palette_choice == "custom"){
        # create own colour palette (based on Cooney et al (2022) Fig. 2)
        cols <- c("#3e9eb5ff", "#eacc2cff", "#f82202ff")
        srpal <- colorRampPalette(cols)(nquants_sr)
      }
    }
    

  }
  
  # Set up basic legend title
  if(metric == "centr-dist"){
    lgd_metric <- paste0(avg_par, " distance\nto centroid ")
  } else if(metric == "nn-k"){
    lgd_metric <- paste0(avg_par, " distance\nto nearest neighbours ")
  } else if(metric == "nn-count"){
    lgd_metric <- paste0(avg_par, " number\nof nearest neighbours\nin radius r ")
  }
  
  
  # set up png saving parameters
  if(plot_sr == FALSE){
    png_width <- 420
    png_height <- 300
  } else if(plot_sr == TRUE){
    png_width <- 320
    png_height <- 400
  }
  
  # initialise png saving
  if(save_png == TRUE){
    png(
      here::here(
        "2_Patches", "4_OutputPlots", "3_Spatial_mapping", pam_res, space, pam_type,
        png_filename
      ), 
      width = png_width, height = png_height, units = "mm", pointsize = 24, res = 100,
    )
  }
  
  # set up plotting grid and remove species richness layer of raster, if not plotting
  if(plot_sr == FALSE){
    
    sr_index <- which(layer_names == "species_richness")
    new_layer_names <- layer_names[-sr_index]
    div_rast <- div_rast[[new_layer_names]]
    
    par(mfrow = c(2, 1))
    
  } else if(plot_sr == TRUE){
    
    par(mfrow = c(3, 1))
    
  }
  
  # set font size
  cexval <- 0.7
  
  # get updated names of layers, to apply plotting across
  layers <- terra::names(div_rast)
  
  # just use a for loop, for ease (could also do as an apply function)
  for(layer in layers){
    
    # set up colour breaks, palette and sex-specific legend title
    if(layer == "M"){
      breaks <- div_breaks
      colpal <- divpal
      lgd_title <- paste0(lgd_metric, "(males)")
    } else if(layer == "F"){
      breaks <- div_breaks
      colpal <- divpal
      lgd_title <- paste0(lgd_metric, "(females)")
    } else if(layer == "U"){
      breaks <- div_breaks
      colpal <- divpal
      lgd_title <- paste0(lgd_metric, "(unknown sex)")
    } else if(layer == "species_richness"){
      breaks <- sr_breaks
      colpal <- srpal
      lgd_title <- "species richness"
    }
    
    
    # first plot male diversity raster to get legend only
    terra::plot(
      div_rast, layer, 
      breaks = breaks, 
      col = colpal, 
      #  colNA = "grey",
      legend = TRUE,
      frame = FALSE,
      axes = FALSE, 
      plg = list(title = lgd_title, x = "left", cex = cexval),
      mar = c(1, 2, 1, 2)
    )
    # overlay grey world map so that NA parts of world are grey
    terra::plot(
      world_spatvec, 
      col="grey80", 
      border = NA, 
      add = TRUE,
      axes = "n"
    )
    # overlay male diversity raster again, but without legend
    terra::plot(
      div_rast, layer, 
      breaks = breaks, 
      col = colpal, 
      # colNA = "grey",
      add = TRUE,
      legend = FALSE,
      mar = c(1, 2, 1, 2)
    )
    
  }
  
  # end png saving, if initialised
  if(save_png == TRUE){
    dev.off()
  }
  
  
}


# Load data ----

# set diversity data filename
div_data_filename <- here::here(
  "2_Patches", "3_OutputData", "4_Diversity_measures", space,
  paste(clade, "patches", space, sex_match, "diversitymetrics.csv", sep = "_")
  )
# read in data
div_data <- read.csv(div_data_filename, row.names = 1, check.names = FALSE)

# load PAM
pam_filename <- paste0("PAM_birds_Behrman", pam_res, "_Pres12_Orig12_Seas12_", pam_type, pam_seas, ".rds")
pam_raw <- readRDS(
  paste(pams_filepath, pam_filename, sep = "/")
)



# Data preparation ----

# filter to top percentile of data, if required
if(sift_div_data == TRUE){
  div_data <- filter_div_data(div_data, metric, separate_sexes = TRUE, pcile = 75)
}

# extract null raster from PAM
null_rast <- extract_null_rast(pam_raw)

# extract matrix values from PAM file
pam <- extract_pam_vals(pam_raw)

# remove raw PAM file (for RAM reasons)
rm(pam_raw)
gc()

# change "Genus species" style to "Genus_species" style in PAM colnames
colnames(pam) <- gsub(" ", "_", colnames(pam))

# get list of species included in diversity data
spp_list <- get_unique_spp(div_data)

# subset pam to only species in diversity data
pam <- subset_pam(pam, spp_list)

# subset diversity data to only species present in PAM
div_data <- subset_data_by_spp(div_data, colnames(pam))

# generate individual list of sex-specific species values of metric of interest (keep or discard sexes
# depending on value of sex_match)
# first set sex list to iterate over
sexes <- set_sex_list(sex_interest)
# now apply a function to extract named vector lists of diversity values for each sex individually
sexed_metric_list <- lapply(sexes, extract_sex_vals, div_data = div_data, metric = metric)

# name the elements of the list
names(sexed_metric_list) <- sexes

# make raster of average diversity in each grid cell - for each sex individually
sexed_div_rasters <- lapply(
  sexes, make_diversity_raster, 
  # hyperparameters
  sexed_metric_vals = sexed_metric_list,
  pam = pam,
  null_rast = null_rast,
  avg_type = avg_par
  )
names(sexed_div_rasters) <- sexes


# get species richness by grid cell
species_richness <- calc_species_richness(pam)

# create a mask for cells with species richness >= threshold
sr_mask <- make_sr_mask(species_richness, sr_threshold)

# save mask for later use
sr_mask_filename <- paste("species_richness_mask_threshold", sr_threshold, sub(".rds", "", pam_filename, fixed = TRUE), clade, sex_match, sep = "_")
saveRDS(
  sr_mask, file = here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", pam_res, space, "sr_masks",
    paste0(sr_mask_filename, ".rds")
  )
)

# create raster of species richness (for plotting later)
sr_raster <- make_sr_raster(species_richness, null_rast)

# first add SR raster to the list of diversity rasters
sexed_div_rasters[["species_richness"]] <- sr_raster

# combine the separate sexed rasters and the species richness raster into one multilayer raster with a layer for each sex
diversity_raster <- combine_raster_list(sexed_div_rasters)

# save raster
div_raster_filename <- paste(clade, "patches", space, sex_match, metric, avg_par, pam_res, "Behrman", pam_type, pam_seas, "raster.tif", sep = ".")
terra::writeRaster(
  diversity_raster,
  filename = here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", pam_res, space,
    div_raster_filename
  )
)

# Plotting ----

# clear environment except variables used to generate raster filename
rm(list=setdiff(ls(), c("clade", "space", "sex_match", "metric", "avg_par", "pam_res", "pam_type", "pam_seas", "sr_threshold", "sift_rast_data", "col_scale_type", "nquants", "palette_choice", "plot_sr",  "apply_sr_threshold", "filter_div_raster", "get_world", "make_colour_breaks", "match_extents", "plot_div_raster")))

# load diversity raster
div_raster_filename <- paste(clade, "patches", space, sex_match, metric, avg_par, pam_res, "Behrman", pam_type, pam_seas, "raster.tif", sep = ".")
diversity_raster <- terra::rast(
  here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", pam_res, space,
    div_raster_filename
  )
)

# load species richness mask
sr_mask_filename <- paste0("species_richness_mask_threshold_", sr_threshold, "_PAM_birds_Behrman", pam_res, "_Pres12_Orig12_Seas12_", pam_type, pam_seas, "_", clade, "_", sex_match, ".rds")
sr_mask <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", pam_res, space, "sr_masks",
    sr_mask_filename
  )
)

# apply species richness threshold to diversity raster (diversity metric layers only)
diversity_raster_masked <- apply_sr_threshold(diversity_raster, sr_mask, exclude_sr_lyr = TRUE)

# sift data to only the upper quartile of each layer
if(sift_rast_data == TRUE){
  diversity_raster_masked <- filter_div_raster(diversity_raster_masked)
}

# get world map and transform to same CRS as diversity raster (for plotting)
world <- get_world(new_crs = terra::crs(diversity_raster), spatvec = TRUE)

# match extent of diversity raster to world map extent (need to do this to plot in base R)
diversity_raster_masked <- match_extents(div_rast = diversity_raster_masked, world_spatvec = world)

# set up colour scale breaks as quantiles
colour_breaks <- make_colour_breaks(diversity_raster_masked, nquants, plot_sr = plot_sr)

# plot

# set png filename, if saving
png_filename <- paste(clade, "patches", space, sex_match, metric, "sr_thresh", sr_threshold, col_scale_type, palette_choice, "colscale", pam_res, "Behrman", pam_type, pam_seas, "png", sep = ".")

plot_div_raster(diversity_raster_masked, 
                world_spatvec = world, 
                div_metric = metric,
                scale_type = col_scale_type, pal_choice = palette_choice,
                plot_sr = TRUE,
                div_breaks = colour_breaks$div_breaks, sr_breaks = colour_breaks$sr_breaks,
                save_png = TRUE,
                png_filename = png_filename)

