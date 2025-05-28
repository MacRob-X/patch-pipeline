# Assemblage-level spatial mapping
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


## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neoaves"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
# N.B. will need to change the date in the pca_all filename if using usmldbl or usmldblr
# (from 240603 to 240806)
space <- "ab"
# Select number of PC axes to retain ("all" or a number)
axes <- "all"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# select sex of interest ("all", "male_female", "male_only", "female_only", "unknown_only")
sex_interest <- "male_female"
# select metric ("centr-dist", "nn-k", "nn-count", "sum.variances", "sum.ranges", "convhull.volume")
metric <- "centr-dist"
# select type of averaging to use ("mean" or "median")
avg_par <- "median"
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
nquants <- 12
# Also plot species richness map?
plot_sr <- TRUE
## Use viridis_c, viridis turbo palette or custom colour palette?
palette_choice <- "viridis"

## END EDITABLE CODE ##




# Load data ----

# Load PCA data
pca_filename <- paste(clade, sex_match, "patches.231030.PCAcolspaces", "rds", sep = ".")
pca_dat <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]][["x"]]

# load PAM
pam_filename <- paste0("PAM_birds_Behrman", pam_res, "_Pres12_Orig12_Seas12_", pam_type, pam_seas, ".rds")
pam_raw <- readRDS(
  paste(pams_filepath, pam_filename, sep = "/")
)

# Data preparation ----

# restrict PCA to requested number of axes
if(axes != "all"){
  pca_dat <- pca_dat[, 1:axes]
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

# get list of species included in PCA data
spp_list <- get_unique_spp(pca_dat)

# subset pam to only species in diversity data
pam <- subset_pam(pam, spp_list)

# subset PCA data to only species present in PAM
pca_dat <- subset_data_by_spp(pca_dat, colnames(pam))

# generate individual list of sex-specific species values of metric of interest (keep or discard sexes
# depending on value of sex_match)
# first set sex list to iterate over
sexes <- set_sex_list(sex_interest)
# now apply a function to extract named vector lists of diversity values for each sex individually
sexed_metric_list <- lapply(sexes, extract_sex_vals, div_data = pca_dat, metric = metric, assemblage_level = TRUE)

# name the elements of the list
names(sexed_metric_list) <- sexes

# calculate value of chosen metric for each occupied grid cell in PAM (i.e., each non-NA row)
# apply across all sexes

div_values <- lapply(
  sexed_metric_list,
  function(x) apply(pam, 1, calc_metric_gridcell, pca_data = x, metric = metric, avg_par = avg_par, min_species = sr_threshold)
)

# div_values <- apply(pam, 1, calc_metric_gridcell, pca_data = sexed_metric_list$M, metric = metric, avg_par = avg_par, min_species = sr_threshold)
# 
# for(i in 2000:2010){(calc_metric_gridcell(pam[i,], pca_data =  sexed_metric_list$M, metric = metric, avg_par = avg_par, min_species = sr_threshold))}

# make a list of individual rasters of assemblage-level diversity values
div_raster <- lapply(div_values, make_assdiv_raster, null_rast = null_rast)

# get species richness by grid cell
species_richness <- calc_species_richness(pam)

# create a mask for cells with species richness >= threshold
sr_mask <- make_sr_mask(species_richness, sr_threshold)

# save mask for later use
sr_mask_filename <- paste("species_richness_mask_threshold", sr_threshold, sub(".rds", "", pam_filename, fixed = TRUE), clade, sex_match, sep = "_")
space_data_path <-  here::here(
  "2_Patches", "3_OutputData", "6_Spatial_mapping", "3_Assemblage_level_mapping", pam_res, space, "sr_masks")
if(!dir.exists(space_data_path)){
  dir.create(space_data_path, recursive = TRUE)
}
saveRDS(
  sr_mask, file = here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", "3_Assemblage_level_mapping", pam_res, space, "sr_masks",
    paste0(sr_mask_filename, ".rds")
  )
)

# create raster of species richness (for plotting later)
sr_raster <- make_sr_raster(species_richness, null_rast)

# first add SR raster to the list of diversity rasters
div_raster[["species_richness"]] <- sr_raster

# combine into a single multilayer raster
div_raster <- combine_raster_list(div_raster)

# save raster
div_raster_filename <- paste("assemblvl", clade, "patches", space, sex_match, avg_par, metric, "sr_thresh", sr_threshold, avg_par, pam_res, "Behrman", pam_type, pam_seas, "raster.tif", sep = ".")
terra::writeRaster(
  div_raster,
  filename = here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", "3_Assemblage_level_mapping", pam_res, space,
    div_raster_filename
  ),
)

# Plotting ----

# clear environment except variables used to generate raster filename
rm(list=setdiff(ls(), c("clade", "space", "sex_match", "metric", "avg_par", "pam_res", "pam_type", "pam_seas", "sr_threshold", "sift_rast_data", "col_scale_type", "nquants", "palette_choice", "plot_sr",  "apply_sr_threshold", "filter_div_raster", "get_world", "make_colour_breaks", "match_extents", "plot_div_raster")))

# load diversity raster
div_raster_filename <- paste("assemblvl", clade, "patches", space, sex_match, avg_par, metric, "sr_thresh", sr_threshold, avg_par, pam_res, "Behrman", pam_type, pam_seas, "raster.tif", sep = ".")
div_raster <- terra::rast(
  here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", "3_Assemblage_level_mapping", pam_res, space,
    div_raster_filename
  )
)

# load species richness mask
sr_mask_filename <- paste0("species_richness_mask_threshold_", sr_threshold, "_PAM_birds_Behrman", pam_res, "_Pres12_Orig12_Seas12_", pam_type, pam_seas, "_", clade, "_", sex_match, ".rds")
sr_mask <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", "3_Assemblage_level_mapping", pam_res, space, "sr_masks",
    sr_mask_filename
  )
)

# apply species richness threshold to diversity raster (diversity metric layers only)
div_raster_masked <- apply_sr_threshold(div_raster, sr_mask, exclude_sr_lyr = TRUE)

# sift data to only the upper quartile of each layer
if(sift_rast_data == TRUE){
  div_raster_masked <- filter_div_raster(div_raster_masked)
}

# get world map and transform to same CRS as diversity raster (for plotting)
world <- get_world(new_crs = terra::crs(div_raster), spatvec = TRUE)

# match extent of diversity raster to world map extent (need to do this to plot in base R)
div_raster_masked <- match_extents(div_rast = div_raster_masked, world_spatvec = world)

# set up colour scale breaks as quantiles
colour_breaks <- make_colour_breaks(div_raster_masked, nquants, plot_sr = plot_sr)

# plot

# set png filename, if saving
png_filename <- paste(clade, "patches", sex_match, avg_par, metric, "sr_thresh", sr_threshold, col_scale_type, palette_choice, "colscale", "Behrman", pam_seas, "png", sep = ".")

space_plot_path <-  here::here(
  "2_Patches", "4_OutputPlots", "3_Spatial_mapping", "3_Assemblage_level_mapping", pam_res, space, pam_type)
if(!dir.exists(space_plot_path)){
  dir.create(space_plot_path, recursive = TRUE)
}

plot_div_raster(div_raster_masked, 
                world_spatvec = world, 
                div_metric = metric,
                scale_type = col_scale_type, pal_choice = palette_choice,
                plot_sr = plot_sr,
                col_breaks = colour_breaks,
                save_png = TRUE,
                png_filename = png_filename,
                assemb_level = TRUE)



# dev ----

plot(log(div_raster))

# create df for ggplot of latitudinal gradient
div_df <- terra::as.data.frame(div_raster[[(c("M", "F"))]], xy = TRUE, wide = FALSE)
# div_df <- tidyr::pivot_longer(div_df, cols = c("M", "F"), names_to = "Layer", values_to = "Centroid_distance")
colnames(div_df) <- c("longitude", "latitude", "sex", "centroid_distance")

# plot latitudinal gradient - all values
div_df %>% 
  ggplot(aes(x = latitude, y = centroid_distance, colour = sex)) + 
  geom_point() + 
  facet_wrap(~ sex) + 
  theme_minimal()

# plot latitudinal gradient - medians
div_df %>% 
  group_by(latitude, sex) %>% 
  summarise(median_value = median(centroid_distance, na.rm = TRUE), .groups = "drop") %>% 
  ggplot(aes(y = median_value, x = latitude, colour = sex)) + 
 # geom_line() + 
  geom_point() + 
  facet_wrap(~ sex, scales = "free_x") + 
  theme_minimal()

# ridgeline plots - multiple histograms (one for each latitude value)
library(ggridges)

# first bin into latitude bins as too many to do full latitude plot
div_df_binned <- div_df %>%
  mutate(lat_bin = round(latitude / 300) * 300)  # Bin to nearest 300 units

# plot ridgeline plot
p_binned <- ggplot(div_df_binned, aes(x = centroid_distance, y = factor(lat_bin), fill = sex)) +
  geom_density_ridges(alpha = 0.7, scale = 2, rel_min_height = 0.01) +
  facet_wrap(~ sex, ncol = 2) +
  labs(
    title = paste0("Distribution of Centroid Distance (", space, ") by Latitude (Binned)"),
    subtitle = "Ridgeline plot with latitude binned to nearest 300 units",
    x = "Centroid Distance",
    y = "Latitude (Binned)",
    fill = "Sex"
  ) +
  theme_ridges() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "none"
  )

# Display the binned version
print(p_binned)
