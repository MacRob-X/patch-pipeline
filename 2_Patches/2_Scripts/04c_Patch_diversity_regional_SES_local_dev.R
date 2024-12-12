## Calculate SES for diversity change in each bioregion separately
## Robert MacDonald
## 15th November 2024
## 
## NON-FUNCTIONALISED VERSION (or minimally functionalised version)
## LOCAL DIVERSITY LOSS ONLY


# clear environment
rm(list=ls())

# Load libraries
library(raster)
library(parallel)
library(dispRity)


# Load shared functions
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

# temporary location for the SES calculating functions
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "regional_ses.R"
  )
)

## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neoaves"
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
# select whether to calculate local or global diversity loss (i.e., mean distance to local or global centroid)
div_loss_type <- "local"
# select number of null distributions to generate
n_sims <- 1000
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
pam_res <- "100km"
# enter PAM files location
pams_filepath <- "X:/cooney_lab/Shared/Rob-MacDonald/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/PAMs"
#pams_filepath <- "X:/cooney_lab/Shared/Rob-MacDonald/SpatialData/BirdLife/BirdLife_Shapefiles_v9/PAMs/100km/Behrmann_cea/"
# # use ecoregions or biomes?
regions <- "biomes"

## Plotting parameters
# Select whether to use binned (based on quantiles) or continuous colour scale
col_scale_type <- "binned"
# If binned, choose the number of quantiles to use
nquants <- 10
# Also plot species richness map?
plot_sr <- TRUE
## Use viridis_c, viridis turbo palette or custom colour palette?
palette_choice <- "viridis"


# Functions ----
# make wrapper function for averaging (allows selection of average type e.g. mean, median)
avg <- function(vals, avg_type, na.rm = TRUE){
  return(get(avg_type)(vals, na.rm = na.rm))
}


# Load data ----

# set PCA filename
pca_filename <- paste(clade, "patches.231030.PCAcolspaces", "rds", sep = ".")
# Load PCA (values only)
pca_dat <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]]$x

# load IUCN Red List data
iucn_filename <- paste0("iucn_2024_", iucn_type, ".csv")
iucn <- read.csv(
  here::here(
    "4_SharedInputData", iucn_filename
  )
)[, c("species_birdtree", "iucn_cat")]

# load region data (either Dinerstein et al 2017 ecoregion data or biome data derived from same)
region_shapes <- load_regions(regions)

# load PAM
pam_filename <- paste0("PAM_birds_Behrman", pam_res, "_Pres12_Orig12_Seas12_", pam_type, pam_seas, ".rds")
pam_raw <- readRDS(
  paste(pams_filepath, pam_filename, sep = "/")
)

# Analysis ----


# extract null raster from PAM file
null_rast <- extract_null_rast(pam_raw)

# extract matrix values from PAM file
pam <- extract_pam_vals(pam_raw)

# remove raw PAM file (for RAM reasons)
rm(pam_raw)

# change "Genus species" style to "Genus_species" style in PAM colnames
colnames(pam) <- gsub(" ", "_", colnames(pam))

# add species and sex columns (from rownames)
# only necessary if we're subsetting in some way
if(sex_match != "all"){
  pca_dat <- pca_spp_sex(pca_dat)
}

# if required, get species for which we have matched male/female data
if(sex_match == "matchedsex"){
  
  # get species to keep
  spp_to_keep <- sex_match_fun(pca_dat)
  
  # subset data to species of interest
  pca_dat <- subset_sex_match(pca_dat, spec_to_keep = spp_to_keep)
  
}

# get list of species included in diversity/PCA data
spp_list <- get_unique_spp(pca_dat)

# subset pam to only species in diversity data
pam <- subset_pam(pam, spp_list)

# subset diversity pca_dat to only species present in PAM
pca_dat <- subset_data_by_spp(pca_dat, colnames(pam))

# convert ecoregions to Behrman CEA CRS
region_shapes <- sf::st_transform(region_shapes, crs = terra::crs(null_rast, proj = TRUE))

# generate individual list of sex-specific species values of metric of interest (keep or discard sexes
# depending on value of sex_match)
# first set sex list to iterate over
sexes <- set_sex_list(sex_interest)

# now apply a function to extract named vector lists of values for each sex individually
sexed_pca_list <- lapply(sexes, extract_sex_vals, div_data = pca_dat)
names(sexed_pca_list) <- sexes




# Calculate SESs for each region ----

# Attempted parallelisation ----------------------------------------------------------

# Note that parallelisation doesn't work here because of RAM limitations - need to export the full global 
# environment to each worker node (could perhaps be overcome by using only two nodes?)
# for parallelisation, set up a cluster using the number of cores (4 less than total number of laptop cores)
# Note that it might be more sensible to initialise the cluster within the calc_region_ses function itself
# - this will increase the overhead time, but maybe will be less heavy in terms of RAM usage
# - the above is incorrect, it's just as RAM-hungry. This might be a problem with the full dataset

# set up parallel cluster (if using)
# no_cores <- parallel::detectCores() - 4
no_cores <- 2
cl <- parallel::makeCluster(no_cores)

wrapper_fn <- function(sex){
   
   pca_sexed <- sexed_pca_list[[sex]]
   
   # apply function across each region (row) in sf object
   sexed_tmp_results <- lapply(1:nrow(region_shapes), function(region_number) {
     
     region_sf <- region_shapes[region_number, ]
     
     calc_region_ses(
       region_sf = region_sf,
       pca_data = pca_sexed,
       pam = pam,
       null_raster = null_rast,
       iucn = iucn,
       metric = metric,
       n_sims = n_sims,
       parallel_run = FALSE,
       append_sf = TRUE,
       # cluster = cl,
       regions = regions
     )
     
     
   })
   
   return(sexed_tmp_results)
   
 }
 
 
 parallel::clusterExport(cl, c("sexes", "wrapper_fn", "sexed_pca_list", "region_shapes", "calc_region_ses", "regions", "metric", "null_rast", "iucn", "n_sims", "dispRity", "calc_iucn_ses"), envir = .GlobalEnv)
 
 parallel::clusterEvalQ(cl, {
   library(Rcpp)
   library(sf)
   library(dispRity)
   # If other packages with Rcpp dependencies are used, load them here as well
 })
 
 results <- parallel::parLapply(cl, sexes, function(sex) {
   try(wrapper_fn(sex))
 })
 results <- parallel::parLapply(cl = cl, sexes, wrapper_fn)
 
 parallel::stopCluster(cl)
 
 # End attempted parallelisation --------------------------------------------------
# ------------------------------------------------------------

# Non-parallelised version (working)
# apply function across each sexed dataset
results <- lapply(sexes, function(sex){
  
  pca_sexed <- sexed_pca_list[[sex]]
  
  # apply function across each region (row) in sf object
  lapply(1:nrow(region_shapes), function(region_number) {
    
    region_sf <- region_shapes[region_number, ]
    
    calc_region_ses(
      region_sf = region_sf,
      pca_data = pca_sexed,
      pam = pam,
      null_raster = null_rast,
      iucn = iucn,
      metric = metric,
      n_sims = n_sims,
      parallel_run = FALSE,
      append_sf = TRUE,
      # cluster = cl,
      regions = regions
    )
  })
  
})

# name elements of the results list
names(results) <- sexes


# for each sex, process results into dfs to bind to sf data
results_dfs <- lapply(sexes, function(sex) {
  
  # get individual sexed result
  sexed_results <- results[[sex]]
  
  # bind results together into a single dataframe
  sexed_results <- as.data.frame(do.call(rbind, sexed_results))
  
  # adjust column names to specify sex
  colnames(sexed_results) <- paste(sex, colnames(sexed_results), sep = "_")
  
  # convert NaN values to NA (where there are e.g. no CR species to be lost)
  sexed_results <- data.frame(lapply(sexed_results, gsub, pattern = NaN, replacement = NA, fixed = TRUE))
  
})


# identify rows to be removed (regions which contain no species or only a single species in our dataset)
rows_to_keep <- which(!(results_dfs[[1]][, 1] == "nospec_region" | results_dfs[[1]][, 1] == "singlespec_region"))

# bind results to sf data
region_results <- cbind(region_shapes, do.call(cbind, results_dfs))

# remove regions for which there are no species
region_results <- region_results[rows_to_keep, ]

# rename results df elements (so the plotting works)
names(results_dfs) <- sexes



# plot results ----

# load results (for plotting)
results_filename <- paste0(paste(clade, space, "regionalSESplotdata", metric, avg_par, regions, paste0(n_sims, "sims"), paste0(iucn_type, "iucn"), sex_match, pam_res, pam_type, pam_seas, sep = "_"), ".rds")
results_dfs <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", "2_Regional_diversity_SES", pam_res, space,
    results_filename
  )
)

# remove names (for processing)
names(results_dfs) <- NULL

# identify rows to be removed (regions which contain no species or only a single species in our dataset)
rows_to_keep <- which(!(results_dfs[[1]][, 1] == "nospec_region" | results_dfs[[1]][, 1] == "singlespec_region"))

# bind results to sf data
region_results <- cbind(region_shapes, do.call(cbind, results_dfs))

# remove regions for which there are no species
region_results <- region_results[rows_to_keep, ]

# rename results df elements (so the plotting works)
sexes <- set_sex_list(sex_interest)
names(results_dfs) <- sexes

# convert to numerics
column_names <- unlist(lapply(sexes, function(sex){
  return(colnames(results_dfs[[sex]]))
}))
region_results[, column_names] <- sapply(sf::st_drop_geometry(region_results[, column_names]), as.numeric)

## PRODUCE RASTER OF RESULTS - layer for each attribute
# produce null raster with a layer for each attribute (e.g., M_ses_CR, M_ses_EN etc.)
# get null raster from raw pam
# extract null raster from PAM file
null_rast <- extract_null_rast(pam_raw, rast_type = "raster")
results_raster <- fasterize::fasterize(region_results, null_rast, field = c("M_ses_CR"))
results_raster <- lapply(column_names, fasterize::fasterize, sf = region_results, raster = null_rast)
# convert elements to spatrasters
results_raster <- lapply(results_raster, terra::rast)
results_raster <- terra::rast(results_raster)
names(results_raster) <- column_names
# remove the SES of all species (which are always 0 by definition)
results_raster <- subset(results_raster, paste0(sexes, "_ses_all"), negate = TRUE)
# remove raw PAM and collect garbage
rm(pam_raw)
gc()


# produce thresholded raster where all values < -2 are set to -2
thresh_rast <- clamp(results_raster, lower = -2, values = TRUE)
# extract values
new_vals <- terra::values(thresh_rast)
# set all values which aren't -2 to NA (so we can )
new_vals[new_vals != -2] <- NA
# now set as raster values
terra::values(thresh_rast) <- new_vals


plot(results_raster[[1]])
terra::polys(terra::as.polygons(thresh_rast[[1]]), border = "black", lwd = 1.5)

# plot all layers - THIS CURRENTLY WORKS but with individual legends
layers <- names(results_raster)
n_cols <- length(sexes)
n_rows = length(layers) %/% n_cols
par(mfcol = c(n_rows, n_cols))
for(layer in layers){
  
  # plot raster of values
  plot(results_raster[[layer]], main = layer)
  
  # plot polygon of values < 2
  terra::polys(terra::as.polygons(thresh_rast[[layer]]), border = "black", lwd = 1.5)
  
}

# same, but with layout matrix
png("test_image.png", width = 1450, height = 1550, res = 72, pointsize = 20)
# First define range for legend and colour scale (make it very slightly higher/lower
# than true range to stop the max/min regions being whited out)
leg_range <- c(min(terra::values(results_raster), na.rm = TRUE), 
              max(terra::values(results_raster), na.rm = TRUE))
leg_increm <- (max(leg_range) - min(leg_range)) / 100
leg_range <- c(min(leg_range) - leg_increm,
               max(leg_range) + leg_increm)
# define colour scale
leg_fill <- viridis::viridis(100, direction = -1)
# define legend title
leg_title <- paste0("SES value (", metric, ")")

layers <- names(results_raster)
layout_matrix <- matrix(
  c(1, 5, 
    2, 6, 
    3, 7, 
    4, 8, 
    9, 9), 
  nrow = 5, ncol = 2, byrow = TRUE
  )
layout(
  heights = c(1, 1, 1, 1, 0.5),
  mat = layout_matrix
  )

# generate A-H labels for plots
labels <- LETTERS[1:8]
names(labels) <- layers

# add map plots (no legends)
for(layer in layers){
  
  # plot raster of values
  plot(
         results_raster[[layer]],         # layer to plot
         main = labels[layer], loc.main = c(-16500,7500), # title and position
         legend = FALSE,                  # no legend
         col = leg_fill,                  # colour scale
         range = leg_range,               # colour scale range
         mar = c(1.1, 1.2, 2.1, 1.2),
         pax = list(tick = 0, lab = 0)    # no ticks or labels on axes
       )
  
  # plot polygon of values < 2
  terra::polys(terra::as.polygons(thresh_rast[[layer]]), 
               col = "#f28383", density = 10, angle = 45,
               border = "#f28383", lwd = 1)
  
}
# add combined legend
par(mar = c(8, 12, 1, 12))  # Reduce margins for the legend
image(x = seq(leg_range[1], leg_range[2], length.out = 100), 
      y = 1, 
      z = matrix(seq(leg_range[1], leg_range[2], length.out = 100), ncol = 1), 
      col = leg_fill, 
      axes = FALSE, 
      xlab = "", 
      ylab = "")
axis(1, at = seq(leg_range[1], leg_range[2], length.out = 5), 
     labels = round(seq(leg_range[1], leg_range[2], length.out = 5), 2),
     cex.axis = 1.4)
mtext(leg_title, side = 1, line = 3, cex = 1)

dev.off()




plot(1:10, 1:10, type = "n", axes = FALSE, xlab = "", ylab = "")  # Placeholder for legend plot
image(1, seq(min(leg_range), max(leg_range), length.out = 10), 
      t(seq_along(leg_range)), col = leg_fill, axes = FALSE)
axis(1, at = seq(min(leg_range), max(leg_range), length.out = 5), 
     labels = round(seq(min(leg_range), max(leg_range), length.out = 5), 2))
mtext(leg_title, side = 1, line = 2)
image(z = leg_vals, t(seq_along(leg_range)), col = leg_fill, axes = FALSE)



plot(results_raster[[1]], axes = FALSE, legend = TRUE, plg = c(x = "top"))
legend("bottom", legend = paste0("SES value (", metric, ")"),)
add_legend
plot(results_raster[[1]], type = "n", axes = FALSE, legend = TRUE,  size = c(1, 1))



# plot combined legend
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
terra::add_legend("bottom", legend = "SES")
legend("bottom", legend = "SES")
legend(x = "top",inset = 0,
       legend = "SES", 
       col=plot_colors, lwd=5, cex=.5, horiz = TRUE)


# produce polygon of regions for which SES <= -2
results_trunc <- results_raster
values(results_trunc) <- ifelse(values(results_raster) <= -2, values(results_raster), NA)
loss_polys <- lapply(names(results_trunc), function(layer){
  
  polygon <- terra::as.polygons(results_trunc[[layer]])
  
})

# produce tidy df of boundaries - NOTE: THIS DOESN'T CURRENTLY WORK FOR WAHT i WANT IT TO DO (create tidy df of polygons)
polys_df <- lapply(1:length(loss_polys), function(polygon){
  polygon_df <- as.data.frame(loss_polys[[polygon]])
  polygon_df$lyr <- names(results_raster)[polygon]
  return(polygon_df)
})
polys_df <- dplyr::bind_rows(polys_df)

# produce tidy df of SES values
results_tidy_df <- as.data.frame(results_raster, xy = TRUE)
results_tidy_df <- tidyr::pivot_longer(results_tidy_df, cols = names(results_raster), names_to = "lyr", values_to = "value")

# plot with boundaries of SES <= 2
ggplot() + 
  geom_raster(data = results_tidy_df, aes(x = x, y = y, fill = value)) +
 # geom_path(data = polys_df, aes(x = x, y = y, group = group), color = "white") +
  facet_wrap(~lyr) +
  scale_fill_viridis_c()

# plot raster
library(ggplot2)
library(tidyterra)

# THIS CURRENTLY ALMOST WORKS but plots the threshold polygon of the first layer multiple times,
# rather than that of each layer
ggplot() + 
  geom_spatraster(data = results_raster) + 
  geom_spatvector(data = terra::as.polygons(thresh_rast), colour = "white", fill = NA) + 
 # geom_sf(data = loss_polys[[1]], colour = "white", fill = NA) + 
  facet_wrap(~lyr, ncol = 2, dir = "v") + 
  scale_fill_viridis_c(direction = -1, na.value = "transparent")

# create individual polygons of each layer
layers <- names(results_raster)
thresh_polygons <- lapply(layers, function(layer){
  
  polygon <- terra::as.polygons(thresh_rast[[layer]])
  
  # ensure correct CRS
  polygon <- terra::project(polygon, crs(results_raster))
  
  # ensure correct extent
  polygon <- terra::crop(polygon, terra::ext(results_raster))
  
  return(polygon)
  
})
names(thresh_polygons) <- layers
# combine into spatvectorcollection
thresh_polygons_svc <- terra::svc(thresh_polygons)

# plot with base R - THIS ALMOST WORKS - plotting each individually works, but when you set add = TRUE
# the overlay fails for some reason



terra::plot(results_raster)
par(new = TRUE, bg = "transparent", mfrow = c(n_rows, n_cols))
terra::plot(thresh_polygons_svc, box = FALSE, axes = FALSE, main = "")
terra::plot(thresh_polygons_svc, box = FALSE, axes = FALSE, add = TRUE)

########## ALMOST WORKING - 2024-12-10

# TRY SAVING AND OVERLAYING
# set number of columns and rows
n_cols <- length(sexes)
n_rows = length(layers) %/% n_cols

# set legend range
leg_range = c(max(terra::values(results_raster), na.rm = TRUE), min(terra::values(results_raster), na.rm = TRUE))

# save rasters png
png("rasters.png", width = 2000, height = 1000)
par(mfrow = c(n_rows, n_cols))
terra::plot(results_raster, mar = c(0,0,0,0), buffer = FALSE)
dev.off()
# save polygons png with transparent background
png("polygons.png", bg = "transparent", width = 2000, height = 1000)
par(mfrow = c(n_rows, n_cols))
terra::plot(thresh_polygons_svc, box = TRUE, axes = FALSE, main = "", buffer = FALSE, mar = c(0,0,0,0))
dev.off()
# overlay polygons png on rasters png
terrainr::combine_overlays(
  "rasters.png",
  "polygons.png",
  output_file = "rasters_and_polygons.png"
  
)

ggplot() + 
  geom_spatvector(data = thresh_polygons, colour = "black", fill = NA) + 
  facet_wrap(~lyr, ncol = 2, dir = "v")

# wrapper function to feed into apply
plot_layer <- function(lyr){
  
  p <- ggplot() + 
    geom_sf(data = region_results, aes(fill = as.numeric(get(lyr)), colour = as.numeric(get(lyr)))) + 
    scale_fill_viridis_c(name = lyr, option = "viridis", direction = -1) + 
    scale_colour_viridis_c(name = lyr, option = "viridis", direction = -1)
  p + labs(fill = lyr)
  return(p)
}
layers <- colnames(results_dfs[["F"]])[2:5]

plots <- lapply(layers, plot_layer)

png(paste(clade, regions, "ses_F.png", sep = "_"), width = 2000, height = 1000)
gridExtra::grid.arrange(grobs = plots, nrow = 2)
dev.off()


# Format results

# save intermediate results
results_filename <- paste0(paste(clade, space, "regionalSESplotdata", metric, avg_par, regions, paste0(n_sims, "sims"), paste0(iucn_type, "iucn"), sex_match, pam_res, pam_type, pam_seas, sep = "_"), ".rds")
saveRDS(
  results_dfs, 
  here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", "2_Regional_diversity_SES", pam_res, space,
    results_filename
    )
  )

