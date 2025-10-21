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

# temporary location for the SES calculating functions
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "regional_ses.R"
  )
)

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
# select metric ("centr-dist", "nn-k", "nn-count", "sum.variances", "sum.ranges")
metric <- "centr-dist"
# select type of averaging to use ("mean" or "median")
avg_par <- "mean"
# select whether to calculate local or global diversity loss (i.e., mean distance to local or global centroid)
div_loss_type <- "local"
# select number of null distributions to generate
n_sims <- 2
# select whether to exclude grid cells with species richness below a certain threshold (e.g. 5)
# set as 0 if no threshold wanted
sr_threshold <- 5
# select whether the full results (e.g. SES, SR, raw diversity) or only SES results (for
# plotting) are to be returned
ses_only <- FALSE
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
pam_res <- "200km"
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
pca_filename <- paste(clade, sex_match, "patches.250716.PCAcolspaces", "rds", sep = ".")
# Load PCA (values only)
pca_dat <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", clade, "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]]$x

# load IUCN Red List data
iucn_filename <- paste0("neognath_iucn_2024_", iucn_type, ".csv")
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
if(sex_match != "allspecimens"){
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
sexed_pca_list <- lapply(sexes, extract_sex_vals, div_data = pca_dat, assemblage_level = TRUE)
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
# no_cores <- 2
# cl <- parallel::makeCluster(no_cores)
# 
# wrapper_fn <- function(sex){
#    
#    pca_sexed <- sexed_pca_list[[sex]]
#    
#    # apply function across each region (row) in sf object
#    sexed_tmp_results <- lapply(1:nrow(region_shapes), function(region_number) {
#      
#      region_sf <- region_shapes[region_number, ]
#      
#      calc_region_ses(
#        region_sf = region_sf,
#        pca_data = pca_sexed,
#        pam = pam,
#        null_raster = null_rast,
#        iucn = iucn,
#        metric = metric,
#        n_sims = n_sims,
#        parallel_run = FALSE,
#        append_sf = TRUE,
#        # cluster = cl,
#        regions = regions
#      )
#      
#      
#    })
#    
#    return(sexed_tmp_results)
#    
#  }
#  
#  
#  parallel::clusterExport(cl, c("sexes", "wrapper_fn", "sexed_pca_list", "region_shapes", "calc_region_ses", "regions", "metric", "null_rast", "iucn", "n_sims", "dispRity", "calc_iucn_ses"), envir = .GlobalEnv)
#  
#  parallel::clusterEvalQ(cl, {
#    library(Rcpp)
#    library(sf)
#    library(dispRity)
#    # If other packages with Rcpp dependencies are used, load them here as well
#  })
#  
#  results <- parallel::parLapply(cl, sexes, function(sex) {
#    try(wrapper_fn(sex))
#  })
#  results <- parallel::parLapply(cl = cl, sexes, wrapper_fn)
#  
#  parallel::stopCluster(cl)
 
 # End attempted parallelisation --------------------------------------------------
# ------------------------------------------------------------

# Non-parallelised version (working)

# get number of regions to apply across (this method is much faster than nrow())
n_regions <- length(attr(region_shapes, "row.names"))

# apply function across each sexed dataset
results <- lapply(sexes, function(sex){
  
  pca_sexed <- sexed_pca_list[[sex]]
  
  # apply function across each region (row) in sf object
  lapply(1:n_regions, function(region_number) {
    
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
      append_sf = ses_only,
      # cluster = cl,
      regions = regions,
      centroid = div_loss_type
    )
  })
  
})

# name elements of the results list
names(results) <- sexes


# for each sex, process results into dfs to bind to sf data
results_dfs_full <- lapply(sexes, function(sex) {
  
  # get individual sexed result
  sexed_results <- results[[sex]]
  
  # bind results together into a single dataframe
  sexed_results <- as.data.frame(do.call(rbind, sexed_results))
  
  # add column for sex
  sexed_results$sex <- rep(sex, nrow(sexed_results))
  
  # convert NaN values to NA (where there are e.g. no CR species to be lost)
  sexed_results <- data.frame(lapply(sexed_results, gsub, pattern = NaN, replacement = NA, fixed = TRUE))
  
})

# identify rows to be removed (regions which contain no species or only a single species in our dataset)
rows_to_keep <- which(!(results_dfs_full[[1]][, 1] == "nospec_region" | results_dfs_full[[1]][, 1] == "singlespec_region"))

# rename results df elements (so the plotting works)
names(results_dfs_full) <- sexes

# combine into single df
results_dfs_full <- do.call(rbind, results_dfs_full)

# create ses-only results (to be appended to region shapefiles)
results_dfs_ses <- lapply(sexes, function(sex){
  
  ses_results_sexed <- lapply(results[[sex]], function(element){
    ses_only <- element[, "ses"]
  })
  
  # bind results together into a single dataframe
  ses_results_sexed <- as.data.frame(do.call(rbind, ses_results_sexed))
  
  colnames(ses_results_sexed) <- paste("ses", c("all", "CR", "EN", "VU", "NT"), sep = "_")
  
  # adjust column names to specify sex
  colnames(ses_results_sexed) <- paste(sex, colnames(ses_results_sexed), sep = "_")
  
  # convert NaN values to NA (where there are e.g. no CR species to be lost)
  ses_results_sexed <- data.frame(lapply(ses_results_sexed, gsub, pattern = NaN, replacement = NA, fixed = TRUE))
  
})
# identify rows to be removed (regions which contain no species or only a single species in our dataset)
rows_to_keep <- which(!(results_dfs_ses[[1]][, 1] == "nospec_region" | results_dfs_ses[[1]][, 1] == "singlespec_region"))

# rename results df elements (so the plotting works)
names(results_dfs_ses) <- sexes


# bind ses results to sf data
region_results <- cbind(region_shapes, do.call(cbind, results_dfs_ses))
# remove regions for which there are no species
region_results <- region_results[rows_to_keep, ]



# save intermediate results

# as RDS for plot data
results_filename <- paste0(paste(clade, space, "regionalSESplotdata", div_loss_type, metric, avg_par, regions, paste0(n_sims, "sims"), paste0(iucn_type, "iucn"), sex_match, pam_res, pam_type, pam_seas, sep = "_"), ".rds")
saveRDS(
  results_dfs_ses, 
  here::here(
    "2_Patches", "3_OutputData", clade, "6_Spatial_mapping", "2_Regional_diversity_SES", pam_res, space,
    results_filename
  )
)

# save as csv for full data
results_filename <- paste0(paste(clade, space, "regionalSESfulldata", div_loss_type, metric, avg_par, regions, paste0(n_sims, "sims"), paste0(iucn_type, "iucn"), sex_match, pam_res, pam_type, pam_seas, sep = "_"), ".csv")
write.csv(
  results_dfs_full,
  here::here(
    "2_Patches", "3_OutputData", clade, "6_Spatial_mapping", "2_Regional_diversity_SES", pam_res, space,
    results_filename
  )
)


# plot results ----

# load results (for plotting)
results_filename <- paste0(paste(clade, space, "regionalSESplotdata", div_loss_type, metric, avg_par, regions, paste0(n_sims, "sims"), paste0(iucn_type, "iucn"), sex_match, pam_res, pam_type, pam_seas, sep = "_"), ".rds")
results_dfs <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", clade, "6_Spatial_mapping", "2_Regional_diversity_SES", pam_res, space,
    results_filename
  )
)

# remove names (for processing)
names(results_dfs) <- NULL

# identify rows to be removed (regions which contain no species or only a single species in our dataset)
rows_to_keep <- which(!(results_dfs[[1]][, 1] == "nospec_region" | results_dfs[[1]][, 1] == "singlespec_region"))

# load region data (either Dinerstein et al 2017 ecoregion data or biome data derived from same)
region_shapes <- load_regions(regions)

# bind results to sf data
region_results <- cbind(region_shapes, do.call(cbind, results_dfs))

# save region results as shapefiles
# i.e. ecoregion data with diversity values attached
results_sh_filename <- paste0(paste(clade, space, "regionalSESplotdata", div_loss_type, metric, avg_par, regions, paste0(n_sims, "sims"), paste0(iucn_type, "iucn"), sex_match, pam_res, pam_type, pam_seas, sep = "_"), ".shp")
output_folder <- here::here("2_Patches", "3_OutputData", clade, "6_Spatial_mapping", "2_Regional_diversity_SES", pam_res, space, "shapefiles")
if(!dir.exists(output_folder)){
  dir.create(output_folder, recursive = TRUE)
}
sf::st_write(region_shapes, paste(output_folder, results_sh_filename, sep = "/"),append = FALSE)

# save as dfs with biome info but no geometry
region_results_dfs_nogeom <- sf::st_drop_geometry(region_results)
results_nogeom_filename <-  paste0(paste(clade, space, "regionalSESplotdata", div_loss_type, metric, avg_par, regions, paste0(n_sims, "sims"), paste0(iucn_type, "iucn"), sex_match, pam_res, pam_type, pam_seas, sep = "_"), ".csv")
write.csv(
  region_results_dfs_nogeom,
  here::here(
    "2_Patches", "3_OutputData", clade, "6_Spatial_mapping", "2_Regional_diversity_SES", pam_res, space,
    results_nogeom_filename
  )
)

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
# load PAM
pam_filename <- paste0("PAM_birds_Behrman", pam_res, "_Pres12_Orig12_Seas12_", pam_type, pam_seas, ".rds")
pam_raw <- readRDS(
  paste(pams_filepath, pam_filename, sep = "/")
)
# extract null raster from PAM file
null_rast <- extract_null_rast(pam_raw, rast_type = "raster")
# remove raw PAM (for RAM)
rm(pam_raw)
# produce list of rasterlayers of each results - each SES size (e.g. -CR, -EN etc.) has a raster, for
# each sex
results_raster <- lapply(column_names, fasterize::fasterize, sf = region_results, raster = null_rast)
# convert elements to spatrasters
results_raster <- lapply(results_raster, terra::rast)
# combine list elements into a single multilayer raster, with a layer for each element
results_raster <- terra::rast(results_raster)
names(results_raster) <- column_names
# remove the SES of all species - i.e., pre-removing IUCN catefories - which are always 0 by definition
results_raster <- subset(results_raster, paste0(sexes, "_ses_all"), negate = TRUE)
# collect garbage
gc()


# produce thresholded raster where all values < -2 are set to -2
thresh_rast <- clamp(results_raster, lower = -2, values = TRUE)
# extract values
new_vals <- terra::values(thresh_rast)
# set all values which aren't -2 to NA (so we can )
new_vals[new_vals != -2] <- NA
# now set as raster values
terra::values(thresh_rast) <- new_vals




## Plot ---
## all layers in one plot, with regions with SES < -2 outlined (if requested)
## Set Outlining parameter
outline <- TRUE

# First define range for legend and colour scale (make it very slightly higher/lower
# than true range to stop the max/min regions being whited out)
leg_range <- c(min(terra::values(results_raster), na.rm = TRUE), 
               max(terra::values(results_raster), na.rm = TRUE))
leg_increm <- (max(leg_range) - min(leg_range)) / 100
leg_range <- c(min(leg_range) - leg_increm,
               max(leg_range) + leg_increm)
# define colour scale
if(palette_choice == "viridis"){
  leg_fill <- viridis::viridis(100, direction = -1)
} else if(palette_choice == "custom"){
  # custom option 1: truncated viridis (so it's a one-colour gradient)
 # leg_fill <- viridis::plasma(125, direction = -1)[26:125]
 # custom option 2: red-white gradient
 leg_fill <- rgb(colorRamp(c("red", "white"), space = "Lab", interpolate = "linear")(seq(from = 0, to = 1, by = 0.01)), maxColorValue = 255)
}

# define legend title
leg_title <- paste0("SES value (", metric, ")")

# get raster layers
layers <- names(results_raster)

# generate alphabetical labels for plots
labels <- LETTERS[1:length(layers)]
names(labels) <- layers

# set number of columns and rows
n_cols <- length(sexes)
n_rows = (length(layers) %/% n_cols) + 1
# Define layout matrix of plots within figure
layout_matrix <- matrix(
  c(1, 5, 
    2, 6, 
    3, 7, 
    4, 8, 
    9, 9), # legend across two rows at the bottom
  nrow = n_rows, ncol = n_cols, byrow = TRUE
)

# set png filename
png_filename <- paste0(paste(clade, space, "regionalSES", div_loss_type, metric, avg_par, regions, paste0(n_sims, "sims"), paste0(iucn_type, "iucn"), sex_match, pam_res, pam_type, pam_seas, sep = "_"), ".png")
# initialise png saving
png(
  here::here(
    "2_Patches", "4_OutputPlots", clade, "3_Spatial_mapping", "2_Regional_diversity_SES", pam_res, space,
    png_filename
  ), 
  width = 1450, height = 1550, res = 72, pointsize = 20
  )

# set layout of plots
layout(
  heights = c(1, 1, 1, 1, 0.5),
  mat = layout_matrix
  )


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
  
  # plot polygon of values < 2 (if requested)
  if(outline == TRUE){
    
    if(!all(is.na(terra::values(thresh_rast[[layer]])))){
      
      terra::polys(terra::as.polygons(thresh_rast[[layer]]), 
                col = "red", density = 10, angle = 45,
                border = "red", lwd = 1)
      }
    
    }
  
}
# add combined legend
par(mar = c(8, 12, 1, 12))  # set margins for the legend
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

# end png device
dev.off()





