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
library(ggplot2)

# Load shared functions
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "mapping.R"
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
# select whether to calculate local or global diversity loss (i.e., mean distance to local or global centroid)
div_loss_type <- "local"
# select number of null distributions to generate
n_sims <- 10
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
pam_res <- "200km"
# enter PAM files location
pams_filepath <- "X:/cooney_lab/Shared/Rob-MacDonald/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/PAMs"
#pams_filepath <- "X:/cooney_lab/Shared/Rob-MacDonald/SpatialData/BirdLife/BirdLife_Shapefiles_v9/PAMs/100km/Behrmann_cea/"
# # use ecoregions or biomes?
regions <- "ecoregions"

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

# for parallelisation, set up a cluster using the number of cores (4 less than total number of laptop cores)
# Note that it might be more sensible to initialise the cluster within the calc_region_ses function itself
# - this will increase the overhead time, but maybe will be less heavy in terms of RAM usage
# - the above is incorrect, it's just as RAM-hungry. This might be a problem with the full dataset

# set up parallel cluster (if using)
# no_cores <- parallel::detectCores() - 4
# cl <- parallel::makeCluster(no_cores)
results <- lapply(1:nrow(region_shapes), function(region_number) {
  
  region_sf <- region_shapes[region_number, ]
  
  calc_region_ses(
    region_sf = region_sf,
    pca_data = sexed_pca_list$M,
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
# stop the cluster
# stopCluster(cl)


# bind results together into a single matrix
results_all <- do.call(rbind, results)

# convert NaN values to NA (where there are e.g. no CR species to be lost)
results_all <- sub(NaN, NA, x = results_all)

# bind to sf data
region_results <- cbind(region_shapes, results_all)

# remove regions for which there are no species or only a single species
region_results <- region_results[region_results$ses_all != "nospec_region" & region_results$ses_all != "singlespec_region", ]

# plot results
# wrapper function to feed into apply
plot_layer <- function(lyr){
  
  p <- ggplot() + 
    geom_sf(data = region_results, aes(fill = as.numeric(get(lyr)), colour = as.numeric(get(lyr)))) + 
    scale_fill_viridis_c() + 
    scale_colour_viridis_c()
  p + labs(fill = lyr)
  return(p)
}
layers <- colnames(results_all)[2:5]

plots <- lapply(layers, plot_layer)

png("ecoregion_ses_M.png", width = 2000, height = 1000)
gridExtra::grid.arrange(grobs = plots, nrow = 2)
dev.off()


