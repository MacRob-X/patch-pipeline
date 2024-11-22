## Calculate SES for diversity change in each bioregion separately
## Robert MacDonald
## 15th November 2024

# clear environment
rm(list=ls())

# Load libraries
library(raster)
library(parallel)

# Load shared functions
source(
  here::here(
    "2_Patches", "2_Scripts", "2_BetaVersions", "04_shared_mapping_functions.R"
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
div_loss_type <- "global"
# select number of null distributions to generate
n_sims <- 1000
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
pam_res <- "200km"
# enter PAM files location
pams_filepath <- "X:/cooney_lab/Shared/Rob-MacDonald/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/PAMs"
#pams_filepath <- "X:/cooney_lab/Shared/Rob-MacDonald/SpatialData/BirdLife/BirdLife_Shapefiles_v9/PAMs/100km/Behrmann_cea/"
# # use ecoregions or biomes?
regions <- "biomes"

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


# Functions ----

# attach IUCN data to diversity data
attach_iucn <- function(data, iucn_data){
  
  # get species and sex columns in data
  data <- cbind(data, get_spp_sex(data))
  
  # check if IUCN data has a 'species' column
  if("species" %in% colnames(iucn)){
    
    # join IUCN to diversity data
    data <- dplyr::left_join(data, iucn_data, by = "species")
    
  } else if("species_birdtree" %in% colnames(iucn)){
    
    # join IUCN data to diversity data
    data <- dplyr::left_join(data, iucn_data, by = dplyr::join_by(species == species_birdtree))
    
  }
  
  return(data)
  
}

# calculate SES for loss of individual level (cumulative) of IUCN
calc_iucn_ses <- function(data, iucn_cat, n_sims = 1000, level = c("local", "global")){
  
  # initialise output dataframe
  results <- matrix(NA, nrow = 1, ncol = 6)
  
  # if working on global level, just use the centroid distances already calculated (i.e., use the global
  # centroid to calculate all SES)
  if(level == "global"){
    
    
    
  }
  
}

# calculate SES of metric given for a given area with sequential removal of threatened species
#' Title
#'
#' @param pam species presence-absence matrix
#' @param pca_data PCA values for each specimen (only needed if level = "local")
#' @param iucn_data IUCN category data associated with each species
#' @param ecoreg ecoregion data (as a single row of an sf object)
#' @param sexed_metric_vals values of diversity metric, as a list of named vectors with one vector for each sex
#' @param null_terrast nulla terra raster
#' @param avg type of average to use (mean or median)
#' @param level level on which to calculate metric (global - uses predefined metric values - or local - calculates metric for subset of species in each bioregion)
#'
#' @return
#' @export
#'
#' @examples
calc_region_ses <- function(
    pam, pca_data = NULL, 
    iucn_data, 
    ecoreg, 
    sexed_metric_vals = NULL, 
    null_terrast, 
    avg = c("mean", "median"), 
    level = c("local", "global"),
    n_sims = 1000
    ){
  
  # set number of grid cells in PAM
  ncells <- nrow(pam)
  
  # create null raster, subsetted to extent of ecoregion
  # first assign values to null_rast to keep track of grid cells
  null_terrast_sub <- null_terrast
  terra::values(null_terrast_sub) <- 1:ncells
  null_terrast_sub <- terra::mask(terra::crop(null_terrast_sub, ecoreg), ecoreg)
  
  # use values of subsetted raster to subset PAM to ROI
  pam_sub <- pam[terra::values(null_terrast_sub), ]
  
  # get species which are present in subsetted PAM
  # use different method if ecoregion occupies 1 or fewer grid cells
  if(length(terra::values(null_terrast_sub)) < 2){
    
    species_pres <- names(which(!is.na(pam_sub)))
    
  } else {
    
    species_pres <- names(which(apply(pam_sub, 2, function(x) any(!(is.na(x))))))
    
  }
  
  # check if there are as many or more species in ROI than the threshold SR
  if(length(species_pres) >= sr_threshold){
    
    # if there are enough species, calculate SES for males and females 
    
    # set IUCN levels
    iucn_levels <- c("CR", "EN", "VU", "NT", "LC")
    
    # calculate SES for cumulative loss of each IUCN category
    res_region <- calc_iucn_ses()
    
    
    
    avg_metric_M <- get(avg)(metric_M[species_pres], na.rm = TRUE)
    avg_metric_F <- get(avg)(metric_F[species_pres], na.rm = TRUE) 
    
  } else {
    
    # if there are fewer, set values = NA
    avg_metric_M <- NA
    avg_metric_F <- NA
  }
  
  
}

# Load data ----

# load diversity metric data 
div_data_filename <- here::here(
  "2_Patches", "3_OutputData", "4_Diversity_measures", space,
  paste(clade, "patches", space, sex_match, "diversitymetrics.csv", sep = "_")
)
div_data <- read.csv(div_data_filename, row.names = 1, check.names = FALSE)

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

# get list of species included in diversity data
spp_list <- get_unique_spp(div_data)

# subset pam to only species in diversity data
pam <- subset_pam(pam, spp_list)

# subset diversity data to only species present in PAM
div_data <- subset_data_by_spp(div_data, colnames(pam))

# attach IUCN categories to diversity data
div_data <- attach_iucn(div_data, iucn)

# clear RAM
gc()

# convert ecoregions to Behrman CEA CRS
ecoregions <- sf::st_transform(ecoregions, crs = terra::crs(null_rast, proj = TRUE))

# generate individual list of sex-specific species values of metric of interest (keep or discard sexes
# depending on value of sex_match)
# first set sex list to iterate over
sexes <- set_sex_list(sex_interest)
# now apply a function to extract named vector lists of diversity values for each sex individually
sexed_metric_list <- lapply(sexes, extract_sex_vals, div_data = div_data, metric = metric)

# calculate SES for each ecoregion

# test
ses_biome_1 <- calc_region_ses(pam, iucn_data = iucn, ecoreg = region_shapes[1, ], sexed_metric_vals = sexed_metric_list, null_terrast = null_rast, avg = avg_par, level = div_loss_type)
