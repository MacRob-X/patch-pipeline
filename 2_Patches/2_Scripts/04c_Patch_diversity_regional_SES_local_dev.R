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

# subset PCA data to only species for which we 
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



## Functionalise this later

# set dispRity metric
if(metric == "centr-dist"){
  metric_get <- "centroids"
} else if (metric == "nn-k"){
  metric_get <- "mean.nn.dist"
}else if (metric == "nn-count"){
  metric_get <- "count.neighbours"
}

# set up place to store simulations
a <- n_sims # change depending on number of simulations  
sims <- matrix(NA, nrow=1, ncol = a) # place to store simulations for SES

# set up results matrix
results <- matrix(NA, nrow=5, ncol = 6)
colnames(results) <- c("iucn", "species_richness", "raw", "null_mean", "null_sd", "null_se")

# get male pca data
pca_M <- sexed_pca_list$M

# set region of interest
roi <- region_shapes[5, ]

# create null rast, subsetted to extent of ecoregion
# first assign values to null_rast to keep track of grid cells
null_rast_sub <- null_rast
# get number of grid cells in full raster
ncells <- length(terra::values(null_rast))
terra::values(null_rast_sub) <- 1:ncells
# crop raster to ecoregion extent
null_rast_sub <- terra::crop(null_rast_sub, roi)
# mask raster to ecoregion extent
null_rast_sub <- terra::mask(null_rast_sub, roi)


# use values of subsetted raster to subset PAM to ROI
pam_sub <- pam[terra::values(null_rast_sub), ]

# get species which are present in subsetted PAM
if(length(terra::values(null_rast_sub)) < 2){
  # if ROI consists of 1 or fewer grid cells
  species_pres <- names(which(!is.na(pam_sub)))
} else {
  # if ROI consists of more than 1 grid cell
  species_pres <- names(which(apply(pam_sub, 2, function(x) any(!(is.na(x))))))
}

# get collections of species within region with threat levels sequentially trimmed
all <- species_pres
no_CR <- iucn[iucn$species_birdtree %in% species_pres & iucn$iucn_cat != "CR", "species_birdtree"]
no_EN <- iucn[iucn$species_birdtree %in% species_pres & iucn$iucn_cat != "CR" & iucn$iucn_cat != "EN", "species_birdtree"]
no_VU <- iucn[iucn$species_birdtree %in% species_pres & iucn$iucn_cat != "CR" & iucn$iucn_cat != "EN" & iucn$iucn_cat != "VU", "species_birdtree"]
no_NT <- iucn[iucn$species_birdtree %in% species_pres & iucn$iucn_cat != "CR" & iucn$iucn_cat != "EN" & iucn$iucn_cat != "VU" & iucn$iucn_cat != "NT", "species_birdtree"]



# clip PCA data to only species present in regions (all IUCN levels)
pca_region <- pca_M[all, ]

# calculate average of metric of all species in ROI
roi_metric_avg <- avg(dispRity(pca_region, metric = get(metric_get))$disparity[[1]][[1]], avg_type = avg_par, na.rm = TRUE)

# add species richness and average metric to results table
results[1, "species_richness"] <- length(all)
results[1, "raw"] <- roi_metric_avg

# Now do the same but for CR species only, comparing to n simulations of random species loss to calculate
# SES

# clip PCA data to only species present in regions (no CR species)
pca_subset <- pca_region[no_CR, ]

# calculate mean distance to centroid of all species in ROI in IUCN cats of interest
roi_metric_avg <- avg(dispRity(pca_subset, metric = get(metric_get))$disparity[[1]][[1]], avg_type = avg_par, na.rm = TRUE)

# add species richness and average metric to results table
results[2, "species_richness"] <- length(no_CR)
results[2, "raw"] <- roi_metric_avg

# for parallelisation, set up a cluster using the number of cores (2 less than total number of laptop cores)
no_cores <- parallel::detectCores() - 2
cl <- parallel::makeCluster(no_cores)
# export necessary objects to the cluster
parallel::clusterExport(cl, c("all", "no_CR", "no_EN", "no_VU", "no_NT", "pca_region", "metric_get", "dispRity", "avg", "avg_par", metric_get))

# Calculate SES compared to simulations of random loss of species
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      coms <- sample(x = all, size = length(no_CR), replace = FALSE)
      trait_vec <- matrix(pca_region[coms, ], ncol = dim(pca_region)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <- avg(dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]], avg_type = avg_par, na.rm = TRUE)
      
      return(disparity_value)
    }
  )
)

# add to results table
results[2, "null_mean"] <- mean(sims)
results[2, "null_sd"] <- sd(sims)
results[2, "null_se"] <- sd(sims)/sqrt(length(sims))


# check if there are as many or more species in ROI than the threshold SR
if(length(species_pres) >= sr_threshold){
  
  # if there are enough species, calculate mean distance to centroid of all species in space
  
  
  
} else {
  
  return(NA)
  
}

