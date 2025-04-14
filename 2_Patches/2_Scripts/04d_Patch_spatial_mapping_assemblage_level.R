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
# select metric ("centr-dist", "nn-k", "nn-count", "sum.variances", "sum.ranges", "convhull.volume")
metric <- "centr-dist"
# select type of averaging to use ("mean" or "median")
avg_par <- "mean"
# select whether to exclude grid cells with species richness below a certain threshold (e.g. 5)
# set as 0 if no threshold wanted
sr_threshold <- 10
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
div_values <- apply(pam, 1, calc_metric_gridcell, pca_data = sexed_metric_list$M, metric = metric, avg_par = avg_par, min_species = sr_threshold)

# make a raster of assemblage-level diversity values
div_raster <- make_assdiv_raster(div_values, null_rast)
