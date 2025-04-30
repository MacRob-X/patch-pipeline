# Testing assemblage-level convex hull volume calculation


# make wrapper function for averaging (allows selection of average type e.g. mean, median)
avg <- function(vals, avg_type, na.rm = TRUE){
  return(get(avg_type)(vals, na.rm = na.rm))
}

# Calculate assemblage-level disparity metric
# calculate value of chosen metric for each occupied grid cell in PAM (i.e., each non-NA row)
calc_metric_gridcell <- function(pam_row, pca_data, metric, avg_par, min_species = 0, ...){
  
  # set dispRity metric
  if(metric == "centr-dist"){
    metric_get <- get("centroids")
  } else if (metric == "nn-k"){
    metric_get <- get("mean.nn.dist")
  }else if (metric == "nn-count"){
    metric_get <- get("count.neighbours")
  } else if(metric == "convhull.volume"){
    metric_get <- get("convhull.volume")
  } else if(metric == "displacements"){
    metric_get <- get("displacements")
  } else if(metric == "sum.variances"){
    metric_get <- list(sum, variances)
  } else if(metric == "sum.ranges"){
    metric_get <- list(sum, ranges)
  } 
  
  # get species present in grid cell
  species_pres <- names(pam_row)[which(!is.na(pam_row))]
  
  # if no species are present, one species is present, or if there are fewer species present than 
  # cutoff value, return NA
  if(length(species_pres) <= 1 | length(species_pres) < min_species){
    return(NA)
  }
  
  # subset PCA data to these species only
  pca_data <- pca_data[which(rownames(pca_data) %in% species_pres), ]
  
  # calculate value of metric for these species only
  div_values <- dispRity(
    pca_data,
    metric = metric_get,
    ...
  )$disparity[[1]][[1]]
  
  # calculate average if output is dimension-level 2
  metric_type <- dispRity::make.metric(metric_get, silent = TRUE)[["type"]]
  if(metric_type == "level2"){
    
    # calculate average (using wrapper function)
    avg_div_val <- avg(div_values, avg_par, na.rm = TRUE)
    
    # return average
    return(avg_div_val)
    
  } else {
    
    # else return single level1 value
    return(div_values)
    
  }
  
  
  
  
}

str(pam)
str(sexed_metric_list$M)
str(metric)
str(avg_par)
str(sr_threshold)


div_values <- apply(pam, 1, calc_metric_gridcell, pca_data = sexed_metric_list$M, metric = metric, avg_par = avg_par, min_species = sr_threshold)


# try out hypervolume package
install.packages("hypervolume")
hypervolume::hypervolume(pca_dat)


# 22/04/2025

# save csv of male/female PCA spaces to test 

# Probe PCA space to find best disparity metric
# Largely based on code from T Guillerme's shinyapp
# https://tguillerme.shinyapps.io/moms/


# Clear environment
rm(list=ls())


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
# Select number of PC axes to retain ("all" or a number)
axes <- "all"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# select sex of interest ("all", "male_female", "male_only", "female_only", "unknown_only")
sex_interest <- "male_female"



# Load data ----

# Load PCA data
pca_filename <- paste(clade, sex_match, "patches.231030.PCAcolspaces", "rds", sep = ".")
pca_dat <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]][["x"]]


# Data preparation ----

# restrict PCA to requested number of axes
if(axes != "all"){
  pca_dat <- pca_dat[, 1:axes]
}


# generate individual list of sex-specific species values of metric of interest (keep or discard sexes
# depending on value of sex_match)
# first set sex list to iterate over
sexes <- set_sex_list(sex_interest)
# now apply a function to extract named vector lists of diversity values for each sex individually
sexed_metric_list <- lapply(sexes, extract_sex_vals, div_data = pca_dat, assemblage_level = TRUE)

# name the elements of the list
names(sexed_metric_list) <- sexes

# create function to sample a single PCA space
sample_space <- function(space, sample_size, seed){
  
  # set seed for reproducibility
  set.seed(seed)
  
  # get rownames of random species to retain
  to_retain <- sample(rownames(space), sample_size)
  
  # subset space by rownames and return
  return(space[to_retain, , drop = FALSE])
  
}

# Randomly remove species to limit to 1000 observations
spaces_to_save <- lapply(sexed_metric_list, sample_space, sample_size = 1000, seed = 10)

# save as CSV to upload to moms shinyapp
for(sex in sexes){
  write.csv(spaces_to_save[[sex]], file = paste0(sex, ".", space, ".randomsubset.csv"))
}
lapply(spaces_to_save, function(space) write.csv(space, file = ""))


# 30/04/2025 ----
# View distributions of centroid distances as you remove IUCN categories


# clear environment
rm(list=ls())

# Load libraries
library(dispRity)
library(dplyr)
library(ggplot2)
library(extrafont)


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
# select metric ("centr-dist", "nn-k", "nn-count", "nn-all", "sum.variances", "sum.ranges")
metric <- "centr-dist"
# select whether to calculate local or global diversity loss (i.e., mean distance to local or global centroid)
div_loss_type <- "global"


# select whther to use liberal, conservative, or nominate IUCN data
# "liberal" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the highest threat level of the multiple species
# "conservative" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the lowest threat level of the multiple species
# "nominate" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the threat level of the BL species that corresponds to the nominate subspecies
iucn_type <- "nominate"



## Plotting parameters
# Select whether to use binned (based on quantiles) or continuous colour scale
col_scale_type <- "binned"
# If binned, choose the number of quantiles to use
nquants <- 10
## Use viridis_c, viridis turbo palette or custom colour palette?
palette_choice <- "viridis"
## Plot on log scale?
log_choice <- FALSE

# Load data ----

# set PCA filename
pca_filename <- paste(clade, sex_match, "patches.231030.PCAcolspaces", "rds", sep = ".")
# Load PCA (values only)
pca_mat <- readRDS(
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


# Prepare data ----

# generate individual list of sex-specific species values of metric of interest (keep or discard sexes
# depending on value of sex_match)
# first set sex list to iterate over
sexes <- set_sex_list(sex_interest)

# now apply a function to extract named vector lists of values for each sex individually
sexed_pca_list <- lapply(sexes, extract_sex_vals, div_data = pca_mat, assemblage_level = TRUE)
names(sexed_pca_list) <- sexes

# define list of IUCN categories
iucn_cats <- c("All", "CR", "EN", "VU", "NT")

# define function to create lists of trimmed species
metric_dist_iucn <- function(pca_data, iucn_cat, iucn, metric, sex){
  
  # get species in data
  all_species <- rownames(pca_data)
  
  # get list of species within region with threat level of interest trimmed (inclusive - i.e. trims 
  # threat level of interest and any more threatened categories)
  if(iucn_cat == "All"){
    trimmed_species <- all_species
  } else if(iucn_cat == "CR"){
    trimmed_species <- iucn[iucn$species_birdtree %in% all_species & iucn$iucn_cat != "CR", "species_birdtree"]
  } else if(iucn_cat == "EN"){
    trimmed_species <- iucn[iucn$species_birdtree %in% all_species & iucn$iucn_cat != "CR" & iucn$iucn_cat != "EN", "species_birdtree"]
  } else if(iucn_cat == "VU"){
    trimmed_species <- iucn[iucn$species_birdtree %in% all_species & iucn$iucn_cat != "CR" & iucn$iucn_cat != "EN" & iucn$iucn_cat != "VU", "species_birdtree"]
  } else if(iucn_cat == "NT"){
    trimmed_species <- iucn[iucn$species_birdtree %in% all_species & iucn$iucn_cat != "CR" & iucn$iucn_cat != "EN" & iucn$iucn_cat != "VU" & iucn$iucn_cat != "NT", "species_birdtree"]
  }
  
  # check if more than one species in trimmed species list
  # if not, set all diversity metrics to 0 and end function
  if(length(trimmed_species) <= 1){
    
    return("No species in category")

    
  }
  
  # clip PCA data to only species present in regions (no CR species)
  pca_subset <- pca_data[trimmed_species, ]
  
  # set metric
  metric_get <- set_metric(metric)
  
  # calculate metric of all species in IUCN cats of interest
  metric <- dispRity(pca_subset, metric = metric_get)$disparity[[1]][[1]]
  
  # convert to df with columns for species, sex, and IUCN level
  metric_df <- data.frame(species = rownames(pca_subset), iucn_loss = iucn_cat, metric = metric, sex = sex)
  
  # return df
  return(metric_df)
}

metric_dists <- lapply(
  sexes,
  function(sex) {
    # get sexed data
    sexed_pca <- sexed_pca_list[[sex]]
    # get metric distributions
    metric_dists_sex <- lapply(iucn_cats, metric_dist_iucn, pca_data = sexed_pca, iucn = iucn, metric = metric, sex = sex)
    # bind into one df
    metric_dists_sex <- do.call(rbind, metric_dists_sex)
    return(metric_dists_sex)
  }
)

# bind into a single df
metric_dists <- do.call(rbind, metric_dists)

# convert IUCN category to factor and set levels
metric_dists$iucn_loss <- factor(metric_dists$iucn_loss, levels = iucn_cats)

# plot overlain density plots
metric_dists %>% 
  ggplot(aes(x = metric, fill = iucn_loss)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~iucn_loss + sex, nrow = 5)
