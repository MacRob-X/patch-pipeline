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



# 01/05/2025 ----
# SES of aesthetic attractiveness loss
# using iratebirds data

# Calculate Standard Effect Size (SES) of change in diversity metric as species at 
# different IUCN threat levels are sequentially lost
# Robert MacDonald
# Edited 29/08/2024

# Note: this code is very heavily based on code from Hughes et al. (2022) Curr Biol

# attach libraries
library(dplyr)
library(dispRity)
library(ggplot2)
library(parallel)
library(extrafont)

# clear environment
rm(list=ls())

# custom dispRity metrics
source(
  here::here(
    "3_SharedScripts", "dispRity_metric_functions.R"
  )
)

## EDITABLE CODE ## ----

# select number of null distributions to generate
n_sims <- 1000
# select whther to use liberal, conservative, or nominate IUCN data
# "liberal" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the highest threat level of the multiple species
# "conservative" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the lowest threat level of the multiple species
# "nominate" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the threat level of the BL species that corresponds to the nominate subspecies
iucn_type <- "nominate"
## END EDITABLE CODE ##

# Load data ----

# load iratebirds data (PCA of whichever colourspace - generated in 02_Patch_Analyse_features.R)
# Filter to Neoaves
irate_master <- read.csv(
  here::here(
    "4_SharedInputData", "iratebirds_data", "iratebirds_final_predictions_average_fullmodel_subsetmodel_151122.csv"
  ),sep = ";" 
) %>% 
  filter(
    order != "Struthioniformes", 
    order != "Rheiformes", 
    order != "Casuariiformes", 
    order != "Apterygiformes", 
    order != "Tinamiformes", 
    order != "Galliformes", 
    order != "Anseriformes"
  )

# load IUCN Red List data
iucn_filename <- paste0("IUCN_RedList_data_130324", ".rds")
iucn <- readRDS(
  here::here(
    "4_SharedInputData", iucn_filename
  )
) %>% 
  select(scientific_name, category) %>% 
  mutate(scientific_name = gsub(" ", "_", scientific_name)) %>% 
  rename(species_birdlife = scientific_name, iucn_cat = category)


# join iratebirds data to IUCN data
# NOTE that these are not properly matched, so we're going to lose lots of species and possibly
# bias the results - this is messing around only


####################    
### 2: GLOBAL ANALYSIS
####################   
## a) Morphological diversity
## NB: To generate Table S1, the code below can be ran on individual PCs.
## E.g. trait <- trait[,1]



# Matrix containing attractiveness data (remove duplicates):
trait <- irate_master %>% 
  as.data.frame() %>% 
  rename(attract = predicted_attractiveness_full_model) %>% 
  mutate(
    species = sci_name,
    attract = as.numeric(gsub(",", ".", attract))
  ) %>% 
  select(species, attract) %>% 
  mutate(species = gsub(" ", "_", species)) %>% 
  # Filter out subspecies
  filter(stringr::str_count(species, "_") == 1) %>% 
  distinct()

# Species matched to IUCN categories:
species <- iucn %>% 
  rename(
    species = species_birdlife
  ) %>% 
  select(
    species, iucn_cat
  ) %>% 
  filter(
    species %in% trait[, "species"]
  )

# remove species not in IUCN data
trait <- trait[trait[, "species"] %in% species[, "species"], ]

# check species are identical in both matrices
spec_iucn <- sort(species[, "species"])
spec_trait <- sort(trait[, "species"])
identical(spec_iucn, spec_trait)
rm(spec_iucn, spec_trait)

# add rownames
rownames(species) <- species$species
rownames(trait) <- trait$species
species <- species[rownames(species) %in% rownames(trait), ]
# reorder to match trait data
species <- species[match(rownames(trait), rownames(species)), ]

# remove species column from trait data
trait <- trait %>% 
  select(attract)


### OUTPUT DATA FRAME
res <- matrix(NA, nrow=5, ncol = 6)

a <- n_sims #change depending on number of simulations  
sims <- matrix(NA, nrow=1, ncol = a) #place to store simulations for SES


all <- rownames(species)
# extract trait for species in the random community
trait_vec <- matrix(trait[all,], ncol=dim(trait)[2])
# calculate SR & median aesthetic attractiveness
res[1,2] <- length(all)
res[1,3] <- median(trait_vec)

# get collections of species with threat levels sequentially trimmed
noCR <- rownames(species)[species$iucn_cat != "CR"]
noEN <- rownames(species)[species$iucn_cat != "CR" & species$iucn_cat != "EN"]
noVU <- rownames(species)[species$iucn_cat != "CR" & species$iucn_cat != "EN" & species$iucn_cat != "VU"]
noNT <- rownames(species)[species$iucn_cat != "CR" & species$iucn_cat != "EN" & species$iucn_cat != "VU" & species$iucn_cat != "NT"]

# for parallelisation, set up a cluster using the number of cores (2 less than total number of laptop cores)
no_cores <- parallel::detectCores() - 2
cl <- parallel::makeCluster(no_cores)
# export necessary objects to the cluster
parallel::clusterExport(cl, c("all", "noCR", "noEN", "noVU", "noNT", "trait"))

# extract trait for species in the random community
trait_vec <- matrix(trait[noCR,], ncol=dim(trait)[2])
# calculate SR & Mean distance to centroid
res[2,2] <- length(noCR)
res[2,3] <- median(trait_vec)

# Use parLapply to replace the for loop with a parallelised apply function for speed
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      coms <- sample(x = all, size = length(noCR), replace = FALSE)
      trait_vec <- matrix(trait[coms, ], ncol = dim(trait)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <-median(trait_vec)
      
      return(disparity_value)
    }
  )
)


# add results to table
res[2,4] <- mean(sims)
res[2,5] <- sd(sims)
res[2,6] <- sd(sims)/sqrt(length(sims))


# extract trait for species in the random community
trait_vec <- matrix(trait[noEN,], ncol=dim(trait)[2])
# calculate SR & Mean distance to centroid
res[3,2] <- length(noEN)
res[3,3] <- median(trait_vec)

# Use parLapply to replace the for loop with a parallelised apply function for speed
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      coms <- sample(x = all, size = length(noEN), replace = FALSE)
      trait_vec <- matrix(trait[coms, ], ncol = dim(trait)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <- median(trait_vec)
      
      return(disparity_value)
    }
  )
)

# add results to table
res[3,4] <- mean(sims)
res[3,5] <- sd(sims)
res[3,6] <- sd(sims)/sqrt(length(sims))

# extract trait for species in the random community
trait_vec <- matrix(trait[noVU,], ncol=dim(trait)[2])
# calculate SR & Mean distance to centroid
res[4,2] <- length(noVU)
res[4,3] <- median(trait_vec)

# Use parLapply to replace the for loop with a parallelised apply function for speed
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      coms <- sample(x = all, size = length(noVU), replace = FALSE)
      trait_vec <- matrix(trait[coms, ], ncol = dim(trait)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <- median(trait_vec)
      
      return(disparity_value)
    }
  )
)


res[4,4] <- mean(sims)
res[4,5] <- sd(sims)
res[4,6] <- sd(sims)/sqrt(length(sims))

# extract trait for species in the random community
trait_vec <- matrix(trait[noNT,], ncol=dim(trait)[2])
# calculate SR & Mean distance to centroid
res[5,2] <- length(noNT)
res[5,3] <- median(trait_vec)

# Use parLapply to replace the for loop with a parallelised apply function for speed
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      coms <- sample(x = all, size = length(noNT), replace = FALSE)
      trait_vec <- matrix(trait[coms, ], ncol = dim(trait)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <- median(trait_vec)
      
      return(disparity_value)
    }
  )
)

# stop the cluster
stopCluster(cl)

res[5,4] <- mean(sims)
res[5,5] <- sd(sims)
res[5,6] <- sd(sims)/sqrt(length(sims))

res <- as.data.frame(res)

colnames(res) <- c("IUCN", "SR", "Raw", "NULL_MEAN", "NULL_SD", "NULL_SE")

res$IUCN[1] <- "All Species Retained"
res$IUCN[2] <- "CR Lost"
res$IUCN[3] <- "EN Lost"
res$IUCN[4] <- "VU Lost"
res$IUCN[5] <- "NT Lost (LC Only Retained)"

res$SR <- as.numeric(res$SR)
res$Raw <- as.numeric(res$Raw)
res$NULL_MEAN <- as.numeric(res$NULL_MEAN)
res$NULL_SD <- as.numeric(res$NULL_SD)
res$NULL_SE <- as.numeric(res$NULL_SE)

res$SES <- (res$Raw - res$NULL_MEAN)/res$NULL_SD
str(res)
res$SES[1] <- 0

res$Metric <- "Trait Diversity"


# Save results as csv
# set filename
filename <- paste0( "iratebirds_", "SES_", "nsims", n_sims, "_", "attractiveness", "_", iucn_type, "-iucn", ".csv")
write.csv(
  res, 
  here::here(
    "2_Patches", "3_OutputData", "4_Diversity_measures", filename
  ), 
  row.names = FALSE
)


# Plot results ----

# clear environment (except chosen variables)
rm(list=setdiff(ls(), c("clade", "space", "metric", "n_sims", "iucn_type")))

# load in results and get species richness as a proportion
filename <- paste0( "iratebirds_", "SES_", "nsims", n_sims, "_", "attractiveness", "_", iucn_type, "-iucn", ".csv")
res <- read.csv(
  here::here(
    "2_Patches", "3_OutputData", "4_Diversity_measures", filename
  )
) %>% 
  mutate(
    sr_prop = SR / max(SR)
  )

# convert IUCN cat lost to factor and adjust levels
res$IUCN <- factor(res$IUCN, levels = c("All Species Retained", "CR Lost", "EN Lost", "VU Lost", "NT Lost (LC Only Retained)"))


# plot (subset by sex if required)
# note that there's no point in adding error bars because the error's so small you can't see it
p <- res %>% 
  ggplot() + 
  geom_point(aes(x = sr_prop, y = SES, colour = IUCN), size = 2.5, shape = 17) + 
  geom_line(aes(x = sr_prop, y = SES), alpha = 0.6) +
  #  geom_errorbar(aes(ymin = SES - NULL_SE, ymax = SES + NULL_SE), width=0.005) +
  scale_x_reverse() + 
  #  scale_color_discrete(type = viridisLite::rocket(n = nlevels(res$IUCN))) + 
  scale_color_discrete(type = hcl.colors(n = nlevels(res$IUCN)+1, palette = "Terrain")) + 
  xlab("Remaining proportional species richness") + ylab("Standard Effect Size") + 
  labs(colour = "IUCN Status Lost") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = 2, linetype = "dashed") +
  theme_light() + 
  #  ggtitle(space) + 
  theme_bw() + 
  theme(text=element_text(size=12,  family="Century Gothic"))

p

# Save as png
png_filename <- paste0("iratebirds_attractiveness", "_SES_", "nsims", n_sims, "_", iucn_type, "-iucn", ".png")
png(
  here::here(
    "2_Patches", "4_OutputPlots", "2_Diversity_measures", "1_Diversity_extinction_risk", 
    png_filename
  ),
  width = 1500, height = 2000/3, res = 150
)
p
dev.off()




############ fixing bug in SES diversity code (completely separate from above code)

if(sex == "All"){
  spp_sample <- sample(x = unique_species, size = length(noCR)/2, replace = FALSE)
  coms <- c(paste0(spp_sample, "-M"), paste0(spp_sample, "-F"))
}
coms <- sample(x = all, size = length(noCR), replace = FALSE)

# 08/05/2025 ----

# Creating interactive 2D Neoaves UMAP plot (to get species locations)

# load space
load_umap <- TRUE
umap_filepath <- here::here(
  "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
  paste(clade, sex_match, "patches.nn.25.mindist.0.1", space, "pca1-12", "UMAP.rds", sep = ".")
)

plot_space <- readr::read_rds(umap_filepath) %>% 
  magrittr::extract2("layout") %>% 
  as.data.frame() %>% 
  rename(
    UMAP1 = V1, UMAP2 = V2
  ) %>% 
  mutate(
    UMAP1 = -UMAP1,
    UMAP2 = -UMAP2,
    spp_sex = rownames(.)
  )
x_axis <- "UMAP2"
y_axis <- "UMAP1"


# create plot
p <- plotly::plot_ly(
  plot_space,
  x = ~ UMAP2, y = ~UMAP1, 
  text = ~spp_sex,
  hoverinfo = "text"
)

# save 
htmlwidgets::saveWidget(p, paste0("neoaves_umap_2d.html"))


# 09/05/2025 ----
# Get list of species present in UK
# clear environment
rm(list=ls())

# Load libraries
library(raster)
library(dplyr)

# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neoaves"
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

# Load shared functions
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "mapping.R"
  )
)

# Load PCA data
pca_filename <- paste(clade, sex_match, "patches.231030.PCAcolspaces", "rds", sep = ".")
pca_dat <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]][["x"]]


# load PAM
# select whether to use liberal, conservative, or specified PAM
pam_type <- "conservative"
# clip PAM to land only? ("_clipped" or "")
pam_seas <- "_clipped"
# select PAM grid cell resolution
pam_res <- "50km"
# enter PAM files location
pams_filepath <- "X:/cooney_lab/Shared/Rob-MacDonald/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/PAMs"

pam_filename <- paste0("PAM_birds_Behrman", pam_res, "_Pres12_Orig12_Seas12_", pam_type, pam_seas, ".rds")
pam_raw <- readRDS(
  paste(pams_filepath, pam_filename, sep = "/")
)

# extract null raster from PAM file
null_rast <- extract_null_rast(pam_raw)

# extract matrix values from PAM file
pam <- extract_pam_vals(pam_raw)

# remove raw PAM file (for RAM reasons)
rm(pam_raw)

# change "Genus species" style to "Genus_species" style in PAM colnames
colnames(pam) <- gsub(" ", "_", colnames(pam))

# get list of species included in diversity data
spp_list <- get_unique_spp(pca_dat)

# subset pam to only species in diversity data
pam <- subset_pam(pam, spp_list)

# subset diversity data to only species present in PAM
pca_dat <- subset_data_by_spp(pca_dat, colnames(pam))

gc()

# get map of UK
# get world map for plotting
get_country <- function(new_crs, spatvec, country){
  
  # get country map from rnaturalearth
  country <- rnaturalearth::ne_countries(country = country, returnclass = "sf")
  
  # transform to correct crs
  country <- sf::st_transform(country, crs = new_crs)
  
  # transform to spatvector, if required (need to do this if using base R plot)
  if(spatvec == TRUE){
    
    # transform to spatvector
    country <- terra::vect(country)
    
  }
  
  return(country)
}
# get uk map and transform to same CRS as diversity raster (for plotting)
uk <- get_country(new_crs = terra::crs(null_rast), spatvec = TRUE, country = "united kingdom")
# get papua new guinea map and transform to same CRS as diversity raster (for plotting)
pap_new_guinea <- get_country(new_crs = terra::crs(null_rast), spatvec = TRUE, country = "papua new guinea")

# create null raster, subsetted to extent of ecoregion
# set number of grid cells in PAM
ncells <- nrow(pam)
# first assign values to null_rast to keep track of grid cells
null_terrast_uk <- null_rast
terra::values(null_terrast_uk) <- 1:ncells
null_terrast_uk <- terra::mask(terra::crop(null_terrast_uk, uk), uk)

# use values of subsetted raster to subset PAM to ROI
pam_uk <- pam[terra::values(null_terrast_uk), ]

# get species which are present in subsetted PAM
# use different method if ecoregion occupies 1 or fewer grid cells
if(length(terra::values(null_terrast_uk)) < 2){
  
  uk_species_pres <- names(which(!is.na(pam_uk)))
  
} else {
  
  uk_species_pres <- names(which(apply(pam_uk, 2, function(x) any(!(is.na(x))))))
  
}

# same for PGN
# first assign values to null_rast to keep track of grid cells
null_terrast_pgn <- null_rast
terra::values(null_terrast_pgn) <- 1:ncells
null_terrast_pgn <- terra::mask(terra::crop(null_terrast_pgn, pap_new_guinea), pap_new_guinea)

# use values of subsetted raster to subset PAM to ROI
pam_pgn <- pam[terra::values(null_terrast_pgn), ]

# get species which are present in subsetted PAM
# use different method if ecoregion occupies 1 or fewer grid cells
if(length(terra::values(null_terrast_pgn)) < 2){
  
  pgn_species_pres <- names(which(!is.na(pam_pgn)))
  
} else {
  
  pgn_species_pres <- names(which(apply(pam_pgn, 2, function(x) any(!(is.na(x))))))
  
}



# Load UMAP

umap <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
    paste(clade, sex_match, "patches.nn.25.mindist.0.1", space, "pca1-12", "UMAP.rds", sep = ".")
  )
) %>% 
  magrittr::extract2("layout") %>% 
  as.data.frame() %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  )

# subset to UK species only
uk_umap <- umap %>% 
  filter(
    species %in% uk_species_pres
  )

# subset to Papua New Guinea species only
pgn_umap <- umap %>% 
  filter(
    species %in% pgn_species_pres
  )

# save
saveRDS(uk_umap, 
        here::here(
          "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
          "UK_Neoaves_UMAP.rds"
        ))
saveRDS(pgn_umap, 
        here::here(
          "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
          "PapuaNewGuinea_Neoaves_UMAP.rds"
        ))


# Create subsetted PCA data to get centroid distances
uk_pca <- pca_dat %>% 
  as.data.frame() %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  ) %>% 
  filter(
    species %in% uk_species_pres,
    sex == "M"
  )
pgn_pca <- pca_dat %>% 
  as.data.frame() %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  ) %>% 
  filter(
    species %in% pgn_species_pres,
    sex == "M"
  )

library(dispRity)

uk_centroids <- dispRity::dispRity(as.matrix(uk_pca[, 1:30]), metric = centroids, centroid = 0)$disparity[[1]][[1]]
rownames(uk_centroids) <- rownames(uk_pca)
pgn_centroids <- dispRity::dispRity(as.matrix(pgn_pca[, 1:30]), metric = centroids, centroid = 0)$disparity[[1]][[1]]
rownames(pgn_centroids) <- rownames(pgn_pca)

# interactive umap
uk_umap %>% 
  filter(sex == "M") %>% 
plotly::plot_ly(x = ~-V2,
                y = ~-V1,
                text = ~species,
                hoverinfo = "text")

# get mean centroid distance
mean(uk_centroids)
mean(pgn_centroids)


# plot UMAP with arrows going from points to centroid
uk_umap_M <- uk_umap %>%
  filter(sex == "M")
pgn_umap_M <- pgn_umap %>% 
  filter(sex == "M")
# plot arrows only - UK
svg("uk_arrows.svg")
plot(-uk_umap_M$V2, -uk_umap_M$V1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "", asp = T, lwd = 0.5)
# Add arrows from each point to the origin
arrows(x0 = 0, , y0 = 0,  , x1 = -uk_umap_M$V2, y1 = -uk_umap_M$V1, col = "#FA2481", length = 0.05)
dev.off()

# plot arrows only - PNG
svg("pgn_arrows.svg")
plot(-pgn_umap_M$V2, -pgn_umap_M$V1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "", asp = T, lwd = 0.5)
# Add arrows from each point to the origin
arrows(x0 = 0, , y0 = 0,  , x1 = -pgn_umap_M$V2, y1 = -pgn_umap_M$V1, col = "#FA2481", length = 0.05)
dev.off()


# 12/05/2025 ----
# Plot SES diversity change for Pint of Science
# first run editable code in 03d_Patch_diversity_SES.R

# clear environment (except chosen variables)
rm(list=setdiff(ls(), c("clade", "space", "metric", "n_sims", "iucn_type", "avg_par")))

# load in results and get species richness as a proportion
# create dataframe to populate
res <- data.frame(matrix(NA, nrow = 0, ncol = 9))
colnames(res) <- c("IUCN", "SR", "Raw", "NULL_MEAN", "NULL_SD", "NULL_SE", "SES", "Metric", "SR_prop")
for (sex in c("M", "F", "All")){
  # set file path
  filename <- paste0(clade, "_patch_", space,  "_SES_", "nsims", n_sims, "_", avg_par, "_", metric, "_", sex, "_", iucn_type, "-iucn", ".csv")
  # load file and calculate proportional species richness
  temp_res <- read.csv(
    here::here(
      "2_Patches", "3_OutputData", "4_Diversity_measures", filename
    )
  ) %>% 
    mutate(
      sex = get("sex"),
      sr_prop = SR / max(SR)
    )
  # add to dataframe
  res <- rbind(res, temp_res)
}

# Adjust IUCN category names
res$IUCN <- gsub("All Species Retained", "No species extinct", res$IUCN)
res$IUCN <- gsub("Lost", "extinct", res$IUCN)
res$IUCN <- gsub("Only", "only", res$IUCN)
res$IUCN <- gsub("Retained", "surviving", res$IUCN)


# convert IUCN cat lost to factor and adjust levels
res$IUCN <- factor(res$IUCN, levels = c("No species extinct", "CR extinct", "EN extinct", "VU extinct", "NT extinct (LC only surviving)"))
# same for sex (to make plots appear in All-M-F order)
res$sex <- factor(res$sex, levels = c("All", "M", "F"))

# convert sex label list and labelling function to change facet labels
sex_vals <- list(
  'All' = "All specimens",
  'M' = "Males",
  'F' = "Females"
)
sex_labeller <- function(variable, value){
  return(sex_vals[value])
}

# plot (subset by sex if required)
# note that there's no point in adding error bars because the error's so small you can't see it
p <- res %>% 
  filter(
    sex == "All"
  ) %>%
  ggplot() + 
  geom_point(aes(x = sr_prop*100, y = SES, colour = IUCN), size = 3, shape = 17) + 
  geom_line(aes(x = sr_prop*100, y = SES), alpha = 0.6) +
  #  geom_errorbar(aes(ymin = SES - NULL_SE, ymax = SES + NULL_SE), width=0.005) +
  scale_x_reverse() + 
  #  scale_color_discrete(type = viridisLite::rocket(n = nlevels(res$IUCN))) + 
  scale_color_discrete(type = hcl.colors(n = nlevels(res$IUCN)+1, palette = "Terrain")) + 
  xlab("% of surviving species") + ylab("Change in colour diversity") + 
  labs(colour = "Species lost") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = -2, linetype = "dashed") +
  theme_light() + 
  # facet_wrap(~ sex, ncol = 3, labeller = sex_labeller) + 
  #  ggtitle(space) + 
  theme_bw() + 
  theme(text=element_text(size=12,  family="Century Gothic"))

p


# Save as svg
svg_filename <- paste0("G:/My Drive/Outreach/Pint of Science 2025/Resources/Images/Diversity loss/Neoaves_SES_diversity_loss.svg")
svg(
  svg_filename,
  width = 6, height = 6
)
p
dev.off()


# Find out which pheontypes are being lost ----

# Load UMAP
umap <- readRDS(
  paste0("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/2_UMAP/",
         "Neoaves.matchedsex.patches.nn.25.mindist.0.1.lab.pca1-12.UMAP.rds"
         )
) %>% 
  magrittr::extract2("layout") %>% 
  as.data.frame() %>% 
  rename(
    UMAP1 = V1,
    UMAP2 = V2
  ) %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2),
    UMAP1 = -UMAP1,
    UMAP2 = -UMAP2
  )

# load IUCN
iucn <- read.csv(
  "G:/My Drive/patch-pipeline/4_SharedInputData/iucn_2024_nominate.csv"
) %>% 
  rename(
    species = species_birdtree
  )

# join
dat <- umap %>% 
  left_join(
    iucn,
    by = "species"
  )

# filter to only threatened species and plot heatmap
dat %>% 
  filter(
    !(iucn_cat %in% c("LC", "DD") )
  ) %>% 
  ggplot(aes(x = UMAP2, y = UMAP1)) + 
  geom_hex(bins = 50)

# save as RDS to plot as colour grids
saveRDS(dat %>% 
          filter(
            !(iucn_cat %in% c("LC", "DD") )
          ),
        "Neoaves_threatened_only.rds")

# Load PCA to check which species are threatened and colourful
pca <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    "Neoaves.matchedsex.patches.231030.PCAcolspaces.rds"
  )
)[["lab"]] %>% 
  magrittr::extract2("x")

# get centroid distances
dists <- dispRity::dispRity(pca, dispRity::centroids)$disparity[[1]][[1]]
rownames(dists) <- rownames(pca)
# join to IUCN and filter
dists <- dists %>% 
  as.data.frame() %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2),
  ) %>% 
  left_join(
    iucn,
    by = "species"
  ) %>% 
  filter(
    !(iucn_cat %in% c("LC", "DD", "NT") )
  )



# 19/05/2025 ----
# PCA Loading heatmaps (for colour grids)

# clear environment
rm(list=ls())

# Load libraries ----
library(dplyr)
library(ggplot2)
library(terra)

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

axes <- paste0("PC", 1:2)
pixel_type <- "lab"
col_scale_type <- "axis_channel"

# Load data ----

# Load PCA data
pca_filename <- paste(clade, sex_match, "patches.231030.PCAcolspaces", "rds", sep = ".")
pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]]
pca_space <- pca_all

# Function to generate heatmaps
# Make a function to plot the PCA loading heatmaps
# 
# Note - requires dplyr, ggplot, tidyterra and raster packages to be loaded

#' Plot PCA loading heatmaps
#' Heatmaps show how pixel/channel combinations are loaded on each PC axis
#' Note - function requires dplyr, ggplot and raster packages to be loaded
#'
#' @param pca_space A prcomp object containing the PCA colour pattern space
#' @param axes PC axes to plot heatmaps for
#' @param pixel_type What do the pixel values represent? One of "srgb", "vurgb", "vrgb", "urgb", "tcsxyz", "lab", "rawusml", "jndxyzlum", "jndxyzlumr", "tcsxyzlum", "tcsxyzlumr", "rawusmldbl", "rawusmldblr"
#' @param sample_raster A sample raster from which to extract xy co-ordinates of the pixels
#' @param write_path Path and filename to write plot to, if saving as a PNG
#' @param rgb_colours Should the plot be plotted in pseudo-RGB colouring? Only possible for pixrl_type == "srgb"
#' @param const_col_scale Should the plot use a unified colour scale for all PC axes, or should each axis have its own colour scale? Currently only const_col_scale == TRUE is implemented for any pixel type other than "srgb", "vrgb", "urgb", "vurgb".
#'
#' @return If write_path == NULL, returns a ggplot object of the heatmaps
#' @export
#'
#' @examples
loading_heatmap <- function(pca_space, axes, pixel_type, sample_raster, write_path = NULL, rgb_colours = FALSE, col_scale_type = c("constant", "axis", "axis_channel")){
  
  # define number of channels based on pixel type
  if(pixel_type == "srgb" | pixel_type == "vurgb" | pixel_type == "vrgb" | pixel_type == "urgb" |
     pixel_type == "tcsxyz" | pixel_type == "lab"){
    n_channels <- 3
  } else if(pixel_type == "rawusml" | pixel_type == "jndxyzlum" | pixel_type == "jndxyzlumr" | 
            pixel_type == "tcsxyzlum" | pixel_type == "tcsxyzlumr") {
    n_channels <- 4
  } else if(pixel_type == "rawusmldbl" | pixel_type == "rawusmldblr"){
    n_channels <- 5
  }
  
  # define channel names
  if(pixel_type == "srgb" | pixel_type == "vurgb" | pixel_type == "vrgb" | pixel_type == "urgb"){
    channel_names <- c("R", "G", "B")
  } else if(pixel_type == "tcsxyz"){
    channel_names <- c("x", "y", "z")
  } else if(pixel_type == "lab"){
    channel_names <- c("L", "a", "b")
  } else if(pixel_type == "rawusml"){
    channel_names <- c("U", "S", "M", "L")
  } else if(pixel_type == "jndxyzlum"){
    channel_names <- c("jndx", "jndy", "jndz", "jndlum")
  } else if(pixel_type == "jndxyzlumr"){
    channel_names <- c("jndx", "jndy", "jndz", "jndlumr")
  } else if(pixel_type == "tcsxyzlum"){
    channel_names <- c("x", "y", "z", "lum")
  } else if(pixel_type == "tcsxyzlumr"){
    channel_names <- c("x", "y", "z", "lumr")
  } else if(pixel_type == "rawusmldbl"){
    channel_names <- c("U", "S", "M", "L", "DBL")
  } else if(pixel_type == "rawusmldblr"){
    channel_names <- c("U", "S", "M", "L", "DBLr")
  }
  
  # create named lists of axes and channels
  names(axes) <- axes
  names(channel_names) <- channel_names
  
  # extract loadings from PCA space
  loadings <- pca_space %>% 
    magrittr::extract2(
      "rotation"
    ) 
  
  # subset to PCs of interest
  loadings <- as.data.frame(loadings[, axes])
  colnames(loadings) <- axes
  
  # set body part names
  body_parts <- c("thr", "cro", "bre", "nap", "bel", "man", "cov", "rum", "fli", "tai")
  # set numbers as names (for plotting later)
  names(body_parts) <- 1:length(body_parts)
  
  # add columns to specify channel and body part from rownames
  loadings <- loadings %>% 
    mutate(
      channel = sapply(strsplit(rownames(.), split = ".", fixed = TRUE), "[", 1),     # need backslash to escape full stop (as automatically recognised as regular expression). Need second backslash to escape first backslash
      body_part = sapply(strsplit(rownames(.), split = ".", fixed = TRUE), "[", 2)
    ) %>% 
    # set levels of body part (to ensure displayed correctly in grid)
    mutate(
      body_part = factor(body_part, levels = body_parts)
    )
  
  # add xy coordinates (used to transform to raster for easier plotting)
  # first get body part matrix in the correct layout - reference matrix to get xy coords
  bp_mat <- matrix(data = body_parts, nrow = 5, ncol = 2, byrow = TRUE)
  # reverse row order in matrix (otherwise grid will be plotted upside down)
  bp_mat <- apply(bp_mat, 2, rev)
  loadings_xy <- loadings %>% 
    mutate(
      x = sapply(body_part, function(val) which(bp_mat == val, arr.ind = TRUE)[2]),
      y = sapply(body_part, function(val) which(bp_mat == val, arr.ind = TRUE)[1])
    )
  # pivot longer
  loadings_xy_long <- loadings_xy %>% 
    tidyr::pivot_longer(
      cols = tidyr::starts_with("PC"),
      names_to = "PC",
      values_to = "loading"
    )
  
  if(col_scale_type == "constant" | col_scale_type == "axis"){
    ## Plot individual PCs
    ## For use with any pixel type
    
    # create a separate raster for each channel by subsetting
    # first initialise empty list
    n_axes <- length(axes)
    pc_spatrasters <- vector(mode = "list", length = n_axes)
    names(pc_spatrasters) <- axes
    for(axis in axes){
      # get data to use to create raster
      ras_dat <- loadings_xy_long[loadings_xy_long$PC == axis, ]
      # pivot wider to get channels in different layers
      ras_dat <- ras_dat %>% 
        tidyr::pivot_wider(
          names_from = "channel",
          values_from = "loading"
        )
      # order by body part
      ras_dat <- ras_dat[order(ras_dat$body_part), ]
      ras_dat <- ras_dat[, c("x", "y", channel_names)]
      # create raster and convert to spatraster
      pc_raster <- terra::rast(
        raster::rasterFromXYZ(ras_dat)
      )
      # assign to list
      pc_spatrasters[[axis]] <- pc_raster
      
    }
    
    # calculate colour scale
    
    if(col_scale_type == "constant"){
      # use grand min and max for constant colour scale
      col_lims <- c(min(loadings_xy_long[["loading"]]), max(loadings_xy_long[["loading"]]))
    } else if(col_scale_type == "axis"){
      # # use min and max of individual axes for individual colour scales
      # col_lims <- c(min(loadings_xy_long[loadings_xy_long[["PC"]] == axis, "loading"]), max(loadings_xy_long[loadings_xy_long[["PC"]] == axis, "loading"]))
    }
    
    # get limits of x and y coordinates
    x_min <- min(loadings_xy_long$x) - 0.5
    x_max <- max(loadings_xy_long$x) + 0.5
    y_min <- min(loadings_xy_long$y) - 0.5
    y_max <- max(loadings_xy_long$y) + 0.5
    
    # Create plot for each spatraster
    
    # first initialise list to store plots
    pc_plots <- vector(mode = "list", length = n_axes)
    
    # add names of list elements
    names(pc_plots) <- axes
    
    # calculate aspect ratio
    asp <- length(unique(loadings_xy_long$y)) / length(unique(loadings_xy_long$x))
    
    # loop to create each raster
    for(axis in axes){
      
      # create axis plot
      p <- ggplot() + 
        tidyterra::geom_spatraster(data = pc_spatrasters[[axis]]) + 
        xlim(x_min, x_max) + 
        ylim(y_min, y_max) +
        # scale_fill_viridis_c(option = "plasma", limits = col_lims) + 
        facet_grid(cols = vars(lyr)) + 
        theme(aspect.ratio = asp,
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.margin = margin(0, 0, 0, 0))
        
        
      
      
      
      # add colour scale
      if(col_scale_type == "constant"){
        p <- p + 
          scale_fill_viridis_c(option = "plasma", limits = col_lims) + 
          guides(fill = "none") # remove colour legend
      } else {
        p <- p + 
          scale_fill_viridis_c(option = "plasma")
      }
      
      # add to list
      pc_plots[[axis]] <- p
      
    }
    
    if(col_scale_type == "constant"){
      
      # create plot to extract legend
      leg_plot <- ggplot() +
        tidyterra::geom_spatraster(data = pc_spatrasters[[1]]) + 
        facet_wrap(~lyr) + 
        scale_fill_viridis_c(option = "plasma", limits = col_lims) +
        theme(legend.position = "right")
      # extract legend grob to pass to ggarrange
      leg_grob <- ggpubr::get_legend(leg_plot)
      
      # combine individual channel plots into one figure
      loadings_heatmap <- ggpubr::ggarrange(
        plotlist = pc_plots, 
        nrow = n_axes,
        # ncol = n_channels, 
        labels = axes, 
        # label.x = 0.5,
        # label.y = 0.2,
        align = "hv",
        widths = c(1, 1, 1),  # Adjust widths to minimize space
        heights = c(1, 1, 1),
        common.legend = TRUE,
        legend.grob = leg_grob,
        legend = "right"
      )
      
    } else if(col_scale_type == "axis"){
      
      # combine individual channel plots into one figure
      loadings_heatmap <- ggpubr::ggarrange(
        plotlist = pc_plots, 
        nrow = n_axes,
        # ncol = n_channels, 
        labels = axes, 
        # label.x = 0.5,
        # label.y = 0.2,
        align = "hv",
        widths = c(1, 1, 1),  # Adjust widths to minimize space
        heights = c(1, 1, 1),
        common.legend = FALSE,
        legend = "right"
      )
      
    }
    
    return(loadings_heatmap)
  }
  
  # For an individual legend for each channel of each axis
  if(col_scale_type == "axis_channel"){
    
    # Create a template raster
    template_rast <- terra::rast(bp_mat)
    
    # create rasters
    rasters <- lapply(
      axes, function(axis){
        
        lapply(
          channel_names, function(ch){
            
            dat_for_rast <- loadings_xy_long %>% 
              filter(PC == axis, channel == ch) %>% 
              arrange(body_part)
            
            # get values for raster
            vals <- dat_for_rast %>% 
              pull(loading) 
            
            # create raster
            rast <- dat_for_rast %>% 
              select(x, y) %>% 
              as.matrix() %>% 
              terra::rasterize(
                y = template_rast,
                values = rep(NA, times = nrow(.))
              )
            # assign values
            terra::values(rast) <- vals
            
            return(rast)
            
          }
        )
        
      }
    )
    
    # combine into a single raster
    rasters <- unlist(rasters, recursive = FALSE)
    heatmap_raster <- terra::rast(rasters)
    
    # plot
    p <- terra::plot(heatmap_raster, nc = length(channel_names), nr = length(axes), legend = TRUE, col = map.pal("plasma", 100))
    
    return(p)
    
  }


  
  
  
  rasters <- loadings_xy_long %>% 
    filter(PC == "PC1", channel == "L") %>% 
    select(x, y) %>% 
    as.matrix() %>% 
    terra::rasterize(
      y = template_rast,
      values = loadings_xy_long$loading[which(loadings_xy_long$PC == "PC1" & loadings_xy_long$channel == "L")]
    )
  
  
  # pivot longer (i.e. add single column for PC axis number)
  loadings_long <- loadings %>% 
    tidyr::pivot_longer(
      cols = tidyr::starts_with("PC"),
      names_to = "PC",
      values_to = "loading"
    )
  
  # Create a function to map loading values to colors
  colour_map <- function(loading_value, limits) {
    # Scale the value between 0 and 1
    scaled_value <- (loading_value - limits[1]) / (limits[2] - limits[1])
    # Return the corresponding color from the plasma palette
    return(viridis::viridis(100, option = "plasma")[ceiling(scaled_value * 99) + 1])
  }
  
  names(axes) <- axes
  names(channel_names) <- channel_names
  
  # create plot for each channel of each axis
  grobs <- lapply(axes, function(axis){
    
    # set up colour scale
    if(const_col_scale == TRUE){
      # use grand min and max for constant colour scale
      col_lims <- c(min(loadings_long[["loading"]]), max(loadings_long[["loading"]]))
    } else if(const_col_scale == FALSE){
      # use min and max of individual axes for individual colour scales
      col_lims <- c(min(loadings_long[loadings_long[["PC"]] == axis, "loading"]), max(loadings_long[loadings_long[["PC"]] == axis, "loading"]))
    }
    
    col_scale <- viridisLite::plasma(n = 256)
    
    lapply(channel_names, function(ch){
      
      # subset loadings to axis/channel of interest
      axch_lo <- loadings_long %>% 
        filter(
          channel == ch,
          PC == axis
        )
      
      # Set up the plotting area
      # grid::grid.newpage()
      # grid::pushViewport(grid::viewport(layout = grid::grid.layout(5, 2)))
      
      # Plot each color in a 5x2 grid without white space
      # first switch the names and elements of body_parts around, for indexing
      patch_nums <- setNames(as.numeric(names(body_parts)), body_parts)
      rect_grobs <- list()
      for (i in patch_nums) {
        
        # set up plotting params
        row <- ((i - 1) %/% 2) + 1
        col <- ((i - 1) %% 2) + 1
        
        # get loading for this body part
        loading_value <- as.numeric(axch_lo[body_parts == names(patch_nums[i]), "loading"])
        
        # Map loading value to color
        rect_color <- colour_map(loading_value, col_lims)
        
        rect_grobs[[i]] <- grid::rectGrob(x = (col - 1) / 2, y = (5 - row) / 5,
                        width = 0.5, height = 0.2,
                        just = c("left", "bottom"),
                        gp = grid::gpar(fill = rect_color,
                                        col = NA))
      }
      # combine rectangle into a single grob
      axch_grob <- grid::gTree(
        children = do.call(grid::gList, rect_grobs),
        name = paste(axis, ch, sep = "_")
      )
      
      return(axch_grob)
      
    })
      
  })
  
  # plot the grids using gridExtra
  
  # extract grobs into one list
  grobs_unified <- unlist(grobs, recursive = FALSE)
  
  # plot
  gridExtra::grid.arrange(
    grobs = grobs_unified,
    nrow = length(axes),
    ncol = length(channel_names),
    padding = grid::unit(1, "cm")
  )
  
  # Now actually plot the grids
  # Create a new, empty page
  grid::grid.newpage()
  # Create a layout with n_rows as the number of axes and n__cols as the number of channels
  n_rows <- length(axes)
  n_cols <- length(channel_names)
  n_grids <- n_rows * n_cols
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n_rows, n_cols)))
  
  n_grids <- length(axes)
  for(row_index in 1:n_rows){
    
    
    for(col_index in 1:n_cols){
      print(paste0("axis = ", axes[[row_index]], ", channel = ", channel_names[[col_index]], ", grid position = (", row_index, ", ", col_index, ")"))
      # Create a viewport for this position
      grid::pushViewport(grid::viewport(layout.pos.row = row_index, layout.pos.col = col_index))
      
      # Draw the grob
      grid::grid.draw(grobs[[axes[row_index]]][[channel_names[col_index]]])
      
      # Exit this viewport
      grid::popViewport()
    }
    
  }
  
  # END OF COLOUR GRID BIT #
  
  
  # get mapping of pixel coordinates to pixel values from sample raster
  
  # aggregate
  agg.fact <- 1
  sample_raster <- suppressWarnings(aggregate(sample_raster, fact=agg.fact))
  # convert raster to dataframe with xy coordinates
  dat <- as.data.frame(sample_raster, xy=T)
  # subset to complete cases
  dat <- dat[complete.cases(dat),]
  # add pixel number as column
  dat$pixel_number <- as.numeric(rownames(dat))
  
  # create copies of each coordinate pair (one for each input channel)
  coords <- dat[, c("x", "y", "pixel_number")]
  for(i in 2:n_channels){
    coords <- rbind(coords, dat[, c("x", "y", "pixel_number")])
  }
  
  # order according to pixel number
  coords <- coords[order(coords$pixel_number), ]
  
  
  # create df of channel type
  channels <- data.frame(
    channel = rep(
      channel_names, length.out = (nrow(dat) * n_channels)
    )
  )
  
  # join xy coordinates to loadings and add channel specification
  loading_coords <- cbind(coords, channels, loadings)
  
  if(rgb_colours == FALSE){
    
    ## Plot individual channels ----
    ## For use with any pixel type
    
    # create a separate raster for each channel by subsetting
    # first initialise empty list
    ch_spatrasters <- vector(mode = "list", length = n_channels)
    names(ch_spatrasters) <- channel_names
    for(ch in channel_names){
      # get data to use to create raster
      ras_dat <- loading_coords[loading_coords$channel == ch, ]
      # order by row/'pixel' number
      ras_dat <- ras_dat[order(ras_dat$pixel_number), ]
      ras_dat <- ras_dat[, c("x", "y", axes)]
      # create raster and convert to spatraster
      ch_raster <- terra::rast(
        rasterFromXYZ(ras_dat)
      )
      # assign to list
      ch_spatrasters[[ch]] <- ch_raster
      
    }
    
    # calculate aspect ratio
    asp <- length(unique(loading_coords$y)) / length(unique(loading_coords$x))
    
    # calculate colour scale
    if(const_col_scale == TRUE){
      # calculate grand max and min of loadings to set the same colour scale for each raster
      col_lims <- c(min(loadings), max(loadings))
    }
    
    # get limits of x and y coordinates
    x_min <- min(loading_coords$x) - 0.5
    x_max <- max(loading_coords$x) + 0.5
    y_min <- min(loading_coords$y) - 0.5
    y_max <- max(loading_coords$y) + 0.5
    
    # Create plot for each spatraster
    
    # first initialise list to store plots
    ch_plots <- vector(mode = "list", length = n_channels)
    
    # add names of list elements
    names(ch_plots) <- channel_names
    
    # loop to create each raster
    for(ch in channel_names){
      
      # create channel plot
      p <- ggplot() + 
        tidyterra::geom_spatraster(data = ch_spatrasters[[ch]]) + 
        xlim(x_min, x_max) + 
        ylim(y_min, y_max) +
        # scale_fill_viridis_c(option = "plasma", limits = col_lims) + 
        facet_grid(rows = vars(lyr)) + 
        theme(aspect.ratio = asp,
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.margin = margin(0, 0, 0, 0)) + 
        guides(fill = "none") # remove colour legend
      
      # add colour scale
      if(const_col_scale == TRUE){
        p <- p + 
          scale_fill_viridis_c(option = "plasma", limits = col_lims)
      } else {
        p <- p + 
          scale_fill_viridis_c(option = "plasma")
      }
      
      # add to list
      ch_plots[[ch]] <- p
      
    }
    
    # create plot to extract legend
    leg_plot <- ggplot() +
      tidyterra::geom_spatraster(data = ch_spatrasters[[1]]) + 
      facet_wrap(~lyr) + 
      scale_fill_viridis_c(option = "plasma", limits = col_lims) +
      theme(legend.position = "right")
    # extract legend grob to pass to ggarrange
    leg_grob <- ggpubr::get_legend(leg_plot)
    
    # combine individual channel plots into one figure
    loadings_heatmap <- ggpubr::ggarrange(
      plotlist = ch_plots, 
      ncol = n_channels, 
      labels = channel_names, 
      # label.x = 0.5,
      # label.y = 0.2,
      align = "hv",
      widths = c(1, 1, 1),  # Adjust widths to minimize space
      heights = c(1, 1, 1),
      common.legend = TRUE,
      legend.grob = leg_grob,
      legend = "right"
    )
    
    # write to file, if specified
    if(!is.null(write_path)){
      
      # set png width and height
      plot_width = 500*n_channels
      plot_height = 150*length(axes)
      
      # initialise saving
      png(write_path, width = plot_width, height = plot_height, res = 150)
      
      print(loadings_heatmap)
      
      dev.off()
    } else{
      return(loadings_heatmap)
    }
  } else if(rgb_colours == TRUE){
    
    # check that the pixel type is RGB
    if(!identical(channel_names, c("R", "G", "B"))) {
      stop("Pixel type must be RGB of some kind if rgb_colours == TRUE is selected")
    }
    
    # create function to normalise loading values in each axis/channel combination between the minimum and 1
    normalise <- function(x, max = NULL){
      if(is.null(max)){
        x / max(x)
      } else {
        if(max(x) > max){
          x / max(x)
        } else{
          x / max
        }
      }
    }
    # alternative normalisation function - normalises between 0 and 1, which means the lowest 
    # loading gets reset to 0 (could be confusing)
    # # normalise <- function(x, max = NULL){
    #   if(is.null(max)){
    #     (x - min(x)) / (max(x) - min(x))
    #   } else {
    #     if(max(x) > max){
    #       (x - min(x)) / (max(x) - min(x))
    #     } else{
    #       (x - min(x)) / (max - min(x))
    #     }
    #   }
    # }
    
    
    # initialise df to store data to plot in
    plot_dat <- data.frame(x = numeric(), y = numeric(), dir = character(), hex = character(), axis = character())
    
    # initialise for loop to run over each PC axis
    for(axis in axes){
      
      # select the PC of interest
      img_dat <- loading_coords[, c("x", "y", "channel", axis)] 
      colnames(img_dat) <- c(colnames(img_dat)[1:3], "pc_values")
      
      # add column for positive loadings only (set negative loadings to 0)
      img_dat$pos_pc_values <- img_dat$pc_values
      img_dat$pos_pc_values[img_dat$pc_values < 0] <- 0
      # add column for negative loadings only (set positive loadings to 0)
      img_dat$neg_pc_values <- img_dat$pc_values
      img_dat$neg_pc_values[img_dat$pc_values > 0] <- 0
      # make negative loadings positive (to allow for normalisation and conversion to RGB)
      img_dat$neg_pc_values <- -(img_dat$neg_pc_values)
      
      # normalise loadings between 0 and 1 (can either do this using the grand min/max across all PCs
      # of interest or just using the PC of interest. Using grand min will produce a plot with the same 
      # colour scale across all PCs, using PC-specific min/max will produce plot with colour scales specific to 
      # each PC axis)
      if(const_col_scale == TRUE){
        img_dat <- img_dat %>% 
          mutate(
            pos_pc_values = normalise(
              pos_pc_values, 
              max = max(loading_coords[, axes], -(min(loading_coords[, axes])))
            ),
            neg_pc_values = normalise(
              neg_pc_values, 
              max = max(loading_coords[, axes], -(min(loading_coords[, axes])))
            )
          ) 
      } else if(const_col_scale == FALSE){
        img_dat <- img_dat %>% 
          mutate(
            pos_pc_values = normalise(
              pos_pc_values,
              max = max(loading_coords[, axis], -(min(loading_coords[, axis])))
            ),
            neg_pc_values = normalise(
              neg_pc_values,
              max = max(loading_coords[, axis], -(min(loading_coords[, axis])))
            )
          )
      }
      
      # select only positive loadings, pivot wider to get separate columns for each channel and save as new df
      pos_img_dat <- img_dat %>% 
        dplyr::select(
          x, y, channel, pos_pc_values
        ) %>% 
        tidyr::pivot_wider(
          names_from = channel, values_from = pos_pc_values
        )
      # and the same for negative loadings
      neg_img_dat <- img_dat %>% 
        dplyr::select(
          x, y, channel, neg_pc_values
        ) %>% 
        tidyr::pivot_wider(
          names_from = channel, values_from = neg_pc_values
        )
      
      # convert rgb values to hex and add to new df for plotting. Pivot longer to allow faceting in plot
      plot_axis_dat <- pos_img_dat %>% 
        dplyr::select(
          x, y
        ) %>% 
        mutate(
          pos = rgb(pos_img_dat$R, pos_img_dat$G, pos_img_dat$B),
          neg = rgb(neg_img_dat$R, neg_img_dat$G, neg_img_dat$B)
        ) %>% 
        tidyr::pivot_longer(
          cols = c(pos, neg),
          names_to = "dir", 
          values_to = "hex"
        ) %>% 
        mutate(
          axis = axis
        )
      
      # append to dataframe for plotting
      plot_dat <- rbind(plot_dat, plot_axis_dat)
      
    }
    
    # convert axis to factor (so the axes are plotted in the right order)
    plot_dat$axis <- factor(plot_dat$axis, levels = axes)
    
    # calculate aspect ratio
    asp <- length(unique(img_dat$y)) / length(unique(img_dat$x))
    
    
    # Create the ggplot with points colored by the RGB values, faceted by axis and positive/negative loadings
    loadings_heatmap <- ggplot(plot_dat, aes(x = x, y = y)) +
      geom_point(aes(colour = hex), size = 1) +  # Adjust size as needed
      facet_grid(cols = vars(dir), rows = vars(axis)) + # facet by positive/negative
      scale_color_identity() +  # Use the colors as is
      theme_void() +   # Remove all axes, gridlines, etc.
      theme(aspect.ratio = asp, # Keep aspect ratio square if x and y are on the same scale
            plot.background = element_rect(fill = 'white'),
            str)
    
    if(!is.null(write_path)){
      
      # set png width and height
      plot_width = 1000
      plot_height = 137*length(axes)
      
      # save
      png(write_path, width = plot_width, height = plot_height, res = 240)
      print(loadings_heatmap)
      dev.off()
    } else {
      return(loadings_heatmap)
    }
    
  }
  
}

