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
 #   order != "Galliformes", 
 #   order != "Anseriformes"
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
# calculate SR & mean aesthetic attractiveness
res[1,2] <- length(all)
res[1,3] <- mean(trait_vec)

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
res[2,3] <- mean(trait_vec)

# Use parLapply to replace the for loop with a parallelised apply function for speed
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      coms <- sample(x = all, size = length(noCR), replace = FALSE)
      trait_vec <- matrix(trait[coms, ], ncol = dim(trait)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <-mean(trait_vec)
      
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
res[3,3] <- mean(trait_vec)

# Use parLapply to replace the for loop with a parallelised apply function for speed
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      coms <- sample(x = all, size = length(noEN), replace = FALSE)
      trait_vec <- matrix(trait[coms, ], ncol = dim(trait)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <- mean(trait_vec)
      
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
res[4,3] <- mean(trait_vec)

# Use parLapply to replace the for loop with a parallelised apply function for speed
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      coms <- sample(x = all, size = length(noVU), replace = FALSE)
      trait_vec <- matrix(trait[coms, ], ncol = dim(trait)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <- mean(trait_vec)
      
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
res[5,3] <- mean(trait_vec)

# Use parLapply to replace the for loop with a parallelised apply function for speed
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      coms <- sample(x = all, size = length(noNT), replace = FALSE)
      trait_vec <- matrix(trait[coms, ], ncol = dim(trait)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <- mean(trait_vec)
      
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
filename <- paste0(clade, "_iratebirds_", "SES_", "nsims", n_sims, "_", "attractiveness", "_", iucn_type, "-iucn", ".csv")
write.csv(
  res, 
  here::here(
    "2_Patches", "3_OutputData", clade, "4_Diversity_measures", filename
  ), 
  row.names = FALSE
)


# Plot results ----

# clear environment (except chosen variables)
rm(list=setdiff(ls(), c("clade", "space", "metric", "n_sims", "iucn_type")))

# load in results and get species richness as a proportion
filename <- paste0(clade, "_iratebirds_", "SES_", "nsims", n_sims, "_", "attractiveness", "_", iucn_type, "-iucn", ".csv")
res <- read.csv(
  here::here(
    "2_Patches", "3_OutputData",clade, "4_Diversity_measures", filename
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


# Find out which phenotypes are being lost ----

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

# filter to only non-threatened species
lc_umap <- dat %>% 
  filter(
    iucn_cat %in% c("LC", "NT")
  )

# filter to original, but with DD species removed
orig_umap <- dat %>% 
  filter(
    iucn_cat != "DD"
  )



# let's try creating a heatmap of the proportional loss in a gridded colourspace
# get UMAP1 and UMAP2 limits (padded by 0.1)
orig_limits <- list(
  UMAP1 = c(min(orig_umap$UMAP1) - 0.1, max(orig_umap$UMAP1) + 0.1),
  UMAP2 = c(min(orig_umap$UMAP2) - 0.1, max(orig_umap$UMAP2) + 0.1)
)
# define N bins for each axis
n_bins <- 25
bins <- list(
  UMAP1 = seq(from = orig_limits$UMAP1[[1]], to = orig_limits$UMAP1[[2]], length.out = n_bins+1),
  UMAP2 = seq(from = orig_limits$UMAP2[[1]], to = orig_limits$UMAP2[[2]], length.out = n_bins+1)
)
# get bin length
bin_lengths <- c(UMAP1 = bins$UMAP1[[2]] - bins$UMAP1[[1]], UMAP2 = bins$UMAP2[[2]] - bins$UMAP2[[1]])


# initialise a matrix to store results in
bin_counts_orig <- matrix(
  nrow = n_bins,
  ncol = n_bins
)
# count number of species in each bin in original umap
# UMAP1
for(umap1_bin_start in bins$UMAP1){
  
  # get bin indices
  umap1_bin_start_index <- which(bins$UMAP1 == umap1_bin_start)
  umap1_bin_end_index <- which(bins$UMAP1 == umap1_bin_start) + 1
  
  
  # skip iteration if final bin start index
  if(umap1_bin_start_index == n_bins+1){
    break
  }
  
  # define bin end
  umap1_bin_end <- bins$UMAP1[umap1_bin_end_index]
  
  # loop across UMAP2 bins
  for(umap2_bin_start in bins$UMAP2){
    
    # get bin indices
    umap2_bin_start_index <- which(bins$UMAP2 == umap2_bin_start)
    umap2_bin_end_index <- which(bins$UMAP2 == umap2_bin_start) + 1
    
    # skip iteration if final bin start index
    if(umap2_bin_start_index == n_bins+1){
      break
    }
    
    # get UMAP2 bin end
    umap2_bin_end <- bins$UMAP2[umap2_bin_end_index]
    
    # get species/sex within bin
    species_pres <- orig_umap[umap1_bin_start < orig_umap$UMAP1 & orig_umap$UMAP1 < umap1_bin_end & umap2_bin_start < orig_umap$UMAP2 & orig_umap$UMAP2 < umap2_bin_end, c("species", "sex")]
    
    # count number of species in bin
    species_count <- nrow(species_pres)
    
    # add species count to matrix
    bin_counts_orig[umap2_bin_start_index, umap1_bin_start_index] <- species_count
    

  }
  
}
bin_counts_lc <- matrix(
  nrow = n_bins,
  ncol = n_bins
)
# count number of species in each bin in lcinal umap
# UMAP1
for(umap1_bin_start in bins$UMAP1){
  
  # get bin indices
  umap1_bin_start_index <- which(bins$UMAP1 == umap1_bin_start)
  umap1_bin_end_index <- which(bins$UMAP1 == umap1_bin_start) + 1
  
  
  # skip iteration if final bin start index
  if(umap1_bin_start_index == n_bins+1){
    break
  }
  
  # define bin end
  umap1_bin_end <- bins$UMAP1[umap1_bin_end_index]
  
  # loop across UMAP2 bins
  for(umap2_bin_start in bins$UMAP2){
    
    # get bin indices
    umap2_bin_start_index <- which(bins$UMAP2 == umap2_bin_start)
    umap2_bin_end_index <- which(bins$UMAP2 == umap2_bin_start) + 1
    
    # skip iteration if final bin start index
    if(umap2_bin_start_index == n_bins+1){
      break
    }
    
    # get UMAP2 bin end
    umap2_bin_end <- bins$UMAP2[umap2_bin_end_index]
    
    # get species/sex within bin
    species_pres <- lc_umap[umap1_bin_start < lc_umap$UMAP1 & lc_umap$UMAP1 < umap1_bin_end & umap2_bin_start < lc_umap$UMAP2 & lc_umap$UMAP2 < umap2_bin_end, c("species", "sex")]
    
    # count number of species in bin
    species_count <- nrow(species_pres)
    
    # add species count to matrix
    bin_counts_lc[umap2_bin_start_index, umap1_bin_start_index] <- species_count
    
    
  }
  
}

# Now get proportional reduction matrix
bin_counts_prop <- 1 - (bin_counts_lc / bin_counts_orig)
# check if any cells have value 0 - if not, we can use this to replace NaN
which(bin_counts_prop == 0)
# there are some - let's set NaN as NA instead and use a heatmap function that can handle it
bin_counts_prop[is.nan(bin_counts_prop)] <- NA

heatmap(bin_counts_prop, Rowv = NA, Colv = NA, col = viridisLite::plasma(256), )



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
# PCA Loading heatmaps (for colour grids) ----

# clear environment
rm(list=ls())

# Load libraries ----
library(dplyr)
library(ggplot2)
library(terra)
library(tidyterra)

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

axes <- paste0("PC", 1:8)
pixel_type <- "lab"
col_scale_type <- "constant"

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
#' @param write_path Path and filename to write plot to, if saving as a PNG
#' @param rgb_colours Should the plot be plotted in pseudo-RGB colouring? Only possible for pixrl_type == "srgb"
#' @param col_scale_type Should the plot use a unified colour scale for all PC axes, or should each axis have its own colour scale? Individual colour scales for each channel of each axis are also supported
#' @param stack_channels If TRUE, take absolute of each patch loadings and stack channels to give absolute measure of patch important (not compatible with col_scale_type = "axis_channel")
#'
#' @return If write_path == NULL, returns a ggplot object of the heatmaps
#' @export
#'
#' @examples
patch_loading_heatmap <- function(pca_space, axes, pixel_type, write_path = NULL, rgb_colours = FALSE, col_scale_type = c("constant", "axis", "axis_channel"), stack_channels = FALSE){
  
  # check col_scale_type and stack_channels parameters are compatible
  if(col_scale_type == "axis_channel" & stack_channels == TRUE){
    
    col_scale_type <- "axis"
    
    warning("'axis_channel' colour scale type not compatible with stacked channel heatmap. Setting col_scale_type = 'axis'")
    
  }
  
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
  
  # for non-stacked (individual channel) plot
  if(stack_channels == FALSE){
    
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
      # p <- ggplot() + 
      #   tidyterra::geom_spatraster(data = heatmap_raster) + 
      #   facet_wrap(~lyr, ncol = n_channels) + 
      #   scale_fill_viridis_c(option = "plasma") + 
      #   theme_minimal()
      p <- terra::plot(heatmap_raster, nc = length(channel_names), nr = length(axes), legend = TRUE, col = map.pal("plasma", 100))
      
      return(p)
      
    }
    
  } else if(stack_channels == TRUE){
    
    # for summation of absolute patch loadings across channels
    
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
      # take the absolute and sum the individual absolute channel loadings
    ras_dat <- ras_dat %>% 
        mutate(
          absolute_loading = rowSums(abs(select(., all_of(unname(channel_names)))))
          ) %>% 
      select(
        - all_of(unname(channel_names))
      )
      # order by body part
      ras_dat <- ras_dat[order(ras_dat$body_part), ]
      ras_dat <- ras_dat[, c("x", "y", "absolute_loading")]
      # create raster and convert to spatraster
      pc_raster <- terra::rast(
        raster::rasterFromXYZ(ras_dat)
      )
      # assign to list
      pc_spatrasters[[axis]] <- pc_raster
      
    }
    
    # calculate colour scale
    
    if(col_scale_type == "constant"){
      all_vals <- lapply(pc_spatrasters, terra::values) %>% unlist()
      # use grand min and max for constant colour scale
      col_lims <- c(
        min(all_vals),
        max(all_vals)
      )
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
  
  

}

patch_loading_heatmap(pca_all, axes, pixel_type, col_scale_type = col_scale_type)

patch_loading_heatmap(pca_all, axes, pixel_type, col_scale_type = "constant", stack_channels = TRUE)

p <- patch_loading_heatmap(pca_all, axes, pixel_type, col_scale_type = col_scale_type, stack_channels = TRUE)

folder_path <- here::here(
  "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", space, "loading_heatmaps"
)
if(!dir.exists(folder_path)){
  dir.create(folder_path, recursive = TRUE)
}
filename <- paste(clade, sex_match, sex_interest, "patches", paste0(axes[1], "-", tail(axes, 1)), "stackedchannels", col_scale_type, "loading_heatmaps.png", sep = "_")
png(paste(folder_path, filename, sep = "/"), width = 4, height = 30, units = "in", res = 600)
p
dev.off()

# 21/05/2025 ----
# Plot centroid distance on UMAP


# clear environment
rm(list=ls())

# Load libraries ----
library(dplyr)
library(ggplot2)

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

# Load data ----

# Load PCA data
pca_filename <- paste(clade, sex_match, "patches.231030.PCAcolspaces", "rds", sep = ".")
pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]][["x"]]

# Load UMAP
umap_filename <- paste(clade, sex_match, "patches", space, "pca.canonUMAP", "rds", sep = ".")
umap_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
    umap_filename
  )
)[["layout"]]

# Get centroid distances
centr_dists <- dispRity::dispRity(pca_all, dispRity::centroids)[["disparity"]][[1]][["elements"]]
rownames(centr_dists) <- rownames(pca_all)

# add spp/sex column (to join to UMAP data)
centr_dists <- centr_dists %>% 
  as.data.frame() %>% 
  mutate(spp_sex = rownames(.)) %>% 
  rename(centr_dist = V1)
umap_all <- umap_all %>% 
  as.data.frame() %>% 
  mutate(spp_sex = rownames(.)) %>% 
  rename(
    UMAP1 = V1,
    UMAP2 = V2
  )

# join by spp/sex
plot_data <- umap_all %>% 
  inner_join(
    centr_dists,
    by = "spp_sex"
  )

# plot, with point colour as distance to centroid
plot_data %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, colour = centr_dist)) + 
  geom_point() + 
  scale_colour_viridis_c()

plotly::plot_ly(
  data = plot_data,
  x = ~UMAP1, y = ~UMAP2,
  color = ~log(centr_dist),
  type = "scatter",
  text = ~spp_sex,
  colors = viridisLite::turbo(100)
)


# 04/06/2025 ----

# Meta-visualisation of colour space
# Ma et al. 2023 Nat. Comm.
# https://doi.org/10.1038/s41467-023-36492-2
# https://github.com/rongstat/meta-visualization/tree/main
# https://rongstat.github.io/metaviz_guide.io/user_guide.html

# clear environment
rm(list=ls())

# Load libraries
library(dplyr)
library(ggplot2)

# Load shared functions ----
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "meta_visualisation.R"
  )
)

## EDITABLE CODE ## ----
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Passeriformes"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# select sex of interest ("all", "male_female", "male_only", "female_only", "unknown_only")
sex_interest <- "male_female"
# choose colour space to use
space <- "lab"

# Load in PCA colourspace
pca_space <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    paste(clade, sex_match,  "patches.231030.PCAcolspaces.rds", sep = ".")
  )
) %>% 
  magrittr::extract2(space) %>% 
  magrittr::extract2("x")

# test on a small subset first - 1000 species
all_species <- unique(sapply(strsplit(rownames(pca_space), split = "-"), "[", 1))
random_species <- sample(all_species, size = 1000, replace = FALSE)
random_sexed_species <- c(
  paste0(random_species, "-M"),
  paste0(random_species, "-F")
)
random_subset <- pca_space[which(rownames(pca_space) %in% random_sexed_species), ]
info <- rownames(random_subset)

# get candidate visualisations
candidate.out <- candidate.visual(pca_space, methods = c("MDS", "iMDS", "Sammon", "LLE", "HLLE","Isomap",
                                                         "kPCA", "LEIM", "UMAP", "tSNE"),
                                  umap.k = c(15, 20)
                                  )

# plot all visualisations
for(k in 1:length(candidate.out$method_name)){
  data.plot = data.frame(dim1=candidate.out$embed.list[[k]][,1], dim2=candidate.out$embed.list[[k]][,2])
  p <- ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5) +
    scale_shape_manual(values=1:nlevels(data.plot$cluster)) + ggtitle(candidate.out$method[k])
  print(p)
}

# get meta visualisation
ensemble.out <- ensemble.viz(candidate.out$embed.list, candidate.out$method_name)

# assess candidate visualisations - boxplots of eigenscores
n <- dim(pca_space)[1]
data.plot = data.frame(eigen.score = c(ensemble.out$eigenscore), method = rep(candidate.out$method_name, each=n))
ggplot(data.plot, aes(x=reorder(method, eigen.score, FUN=median), y=eigen.score)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("eigenscore") + xlab("method")

# visualise eigenscores
for(k in 1:length(candidate.out$method_name)){
  data.plot = data.frame(dim1=candidate.out[[1]][[k]][,1], dim2=candidate.out[[1]][[k]][,2], 
                         eigenscore = c(ensemble.out[[2]][,k]))
  p <- ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(color=eigenscore)) +
    ggtitle(candidate.out$method[k])
  print(p)
}

# apply UMAP to meta-distance matrix to obtain final meta-visualisation
ensemble.data=umap(as.dist(ensemble.out$ensemble.dist.mat),  n_neighbors = 50)

data.plot = data.frame(dim1=ensemble.data[,1], dim2=ensemble.data[,2])
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5) +
  scale_shape_manual(values=1:nlevels(data.plot$cluster))+ ggtitle("meta-spec")

rownames(data.plot) <- info

# save to plot with colour grids
saveRDS(data.plot, file = "meta_visualisation.rds")

plot_space <- data.plot
x_axis <- "dim1"
y_axis <- "dim2"
xlabel <- "dim1"
ylabel <- "dim2"

# plot coloured by distance to centroid
centr_dists <- dispRity::dispRity(pca_space, metric = dispRity::centroids)$disparity[[1]][[1]]
rownames(centr_dists) <- rownames(pca_space)

plot_space <- merge(plot_space, centr_dists, by=0)

plot_space %>% 
  ggplot(aes(x = dim1, y = dim2, colour = log(centr_dists))) + 
  geom_point() + 
  scale_colour_viridis_c()
# actually doesn't represent this very well


# 18/06/2025 ----

# interactive species richness map ----

# load in PAM as in mapping scripts


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


## EDITABLE CODE ##
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

gc()

# get species richness by grid cell
species_richness <- calc_species_richness(pam)

# create raster of species richness 
sr_raster <- make_sr_raster(species_richness, null_rast)

# from claude

library(terra)
library(leaflet)
library(sf)

# 1. Convert to dataframe to get individual cell coordinates
richness_df <- as.data.frame(sr_raster, xy = TRUE, na.rm = TRUE)
names(richness_df) <- c("x", "y", "richness")  # Keep original names first

# 2. Transform coordinates from projected CRS to WGS84 (lat/lon)
# Create an sf object with the projected coordinates
richness_sf_points <- st_as_sf(richness_df, coords = c("x", "y"), crs = st_crs(sr_raster))

# Transform to WGS84 (EPSG:4326)
richness_sf_lonlat <- st_transform(richness_sf_points, crs = 4326)

# Extract the transformed coordinates
coords_lonlat <- st_coordinates(richness_sf_lonlat)
richness_df$lon <- coords_lonlat[, 1]
richness_df$lat <- coords_lonlat[, 2]

# 3. Get the cell numbers for non-NA cells (same as before)
cell_numbers <- which(!is.na(values(sr_raster)))

# 4. Create species lists for each non-NA cell
species_lists <- sapply(cell_numbers, function(cell_num) {
  if(cell_num <= nrow(pam)) {
    x <- pam[cell_num, ]
    species_present <- colnames(pam)[x == 1 & !is.na(x)]
    if(length(species_present) > 0) {
      paste(species_present, collapse = "<br>")
    } else {
      "No species recorded"
    }
  } else {
    "Cell index out of range"
  }
})

# 5. Add species lists to dataframe
richness_df$species_list <- species_lists

# 6. Create the interactive map
pal <- colorNumeric(palette = "YlOrRd", domain = richness_df$richness)

leaflet(richness_df) %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~lon, 
    lat = ~lat,
    radius = 4,
    color = "white",
    weight = 1,
    fillColor = ~pal(richness),
    fillOpacity = 0.8,
    popup = ~paste(
      "<b>Species Richness:</b>", richness, "<br><br>",
      "<b>Species present:</b><br>", species_list
    )
  ) %>%
  addLegend(
    "bottomright", 
    pal = pal, 
    values = ~richness,
    title = "Species<br>Richness"
  )

# here are some species present in a species-depauperate but highly colour diverse grid cell in
# the Tiris-Zemmour region of Northern Mauritania
tiris_zemmour_spp <- c("Pterocles_senegallus",
                       "Pterocles_coronatus",
                       "Oenanthe_leucopyga",
                       "Eremopterix_nigriceps",
                       "Corvus_ruficollis",
                       "Ammomanes_deserti",
                       "Alaemon_alaudipes",
                       "Lanius_excubitor",
                       "Chlamydotis_undulata")

# let's plot the UMAP space with these species highlighted in red to see how they're distributed
# in colour space
umap_filepath <- here::here(
  "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
  paste(clade, sex_match, "patches.nn.25.mindist.0.1", space, "pca1-12.UMAP.rds", sep = ".")
)
umap_space <- readRDS(umap_filepath) %>% 
  magrittr::extract2("layout")

# plot space without t-z-spp
t_z_mf <- c(
  paste0(tiris_zemmour_spp, "-M"),
  paste0(tiris_zemmour_spp, "-F")
)
umap_space %>% 
  as.data.frame() %>% 
  mutate(
    spp_sex = rownames(.)
  ) %>% 
  mutate(
    t_z = ifelse(spp_sex %in% t_z_mf, 1, 0.5)
  ) %>% 
  ggplot(aes(x = V1, y = V2, colour = t_z, alpha = t_z)) + 
  geom_point() + 
  scale_colour_viridis_c("turbo", direction = -1)

# now the same for South Georgian species
south_georgia_spp <- c(
  "Thalassoica_antarctica",
  "Thalassarche_melanophrys",
  "Thalassarche_chrysostoma",
  "Sterna_vittata",
  "Puffinus_puffinus",
  "Puffinus_griseus",
  "Puffinus_gravis",
  "Pterodroma_mollis",
  "Pterodroma_incerta",
  "Procellaria_cinerea",
  "Procellaria_aequinoctialis",
  "Phoebetria_palpebrata",
  "Phoebetria_fusca",
  "Pelecanoides_urinatrix",
  "Pagodroma_nivea",
  "Pachyptila_turtur",
  "Pachyptila_desolata",
  "Pachyptila_belcheri",
  "Lugensa_brevirostris",
  "Garrodia_nereis",
  "Fulmarus_glacialoides",
  "Diomedea_exulans",
  "Daption_capense",
  "Catharacta_antarctica",
  "Pterodroma_macroptera"
)

s_g_mf <- c(
  paste0(south_georgia_spp, "-M"),
  paste0(south_georgia_spp, "-F")
)
umap_space %>% 
  as.data.frame() %>% 
  mutate(
    spp_sex = rownames(.)
  ) %>% 
  mutate(
    s_g = ifelse(spp_sex %in% s_g_mf, 1, 0.5)
  ) %>% 
  ggplot(aes(x = V1, y = V2, colour = s_g, alpha = s_g)) + 
  geom_point() + 
  scale_colour_viridis_c("turbo", direction = -1)


# 26/06/2025 ----

# Create function to plot colour grids ----
# function is in R/plotting.R

# load data and call function

# Load libraries ---- 
library(dplyr)
library(ggplot2)
library(extrafont)

# clear environment
rm(list=ls())

# load plotting functions
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "plotting.R"
  )
)

# Choose parameters ----
## EDITABLE CODE ##
## Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neognaths"
## Choose colour space mapping
## Note that if usmldbl or usmldblr, will have to manually modify date in filename below
space <- "lab"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# plot all specimens ("all"), males only ("M"), or females only ("F")?
sex_focus <- "all"
## UMAP or PCA space?
## FALSE - use the PCA space
## TRUE - load a UMAP space from a file
## "perform" - load a PCA space from a file, then perform UMAP on it
## Note that if umap == TRUE, user will have to manually set path to umap file
load_umap <- FALSE
umap_filepath <- here::here(
  "2_Patches", "3_OutputData", clade, "2_PCA_ColourPattern_spaces", "2_UMAP",
  paste(clade, sex_match, "patches", "nn.25.mindist.0.1", space, "UMAP.rds", sep = ".")
)
# umap_filepath <- here::here(
#   "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
#   "Neoaves.patches.pca.jndxyzlumr.UMAPs.iterations.240326.rds"
# )
# If load_umap == FALSE, can set perform_umap == TRUE
# to perform UMAP on the loaded PCA and then plot the UMAP axes
perform_umap <- FALSE
## Choose PCs to plot (if plotting PCA)
x_axis <- "PC2"
y_axis <- "PC3"
## Choose aspect ratio 
## -"square" fixes to size 2500x2500px 
## - "wrap" sets the axis with the larger range to 2500px and scales the other accordingly)
## "square" is better if you want to combine multiple plots into one figure, "wrap" 
## is better for a single plot as it wraps the size of the plot to the range sizes
## Note that the aspect ratio of the plot itself is always fixed to 1, this parameter is 
## solely to control the aspect ratio of the output png
asp_ratio <- "wrap"
## Choose font "default" or name font (e.g. "Century Gothic", "AvenirNext LT Pro Regular)
## Use extrafont::fonts() to see available fonts
font_par <- "Century Gothic"

# Load data ----
umap <- readr::read_rds(umap_filepath)
pca_filepath <- here::here(
  "2_Patches", "3_OutputData", clade, "2_PCA_ColourPattern_spaces", "1_Raw_PCA", 
  paste(clade, sex_match, "patches.250716.PCAcolspaces", "rds", sep = ".")
  #    paste0("Neoaves.patches.231030.PCAcolspaces.", space, ".240829.rds")
)
# load data from file
prcomp <- readr::read_rds(pca_filepath)[[space]]

pca_mat <- prcomp[["x"]]

# set path to colour grids folder
grid_path <- "C:/Users/bop23rxm/Documents/colour_grids_repositioned"

# set save path
save_path <- here::here(
  "junk",
  "test_colspace_plotting.png"
)
save_path <- here::here(
  "2_Patches", "4_OutputPlots", clade, "1_Colourspace_visualisation", space,
  paste(clade, sex_match, sex_focus, "patches_UMAPnn25mindist0.1.png", sep = "_")
)

# try plotting space with the function, with different inputs

# UMAP object input
plot_patch_grids(
  umap_obj = umap,
  colour_grid_path = grid_path, 
  asp_ratio = "wrap",
  save_as = "png",
  save_path = save_path
)

# prcomp object input
plot_patch_grids(
  prcomp_obj = prcomp,
  x_axis = "PC25", y_axis = "PC26",
  colour_grid_path = grid_path, 
  asp_ratio = "wrap",
  save_as = "png",
  save_path = save_path
)

# data matrix input (pca matrix)
# filter to males only for speed
pca_mat <- pca_mat %>% 
    as.data.frame() %>% 
    mutate(
      species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
      sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
    ) %>% 
    filter(
      sex == "M"
    ) %>% 
    select(
      -species, -sex
    )

plot_patch_grids(
  data_matrix = pca_mat,
  x_axis = "Comp1", y_axis = "Comp2",
  colour_grid_path = grid_path, 
  asp_ratio = "wrap",
  save_as = "png",
  save_path = save_path
)

# individual column input
plot_patch_grids(
  x_axis = pca_mat$PC1, y_axis = pca_mat$PC5,
  x_label = "PC1", y_label = "PC5",
  row_names = rownames(pca_mat),
  colour_grid_path = grid_path,
  asp_ratio = "wrap",
  save_as = "png",
  save_path = save_path
)



# try plotting four plots together (UMAP and first 6 PC axes)
folder_path <- here::here(
  "2_Patches", "4_OutputPlots", clade, "1_Colourspace_visualisation", space
)
filename <- paste(clade, sex_match, sex_focus, "patches_PC1-PC8.png", sep = "_")
plot_four_cg(
  prcomp, "PC1", "PC2",
  prcomp, "PC3", "PC4",
  prcomp, "PC5", "PC6",
  prcomp, "PC7", "PC8",
  cg_path = grid_path,
  save_type = "png",
  write_folder = folder_path,
#  thin_number = 500,
  file_name = filename
)


# 07/07/2025 ----
# Plot PACA (from KMult)

# load plotting functions
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "plotting.R"
  )
)

# load Kmult
kmult <- readRDS("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/3_Phylogenetic_signal/lab/kMult.M.hackettTrees1.rds")

# set path to colour grids folder
grid_path <- "C:/Users/bop23rxm/Documents/colour_grids_repositioned"

# set save path
save_path <- here::here(
  "junk",
  "test_pacaspace_plotting.png"
)

# extract PACA embeddings
paca_mat <- kmult$PACA$x

# add sex to rownames
rownames(paca_mat) <- paste0(rownames(pca_mat), "-M")

# plot
plot_patch_grids(
  data_matrix = paca_mat,
  x_axis = "Comp1", y_axis = "Comp2",
  colour_grid_path = grid_path, 
  asp_ratio = "wrap",
  save_as = "png",
  save_path = save_path
)

plot_four_cg(
  paca_mat, "Comp1", "Comp2",
  paca_mat, "Comp3", "Comp4",
  paca_mat, "Comp5", "Comp6",
  paca_mat, "Comp7", "Comp8",
  cg_path = grid_path,
  save_type = "png",
  write_folder = here::here("junk"),
  # thin_number = 500,
  file_name = "test_paca_axes1-8.png"
)

# 29/07/2025 ----
# create PCA space with species and distance to centroid for Jas ----

library(dispRity)
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "plotting.R"
  )
)

# load PCA space
spaces <- readRDS("G:/My Drive/patch-pipeline/2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/1_Raw_PCA/Neoaves.allspecimens.patches.231030.PCAcolspaces.rds")
lab <- spaces$lab
usml <- spaces$usml

# extract coords
lab <- lab$x
usml <- usml$x

# get distances ot centroid
lab_dists <- dispRity(lab, centroids)$disparity[[1]][[1]]
rownames(lab_dists) <- rownames(lab)
usml_dists <- dispRity(usml, centroids)$disparity[[1]][[1]]
rownames(usml_dists) <- rownames(usml)

dists <- data.frame(
  species = sapply(strsplit(rownames(lab), split = "-"), "[", 1),
  sex = sapply(strsplit(rownames(lab), split = "-"), "[", 2),
  lab_centr_dist = lab_dists,
  usml_centr_dist = usml_dists
  )

write.csv(dists, file = here::here("junk", "centr_dists_for_jas.csv"))

plot_matrix <- as.matrix(dists[, c("lab_centr_dist", "usml_centr_dist")])
plot_matrix[, "lab_centr_dist"] <- scale(plot_matrix[, "lab_centr_dist"])
plot_matrix[, "usml_centr_dist"] <- scale(plot_matrix[, "usml_centr_dist"])
plot_patch_grids(
  x_axis = "lab_centr_dist",
  y_axis = "usml_centr_dist",
  data_matrix = plot_matrix,
  colour_grid_path = "C:/Users/bop23rxm/Documents/colour_grids_repositioned",
  asp_ratio = "wrap",
  save_as = "png",
  save_path = here::here("junk", "centr_dists_for_jas.png")
    
)

# 02/09/2025 ----
# Count neognath species ----

# libraries
library(dplyr)

# load patch data
px <- readRDS("./2_Patches/1_InputData/patches.250716.rds")


# load taxonomic data
taxo <- read.csv("./4_SharedInputData/BLIOCPhyloMasterTax_2019_10_28.csv")

# bind and remove any palaeognaths
px <- px %>% 
  left_join(
    taxo,
    by = join_by(species == TipLabel)
  ) %>% 
  filter(
    IOCOrder != "STRUTHIONIFORMES", 
    IOCOrder != "RHEIFORMES", 
    IOCOrder != "CASUARIIFORMES", 
    IOCOrder != "APTERYGIFORMES", 
    IOCOrder != "TINAMIFORMES", 
  )

# count unique species in our data
length(unique(px$species))

# filter to only species for which we have male and female data
px_meta <- px %>% 
  select(species, sex) %>% 
  distinct(species, sex)

px_m <- px_meta %>% 
  filter(sex == "M") %>% 
  select(species)
px_f <- px_meta %>% 
  filter(sex == "F") %>% 
  select(species)
px_matched <- px_m %>% 
  inner_join(px_f)
# count species for which we have matched data
n_matched <- nrow(px_matched)

# count total neognath species in taxonomy
n_taxo <- taxo %>% 
  filter(
    IOCOrder != "STRUTHIONIFORMES", 
    IOCOrder != "RHEIFORMES", 
    IOCOrder != "CASUARIIFORMES", 
    IOCOrder != "APTERYGIFORMES", 
    IOCOrder != "TINAMIFORMES", 
  ) %>% 
  nrow()

# get percentage coverage
pc_cov <- n_matched / n_taxo * 100

# get generic coverage
n_gen <- px_matched %>% 
  left_join(taxo, by = join_by(species == TipLabel)) %>% 
  select(GenusName) %>% 
  distinct() %>% 
  nrow()
n_gen_taxo <- taxo %>% 
  filter(
    IOCOrder != "STRUTHIONIFORMES", 
    IOCOrder != "RHEIFORMES", 
    IOCOrder != "CASUARIIFORMES", 
    IOCOrder != "APTERYGIFORMES", 
    IOCOrder != "TINAMIFORMES", 
  ) %>% 
  select(GenusName) %>% 
  distinct() %>% 
  nrow()
pc_gen_cov <- n_gen / n_gen_taxo * 100

# 17/09/2025 ----
# Investigate median SES noCR decline ----

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

# make wrapper function for averaging (allows selection of average type e.g. mean, median)
avg <- function(vals, avg_type, na.rm = TRUE){
  return(get(avg_type)(vals, na.rm = na.rm))
}
# make wrapper function for calculating either std dev (if using mean) or median absolute deviation (if using median)
dev_wrap <- function(vals, avg_type, na.rm = TRUE){
  if(avg_type == "mean"){
    return(sd(vals))
  } else if(avg_type == "median"){
    return(mad(vals))
  }
}

## EDITABLE CODE ## ----
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neognaths"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
space <- "lab"
# select sex ("M", "F", "All")
sex <- "M"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# select metric ("centr-dist", "nn-k", "nn-count")
# note that nn-k is EXTREMELY slow to run - needs parallelisation (but will probably still
# be too slow to run)
metric <- "centr-dist"
# select type of averaging to use ("mean" or "median")
avg_par <- "median"
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

# load patch data (PCA of whichever colourspace - generated in 02_Patch_Analyse_features.R)
pca_filename <- paste(clade, sex_match,  "patches.250716.PCAcolspaces", "rds", sep = ".")
pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", clade, "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]]

# load IUCN Red List data
iucn_filename <- paste0("neognath_iucn_2024_", iucn_type, ".csv")
iucn <- read.csv(
  here::here(
    "4_SharedInputData", iucn_filename
  )
)

####################    
### 2: GLOBAL ANALYSIS
####################   
## a) Morphological diversity
## NB: To generate Table S1, the code below can be ran on individual PCs.
## E.g. traits <- trait[,1]



# Matrix containing PCs:
trait <- pca_all %>% 
  magrittr::extract2(
    "x"
  ) %>% 
  as.data.frame() %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  ) 

# Species matched to IUCN categories:
species <- iucn %>% 
  rename(
    species = species_birdtree
  ) %>% 
  select(
    species, iucn_cat
  ) %>% 
  filter(
    species %in% trait[, "species"]
  )

# remove species not in IUCN data (now unnecessary since I've crossmatched the iucn data)
trait <- trait[trait[, "species"] %in% species[, "species"], ]

# check species are identical in both matrices
spec_iucn <- sort(species[, "species"])
spec_trait <- sort(unique(trait[, "species"]))
identical(spec_iucn, spec_trait)
rm(spec_iucn, spec_trait)

# add extra rows in species for M/F/U
species$sex <- "M"
species_F <- species; species_F$sex <- "F"
species <- rbind(species, species_F)
rm(species_F)
rownames(species) <- paste(species$species, species$sex, sep = "-")
species <- species[rownames(species) %in% rownames(trait), ]
# reorder to match trait data
species <- species[match(rownames(trait), rownames(species)), ]

# Subset to chosen sex
if(sex != "All"){
  trait <- trait[trait[, "sex"] == sex, ]
  species <- species[species[, "sex"] == sex, ]
  unique_species <- NULL # so the parallelisation cluster export doesn't throw an error
} else {
  # get unique species (to subset species/sex pairs later)
  unique_species <- unique(trait[["species"]])
  
}


# Use all PCS to account for 100% variation in traits
### All PCs
traits <- as.matrix(trait[,1:(ncol(trait) - 2)])

### OUTPUT DATA FRAME
res <- matrix(NA, nrow=3, ncol = 6)
colnames(res) <- c("sex", "subset", "mean", "sd", "median", "mad")
res[, "sex"] <- c(sex, sex, sex)
res[, "subset"] <- c("all", "noCR", "CR")


a <- n_sims #change depending on number of simulations  
sims <- matrix(NA, nrow=1, ncol = a) #place to store simulations for SES

# set dispRity metric
if(metric == "centr-dist"){
  metric_get <- "centroids"
} else if (metric == "nn-k"){
  metric_get <- "mean.nn.dist"
}else if (metric == "nn-count"){
  metric_get <- "count.neighbours"
}


all <- rownames(species)
# extract traits for species in the random community
trait_vec <- matrix(traits[all,], ncol=dim(traits)[2])
# calculate SR & Mean distance to centroid
all_distrib <- dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]]

# populate results
res[1, "mean"] <- mean(all_distrib)
res[1, "median"] <- median(all_distrib)
res[1, "sd"] <- sd(all_distrib)
res[1, "mad"] <- mad(all_distrib)

# get collections of species with threat levels sequentially trimmed
noCR <- rownames(species)[species$iucn_cat != "CR"]
CR <- rownames(species)[species$iucn_cat == "CR"]
EN <- rownames(species)[species$iucn_cat == "EN"]
VU <- rownames(species)[species$iucn_cat == "VU"]


# extract traits for species in the random community
trait_vec <- matrix(traits[noCR,], ncol=dim(traits)[2])
# calculate SR & Mean distance to centroid
noCR_distrib <- dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]]

# populate results
res[2, "mean"] <- mean(noCR_distrib)
res[2, "median"] <- median(noCR_distrib)
res[2, "sd"] <- sd(noCR_distrib)
res[2, "mad"] <- mad(noCR_distrib)

# for CR species only
# extract traits for species in the random community
trait_vec <- matrix(traits[CR,], ncol=dim(traits)[2])
# calculate SR & Mean distance to centroid
CR_distrib <- dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]]

# populate results
res[3, "mean"] <- mean(CR_distrib)
res[3, "median"] <- median(CR_distrib)
res[3, "sd"] <- sd(CR_distrib)
res[3, "mad"] <- mad(CR_distrib)

write.csv(res, here::here("junk", paste0("distribs_", sex, ".csv")))


# collate sex-specific stats (after running the above for all sexes)
all_distribs <- rbind(
  read.csv(here::here("junk", "distribs_All.csv")),
  read.csv(here::here("junk", "distribs_M.csv")),
  read.csv(here::here("junk", "distribs_F.csv"))
)
write.csv(all_distribs, here::here("junk", "distribs_collated.csv"))


# 29/09/2025
# Plot space with particular taxonomic groups highlighted ----

# Clear environment
rm(list=ls())

# load plotting functions ----
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "plotting.R"
  )
)


## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neognaths"
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
## Choose aspect ratio 
## -"square" fixes to size 2500x2500px 
## - "wrap" sets the axis with the larger range to 2500px and scales the other accordingly)
## "square" is better if you want to combine multiple plots into one figure, "wrap" 
## is better for a single plot as it wraps the size of the plot to the range sizes
## Note that the aspect ratio of the plot itself is always fixed to 1, this parameter is 
## solely to control the aspect ratio of the output png
asp_ratio <- "wrap"
## Choose font "default" or name font (e.g. "Century Gothic")
## Use extrafont::fonts() to see available fonts
font_par <- "Century Gothic"
# set path to colour grids folder
grid_path <- "C:/Users/bop23rxm/Documents/colour_grids_repositioned"

# Load data

# pca data
pca_filename <- paste(clade, sex_match, "patches.250716.PCAcolspaces", "rds", sep = ".")
pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", clade, "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]]

# UMAP data
umap_filename <- paste(clade, sex_match, "patches", "nn.25.mindist.0.1", space, "UMAP.rds", sep = ".")
umap_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", clade, "2_PCA_ColourPattern_spaces", "2_UMAP",
    umap_filename
  )
)

# load taxonomic data
taxo_master <- read.csv(
  here::here(
    "4_SharedInputData", "BLIOCPhyloMasterTax_2019_10_28.csv"
  )
)


# Function to plot space with particular group highlighted
# Modified version of the main colour grid plotting function
# Plot patches space with colour grids
plot_patch_taxogroups <- function(
    x_axis = NULL, y_axis = NULL,
    data_matrix = NULL,
    umap_obj = NULL, 
    prcomp_obj = NULL,
    colour_grid_path,
    row_names = NULL,
    x_label = NULL, y_label = NULL,
    asp_ratio = c("wrap", "square"),
    save_as = NULL, # currently only "png" and "plot_object" is implemented
    save_path = NULL,
    font_par = NULL,
    thin_number = NULL,
    fourplot = FALSE, fourplot_size = c(10, 10.5),
    taxonomy = NULL, # taxonomy - must contain family, order, or taxon_subgroup column
    taxo_level = c("family", "order", "taxon_subgroup"), # choose whether to highlight a family, order, or taxonomic subgroup (family for non-passerines, order for passerines)
    group_name = NULL
    
){
  
  # load dplyr (I should really rewrite everything to base R)
  require(dplyr)
  
  # function to test if axes are atomic, numeric vectors
  is.anv <- function(axis){
    if(is.atomic(axis) & is.numeric(axis) & is.vector(axis)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  # Check only one data type has been provided - stop if not
  if(is.anv(x_axis) & is.anv(y_axis)){
    if(!is.null(data_matrix) | !is.null(umap_obj) | !is.null(prcomp_obj)){
      stop("More than one data type provided.")
    }
  } else if(is.character(x_axis) & is.character(y_axis)){
    if(all(is.null(data_matrix), is.null(umap_obj), is.null(prcomp_obj))){
      stop("No data provided.")
    } else if((!is.null(data_matrix) + !is.null(umap_obj) + !is.null(prcomp_obj)) >= 2){
      stop("More than one data type provided.")
    }
  }
  
  # If highlighting a particular group, check columns are correct
  if(!is.null(taxonomy)){
    
    # check the relevant column is present in the taxonomy
    if(!(taxo_level %in% colnames(taxonomy))){
      stop("Chosen taxonomic subgroup not present in taxonomy")
    }
    
    # check there's a species column in the taxonomy
    if(!("species" %in% colnames(taxonomy))){
      stop("Taxonomy must contain a 'species' column")
    }
    
  }
  

  # Check what type of object has been provided and set up plot space
  
  # if axes have been provided as data columns
  if(is.anv(x_axis) & is.anv(y_axis)){
    
    # check if rownames have been provided and match length of data columns 
    # need the rownames to get the png colour grid names
    if(is.null(row_names)){
      stop("No rownames provided.")
    }
    
    plot_space <- cbind(x_axis, y_axis)
    rownames(plot_space) <- row_names
    
    x_axis <- "axis_1"; y_axis <- "axis_2"
    colnames(plot_space) <- c(x_axis, y_axis)
    
    # check if axis labels have been provided
    if(any(is.null(x_label), is.null(y_label))){
      warning("Axis labels not set. Setting to default labels.")
      # set to defaults if not provided
      x_label <- "Axis 1"; y_label <- "Axis 2"
    }
    
  } else if(is.character(x_axis) & is.character(y_axis)){
    # if axes have been provided as characters (e.g. "PC1", "PC2", "UMAP1", "UMAP2 etc)
    
    # Check if data has been provided as a prcomp object, a UMAP object, or as a data matrix
    # and extract data if necessary
    if(!is.null(prcomp_obj)){
      
      plot_space <- prcomp_obj %>% 
        magrittr::extract2("x") %>% 
        as.data.frame()
      
      # get proportions of variance for each axis
      pca_variance <- prcomp_obj %>% 
        summary() %>% 
        magrittr::extract2("importance") %>% 
        as.data.frame()
      
      
    } else if(!is.null(umap_obj)){
      
      plot_space <- umap_obj %>% 
        magrittr::extract2("layout") %>% 
        as.data.frame() %>% 
        rename(
          UMAP1 = V1, UMAP2 = V2
        )
      x_axis <- "UMAP1"; y_axis <- "UMAP2"
      
    } else if(!is.null(data_matrix)){
      
      plot_space <- data_matrix
      
      # check that named axes are in data matrix
      if(!(x_axis %in% colnames(plot_space) & y_axis %in% colnames(plot_space))){
        stop("Chosen axes are not columns in data_matrix.")
      }
      
    }
    
  } else if(any(is.null(x_axis), is.null(y_axis))){
    
    if(!is.null(umap_obj)){
      
      plot_space <- umap_obj %>% 
        magrittr::extract2("layout") %>% 
        as.data.frame() %>% 
        rename(
          UMAP1 = V1, UMAP2 = V2
        )
      x_axis <- "UMAP1"; y_axis <- "UMAP2"
      
    }
    
  }
  

  
  # if requested, thin to specific number of random species
  if(!is.null(thin_number)){
    
    all_spp <- sapply(strsplit(rownames(plot_space), split = "-"), "[", 1)
    
    unique_spp <- unique(sapply(strsplit(rownames(plot_space), split = "-"), "[", 1))
    
    spp_subset <- sample(unique_spp, thin_number, replace = FALSE)
    
    # subset plot space to only random species
    rows_to_keep <- which(all_spp %in% spp_subset)
    plot_space <- plot_space[rows_to_keep, ]
    
  }
  
  # set axis labels if not yet set
  if(!is.null(umap_obj)){
    x_label <- "UMAP1"
    y_label <- "UMAP2"
  } else if(!is.null(prcomp_obj)) {
    x_label <- paste0(x_axis, " (", round(pca_variance[, x_axis][2], 2) * 100, "% of variance)")
    y_label <- paste0(y_axis, " (", round(pca_variance[, y_axis][2], 2) * 100, "% of variance)")
  } else if(!is.null(data_matrix) & any(is.null(x_label), is.null(y_label))){
    warning("Axis labels not set. Setting to axis names.")
    # set to defaults if not provided
    x_label <- x_axis; y_label <- y_axis
  }
  
  
  # set up saving
  if(!is.null(save_as)){
    if(save_as == "png"){
      
      # set up aspect ratio
      if(asp_ratio == "wrap"){
        # determine axis with greatest range (to determine size of png)
        range_x_axis <- diff(range(plot_space[, x_axis]))
        range_y_axis <- diff(range(plot_space[, y_axis]))
        if(range_x_axis > range_y_axis){
          png_width <- 2500
          png_height <- png_width * (range_y_axis / range_x_axis)
        } else if(range_y_axis > range_x_axis){
          png_height <- 2500
          png_width <- png_height * (range_x_axis / range_y_axis)
        }
      } else if (asp_ratio == "square"){
        png_width <- 2500
        png_height <- 2500
      }
      
      # set up device
      if(is.null(font_par)){
        png(
          save_path,
          width = png_width, height = png_height,
          units = "px",
          pointsize = 40
        )
      } else {
        png(
          save_path,
          width = png_width, height = png_height,
          units = "px",
          pointsize = 40,
          family = font_par
        )
      }
      
    }
  }
  
  # fix aspect ratio as square, if requested
  if(asp_ratio == "square"){
    par(pty = "s")
  }
  
  dev.hold()
  
  plot(plot_space[, y_axis] ~ plot_space[, x_axis], 
       asp = T, 
       type = "n", 
       xlab = "", ylab = "", las = 1,
       # xaxt = "n", yaxt = "n"
       #     xlim = xrange,  ylim = yrange
  )
  title(xlab = x_label, ylab = y_label)
  
  box(lwd = 2)
  # if plotting as a single plot, define colour grid size relative to plot range
  if(fourplot == FALSE){
    rw <- diff(range(plot_space[,1]))/12
    rw_x <- rw / 15
    rw_y <- rw / 9
  } else if(fourplot == TRUE){
    # if plotting as part of a four-panel plot, define grid size as absolute relative to overall plot size
    x_range <- diff(range(plot_space[, x_axis]))
    grid_width_inches <- (fourplot_size[1] / 100) * 1.8
    grid_size_par <- (grid_width_inches / fourplot_size[1]) * x_range
    rw_x <- grid_size_par / 2
    rw_y <- rw_x * (5/3)
  }
  
  
  for(i in 1:nrow(plot_space)){
    fname <- rownames(plot_space)[i]
    fpng <- png::readPNG(paste0(
      colour_grid_path, "/",
      paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
    
    # check if species is in highlighted subgroup - make transparent if not
    image_spp <- strsplit(fname, split="-")[[1]][1]
    if(taxonomy[taxonomy$species == image_spp, taxo_level] != group_name){
      fpng_a <- array(1, dim = c(dim(fpng)[1], dim(fpng)[2], 4))
      fpng_a[, , 1:3] <- fpng
      fpng_a[, , 4] <- 0.1
      fpng <- fpng_a
    }
    
    rasterImage(fpng, 
                xleft = plot_space[i, x_axis] - rw_x, 
                ybottom = plot_space[i, y_axis] - rw_y, 
                xright = plot_space[i, x_axis] + rw_x, 
                ytop = plot_space[i, y_axis] + rw_y)
    cat(paste0("\r", i, " of ", nrow(plot_space), " processed"))
    
    if(isTRUE(all.equal(i/500, as.integer(i/500)))){
      # Force an update to the graphics device after each 500 images
      # This helps prevent memory issues
      if(interactive()) {
        Sys.sleep(0.01)  # Small delay to allow graphics to update
      }
      
      # Free up memory
      gc()
    }
    
  }
  
  dev.flush()
  
  if(!is.null(save_as)){
    dev.off()
  }
  
  
}


# test function

# clean taxonomy column names
taxo_clean <- taxo_master
colnames(taxo_clean) <- sub("TipLabel", "species", colnames(taxo_clean))
colnames(taxo_clean) <- sub("IOCOrder", "order", colnames(taxo_clean))
colnames(taxo_clean) <- sub("BLFamilyLatin", "family", colnames(taxo_clean))
colnames(taxo_clean) <- sub("Taxon_subgroup", "taxon_subgroup", colnames(taxo_clean))

# set taxon subgroup to highlight
hl_group <- "CORACIIFORMES"

# set plot name
plot_path <- here::here("junk", paste0("plot_highlighted_", hl_group, ".png"))

# plot, highlighting a particular taxonomic subgroup
plot_patch_taxogroups(umap_obj = umap_all, colour_grid_path = grid_path, asp_ratio = asp_ratio, save_as = "png", save_path = plot_path, taxonomy = taxo_clean, taxo_level = "taxon_subgroup", group_name = hl_group)

