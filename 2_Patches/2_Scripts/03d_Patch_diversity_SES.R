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

# make wrapper function for averaging (allows selection of average type e.g. mean, median)
avg <- function(vals, avg_type, na.rm = TRUE){
  return(get(avg_type)(vals, na.rm = na.rm))
}

## EDITABLE CODE ## ----
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neoaves"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
space <- "lab"
# select sex ("M", "F", "All")
sex <- "F"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# select metric ("centr-dist", "nn-k", "nn-count")
# note that nn-k is EXTREMELY slow to run - needs parallelisation (but will probably still
# be too slow to run)
metric <- "centr-dist"
# select type of averaging to use ("mean" or "median")
avg_par <- "mean"
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
pca_filename <- paste(clade, sex_match,  "patches.231030.PCAcolspaces", "rds", sep = ".")
pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]]

# load IUCN Red List data
iucn_filename <- paste0("iucn_2024_", iucn_type, ".csv")
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
head(traits)

### OUTPUT DATA FRAME
res <- matrix(NA, nrow=5, ncol = 6)

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
res[1,2] <- length(all)
res[1,3] <- avg(dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]], avg_par)

# get collections of species with threat levels sequentially trimmed
noCR <- rownames(species)[species$iucn_cat != "CR"]
noEN <- rownames(species)[species$iucn_cat != "CR" & species$iucn_cat != "EN"]
noVU <- rownames(species)[species$iucn_cat != "CR" & species$iucn_cat != "EN" & species$iucn_cat != "VU"]
noNT <- rownames(species)[species$iucn_cat != "CR" & species$iucn_cat != "EN" & species$iucn_cat != "VU" & species$iucn_cat != "NT"]

# for parallelisation, set up a cluster using the number of cores (2 less than total number of laptop cores)
no_cores <- parallel::detectCores() - 2
cl <- parallel::makeCluster(no_cores)
# export necessary objects to the cluster
parallel::clusterExport(cl, c("all", "noCR", "noEN", "noVU", "noNT", "traits", "unique_species", "sex", "metric_get", "dispRity", "avg", "avg_par", metric_get))

# extract traits for species in the random community
trait_vec <- matrix(traits[noCR,], ncol=dim(traits)[2])
# calculate SR & median distance to centroid
res[2,2] <- length(noCR)
res[2,3] <- avg(dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]], avg_par)

# Use parLapply to replace the for loop with a parallelised apply function for speed
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      # get random species to remove
      if(sex == "All"){
        spp_sample <- sample(x = unique_species, size = length(noCR)/2, replace = FALSE)
        coms <- c(paste0(spp_sample, "-M"), paste0(spp_sample, "-F"))
      } else {
        coms <- sample(x = all, size = length(noCR), replace = FALSE)
      }
  
      trait_vec <- matrix(traits[coms, ], ncol = dim(traits)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <- avg(dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]], avg_par)
      
      return(disparity_value)
    }
  )
)


# add results to table
res[2,4] <- mean(sims)
res[2,5] <- sd(sims)
res[2,6] <- sd(sims)/sqrt(length(sims))


# extract traits for species in the random community
trait_vec <- matrix(traits[noEN,], ncol=dim(traits)[2])
# calculate SR & Mean distance to centroid
res[3,2] <- length(noEN)
res[3,3] <- avg(dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]], avg_par)

# Use parLapply to replace the for loop with a parallelised apply function for speed
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      # get random species to remove
      if(sex == "All"){
        spp_sample <- sample(x = unique_species, size = length(noEN)/2, replace = FALSE)
        coms <- c(paste0(spp_sample, "-M"), paste0(spp_sample, "-F"))
      } else {
        coms <- sample(x = all, size = length(noEN), replace = FALSE)
      }
      
      trait_vec <- matrix(traits[coms, ], ncol = dim(traits)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <- avg(dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]], avg_par)
      
      return(disparity_value)
    }
  )
)

# add results to table
res[3,4] <- mean(sims)
res[3,5] <- sd(sims)
res[3,6] <- sd(sims)/sqrt(length(sims))

# extract traits for species in the random community
trait_vec <- matrix(traits[noVU,], ncol=dim(traits)[2])
# calculate SR & Mean distance to centroid
res[4,2] <- length(noVU)
res[4,3] <- avg(dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]], avg_par)

# Use parLapply to replace the for loop with a parallelised apply function for speed
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      # get random species to remove
      if(sex == "All"){
        spp_sample <- sample(x = unique_species, size = length(noVU)/2, replace = FALSE)
        coms <- c(paste0(spp_sample, "-M"), paste0(spp_sample, "-F"))
      } else {
        coms <- sample(x = all, size = length(noVU), replace = FALSE)
      }
      
      trait_vec <- matrix(traits[coms, ], ncol = dim(traits)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <- avg(dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]], avg_par)
      
      return(disparity_value)
    }
  )
)


res[4,4] <- mean(sims)
res[4,5] <- sd(sims)
res[4,6] <- sd(sims)/sqrt(length(sims))

# extract traits for species in the random community
trait_vec <- matrix(traits[noNT,], ncol=dim(traits)[2])
# calculate SR & Mean distance to centroid
res[5,2] <- length(noNT)
res[5,3] <- avg(dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]], avg_par)

# Use parLapply to replace the for loop with a parallelised apply function for speed
sims[1, ] <- unlist(
  parallel::parLapply(
    cl = cl,
    1:ncol(sims), 
    function(j) {
      
      # get random species to remove
      if(sex == "All"){
        spp_sample <- sample(x = unique_species, size = length(noNT)/2, replace = FALSE)
        coms <- c(paste0(spp_sample, "-M"), paste0(spp_sample, "-F"))
      } else {
        coms <- sample(x = all, size = length(noNT), replace = FALSE)
      }
      
      trait_vec <- matrix(traits[coms, ], ncol = dim(traits)[2])
      
      # Calculate the disparity value for the current simulation (column j)
      disparity_value <- avg(dispRity(trait_vec, metric = get(metric_get))$disparity[[1]][[1]], avg_par)
      
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
filename <- paste0(clade, "_patch_", space,  "_SES_", "nsims", n_sims, "_", avg_par, "_", metric, "_", sex, "_", iucn_type, "-iucn", ".csv")
write.csv(
  res, 
  here::here(
    "2_Patches", "3_OutputData", "4_Diversity_measures", filename
    ), 
  row.names = FALSE
  )


# Plot results ----

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

# convert IUCN cat lost to factor and adjust levels
res$IUCN <- factor(res$IUCN, levels = c("All Species Retained", "CR Lost", "EN Lost", "VU Lost", "NT Lost (LC Only Retained)"))
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
  # filter(
  #   sex != "All"
  # ) %>%
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
  geom_hline(yintercept = -2, linetype = "dashed") +
  theme_light() + 
  facet_wrap(~ sex, ncol = 3, labeller = sex_labeller) + 
#  ggtitle(space) + 
  theme_bw() + 
  theme(text=element_text(size=12,  family="Century Gothic"))

p

# Save as png
png_filename <- paste0(clade, "_patch_", space,  "_SES_", "nsims", n_sims, "_", avg_par, "_", metric, "_", iucn_type, "-iucn", ".png")
png(
  here::here(
    "2_Patches", "4_OutputPlots", "2_Diversity_measures", "1_Diversity_extinction_risk", 
    png_filename
  ),
  width = 1500, height = 2000/3, res = 150
)
p
dev.off()



