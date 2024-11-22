

# attach libraries
library(dplyr)
library(dispRity)
library(ggplot2)

# clear environment
rm(list=ls())


## EDITABLE CODE ##
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr", "lab")
# N.B. will need to change the date in the pca_all filename if using usmldbl or usmldblr
# (from 240603 to 240806)
space <- "usmldbl"
# select sex ("M", "F", "All")
sex <- "All"
# select metric ("centr-dist", "nn-k")
metric <- "centr-dist"
# select number of null distributions to generate
n_sims <- 1000
## END EDITABLE CODE ##

# Load data ----

# load patch data (PCA of whichever colourspace - generated in 02_Patch_Analyse_features.R)

pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    paste0("Neoaves.patches.231030.PCAcolspaces.", space, ".240806.rds")
  )
)

# load IUCN Red List data
iucn <- readr::read_rds(
  here::here(
    "4_SharedInputData", "IUCN_RedList_data_130324.rds"
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
  mutate(
    species = gsub(" ", "_", scientific_name)
  ) %>% 
  select(
    species, category
  ) %>% 
  filter(
    species %in% trait[, "species"]
  )

# remove species not in IUCN data
trait <- trait[trait[, "species"] %in% species[, "species"], ]

# check species are identical in both matrices
spec_iucn <- sort(species[, "species"])
spec_trait <- sort(unique(trait[, "species"]))
identical(spec_iucn, spec_trait)
rm(spec_iucn, spec_trait)

# add extra rows in species for M/F/U
species$sex <- "M"
species_F <- species; species_F$sex <- "F"
species_U <- species; species_U$sex <- "U"
species <- rbind(species, species_F, species_U)
rm(species_F, species_U)
rownames(species) <- paste(species$species, species$sex, sep = "-")
species <- species[rownames(species) %in% rownames(trait), ]
# reorder to match trait data
species <- species[match(rownames(trait), rownames(species)), ]

# Subset to chosen sex
if(sex != "All"){
  trait <- trait[trait[, "sex"] == sex, ]
  species <- species[species[, "sex"] == sex, ]
}


# Use all PCS to account for 100% variation in traits
### All PCs
traits <- as.matrix(trait[,1:(ncol(trait) - 2)])
head(traits)

### OUTPUT DATA FRAME
res <- matrix(NA, nrow=5, ncol = 6)

a <- n_sims #change depending on number of simulations  
sims <- matrix(NA, nrow=1, ncol = a) #place to store simulations for SES


all <- rownames(species)
# extract traits for species in the random community
trait_vec <- matrix(traits[all,], ncol=dim(traits)[2])
# calculate SR & Mean distance to centroid
res[1,2] <- length(all)
res[1,3] <- mean(dispRity(trait_vec, metric = centroids)$disparity[[1]][[1]])

noCR <- rownames(species)[species$category != "CR"]
# extract traits for species in the random community
trait_vec <- matrix(traits[noCR,], ncol=dim(traits)[2])
# calculate SR & Mean distance to centroid
res[2,2] <- length(noCR)
res[2,3] <- mean(dispRity(trait_vec, metric= centroids)$disparity[[1]][[1]])

for (j in 1:ncol(sims)) {
  
  cat("\r", j, "of", ncol(sims))
  
  coms <- sample(x = all, size = length(noCR), replace = FALSE)
  trait_vec <- matrix(traits[coms,], ncol=dim(traits)[2])
  
  sims[1,j] <- mean(dispRity(trait_vec, metric= centroids)$disparity[[1]][[1]])
  
}

res[2,4] <- mean(sims)
res[2,5] <- sd(sims)
res[2,6] <- sd(sims)/sqrt(length(sims))

noEN <- rownames(species)[species$category != "CR" & species$category != "EN"]
# extract traits for species in the random community
trait_vec <- matrix(traits[noEN,], ncol=dim(traits)[2])
# calculate SR & Mean distance to centroid
res[3,2] <- length(noEN)
res[3,3] <- mean(dispRity(trait_vec, metric= centroids)$disparity[[1]][[1]])


for (j in 1:ncol(sims)) {
  
  cat("\r", j, "of", ncol(sims))
  
  coms <- sample(x = all, size = length(noEN), replace = FALSE)
  trait_vec <- matrix(traits[coms,], ncol=dim(traits)[2])
  
  sims[1,j] <- mean(dispRity(trait_vec, metric= centroids)$disparity[[1]][[1]])
  
}

res[3,4] <- mean(sims)
res[3,5] <- sd(sims)
res[3,6] <- sd(sims)/sqrt(length(sims))

noVU <- rownames(species)[species$category != "CR" & species$category != "EN" & species$category != "VU"]
# extract traits for species in the random community
trait_vec <- matrix(traits[noVU,], ncol=dim(traits)[2])
# calculate SR & Mean distance to centroid
res[4,2] <- length(noVU)
res[4,3] <- mean(dispRity(trait_vec, metric= centroids)$disparity[[1]][[1]])

for (j in 1:ncol(sims)) {
  
  cat("\r", j, "of", ncol(sims))
  
  coms <- sample(x = all, size = length(noVU), replace = FALSE)
  trait_vec <- matrix(traits[coms,], ncol=dim(traits)[2])
  
  sims[1,j] <- mean(dispRity(trait_vec, metric= centroids)$disparity[[1]][[1]])
  
}

res[4,4] <- mean(sims)
res[4,5] <- sd(sims)
res[4,6] <- sd(sims)/sqrt(length(sims))

noNT <- rownames(species)[species$category != "CR" & species$category != "EN" & species$category != "VU" & species$category != "NT"]
# extract traits for species in the random community
trait_vec <- matrix(traits[noNT,], ncol=dim(traits)[2])
# calculate SR & Mean distance to centroid
res[5,2] <- length(noNT)
res[5,3] <- mean(dispRity(trait_vec, metric= centroids)$disparity[[1]][[1]])

for (j in 1:ncol(sims)) {
  
  cat("\r", j, "of", ncol(sims))
  
  coms <- sample(x = all, size = length(noNT), replace = FALSE)
  trait_vec <- matrix(traits[coms,], ncol=dim(traits)[2])
  
  sims[1,j] <- mean(dispRity(trait_vec, metric= centroids)$disparity[[1]][[1]])
  
}

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
res$NULL_SE <- NA

res$SES <- (res$Raw - res$NULL_MEAN)/res$NULL_SD
str(res)
res$SES[1] <- 0

res$Metric <- "Trait Diversity"


# Save results as csv
# set filename
filename <- paste0("patch_", space,  "_SES_", "nsims", n_sims, "_", metric, "_", sex, ".csv")
write.csv(
  res, 
  here::here(
    "2_Patches", "3_OutputData", "4_Diversity_measures", filename
    ), 
  row.names = FALSE
  )


# Plot results ----

# clear environment (except chosen variables)
rm(list=setdiff(ls(), c("space", "metric", "n_sims")))

# load in results and get species richness as a proportion
# create dataframe to populate
res <- data.frame(matrix(NA, nrow = 0, ncol = 9))
colnames(res) <- c("IUCN", "SR", "Raw", "NULL_MEAN", "NULL_SD", "NULL_SE", "SES", "Metric", "SR_prop")
for (sex in c("M", "F", "All")){
  # set file path
  filename <- paste0("patch_", space,  "_SES_", "nsims", n_sims, "_", metric, "_", sex, ".csv")
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
res$IUCN <- as.factor(res$IUCN)
levels(res$IUCN) <- c("All Species Retained", "CR Lost", "EN Lost", "VU Lost", "NT Lost (LC Only Retained)")

# plot (subset by sex if required)
res %>% 
  # filter(
  #   sex != "All"
  # ) %>%
  ggplot(
    aes(x = sr_prop, y = SES, colour = IUCN)
    ) + 
  geom_point() + 
  scale_x_reverse() + 
  xlab("Remaining species richness") + ylab("Standard Effect Size") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = -2, linetype = "dashed") +
  theme_light() + 
  facet_wrap(~ sex, ncol = 1) + 
  ggtitle(space)

