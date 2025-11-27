# SES of aesthetic attractiveness change with extinctions
# Using data from Santangeli et al 2023
# Robert MacDonald
# 2025-09-25


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

# select clade to limit analysis to
clade <- "Neognaths"
# select sex ("M", "F", "All")
sex <- "All"
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
# Filter to clade

if(clade == "Neognaths"){
  irate_master <- read.csv(
    here::here(
      "4_SharedInputData", "iratebirds_data", "Santangeli_et_al_data_and_scripts", "monsterALL2_2_2023.csv"
    )
  ) %>% 
    filter(
      order_ebird2021 != "Struthioniformes", 
      order_ebird2021 != "Rheiformes", 
      order_ebird2021 != "Casuariiformes", 
      order_ebird2021 != "Apterygiformes", 
      order_ebird2021 != "Tinamiformes", 
      #   order_ebird2021 != "Galliformes", 
      #   order_ebird2021 != "Anseriformes"
    )
} else if(clade == "Neoaves"){
  irate_master <- read.csv(
    here::here(
      "4_SharedInputData", "iratebirds_data", "Santangeli_et_al_data_and_scripts", "monsterALL2_2_2023.csv"
    )
  ) %>% 
    filter(
      order_ebird2021 != "Struthioniformes", 
      order_ebird2021 != "Rheiformes", 
      order_ebird2021 != "Casuariiformes", 
      order_ebird2021 != "Apterygiformes", 
      order_ebird2021 != "Tinamiformes", 
      order_ebird2021 != "Galliformes", 
      order_ebird2021 != "Anseriformes"
    )
} else if(clade == "Passeriformes"){
  irate_master <- read.csv(
    here::here(
      "4_SharedInputData", "iratebirds_data", "Santangeli_et_al_data_and_scripts", "monsterALL2_2_2023.csv"
    )
  ) %>% 
    filter(
    order_ebird2021 == "Passeriformes", 
  )
}


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
  rename(
    attract = attractiveness_mean,
    iucn_cat = redlistCategory,
    species = sciname_ebird2019
  ) %>% 
  mutate(
    species = gsub(" ", "_", species),
    attract = as.numeric(gsub(",", ".", attract)),
    sex = ifelse(sex == "male", "M", sex),
    sex = ifelse(sex == "female", "F", sex)
  ) %>% 
  select(species, sex, attract, iucn_cat) #%>% 
  # Filter out subspecies
  # filter(stringr::str_count(species, "_") == 1) %>% 
  # distinct()

# rename IUCN categories and filter out DD, EW, EX species and species for which we don't have a category
trait <- trait %>% 
  mutate(
    iucn_cat = ifelse(iucn_cat == "Critically Endangered", "CR", iucn_cat),
    iucn_cat = ifelse(iucn_cat == "Endangered", "EN", iucn_cat),
    iucn_cat = ifelse(iucn_cat == "Vulnerable", "VU", iucn_cat),
    iucn_cat = ifelse(iucn_cat == "Near Threatened", "NT", iucn_cat),
    iucn_cat = ifelse(iucn_cat == "Least Concern", "LC", iucn_cat),
  ) %>% 
  filter(
    iucn_cat %in% c("CR", "EN", "VU", "NT", "LC")
  )


# add extra (duplicate) rows for undefined sex (largely monomorphic) species
mon_trait <- trait %>% 
  filter(
    sex == "average"
  )
mon_trait_m <- mon_trait %>% 
  mutate(
    sex = "M"
  )
mon_trait_f <- mon_trait %>% 
  mutate(
    sex = "F"
  )
mon_trait <- mon_trait_m %>% 
  bind_rows(mon_trait_f)
trait <- trait %>% 
  filter(
    sex != "average"
  ) 
trait <- trait %>% 
  bind_rows(mon_trait)
rm(mon_trait, mon_trait_f, mon_trait_m)
# keep only spp for which we have both male and female data
trait <- trait %>% 
  group_by(species) %>% 
  filter(n() > 1)

# filter iucn data to only these species
# Species matched to IUCN categories:
species <- trait %>% 
  select(
    species, iucn_cat
  ) %>% 
  distinct(
    species, .keep_all = T
  )


# add extra rows in species for M/F
species$sex <- "M"
species_F <- species; species_F$sex <- "F"
species <- rbind(species, species_F)
# convert to dfs
species <- as.data.frame(species)
trait <- as.data.frame(trait)
# add rownames
rownames(trait) <- paste(trait$species, trait$sex, sep = "-")
rownames(species) <- paste(species$species, species$sex, sep = "-")
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
filename <- paste0(clade, "_santangiratebirds_", "SES_", "nsims", n_sims, "_", "attractiveness", "_", sex, "_", iucn_type, "-iucn", ".csv")
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
# create dataframe to populate
res <- data.frame(matrix(NA, nrow = 0, ncol = 9))
colnames(res) <- c("IUCN", "SR", "Raw", "NULL_MEAN", "NULL_SD", "NULL_SE", "SES", "Metric", "SR_prop")
for (sex in c("M", "F", "All")){
  # set file path
  filename <- paste0(clade, "_santangiratebirds_", "SES_", "nsims", n_sims, "_", "attractiveness", "_", sex, "_", iucn_type, "-iucn", ".csv")
  # load file and calculate proportional species richness
  temp_res <- read.csv(
    here::here(
      "2_Patches", "3_OutputData", clade, "4_Diversity_measures", filename
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
  facet_wrap(~ sex, ncol = 3, labeller = sex_labeller) + 
  #  ggtitle(space) + 
  theme_bw() + 
  theme(text=element_text(size=12,  family="Century Gothic"))

p

# Save as png
png_filename <- paste0(clade, "_santangiratebirds_attractiveness", "_SES_", "nsims", n_sims, "_", iucn_type, "-iucn", ".png")
png(
  here::here(
    "2_Patches", "4_OutputPlots", clade, "2_Diversity_measures", "1_Diversity_extinction_risk", 
    png_filename
  ),
  width = 1500, height = 2000/3, res = 150
)
p
dev.off()
