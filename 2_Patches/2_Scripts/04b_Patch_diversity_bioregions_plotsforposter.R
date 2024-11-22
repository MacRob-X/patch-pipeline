## Calculate diversity metrics averages for ecoregions, biomes
## Robert MacDonald
## 17th October 2024


# clear environment
rm(list=ls())

## Load libraries ----
library(ggplot2)
library(dplyr)
#library(terra)
library(tidyterra)
library(dispRity)

## EDITABLE CODE ## ----
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Passeriformes"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
# N.B. will need to change the date in the pca_all filename if using usmldbl or usmldblr
# (from 240603 to 240806)
space <- "lab"
# select sex (""allspecimens" or "matchedsex")
# "matchedsex" will subset to only species for which we have at least one
# male and female specimen
sex_match <- "matchedsex"
# select metric ("centr-dist", "nn-k", "nn-count")
# note that nn-k is EXTREMELY slow to run - needs parallelisation (but will probably still
# be too slow to run)
metric <- "centr-dist"
# select averaging type ("mean", "median", "mode")
avg <- "mean"
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
pam_res <- "50km"
# enter PAM files location
pams_filepath <- "X:/cooney_lab/Shared/Rob-MacDonald/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/PAMs"
#pams_filepath <- "X:/cooney_lab/Shared/Rob-MacDonald/SpatialData/BirdLife/BirdLife_Shapefiles_v9/PAMs/100km/Behrmann_cea/"
## END EDITABLE CODE ##

# set dispRity metric
if(metric == "centr-dist"){
  metric_get <- "centroids"
} else if (metric == "nn-k"){
  metric_get <- "mean.nn.dist"
}else if (metric == "nn-count"){
  metric_get <- "count.neighbours"
}

## Load data ----

# load patch data (PCA of whichever colourspace - generated in 02_Patch_Analyse_features.R)
pca_filename <- paste(clade, "patches.231030.PCAcolspaces", space, "240925", "rds", sep = ".")
pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)

# load PAM
pam_filename <- paste0("PAM_birds_Behrman", pam_res, "_Pres12_Orig12_Seas12_", pam_type, pam_seas, ".rds")
# or
#pam_filename <- paste0(pam_res, "_Behrmann_cea_PAM_combined.rds")
pam_raw <- readRDS(
  paste(pams_filepath, pam_filename, sep = "/")
)

# load Dinerstein et al 2017 ecoregion shape files
ecoregions <- sf::st_read(
  here::here(
    "4_SharedInputData", "Ecoregions2017_accessed2024-10-07",
    "Ecoregions2017.shp"
  )
)


## Analysis ----

# Set spatial parameters
null_rast <- pam_raw[[2]]
pam_allspec <- pam_raw[[1]]
colnames(pam_allspec) <- gsub(" ", "_", colnames(pam_allspec))
crs <- as.character(raster::crs(null_rast)) # Behrmann cylindrical equal area projection (standard parallels at 30 deg)

# convert ecoregions to Behrman CEA CRS
ecoregions <- sf::st_transform(ecoregions, crs = terra::crs(null_rast, proj = TRUE))


# Extract PCA co-ordinates data and add species and sex columns
pca_dat <- pca_all$x %>% 
  as.data.frame() %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  )# %>% 
# filter(
#   sex == sex
# )

# if specified, clip to only species for which we have at least one male and female specimen
if(sex_match == "matchedsex"){
  
  # get list of species to keep
  species_list <- pca_dat[, c("species", "sex")]
  species_m <- species_list[species_list$sex == "M", "species"]
  species_f <- species_list[species_list$sex == "F", "species"]
  spec_to_keep <- intersect(species_m, species_f)
  
  # subset data to these species
  pca_dat <- pca_dat[pca_dat$species %in% spec_to_keep, ]
  
  # remove variables
  rm(species_list, species_m, species_f, spec_to_keep)
  
  # set list of sexes (for looping later)
  sexes <- c("M", "F")
  
} else if(sex_match == "allspecimens"){
  sexes <- c("M", "F", "U")
}

# subset pam so it includes only species in patch data
spec_list <- unique(pca_dat$species)
spec_keep <- colnames(pam_allspec)[colnames(pam_allspec) %in% spec_list]
pam <- pam_allspec[, spec_keep]

#####################################################
#####################################################
##### NEED TO APPLY SPECIES RICHNESS FILTER #########
##### TO PAM IN SOME WAY                    #########
#####################################################
#####################################################
# actually do I? I'm not sure I do
# probably the right thing to do is to set values of mean
# diversity for ecoregions with <5 species equal to NA

# remove raw PAMs for RAM reasons
rm(pam_allspec, pam_raw)
gc()

# calculate distance to the centroid for species
ctrd_dists <- as.data.frame(
  dispRity(
    pca_dat %>% dplyr::select(-species, -sex) %>% as.matrix(), 
    metric = get(metric_get))$disparity[[1]][[1]]
) %>% 
  rename(
    centroid_distance = V1
  )
rownames(ctrd_dists) <- rownames(pca_dat)

# filter to individual sexes
ctrd_dists_M <- ctrd_dists[grepl("-M", rownames(ctrd_dists)), , drop = FALSE]
rownames(ctrd_dists_M) <- sapply(strsplit(rownames(ctrd_dists_M), split = "-"), "[", 1)
ctrd_dists_F <- ctrd_dists[grepl("-F", rownames(ctrd_dists)), , drop = FALSE]
rownames(ctrd_dists_F) <- sapply(strsplit(rownames(ctrd_dists_F), split = "-"), "[", 1)
ctrd_dists_U <- ctrd_dists[grepl("-U", rownames(ctrd_dists)), , drop = FALSE]
rownames(ctrd_dists_U) <- sapply(strsplit(rownames(ctrd_dists_U), split = "-"), "[", 1)

# convert to named vectors
M_names <- rownames(ctrd_dists_M)
ctrd_dists_M <- ctrd_dists_M$centroid_distance
names(ctrd_dists_M) <- M_names
F_names <- rownames(ctrd_dists_F)
ctrd_dists_F <- ctrd_dists_F$centroid_distance
names(ctrd_dists_F) <- F_names
U_names <- rownames(ctrd_dists_U)
ctrd_dists_U <- ctrd_dists_U$centroid_distance
names(ctrd_dists_U) <- U_names
rm(M_names, F_names, U_names)

# loop to calculate average diversity metric for all species present within each ecoregion

# # first add extra columns for male and female mean metrics
# tmp_df <- data.frame(rep(NA, times = nrow(ecoregions)), rep(NA, times = nrow(ecoregions)))
# colnames(tmp_df) <- c(
#   paste("mean", (sub("-", "_", metric)), "M", sep = "_"),
#   paste("mean", (sub("-", "_", metric)), "F", sep = "_")
# )
# ecoregions <- cbind(ecoregions, tmp_df)
# rm(tmp_df)

# get null terra raster
null_terrast <- terra::rast(null_rast)


# for(i in nrow(ecoregions)){
#   
#   # get ecoregion of interest
#   ecoreg_i <- ecoregions[i, ]
#   
# 
#   # create null rast, subsetted to extent of ecoregion
#   # first assign values to null_rast to keep track of grid cells
#   null_terrast_sub <- null_terrast
#   terra::values(null_terrast_sub) <- 1:length(terra::values(null_terrast_sub))
#   null_terrast_sub <- null_terrast_sub %>% 
#     crop(ecoreg_i) %>% 
#     mask(ecoreg_i)
#   
#   # use values of subsetted raster to subset PAM to ROI
#   pam_sub <- pam[terra::values(null_terrast_sub), ]
# 
#   # get species which are present in subsetted PAM
#   species_pres <- names(which(apply(pam_sub, 2, function(x) any(!(is.na(x))))))
#   
#   # get average distance to centroids of present species for males and females
#   ecoregions$avg_centr_dist_M <- get(avg)(ctrd_dists_M[species_pres])
#   ecoregions$avg_centr_dist_F <- get(avg)(ctrd_dists_F[species_pres])
#   
#   
# }

# function to calculate average value of metric for males and females for each ecoregion
# returns [1, 2] matrix that can be inputted to ecoregion data
calc_avg_metric <- function(ecoreg, metric_M, metric_F, null_terrast, ncells, avg){
  
  # for debugging purposes
  ecoreg_num <- ecoreg$OBJECTID
  ecoreg_name <- ecoreg$ECO_NAME
  
  # create null rast, subsetted to extent of ecoregion
  # first assign values to null_rast to keep track of grid cells
  null_terrast_sub <- null_terrast
  terra::values(null_terrast_sub) <- 1:ncells
  null_terrast_sub <- null_terrast_sub %>% 
    terra::crop(ecoreg) %>% 
    terra::mask(ecoreg)
  
  # use values of subsetted raster to subset PAM to ROI
  pam_sub <- pam[terra::values(null_terrast_sub), ]
  
  # check if 1 or fewer grid cells
  if(length(terra::values(null_terrast_sub)) < 2){
    # get species which are present in subsetted PAM
    species_pres <- names(which(!is.na(pam_sub)))
  } else {
    # get species which are present in subsetted PAM
    species_pres <- names(which(apply(pam_sub, 2, function(x) any(!(is.na(x))))))
  }
  
  # check if there are as many or more species in ROI than the threshold SR
  if(length(species_pres) >= sr_threshold){
    
    # if there are enough species, get average distance to centroids of 
    # present species for males and females
    avg_metric_M <- get(avg)(metric_M[species_pres], na.rm = TRUE)
    avg_metric_F <- get(avg)(metric_F[species_pres], na.rm = TRUE) 
  } else {
    
    # if there are fewer, set values = NA
    avg_metric_M <- NA
    avg_metric_F <- NA
  }
  
  # set species richness value
  sr <- length(species_pres)
  
  avgs <- c(avg_metric_M, avg_metric_F, sr)
  
  ## SET VALUES FOR ECOREGION OBJECTID 207 - "Rock and Ice" - BIOM_NUM 11 - to 
  ## NA - it seems, from looking at Dinerstein et al (2017) and their interactive map
  ## https://ecoregions.appspot.com/ that this should be possibly excluded as it's
  ## basically uninhabited - it has an artificially extremely high species richness
  ## because it spans all of Antarctica, Greenland, parts of the Himalaya, Iceland, 
  ## and Northwestern North America
  if(ecoreg$OBJECTID == 207){
    avgs <- rep(NA, times = 3)
  }
  
  return(avgs)
  
}


# apply function across rows of ecoregions
avg_vals <- lapply(1:nrow(ecoregions), function(i) {
  x <- ecoregions[i, ]
  calc_avg_metric(x, metric_M = ctrd_dists_M, metric_F = ctrd_dists_F, null_terrast = null_terrast, ncells = length(terra::values(null_terrast)), avg = avg)
    })
# convert output to 2-column matrix
avg_vals <- do.call(rbind, avg_vals)
colnames(avg_vals) <- c(
  paste(avg, (sub("-", "_", metric)), "M", sep = "_"),
  paste(avg, (sub("-", "_", metric)), "F", sep = "_"),
  "species_richness"
)

# bind to ecoregion data
ecoregions <- cbind(ecoregions, avg_vals)




# plot ecoregions coloured by species richness and by average metric
plots <- list()

# male plot
pm <- ggplot() + 
  geom_sf(data = ecoregions, aes(fill = get(paste0(avg, "_centr_dist_M")), colour = get(paste0(avg, "_centr_dist_M")))) +
  scale_fill_viridis_c() + 
  scale_colour_viridis_c() + 
  coord_sf(expand = FALSE) + 
  labs(fill = "Mean distance\nto centroid (males)", colour = "Mean distance\nto centroid (males)") + 
#  labs(fill = paste0(avg, "_centr_dist_M"), colour = paste0(avg, "_centr_dist_M")) + 
  theme_minimal() + 
  theme(text=element_text(size=12,  family="Century Gothic"), 
        axis.text.x=element_blank(),  # Remove x-axis tick labels
        axis.text.y=element_blank(),  # Remove y-axis tick labels
        axis.ticks=element_blank(),   # Remove ticks
        axis.title.x=element_blank(), # Remove x-axis title
        axis.title.y=element_blank())


plots[["male_metric"]] <- pm

# female plot
pf <- ggplot() + 
  geom_sf(data = ecoregions, aes(fill = get(paste0(avg, "_centr_dist_F")), colour = get(paste0(avg, "_centr_dist_F")))) +
  scale_fill_viridis_c() + 
  scale_colour_viridis_c() + 
  coord_sf(expand = FALSE) + 
  labs(fill = "Mean distance\nto centroid (females)", colour = "Mean distance\nto centroid (females)") +
#  labs(fill = paste0(avg, "_centr_dist_F"), colour = paste0(avg, "_centr_dist_F")) + 
  theme_minimal() + 
  theme(text=element_text(size=12,  family="Century Gothic"), 
        axis.text.x=element_blank(),  # Remove x-axis tick labels
        axis.text.y=element_blank(),  # Remove y-axis tick labels
        axis.ticks=element_blank(),   # Remove ticks
        axis.title.x=element_blank(), # Remove x-axis title
        axis.title.y=element_blank())

plots[["female_metric"]] <- pf

# sr plot
psr <- ggplot() + 
  geom_sf(data = ecoregions, aes(fill = species_richness, colour = species_richness)) +
  scale_fill_viridis_c() + 
  scale_colour_viridis_c() + 
  coord_sf(expand = FALSE) + 
  labs(fill = "Species richness", colour = "Species richness") +
#  labs(fill = "species_richness", colour = "species_richness") + 
  theme_minimal() + 
  theme(text=element_text(size=12,  family="Century Gothic"), 
        axis.text.x=element_blank(),  # Remove x-axis tick labels
        axis.text.y=element_blank(),  # Remove y-axis tick labels
        axis.ticks=element_blank(),   # Remove ticks
        axis.title.x=element_blank(), # Remove x-axis title
        axis.title.y=element_blank())

plots[["species_richness"]] <- psr

# save plot
png_filename <- paste(clade, "patches", "ecoregions", space, sex_match, avg, metric, "sr_thresh", sr_threshold, pam_res, "Behrman", pam_type, pam_seas, "CG", "png", sep = ".")
png(
  here::here(
    "2_Patches", "4_OutputPlots", "3_Spatial_mapping", pam_res, space,
    png_filename
  ), 
  width = 387, height = 400, units = "mm", pointsize = 30, res = 100
)

gridExtra::grid.arrange(grobs = plots, nrow = 3)

dev.off()


# inspect most species-rich ecoregions
ecoregions %>% 
  filter(species_richness > 0) %>% 
  arrange(desc(mean_centr_dist_F), ) %>% 
  as.data.frame() %>% 
  head(n = 20) 
