## Spatial mapping using Presence-Absence Matrices (PAMs) provided by Chris Cooney & Gavin Thomas
## Robert MacDonald
## 2nd October 2024


# clear environment
rm(list=ls())

## Load libraries ----
library(ggplot2)
library(dplyr)
#library(terra)
library(tidyterra)
library(dispRity)

## EDITABLE CODE ##
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


## Analysis - Calculate mean of metric for each grid cell ----

# Set spatial parameters
null_rast <- pam_raw[[2]]
pam_allspec <- pam_raw[[1]]
colnames(pam_allspec) <- gsub(" ", "_", colnames(pam_allspec))
crs <- as.character(raster::crs(null_rast)) # Behrmann cylindrical equal area projection (standard parallels at 30 deg)


# Extract PCA co-ordinates data and add species and sex columns
pca_dat <- pca_all$x %>% 
  as.data.frame() %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
   )

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
# for interest, print species which are in our colour data but not in the range maps
spec_list[!(spec_list %in% colnames(pam_allspec))]

# clip PCA data to only species in range maps
#################################################################
################# IMPORTANT NOTE          #######################
## doing this actually seems to drastically change the output map
## particularly for India, Arabia, Madagascar and the southern
## tip of Africa
pca_dat <- pca_dat[pca_dat$species %in% colnames(pam_allspec), ]

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

# set species richness threshold
species_richness <- rowSums(pam, na.rm = T)
# create a mask for cells with species richness >= threshold
sr_mask <- ifelse(species_richness >= sr_threshold, TRUE, NA)
# save mask for later use
saveRDS(sr_mask,file = here::here(
  "2_Patches", "3_OutputData", "6_Spatial_mapping", pam_res, space,
  paste(clade, pam_type, "species_richness_mask", sr_threshold, sep = "_")))

# create raster of species richness
sr_raster <- terra::rast(null_rast)
species_richness[species_richness == 0] <- NA
terra::values(sr_raster) <- species_richness
# save sr raster for later use
sr_rast_filename <- paste("speciesrichness", clade, "patches", space, sex_match, pam_res, "Behrman", pam_type, pam_seas, "raster.tif", sep = ".")
terra::writeRaster(
  sr_raster,
  filename = here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", pam_res, space,
    sr_rast_filename
  )
)

# NOTE - I should functionalise this and apply it across a list of the species, rather than
# using a for loop - would be better in terms of RAM
for(sex in sexes){
  
  cat("/r", sex)
  
  # get sexed centroid distances
  ctrd_dists_tmp <- get(paste("ctrd_dists", sex, sep = "_"))
  
  # calculate mean distance to centroid per cell
  # N.B. this can't be parallelised because of RAM limitations
  mean_distance_per_cell <- apply(pam, 1, function(presence) {
    species_present <- which(presence == 1)
    if (length(species_present) > 0) {
      return(mean(ctrd_dists_tmp[species_present], na.rm = TRUE))
    } else {
      return(NA)  # If no species are present in the cell
    }
  })
  
  # create null raster to populate with values
  tmp_rast <- terra::rast(null_rast)
  
  # assign values to raster based on null raster
  terra::values(tmp_rast) <- mean_distance_per_cell
  
  # assign to sex-specific raster
  assign(
    paste("mean_ctrd_dist_rast", sex, sep = "_"),
    tmp_rast
  )
  
  # apply richness threshold to raster
  sr_lim_mean_ctrd_dist_rast <- terra::rast(tmp_rast)
  terra::values(sr_lim_mean_ctrd_dist_rast) <- terra::values(tmp_rast) * sr_mask
  
  # assign to sex-specific raster
  assign(
    paste("sr_lim_mean_ctrd_dist_rast", sex, sep = "_"),
    sr_lim_mean_ctrd_dist_rast
  )
  
}

# combine individual sexed rasters into multilayer raster with a layer for each sex
if(sex_match == "matchedsex"){
  mean_ctrd_dist_rast <- c(mean_ctrd_dist_rast_M, mean_ctrd_dist_rast_F)
  names(mean_ctrd_dist_rast) <-  c(paste0("M_", metric), paste0("F_", metric))
} else if(sex_match == "allspecimens"){
  mean_ctrd_dist_rast <- c(mean_ctrd_dist_rast_M, mean_ctrd_dist_rast_F, mean_ctrd_dist_rast_U)
  names(mean_ctrd_dist_rast) <-  c(paste0("M_", metric), paste0("F_", metric), paste0("U_", metric))
}

## NOT RUN - the below code is unnecessary and it changes the LENGTHUNIT parameter to metres
## even though it should be in km
# # explicitly set the CRS to Behrmann Equal Area (Cylindrical Equal Area with a 
# # specific latitude of true scale set at 30°)
# # Use Well Known Text (WKT) obtained from https://spatialreference.org/ref/esri/54017/
# terra::crs(mean_ctrd_dist_rast) <- 'PROJCRS["World_Behrmann",
#     BASEGEOGCRS["WGS 84",
#         DATUM["World Geodetic System 1984",
#             ELLIPSOID["WGS 84",6378137,298.257223563,
#                 LENGTHUNIT["metre",1]]],
#         PRIMEM["Greenwich",0,
#             ANGLEUNIT["Degree",0.0174532925199433]]],
#     CONVERSION["World_Behrmann",
#         METHOD["Lambert Cylindrical Equal Area",
#             ID["EPSG",9835]],
#         PARAMETER["Latitude of 1st standard parallel",30,
#             ANGLEUNIT["Degree",0.0174532925199433],
#             ID["EPSG",8823]],
#         PARAMETER["Longitude of natural origin",0,
#             ANGLEUNIT["Degree",0.0174532925199433],
#             ID["EPSG",8802]],
#         PARAMETER["False easting",0,
#             LENGTHUNIT["metre",1],
#             ID["EPSG",8806]],
#         PARAMETER["False northing",0,
#             LENGTHUNIT["metre",1],
#             ID["EPSG",8807]]],
#     CS[Cartesian,2],
#         AXIS["(E)",east,
#             ORDER[1],
#             LENGTHUNIT["metre",1]],
#         AXIS["(N)",north,
#             ORDER[2],
#             LENGTHUNIT["metre",1]],
#     USAGE[
#         SCOPE["Not known."],
#         AREA["World."],
#         BBOX[-90,-180,90,180]],
#     ID["ESRI",54017]]'

# save the raster (version without SR threshold applied)
rast_filename <- paste(clade, "patches", space, sex_match, metric, pam_res, "Behrman", pam_type, pam_seas, "raster.tif", sep = ".")
terra::writeRaster(
  mean_ctrd_dist_rast,
  filename = here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", pam_res, space,
    rast_filename
  )
)

# Plotting ----

# remove unnecessary variables
rm(list=setdiff(ls(), c("clade", "space", "sex_match", "metric", "pam_res", "pam_type", "pam_seas", "sr_threshold")))
gc()


# load raster 
rast_filename <- paste(clade, "patches", space, sex_match, metric, pam_res, "Behrman", pam_type, pam_seas, "raster.tif", sep = ".")
div_raster <- terra::rast(
  here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", pam_res, space,
    rast_filename
  )
)

# load species richness raster
sr_rast_filename <- paste("speciesrichness", clade, "patches", space, sex_match, pam_res, "Behrman", pam_type, pam_seas, "raster.tif", sep = ".")
sr_raster <- terra::rast(
  here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", pam_res, space,
    sr_rast_filename
  )
)

# get world map (for plotting)
world <- rnaturalearth::ne_countries(returnclass = "sf")
world <- sf::st_transform(world, crs = terra::crs(div_raster))

# apply species richness cutoff, if required
# load mask for cells with species richness >= threshold
sr_mask <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", pam_res, space,
    paste(clade, pam_type, "species_richness_mask", sr_threshold, sep = "_")
    )
)
# apply richness threshold to raster
sr_lim_div_raster <- terra::rast(div_raster)
terra::values(sr_lim_div_raster) <- terra::values(div_raster) * sr_mask

## GGPLOT VERSION

# generate quantiles for colour scale plotting purposes
# first examine the histograms
# males and females
hist(terra::values(sr_lim_div_raster), breaks = 50)
x_lim <- c(min(terra::values(sr_lim_div_raster), na.rm = T), max(terra::values(sr_lim_div_raster), na.rm = T))
# males only
hist(terra::values(sr_lim_div_raster[[1]]), breaks = 50, xlim = x_lim)
# females only
hist(terra::values(sr_lim_div_raster[[2]]), breaks = 50, xlim = x_lim)

# set up breaks as quantiles (specify number of quantiles)
nquants <- 10
quants <- quantile(terra::values(sr_lim_div_raster), na.rm = TRUE, probs = seq(0, 1, by = 1/nquants))
names(quants) <- round(quants, digits = 1)


## Plot with continuous or binned scale?
col_scale <- "continuous"

# plot the raster
p <- ggplot() + 
  geom_spatraster(data = sr_lim_div_raster) + 
  coord_sf(expand = FALSE) + 
  facet_wrap(~ lyr, ncol = 1)

# add colour scale (binned or continuous)
if(col_scale == "continuous"){
  p <- p + scale_fill_viridis_c(option = "turbo")
} else if(col_scale == "binned"){
  p <- p + scale_fill_binned(breaks = quants,
                    type = "viridis")
}

# show in plot viewer
p



# save image as png
png_filename <- paste(clade, "patches", space, sex_match, metric, "sr_thresh", sr_threshold, pam_res, col_scale, "colscale", "Behrman", pam_type, pam_seas, "png", sep = ".")
png(
  here::here(
    "2_Patches", "4_OutputPlots", "3_Spatial_mapping", pam_res, space,
   png_filename
  ), 
  width = 420, height = 420, units = "mm", pointsize = 24, res = 100
)
p
dev.off()


## Base R version (for binned with rainbow colour scale)

## EDITABLE CODE
## plot species richness?
plot_sr <- TRUE
## Choose font "default" or name font (e.g. "Century Gothic", "Avenir")
## Use extrafont::fonts() to see available fonts
font_par <- "default"
font_path <- "C:/Windows/Fonts/AvenirNextLTPro-Regular.ttf"
## Use viridis_c, viridis turbo palette or custom colour palette?
palette_choice <- "viridis"

# Add the custom font by specifying its path or name
if(font_par != "default"){
  sysfonts::font_add(family = font_par, regular = font_path)
}


# get spatvector version of world map for plotting purposes (with terra::plot)
world_vec <- terra::vect(world)
# create copy of div_raster with same extent as world_vec
# FOR PLOTTING PURPOSES ONLY
ext_match_div_raster <- sr_lim_div_raster
terra::ext(ext_match_div_raster) <- terra::ext(world_vec)
# same for species richness raster
ext_match_sr_raster <- sr_raster
terra::ext(ext_match_sr_raster) <- terra::ext(world_vec)

# get unified range
minval <- min(terra::values(sr_lim_div_raster), na.rm = TRUE)
maxval <- max(terra::values(sr_lim_div_raster), na.rm = TRUE)
nquants <- 12
quants <- quantile(terra::values(sr_lim_div_raster), na.rm = TRUE, probs = seq(0, 1, by = 1/nquants))
sr_quants <- quantile(terra::values(sr_raster), na.rm = TRUE, probs = seq(0, 1, by = 1/(nquants)))

if(palette_choice == "viridis"){
  colpal <- viridisLite::viridis(nquants)
} else if(palette_choice == "turbo"){
  colpal <- viridisLite::turbo(nquants)
} else if(palette_choice == "custom"){
  # create own colour palette (based on Cooney et al (2022) Fig. 2)
  cols <- c("#3e9eb5ff", "#eacc2cff", "#f82202ff")
  colpal <- colorRampPalette(cols)(nquants)
}
palname <- paste0("binned", palette_choice)


# plot
png_filename <- paste(clade, "patches", space, sex_match, metric, "sr_thresh", sr_threshold, pam_res, palname, "colscale", "Behrman", pam_type, pam_seas, "png", sep = ".")
if(plot_sr == TRUE){
    png(
    here::here(
      "2_Patches", "4_OutputPlots", "3_Spatial_mapping", pam_res, space,
      png_filename
    ), 
    width = 320, height = 400, units = "mm", pointsize = 24, res = 100,
   # family = font_par
   ## NOTE that I haven't sorted the custom font size out for this plot yet - currently
   ## using Avenir makes the font massive
  )
  par(mfrow = c(3, 1))
  cexval <- 1.1
  } else if(plot_sr == FALSE){
    png(
      here::here(
        "2_Patches", "4_OutputPlots", "3_Spatial_mapping", pam_res, space,
        png_filename
      ), 
      width = 420, height = 300, units = "mm", pointsize = 24, res = 100,
      family = font_par
    )
    par(mfrow = c(2, 1))
    # set font size
    cexval <- 0.7
    # use custom font
    showtext::showtext_begin()
  }

# first plot male diversity raster to get legend
terra::plot(
  ext_match_div_raster, "M_centr-dist", 
  breaks = quants, 
  col = colpal, 
 #  colNA = "grey",
  legend = TRUE,
  frame = FALSE,
  axes = FALSE, 
  plg = list(title = "Mean distance\nto centroid (males)", x = "left", cex = cexval),
  mar = c(1, 2, 1, 2)
)
# overlay grey world map so that NA parts of world are grey
terra::plot(
  world_vec, 
  col="grey80", 
  border = NA, 
  add = TRUE,
  axes = "n"
  )
# overlay male diversity raster again, but without legend
terra::plot(
  ext_match_div_raster, "M_centr-dist", 
  breaks = quants, 
  col = colpal, 
 # colNA = "grey",
  add = TRUE,
  legend = FALSE,
  plg = list(title = "Mean distance\nto centroid"),
  mar = c(1, 2, 1, 2)
  )
# plot female diversity raster (in separate plot)
# 
terra::plot(
  ext_match_div_raster, "F_centr-dist", 
  breaks = quants, 
  col = colpal, 
  # colNA = "grey",
  legend = TRUE,
  plg = list(title = "Mean distance\nto centroid (females)", x = "left", cex = cexval),
  mar = c(1, 2, 1, 2),
  axes = "n"
  # add = TRUE
)
# then overlay grey world a second time
terra::plot(
  world_vec, 
  col="grey80", 
  border = NA, 
  axes = "n",
  add = TRUE
)
# overlay female diversity raster again
# need to do it this way to get axes to match male map
terra::plot(
  ext_match_div_raster, "F_centr-dist", 
  breaks = quants, 
  col = colpal, 
 # colNA = "grey",
  legend = FALSE,
  mar = c(1, 2, 1, 2),
  add = TRUE
  )
# plot species richness (in a separate plot)
if(plot_sr == TRUE){
  terra::plot(
    ext_match_sr_raster,
    breaks = sr_quants,
    col = colpal,
    mar = c(1, 2, 1, 2),
    legend = TRUE,
    frame = FALSE,
    axes = FALSE,
    plg = list(title = "Species richness", x = "left", cex = cexval)
  )
  # overlay grey world
  terra::plot(
    world_vec, 
    col="grey80", 
    border = NA, 
    axes = "n",
    add = TRUE
  )
  # overlay species richness again
  terra::plot(
    ext_match_sr_raster,
    breaks = sr_quants,
    #  col = colpal,
    mar = c(1, 2, 1, 2),
    legend = FALSE,
    add = TRUE,
  )
  # end custom font
  if(font_par != "default"){
    showtext::showtext_end()
    }
}
dev.off() 

# other (non-working) colour scale options
ggplot() + 
  geom_spatraster(data = sr_lim_div_raster) + 
  ## this scale works but it doesn't bin the palette equally
  ## I need to find a way to define my own colour palette to use
  # scale_fill_whitebox_b(
  #   breaks = quants,
  #   labels = names(quants),
  #   palette = "viridi"
  # ) +
  ## This one doesn't work at all - it has no effect
  # binned_scale(
  #   aesthetics = "color",
  #   scale_name = "stepsn",
  #   palette = function(x) viridisLite::turbo(nquants),
  #   labels = names(quants),
  #   limits = c(min(terra::values(sr_lim_div_raster), na.rm = TRUE), max(terra::values(sr_lim_div_raster), na.rm = TRUE)),
  #   breaks = quants,
  #   show.limits = TRUE,
  #   guide = "colorsteps"
  # ) + 
  ## WORKING SCALE _ DISCRETE 
  ## This works but I can't work out how to pass a turbo colour palette to it
  ## Need to pass it to the 'type' argument I think, but it might not actually
  ## be possible
  scale_fill_binned(breaks = quants,
                    type = "viridis") + 
  #      type = getOption("ggplot2.binned.fill")) +
  #  geom_sf(data = world, colour = "steelblue", fill = NA) + 
  coord_sf(expand = FALSE) + 
  facet_wrap(~ lyr, ncol = 1)


# Downstream analysis ----

# read in raster
rast_filename <- paste(clade, "patches", space, sex_match, metric, pam_res, "Behrman", pam_type, pam_seas, "raster.tif", sep = ".")
div_raster <- terra::rast(
  here::here(
    "2_Patches", "3_OutputData", "6_Spatial_mapping", pam_res, space,
    rast_filename
  )
)

# get CRS
cea_crs <- div_raster %>% 
  raster::crs() %>% 
  as.character()

# read in Dinerstein et al 2017 ecoregion shape files
ecoregions <- sf::st_read(
  here::here(
    "4_SharedInputData", "Ecoregions2017_accessed2024-10-07",
    "Ecoregions2017.shp"
  )
)

# convert to Behrman CEA CRS
ecoregions <- sf::st_transform(ecoregions, crs = sf::st_crs(div_raster))

# plot to inspect
# ecoregions only
ggplot() + 
  geom_sf(data = ecoregions)
# frst 20 ecoregions overlaid on diversity raster (to check extents, crs etc match)
ggplot() + 
  geom_spatraster(data = div_raster) + 
  geom_sf(data = ecoregions[1:20, ]) + 
  facet_wrap(~ lyr, ncol = 1) + 
  coord_sf(expand = FALSE)

# add column to ecoregions object to populate with mean and median values of
# diversity metric
tmp_vect <- vector("numeric", length = nrow(ecoregions))
ecoregions <- cbind(ecoregions, tmp_vect)
colnames(ecoregions) <- gsub(
  pattern = "tmp_vect", 
  replacement = paste("mean", (sub("-", "_", metric)), sep = "_"), 
  x = colnames(ecoregions)
  )

# loop to subset raster to each ecoregion polygon
for(i in nrow(ecoregions)){
  
  # get eccoregion of interest
  ecoreg_i <- ecoregions[i, ]
  
  # crop diversity raster to extent of ecoregion polygon, then mask to outline of polygon
  div_rast_ecoreg <- div_raster %>% 
    terra::crop(ecoreg_i) %>% 
    terra::mask(ecoreg_i)
  
  #
  
}