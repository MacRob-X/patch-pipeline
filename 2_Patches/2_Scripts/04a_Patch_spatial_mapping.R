## Development version of 04a_Patch_spatial_mapping.R
## Robert MacDonald
## 2nd October 2024
## 
## Note that this script works with the ranges_10_species file, which contains range maps for 10 species randomly extracted from the 
## BOTW range maps 
## 
## Code used to produce the ranges_10_species file can be found at the bottom of this script
## This includes useful code for munging the all-species data

# clear environment
rm(list=ls())

## Load libraries ----
library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(dispRity)


## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Passeriformes"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
# N.B. will need to change the date in the pca_all filename if using usmldbl or usmldblr
# (from 240603 to 240806)
space <- "lab"
# select sex ("M", "F", "All")
sex <- "M"
# select metric ("centr-dist", "nn-k", "nn-count")
# note that nn-k is EXTREMELY slow to run - needs parallelisation (but will probably still
# be too slow to run)
metric <- "centr-dist"
# select whther to use liberal, conservative, or nominate IUCN data
# "liberal" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the highest threat level of the multiple species
# "conservative" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the lowest threat level of the multiple species
# "nominate" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the threat level of the BL species that corresponds to the nominate subspecies
iucn_type <- "nominate"
# enter range maps location
rangemaps_filepath <- "X:/cooney_lab/Shared/Rob-MacDonald/SpatialData/BirdLife_dist_maps_2022.2"
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

# get world map (for plotting and for CRS)
# we need to get the CRS from this because when I created the ranges_10_species file I used a CRS that doesn't seem to work
# with the rasterize function
world <- rnaturalearth::ne_countries(returnclass = "sf")

# inspect layers of sample 10-species data
st_layers(
  paste(rangemaps_filepath, "ranges_10", "ranges_10_species.shp", sep = "/")
)

# read in sample 10-species data 
ranges_10 <- sf::st_read(
  dsn = paste(rangemaps_filepath, "ranges_10", "ranges_10_species.shp", sep = "/")
) %>% 
  st_transform(crs = crs(world))


## Munging ----
## Note that prior to this step, taxonomy clashes must be resolved

# Extract PCA co-ordinates data and add species and sex columns
pca_dat <- pca_all$x %>% 
  as.data.frame() %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  )

# calculate distance to the centroid for species
ctrd_dists <- as.data.frame(
  dispRity(
    pca_dat %>% select(-species, -sex), 
    metric = get(metric_get))$disparity[[1]][[1]]
) %>% 
  rename(
    centroid_distance = V1
  )
rownames(ctrd_dists) <- rownames(pca_dat)


# let's just look at male specimens, for simplicity
# filter centroid data to males only
ctrd_dists <- ctrd_dists %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  ) %>% 
  filter(
    sex == "M"
  ) %>% 
  select(
    -sex
  )

# left join the centroid distances to the species range data
ranges_10_ctrds <- ranges_10 %>% 
  left_join(
    ctrd_dists,
    by = join_by(
      sci_nam == species
    )
  ) %>% 
  st_as_sf()

# check the output
ranges_10_ctrds

# Check if all range maps are valid and make valid if not
if(!all(st_is_valid(ranges_10_ctrds))) {
  ranges_10_ctrds <- ranges_10_ctrds %>% 
    st_make_valid()
}

## Visualise raw range maps ----
ggplot() + 
  geom_sf(data = world, colour = "steelblue") + 
  geom_sf(data = ranges_10, aes(colour = sci_nam, fill = sci_nam), alpha = 1/3) + 
  coord_sf(expand = FALSE)

# visualise centroid distances (in individual polygons)
ggplot() + 
  geom_sf(data = world, colour = "steelblue") + 
  geom_sf(data = ranges_10_ctrds, aes(colour = centroid_distance, fill = centroid_distance), alpha = 1/3) + 
  coord_sf(expand = FALSE)

## Creation of raster ----

# Create a template raster
template_raster <- rast(
  extent = ext(-180, 180, -90, 90),
  resolution = 0.5,
  crs = crs(ranges_10_ctrds)
)

# 2 methods to create the raster

## Method 1: directly getting mean of centroid distance ----
ctrds_rast_direct <- ranges_10_ctrds %>% 
  vect() %>% 
  rasterize(
    template_raster,
    field = "centroid_distance",
    fun = "mean"
  )

# inspect ouput
ctrds_rast_direct

# and plot
ggplot() + 
  geom_spatraster(data = ctrds_rast_direct) +
  scale_fill_viridis_c(option = "plasma") + 
  geom_sf(data = world, colour = "darkgrey", fill = NA) + 
  coord_sf(expand = FALSE)

## Method 2: create raster for each species map individually, then stack and calculate mean centroid distance ----

# get list of unique species in ranges data
species_10 <- unique(ranges_10$sci_nam)

# initialise empty list to hold individual species rasters
species_rasters <- vector("list", length = length(species_10))
names(species_rasters) <- species_10

# iterate over species
for(i in 1: length(species_10)) {
  # get species name
  species <- species_10[i]
  # get species map
  species_map <- vect(ranges_10_ctrds[ranges_10_ctrds$sci_nam == species, ])
  # create raster for species
  species_raster <- terra::rasterize(
    species_map,
    template_raster,
    field = "centroid_distance",
    fun = "mean",
    background = NA
  )
  # assign to list of rasters
  species_rasters[[i]] <- species_raster
  
  # print how far through the list of species the loop is
  cat("\r", i, " of ", length(species_10), " species rasterised")
}

# combine rasters into a single stack
species_rasters_stack <- rast(species_rasters)

# get the mean centroid distance for each pixel from this raster stack
ctrds_rast_looped <- app(species_rasters_stack, fun = mean, na.rm = TRUE)

# inspect output
ctrds_rast_looped

# plot
ggplot() + 
  geom_spatraster(data = ctrds_rast_looped) +
  scale_fill_viridis_c(option = "plasma") + 
  geom_sf(data = world, colour = "darkgrey", fill = NA) + 
  coord_sf(expand = FALSE)


#----------------------------------------------------------------------------------------------------------------------#
## Code to produce the species_10_ranges file ----

## Load all-species range map data

# Check layers in range map data
st_layers(
  here::here(
    "4_SharedInputData", "BirdLife_dist_maps_2022.2", "data.gdb"
  )
)

# Read in range map layer
birdlife_range_maps <- st_read(
  dsn = here::here(
    "4_SharedInputData", "BirdLife_dist_maps_2022.2", "data.gdb"
  ),
  layer = "All_Species"
)

## Data munging 
## Note that prior to this step, taxonomy clashes must be resolved

# Extract PCA co-ordinates data and add species and sex columns
pca_dat <- pca_all$x %>% 
  as.data.frame() %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  )

# Get species names separated by _ in range data
# Filter range data to  species’ breeding geographic ranges only (seasonal, 1 or 2) and regions where species are known 
# to be native or reintroduced (origin, 1 or 2) and extant or probably extant (presence, 1 or 2) - after
# Cooney et al (2022) - and trim to species for which we have colour data
all_species_ranges <- birdlife_range_maps %>% 
  mutate(
    sci_name = stringr::str_replace_all(sci_name, " ", "_")
  ) %>% 
  filter(
    seasonal == 1 | seasonal == 2,
    origin == 1 | origin == 2,
    presence == 1 | presence == 2,
  ) %>% 
  semi_join(
    pca_dat,
    by = join_by(sci_name == species)
  ) %>% 
  rename(
    geometry = Shape
  )

# Filter PCA data to remove species which don't have range maps
pca_dat <- pca_dat %>% 
  semi_join(
    all_species_ranges,
    by = join_by(species == sci_name)
  )

# Count remaining species (for which we have both colour data and range maps)
pca_dat %>% 
  summarise(
    species = n_distinct(species)
  )

# choose species
species_10 <- c("Accipiter_brevipes", "Larus_bulleri", "Rhynchopsitta_pachyrhyncha", "Chelictinia_riocourii",
                "Cyclopsitta_diophthalma", "Garrulus_glandarius", "Porphyrio_porphyrio", "Haematopus_finschi",
                "Falco_jugger", "Melospiza_lincolnii")

# get 10 species' ranges in WGS 84 CRS
ranges_10 <- all_species_ranges %>% 
  filter(
    sci_name == species_10[1] |
      sci_name == species_10[2] |
      sci_name == species_10[3] |
      sci_name == species_10[4] |
      sci_name == species_10[5] |
      sci_name == species_10[6] |
      sci_name == species_10[7] |
      sci_name == species_10[8] |
      sci_name == species_10[9] |
      sci_name == species_10[10]
  ) %>% 
  select(
    sisid, sci_name, presence, origin, seasonal, dist_comm, generalisd, yrmodified, version, 
    Shape_Length, Shape_Area, geometry
  ) %>% 
  st_as_sf() %>% 
  st_transform(crs = crs(world))

# inspect range data
ranges_10

# save as shape file
sf::st_write(
  ranges_10,
  here::here(
    "4_SharedInputData", "BirdLife_dist_maps_2022.2", "ranges_10", "ranges_10_species.shp"
  )
)

#----------------------------------------------------------------------------------------------------------------------#
## Apply rasterisation to all species in range map data ----
## Note that this will take an extremely long time to run

## Get directory path (on local machine for speed)
# set path
dir_botw <- "C:/Users/bop23rxm/Documents/BirdLife_dist_maps_2023_1/BOTW.gdb"

# Check layers in range map data
st_layers(
  dir_botw
)

# Read in range map layer
birdlife_range_maps <- st_read(
  dsn = dir_botw,
  layer = "All_Species"
)

# get world map (for plotting and for CRS)
world <- rnaturalearth::ne_countries(returnclass = "sf")

# Load PCA data
pca_all <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "Neoaves.patches.231030.PCAcolspaces.jndxyzlumr.240404.rds"
  )
)

## Munging ----
## Note that prior to this step, taxonomy clashes must be resolved

# Extract PCA co-ordinates data and add species and sex columns
pca_dat <- pca_all$x %>% 
  as.data.frame() %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  )

# calculate distance to the centroid for species
ctrd_dists <- as.data.frame(
  dispRity::dispRity(
    pca_dat %>% select(-species, -sex), 
    metric = dispRity::centroids)$disparity[[1]][[1]]
) %>% 
  rename(
    centroid_distance = V1
  )
rownames(ctrd_dists) <- rownames(pca_dat)

# Get species names separated by _ in range data
# Filter range data to  species’ breeding geographic ranges only (seasonal, 1 or 2) and regions where species are known 
# to be native or reintroduced (origin, 1 or 2) and extant or probably extant (presence, 1 or 2) - after
# Cooney et al (2022) - and trim to species for which we have colour data
all_species_ranges <- birdlife_range_maps %>% 
  mutate(
    sci_name = stringr::str_replace_all(sci_name, " ", "_")
  ) %>% 
  filter(
    seasonal == 1 | seasonal == 2,
    origin == 1 | origin == 2,
    presence == 1 | presence == 2,
  ) %>% 
  semi_join(
    pca_dat,
    by = join_by(sci_name == species)
  ) %>% 
  rename(
    geometry = Shape
  )

# filter centroid data to males only
ctrd_dists <- ctrd_dists %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  ) %>% 
  filter(
    sex == "M"
  ) %>% 
  select(
    -sex
  )

# left join the centroid distances to the species range data
all_species_ranges <- all_species_ranges %>% 
  left_join(
    ctrd_dists,
    by = join_by(
      sci_name == species
    )
  ) %>% 
  st_as_sf()

# check the output
all_species_ranges

# Check if all range maps are valid and make valid if not
if(!all(st_is_valid(all_species_ranges))) {
  all_species_ranges <- all_species_ranges %>% 
    st_make_valid()
}

# remove the original file from environment (to free up RAM)
rm(birdlife_range_maps)

## Convert MULTISURFACE-type geometries to MULTIPOLYGON ----
# Some species range maps are of multisurface-type geometry in the BOTW data
# This type of geometry cannot be converted to a SpatVector and so cannot be rasterised using terra::rasterize
# For this reason need to convert these geometries to multipolygons like the other geometries
# Need to use ogr2ogr facility in GDAL (install with conda forge)
# see https://stackoverflow.com/questions/74479973/casting-from-geomoetrycollection-with-multiple-nested-geometries-r-sf

# get only the MULTISURFACE range maps
multisurface_ranges <- all_species_ranges[!grepl("MULTISURFACE", st_geometry_type(all_species_ranges)), ]

# input and output locations
input_loc <- "C:/Users/bop23rxm/Documents/BirdLife_dist_maps_2023_1/all_species_ranges_with_ctrds_male_WGS_84.gpkg"
output_loc <- "C:/Users/bop23rxm/Documents/BirdLife_dist_maps_2023_1/all_species_ranges_with_ctrds_male_WGS_84_multisurf_recast.gpkg"

# save to .gpkg file
# see https://mapping-in-r-workshop.ryanpeek.org/02_import_export_gpkg#:~:text=Write%20TO%20geopackage,-Now%20we%20have&text=Just%20as%20we%20have%20been,dsn%3D%22data%2Fgpkg_in_R_example.

st_write(
  multisurface_ranges,
  dsn = output_loc,
  layer = "multisurface_ranges"
)

##--NOT RUN - need to do this via Miniforge Prompt command line (see below)------------------------#
# cast to polygon using ogr2ogr (via system command)
# ogr2ogr_loc <- "C:/Users/bop23rxm/AppData/Local/miniforge3/pkgs/libgdal-3.9.0-h4f813f3_3/Library/bin/ogr2ogr.exe"

# system(
#   paste0(
#     ogr2ogr_loc, " ",
#     output_loc, " ",
#     input_loc, " ",
#     "-explodecollections -nlt CONVERT_TO_LINEAR"
#   )
# )
##-------------------------------------------------------------------------------------#



## Miniforge Prompt command to recast multisurface geometries (substitute [input_loc] and [output_loc] for the actual 
## locations:
## C:\> ogr2ogr [output_loc] [input_loc] -explodecollections -nlt MULTIPOLYGON
## 
# note that this splits the complex multisurfaces up into many smaller multipolygons (like, thousands of them)
# so I might want to consider actually using "CONVERT_TO_LINEAR" in the command prompt instead of "MULTIPOLYGON" -
# I think this would simplify it a lot without losing info
# See https://gdal.org/programs/ogr2ogr.html


# inspect layers in newly created object
st_layers(
  output_loc
)

# read in newly created file as sf object (rename geom column 'geometry' to match original data (for some reason the ogr2ogr command changes it to 'geom'))
multisurf_ranges_recast <- st_read(
  dsn = output_loc,
  layer = "multisurface_ranges"
) %>% 
  rename(
    geometry = geom
  )

# inspect
multisurf_ranges_recast


# combine new recast ranges with original ranges (avoiding duplicates)
all_species_ranges_recast <- rbind(
  all_species_ranges[!grepl("MULTISURFACE", st_geometry_type(all_species_ranges)), ],
  multisurf_ranges_recast
)



# save this as new .gpkg file
st_write(
  all_species_ranges_recast,
  dsn = "C:/Users/bop23rxm/Documents/BirdLife_dist_maps_2023_1/all_species_ranges_with_ctrds_male_WGS_84_for_analysis.gpkg",
  layer = "multisurface_ranges"
)

# remove variables from environment
rm(multisurface_ranges, multisurf_ranges_recast, command, input_loc, output_loc)

## Creation of raster ----

# Create a template raster
template_raster <- rast(
  extent = ext(-180, 180, -90, 90),
  resolution = 0.5,
  crs = crs(all_species_ranges)
)

## create raster for each species map individually, then stack and calculate mean centroid distance

# get list of unique species in ranges data
species_list <- unique(all_species_ranges$sci_name)

# initialise empty list to hold individual species rasters
species_rasters <- vector("list", length = length(species_list))
names(species_rasters) <- species_list

# iterate over species
for(i in 5545 : length(species_list)) {
  # get species name
  species <- species_list[i]
  # get species map
  species_map <- vect(all_species_ranges_recast[all_species_ranges_recast$sci_name == species, ])
  # create raster for species
  species_raster <- terra::rasterize(
    species_map,
    template_raster,
    field = "centroid_distance",
    fun = "mean",
    background = NA
  )
  
  # save individual species raster (for sanity)
  savepath <- paste0("C:/Users/bop23rxm/Documents/birdlife_centroid_distance_rasters_individual/", 
                     species, "_M_centr_dist_raster.tif")
  terra::writeRaster(species_raster, savepath, overwrite = TRUE, memfrac = 0.85)
  
  # assign to list of rasters
  species_rasters[[species]] <- species_raster
  
  # print how far through the list of species the loop is
  cat("\r", i, " of ", length(species_list), " species rasterised")
}

# read in individual species rasters from files
# (we won't bother naming the elements as we're going to combine in a single stack)
species_rasters <- vector("list", length = length(species_list))
files <- list.files("C:/Users/bop23rxm/Documents/birdlife_centroid_distance_rasters_individual", full.names = TRUE)
for(i in 1:length(files)){
  species_rasters[[i]] <- terra::rast(
    files[i]
  )
  cat("\r", i)
}

# combine rasters into a single stack
species_rasters_stack <- rast(species_rasters)

# get the mean centroid distance for each pixel from this raster stack
mean_ctrds_rast <- app(species_rasters_stack, fun = mean, na.rm = TRUE)

mean_ctrds_rast <- mean(species_rasters_stack, na.rm = TRUE)

# save this final raster file
terra::writeRaster(mean_ctrds_rast, 
                   here::here(
                     "2_Patches", "3_OutputData", "6_Spatial_data", "jndxyzlumr_allspecies_M_mean_centr_dist_raster.tif"
                   ), 
                   overwrite = TRUE)

# Plot ----
ggplot() + 
  geom_spatraster(data = log(mean_ctrds_rast)) +
  scale_fill_viridis_c(option = "plasma") + 
  geom_sf(data = world, colour = "darkgrey", fill = NA) + 
  coord_sf(expand = FALSE)

# try clipping raster to extent of world - may make terrestrial trends easier to see
mean_ctrds_rast_terr <- terra::crop(mean_ctrds_rast, world, mask = TRUE)

# plot
ggplot() + 
  geom_spatraster(data = log(mean_ctrds_rast_terr)) +
  scale_fill_viridis_c(option = "plasma") + 
  geom_sf(data = world, colour = "darkgrey", fill = NA) + 
  coord_sf(expand = FALSE)
