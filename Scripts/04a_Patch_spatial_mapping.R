## Colour space and spatial mapping ----
## 1st May 2024
## Robert MacDonald

# clear environment
rm(list=ls())

## Load libraries ----
library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)

## Load data ----

# Load PCA data
pca_all <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "Neoaves.patches.231030.PCAcolspaces.jndxyzlumr.240404.rds"
  )
)

# Load taxonomic data
taxo <- readr::read_csv(
  here::here(
    "4_SharedInputData", "BLIOCPhyloMasterTax_2019_10_28.csv"
  )
)

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

# Load raster layer with correct projection, extent and resolution to use as template
# from Hughes et al (2022)
rast_template <- readr::read_rds(
  here::here(
    "4_SharedInputData", "Hughes_etal_2022", "Trimmed_50km_Behrmann.rds"
  )
)

## Data munging ----
## 
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


############################# DEV ###############################

## Extract 10 random species range maps to play with
## I'LL SIMPLIFY TO MAKE IT EASIER/FASTER TO MANIPULATE - BUT I PROBABLY WOULDN'T WANT TO DO THIS FOR THE REAL DATA
## Let's also get rid of some of the fields I'm not interested in and transform to the Behrmann equal area
## cylindrical projection

# choose species
species_10 <- c("Accipiter_brevipes", "Larus_bulleri", "Rhynchopsitta_pachyrhyncha", "Chelictinia_riocourii",
             "Cyclopsitta_diophthalma", "Garrulus_glandarius", "Porphyrio_porphyrio", "Haematopus_finschi",
             "Falco_jugger", "Melospiza_lincolnii")
# set CRS
behrmann_cea <- "+proj=cea"

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
  st_transform(crs = behrmann_cea)
ranges_10

# save as shape file
sf::st_write(
  ranges_10,
  here::here(
    "4_SharedInputData", "BirdLife_dist_maps_2022.2", "ranges_10", "ranges_10_species.shp"
    )
  )

# read in again (if necessary)
ranges_10 <- sf::st_read(
  dsn =here::here(
    "4_SharedInputData", "BirdLife_dist_maps_2022.2", "ranges_10", "ranges_10_species.shp"
    )
)

# Check if all range maps are valid and make valid if not
if(!all(st_is_valid(ranges_10))) {
  ranges_10 <- ranges_10 %>% 
    st_make_valid()
}

# Plot on world map
world <- rnaturalearth::ne_countries(returnclass = "sf") %>% 
  st_transform(crs = crs(ranges_10))

ggplot() + 
  geom_sf(data = world, colour = "steelblue") + 
  geom_sf(data = ranges_10, aes(colour = sci_name, fill = sci_name), alpha = 1/3) + 
  coord_sf(expand = FALSE)


# Extract polygon range maps ontp equal-area grid (Behrmann projection) at 0.5° resolution (~50km at the equator)
# the crs = "+proj=cea" specifies a Behrmann cylindrical equal area projection
# should I ise the extent of the ranges as the extent? or use teh bounding box of the range maps as the extent?
# or just set the extent to be the whole world?
grid <- rast(
  # whole world extent
  ext(-180, 180, -90, 90),
  # bounding box extent
  # ext(
  #   st_bbox(ranges_10)
  # ),
  # ranges extent - note that the rasterize() fn later won't work if this extent is used
  #ext(ranges_10),
  res = 0.5,  # 0.5 degree grid cell, I think (approx 55km^2 squares)
  crs = behrmann_cea
)

grid %>% glimpse()

# # get species range maps crs
# ranges_10 %>% 
#   st_crs() %>% 
#   magrittr::extract2(1)
# 
# # convert species range maps polygons to Behrmann projection
# ranges_10 <- ranges_10 %>% 
#   st_transform(
#     crs = st_crs(grid)
#   )
# 
# # also convert world map crs, for plotting purposes
# world <- world %>% 
#   st_transform(
#     crs = st_crs(grid)
#   )



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
ranges_10_ctrds %>% glimpse()

# Check if all range maps are valid and make valid if not
if(!all(st_is_valid(ranges_10_ctrds))) {
  ranges_10_ctrds <- ranges_10_ctrds %>% 
    st_make_valid()
}

# let's just plot this for sanity
ggplot() + 
  geom_sf(data = world, colour = "steelblue", alpha = 2/3) +
  geom_sf(data = ranges_10_ctrds, aes(colour = centroid_distance, fill = centroid_distance)) + 
  coord_sf(expand = FALSE)


# now rasterise to get a raster with mean centroid distance for each grid cell
ctrds_rast <- ranges_10_ctrds %>% 
  vect() %>% 
  rasterize(
    grid,
    field = "centroid_distance",
    fun = "mean"
  )


# check the output
ctrds_rast %>% glimpse()


# plot
ggplot() + 
  geom_sf(data = world, colour = "steelblue") + 
  geom_spatraster(data = ctrds_rast, aes(fill = centroid_distance)) + 
  scale_fill_viridis_c(option = "plasma") + 
  coord_sf(expand = FALSE)


#--------------------------------------------------------------------------------------------------------------------#
## Attempt 2
## based on "00_script_Produce_Renamed_BL_Shapefiles_AngelaChira.R

# Load data and range maps as at top of page

# Load birdlife taxonomy data
birdlife_taxo_checklist <- st_read(
  dsn = here::here(
    "4_SharedInputData", "BirdLife_dist_maps_2022.2", "data.gdb"
  ),
  layer = "Taxonomic_checklist"
)

# try creating shape file for test species (Abeillia abeillei)

# get row with data/polygon
abeillia_abeillei <- birdlife_range_maps %>% 
  filter(
    sci_name == "Abeillia abeillei"
  )

# write as new shapefile
abeillia_abeillei %>% 
  st_write(
    dsn = here::here(
      "4_SharedInputData", "BirdLife_dist_maps_2022.2", "data.gdb"
    ),
    layer = "sp_NAME_id_i",
    driver = "ESRI Shapefile"
  )

# try reading in shapefile
abeillia_abeillei_reread <- st_read(
  here::here(
  "4_SharedInputData", "BirdLife_dist_maps_2022.2", "data.gdb"
    ), 
  layer = "sp_NAME_id_i"
  )

# plot for sanity
world <- rnaturalearth::ne_countries(returnclass = "sf") %>% 
  st_transform(
    crs = crs(abeillia_abeillei_reread)
  )

ggplot() + 
  geom_sf(data = world, colour = "steelblue") + 
  geom_sf(data = abeillia_abeillei) + 
  coord_sf(expand = FALSE)
# looks good

## now create individual shapefiles for each species
list_spp <- as.character(birdlife_range_maps$sci_name)
list_id <- birdlife_range_maps$sisid
dir <- "C:/Users/bop23rxm/Documents/birdlife_allspecies_shapefiles"

# loop over all species
for(i in 1:length(list_spp)){
  # set filename
  name <- paste(gsub(" ", "_", list_spp[i]), as.character(list_id[i]), i, sep = "_")
  # write shapefile
  st_write(
    birdlife_range_maps[birdlife_range_maps$sci_name %in% list_spp[i], ],
    dsn = dir,
    layer = name,
    driver = "ESRI Shapefile",
    quiet = TRUE
  )
  cat("\r", i, "of ", length(list_spp))
}
# i = 12702



## Attempt 3 -----
## Based on ChatGPT input

# Load data (from BOTW 2023.1)
# Load from local machine storage for speed

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

# create spatvector directly on disc
range_maps_spatvec <- vect(
  dir_botw,
  layer = "All_Species"
)

# get world map (for plotting and for CRS)
world <- rnaturalearth::ne_countries(returnclass = "sf")

# read in sample 10-species data 
st_layers(
  here::here(
    "4_SharedInputData", "BirdLife_dist_maps_2022.2", "ranges_10", "ranges_10_species.shp"
  )
)

ranges_10 <- sf::st_read(
  dsn =here::here(
    "4_SharedInputData", "BirdLife_dist_maps_2022.2", "ranges_10", "ranges_10_species.shp"
  )
) %>% 
  st_transform(crs = crs(world))
  

# as spatvector
ranges_10_spatvec <- vect(
  here::here(
    "4_SharedInputData", "BirdLife_dist_maps_2022.2", "ranges_10", "ranges_10_species.shp"
  ),
  layer = "ranges_10_species"
)
# set CRS
crs(ranges_10_spatvec) <- crs(world)

## Data munging ----
## 
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


# choose 10 species to play with
species_10 <- c("Accipiter_brevipes", "Larus_bulleri", "Rhynchopsitta_pachyrhyncha", "Chelictinia_riocourii",
                "Cyclopsitta_diophthalma", "Garrulus_glandarius", "Porphyrio_porphyrio", "Haematopus_finschi",
                "Falco_jugger", "Melospiza_lincolnii")
# set CRS
behrmann_cea <- "+proj=cea"

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
  st_transform(crs = behrmann_cea)
ranges_10

# Plot on world map for sanity
world <- rnaturalearth::ne_countries(returnclass = "sf") %>% 
  st_transform(crs = crs(ranges_10))

ggplot() + 
  geom_sf(data = world, colour = "steelblue") + 
  geom_sf(data = ranges_10, aes(colour = sci_nam, fill = sci_nam), alpha = 1/3) + 
  coord_sf(expand = FALSE)


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

# plot for sanity
ggplot() + 
  geom_sf(data = world, colour = "steelblue") + 
  geom_sf(data = ranges_10_ctrds, aes(colour = centroid_distance, fill = centroid_distance), alpha = 1/3) + 
  coord_sf(expand = FALSE)

## Create a template raster
template_raster <- rast(
  extent = ext(-180, 180, -90, 90),
  resolution = 0.5,
#  crs = behrmann_cea
  crs = crs(ranges_10_ctrds)
)

# plot raster template
ggplot() + 
  geom_spatraster(data = template_raster)

## Rasterise the range maps
# initialise empty list to hold individual species rasters
species_rasters <- vector("list", length = length(species_10))
# note that I should really initialise this with the correct number of elements to improve RAM usage when looping

# iterate over species
for(species in species_10) {
  species_map <- vect(ranges_10_ctrds[ranges_10_ctrds$sci_nam == species, ])
  species_raster <- terra::rasterize(
    species_map,
    template_raster,
    field = "centroid_distance",
    fun = "mean",
    background = NA
  )
  species_rasters[[species]] <- species_raster
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

## or to get a mean raster directly in one go
ctrds_rast <- ranges_10_ctrds %>% 
  vect() %>% 
  rasterize(
    template_raster,
    field = "centroid_distance",
    fun = "mean"
  )

# inspect ouput
ctrds_rast

# and plot
ggplot() + 
  geom_spatraster(data = ctrds_rast) +
  scale_fill_viridis_c(option = "plasma") + 
  geom_sf(data = world, colour = "darkgrey", fill = NA) + 
  coord_sf(expand = FALSE)