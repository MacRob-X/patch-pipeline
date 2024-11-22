## Shared functions used in spatial mapping scripts (04x)
## Robert MacDonald
## 15th November 2024

# Define functions ----

# load Dinerstein et al 2017 ecoregion or load biome shape files derived from same
# will create biome data if no file already present
#' Title
#'
#' @param regions Defines whether to load ecoregion (Dinerstein et al 2017) or biome data
#'
#' @return an sf object containing the regions
#' @export
#'
#' @examples 
load_regions <- function(regions){
  
  if(regions == "ecoregions"){
    ecoregions <- sf::st_read(
      here::here(
        "4_SharedInputData", "Ecoregions2017_accessed2024-10-07",
        "Ecoregions2017.shp"
      )
    )
  } else if(regions == "biomes"){
    
    # check if biomes file is already present
    if(file.exists(here::here(
      "4_SharedInputData", "Ecoregions2017_accessed2024-10-07",
      "Biomes2017.shp"
    ))){
      
      # load file if present
      ecoregions <- sf::st_read(
        here::here(
          "4_SharedInputData", "Ecoregions2017_accessed2024-10-07",
          "Biomes2017.shp"
        )
      )
      
    } else{
      
      # create biomes file from ecoregions and save    
      print("No biomes shapefile present - creating from ecoregions file")
      
      # load ecoregions
      ecoregions <- sf::st_read(
        here::here(
          "4_SharedInputData", "Ecoregions2017_accessed2024-10-07",
          "Ecoregions2017.shp"
        )
      )
      
      # make geometries valid so we can combine them
      ecoregions <- ecoregions %>% 
        mutate(
          geometry = sf::st_make_valid(geometry)
        )
      # group by biome and then combine geometries within each biome
      # note that this takes about 15 minutes to run
      ecoregions <- ecoregions %>% 
        group_by(BIOME_NUM, BIOME_NAME) %>% 
        summarise(
          geometry = sf::st_union(geometry),
        ) %>% 
        ungroup()
      
      # save as shapefile
      ecoregions %>% 
        sf::st_write(
          here::here(
            "4_SharedInputData", "Ecoregions2017_accessed2024-10-07",
            "Biomes2017.shp"
          )
        )
      
    }
    
  }
}

# extract null raster from PAM file
extract_null_rast <- function(pam_raw){
  
  return(terra::rast(pam_raw[[2]]))
  
}

# extract PAM values
extract_pam_vals <- function(pam_raw){
  
  return(pam_raw[[1]])
  
}


# get list of unique species (of any sex) in data
get_unique_spp <- function(div_data){
  
  # extract species
  species_list <- unique(sapply(strsplit(rownames(div_data), split = "-"), "[", 1))
  
  return(species_list)
  
}

# subset PAM by a list of species (must be a vector)
subset_pam <- function(pam, spp_list){
  
  spp_keep <- colnames(pam)[colnames(pam) %in% spp_list]
  
  return(pam[, spp_keep])
  
}

# subset dataframe/matrix (with species/sex as rownames) by a list of species (vector)
subset_data_by_spp <- function(df, spp_list){
  
  # extract list of species (non-unique)
  df_spp <- sapply(strsplit(rownames(df), split = "-"), "[", 1)
  
  # subset
  return(df[df_spp %in% spp_list, ])
}

# set list of sexes to iterate over
set_sex_list <- function(sex_interest = c("all", "male_female", "male", "female", "unknown")){
  
  if(!(sex_interest %in% c("all", "male_female", "male", "female", "unknown"))){
    stop("Please enter a valid sex_interest ('all', 'male_female', 'male', 'female', 'unknown')")
  }
  
  if(sex_interest == "all"){
    sexes <- c("M", "F", "U")
  } else if(sex_interest == "male_female"){
    sexes <- c("M", "F")
  } else if(sex_interest == "male"){
    sexes <- "M"
  } else if(sex_interest == "female"){
    sexes <- "F"
  } else if(sex_interest == "unknown"){
    sexes <- "U"
  }
  
  return(sexes)
  
}

# extract named vectors of diversity metric for an individual sex
extract_sex_vals <- function(div_data, metric, sex){
  
  # generate sex suffix
  sex_suffix <- paste0("-", sex)
  
  # filter to individual sex
  metric_sex <- div_data[grepl(sex_suffix, rownames(div_data)), metric, drop = FALSE]
  rownames(metric_sex) <- sapply(strsplit(rownames(metric_sex), split = "-"), "[", 1)
  
  return(metric_sex)
  
}

# get list of species (with sex) included in diversity data
get_spp_sex <- function(div_data){
  
  # extract species and sex
  spp_sex <- t(sapply(strsplit(rownames(div_data), split = "-"), "[", 1:2))
  
  # set column names
  colnames(spp_sex) <- c("species", "sex")
  
  return(spp_sex)
  
}