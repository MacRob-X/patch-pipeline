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
extract_null_rast <- function(pam_raw, rast_type = "terra"){
  
  # return terra-type raster
  if(rast_type == "terra"){
    return(terra::rast(pam_raw[[2]]))
  } else if(rast_type == "raster"){
    # return raster-type raster
    return(pam_raw[[2]])
  }
  
  
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
extract_sex_vals <- function(div_data, sex, metric = NULL, assemblage_level = FALSE){
  
  # generate sex suffix
  sex_suffix <- paste0("-", sex)
  
  # if working with diversity metric data, rather than calculating metrics at an assemblage level
  if(!is.null(metric) & assemblage_level == FALSE){
    
    # filter to individual sex
    metric_sex <- div_data[grepl(sex_suffix, rownames(div_data)), metric, drop = FALSE]
    rownames(metric_sex) <- sapply(strsplit(rownames(metric_sex), split = "-"), "[", 1)
    
    return(metric_sex)
    
  } else if(assemblage_level == TRUE){ # If metric is null, assume we're working with PCA data (i.e., keep all columns)
    
    # filter to individual sex  
    pca_sex <- div_data[grepl(sex_suffix, rownames(div_data)), , drop = FALSE]
    rownames(pca_sex) <- sapply(strsplit(rownames(pca_sex), split = "-"), "[", 1)
    
    return(pca_sex)
    }
  
}

# get list of species (with sex) included in diversity data
get_spp_sex <- function(div_data){
  
  # extract species and sex
  spp_sex <- t(sapply(strsplit(rownames(div_data), split = "-"), "[", 1:2))
  
  # set column names
  colnames(spp_sex) <- c("species", "sex")
  
  return(spp_sex)
  
}

# Get list of species for which we have matched male and female data
sex_match_fun <- function(data){
  
  # get list of all species
  species_list <- pca_dat[, c("species", "sex")]
  
  # get male species
  species_m <- species_list[species_list[, "sex"] == "M", "species"]
  # get female species
  species_f <- species_list[species_list[, "sex"] == "F", "species"]
  
  # identify species present in both male and female
  spec_to_keep <- intersect(species_m, species_f)
  
  return(spec_to_keep)
  
}

# Subset to only species we want to keep (e.g. species identifed by sex_match())
subset_sex_match <- function(data, spec_to_keep){
  
  data <- data[data[, "species"] %in% spec_to_keep, ]
  
  # remove species and sex columns
  data <- data[, !colnames(data) %in% c("species", "sex")]
  
  # set class to numeric matrix (for use with dispRity)
  class(data) <- "numeric"
  
  return(data)
  
}

# add species and sex columns
pca_spp_sex <- function(pca_vals){
  
  # extract species and sex
  spp_sex <- t(sapply(strsplit(rownames(pca_vals), split = "-"), "[", 1:2))
  
  # add species and sex columns
  pca_vals <- cbind(
    pca_vals,
    species = spp_sex[, 1],
    sex = spp_sex[, 2]
  )
  
  return(pca_vals)
  
}

# combine individual rasters stored in a list into a single multilayer raster
combine_raster_list <- function(raster_list){
  
  require(terra)
  
  # get number of rasters in list
  n_rasters <- length(raster_list)
  
  # extract elements and combine to single raster
  # (unfortunately I think the only way to do this is recursively)
  multi_rast <- raster_list[[1]]
  for(n in 2:n_rasters){
    # combine multilayer raster with new layer
    multi_rast <- c(multi_rast, raster_list[[n]])
  }
  
  # set layer names as element names from list
  if(!is.null(names(raster_list))){
    names(multi_rast) <- names(raster_list)
  }
  
  return(multi_rast)
  
}

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
  
  # calculate average (using wrapper function)
  avg_div_val <- avg(div_values, avg_par, na.rm = TRUE)
  
  # return average
  return(avg_div_val)
  
}

# create assemblage-level diversity raster - takes diversity values per grid cell as input
make_assdiv_raster <- function(div_vals, null_rast){
  
  # create copy of null rast
  div_raster <- null_rast
  
  # set raster values to species richness values
  terra::values(div_raster) <- div_vals
  
  return(div_raster)
  
}