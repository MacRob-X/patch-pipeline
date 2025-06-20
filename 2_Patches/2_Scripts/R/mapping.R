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
subset_sex_match <- function(data, spec_to_keep, convert_to_numeric = TRUE){
  
  data <- data[data[, "species"] %in% spec_to_keep, ]
  
  # remove species and sex columns
  data <- data[, !colnames(data) %in% c("species", "sex")]
  
  # set class to numeric matrix (for use with dispRity)
  if(convert_to_numeric == TRUE){class(data) <- "numeric"}
  
  return(data)
  
}

# add species and sex columns
pca_spp_sex <- function(pca_vals, bind_to_data = TRUE){
  
  # extract species and sex (old, slower version)
 # spp_sex <- t(sapply(strsplit(rownames(pca_vals), split = "-"), "[", 1:2))
  # more efficient version that avoids transposition
  spp_sex <- do.call(rbind, lapply(strsplit(rownames(pca_vals), split = "-"), "[", 1:2))
  
  if(bind_to_data == TRUE){
    
    # add species and sex columns
    pca_vals <- cbind(
      pca_vals,
      species = spp_sex[, 1],
      sex = spp_sex[, 2]
    )
    
    # return data bound to species and sex columns
    return(pca_vals)
    
  } else if(bind_to_data == FALSE){
    
    colnames(spp_sex) <- c("species", "sex")
    
    # only return species and sex columns
    return(spp_sex)
    
  }
  
  
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

# Set dispRity metric
set_metric <- function(metric){
  
  # set dispRity metric
  if(metric == "centr-dist"){
    metric_get <- get("centroids")
  } else if (metric == "nn-k"){
    metric_get <- get("mean.nn.dist")
  }else if (metric == "nn-count"){
    metric_get <- get("count.neighbours")
  } else if(metric == "nn-all"){
    metric_get <- get("neighbours")
  } else if(metric == "convhull.volume"){
    metric_get <- get("convhull.volume")
  } else if(metric == "displacements"){
    metric_get <- get("displacements")
  } else if(metric == "sum.variances"){
    metric_get <- list(sum, variances)
  } else if(metric == "sum.ranges"){
    metric_get <- list(sum, ranges)
  } 
  
  return(metric_get)
  
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
  
  # if calculating convex hull hypervolume, need at least n+1 species present, where n =
  # number of dimensions
  if(metric == "convhull.volume"){
    if(length(species_pres) <= ncol(pca_data)){
      return(NA)
    }
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
  # note that for metrics using more than one function (e.g. sum of variances), we just test the
  # first function in the list as this provides the level
  if(class(metric_get) == "list"){
    metric_type <- dispRity::make.metric(metric_get[[1]], silent = TRUE)[["type"]]
  } else if(class(metric_get) == "function"){
    metric_type <- dispRity::make.metric(metric_get, silent = TRUE)[["type"]]
  }
  
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

# create assemblage-level diversity raster - takes diversity values per grid cell as input
make_assdiv_raster <- function(div_vals, null_rast){
  
  # create copy of null rast
  div_raster <- null_rast
  
  # set raster values to species richness values
  terra::values(div_raster) <- div_vals
  
  return(div_raster)
  
}


# calculate species richness for each grid cell
calc_species_richness <- function(pam_matrix){
  
  return(rowSums(pam_matrix, na.rm = TRUE))
  
}


# create species richness mask that can be applied to PAM
make_sr_mask <- function(sr_vector, sr_threshold){
  
  # create mask
  sr_mask <- ifelse(sr_vector >= sr_threshold, TRUE, NA)
  
  return(sr_mask)
  
}

# create species richness raster
make_sr_raster <- function(sr_vals, null_rast){
  
  # set values == 0 to NA (so they don't get plotted)
  sr_vals[sr_vals == 0] <- NA
  
  # create copy of null rast
  sr_raster <- null_rast
  
  # set raster values to species richness values
  terra::values(sr_raster) <- sr_vals
  
  return(sr_raster)
  
  
}


apply_sr_threshold <- function(diversity_raster, sr_mask, exclude_sr_lyr = TRUE){
  
  # exclude species richness values from masking, if requested
  if(exclude_sr_lyr == TRUE){
    
    # set aside species richness values
    sr_values <- terra::values(diversity_raster[["species_richness"]])
    
    # apply sr threshold mask to whole raster
    div_rast_masked <- diversity_raster
    terra::values(div_rast_masked) <- terra::values(diversity_raster) * sr_mask
    
    # reset sr values to their original
    terra::values(div_rast_masked[["species_richness"]]) <- sr_values
    
  } else {
    
    # apply sr threshold mask to whole raster
    div_rast_masked <- diversity_raster
    terra::values(div_rast_masked) <- terra::values(diversity_raster) * sr_mask
    
  }
  
  return(div_rast_masked)
  
  
}

# get world map for plotting
get_world <- function(new_crs, spatvec){
  
  # get world map from rnaturalearth
  world <- rnaturalearth::ne_countries(returnclass = "sf")
  
  # transform to correct crs
  world <- sf::st_transform(world, crs = new_crs)
  
  # transform to spatvector, if required (need to do this if using base R plot)
  if(spatvec == TRUE){
    
    # transform to spatvector
    world <- terra::vect(world)
    
  }
  
  return(world)
}

make_colour_breaks <- function(raster, nquants, plot_sr = FALSE){
  
  # calculate species richness quantiles, if required
  if(plot_sr == TRUE){
    sr_quants <- quantile(
      terra::values(raster[["species_richness"]]),
      na.rm = TRUE,
      probs = seq(0, 1, by = 1/nquants)
    )
    names(sr_quants) <- round(sr_quants, digits = 1)
  }
  
  # remove the species_richness layer
  layer_names <- terra::names(raster)
  sr_index <- which(layer_names == "species_richness")
  new_layer_names <- layer_names[-sr_index]
  raster <- raster[[new_layer_names]]
  
  # get quantiles from raster diversity values
  div_quants <- quantile(
    terra::values(raster), 
    na.rm = TRUE, 
    probs = seq(0, 1, by = 1/nquants)
  )
  
  # make names of quants a rounded version of the quantiles (for plotting)
  names(div_quants) <- round(div_quants, digits = 1)
  
  if(plot_sr == TRUE){
    quants <- list(div_breaks = div_quants, sr_breaks = sr_quants)
  } else {
    quants <- div_quants
  }
  
  return(quants)
  
}

# match extents of world and diversity rasters (necessary to plot in base R)
match_extents <- function(div_rast, world_spatvec){
  
  # match extent of diversity raster to that of world spatvector
  terra::ext(div_rast) <- terra::ext(world_spatvec)
  
  return(div_rast)
  
}

# plot individual layers of diversity raster on the same plot
plot_div_raster <- function(div_rast, 
                            world_spatvec,
                            div_metric,
                            scale_type = c("binned", "continuous"), 
                            col_breaks = NULL,
                            pal_choice = c("viridis", "turbo", "custom"), 
                            plot_sr = FALSE,
                            save_png = FALSE,
                            png_filename = NULL,
                            assemb_level = FALSE){
  
  # get names of raster layers
  layer_names <- terra::names(div_rast)
  
  # set up colour scale and palette
  if(plot_sr == TRUE){
    div_breaks <- col_breaks$div_breaks
    sr_breaks <- col_breaks$sr_breaks
  } else {
    div_breaks <- col_breaks
  }
  if(scale_type == "binned"){
    
    nquants <- length(div_breaks) - 1
    
    # set colour palette for diversity metric
    if(pal_choice == "viridis"){
      divpal <- viridisLite::viridis(nquants)
    } else if(pal_choice == "turbo"){
      divpal <- viridisLite::turbo(nquants)
    } else if(palette_choice == "inferno"){
      divpal <- viridisLite::inferno(nquants)
    } else if(pal_choice == "custom"){
      # create own colour palette (based on Cooney et al (2022) Fig. 2)
      cols <- c("#3e9eb5ff", "#eacc2cff", "#f82202ff")
      divpal <- colorRampPalette(cols)(nquants)
    }
    
    # reverse colour palette if using nn-count (as low numbers indicate more unusual colours for this metric)
    if(metric == "nn-count"){
      divpal <- rev(divpal)
    }
    
    # set colour palette for species richness, if plotting
    if(plot_sr == TRUE){
      
      nquants_sr <- length(sr_breaks) - 1
      
      if(palette_choice == "viridis"){
        srpal <- viridisLite::viridis(nquants_sr)
      } else if(palette_choice == "turbo"){
        srpal <- viridisLite::turbo(nquants_sr)
      } else if(palette_choice == "inferno"){
        srpal <- viridisLite::inferno(nquants_sr)
      } else if(palette_choice == "custom"){
        # create own colour palette (based on Cooney et al (2022) Fig. 2)
        cols <- c("#3e9eb5ff", "#eacc2cff", "#f82202ff")
        srpal <- colorRampPalette(cols)(nquants_sr)
      }
    }
    
    
  }
  
  # Set up basic legend title
  if(metric == "centr-dist"){
    lgd_metric <- paste0(avg_par, " distance\nto centroid ")
  } else if(metric == "nn-k"){
    lgd_metric <- paste0(avg_par, " distance\nto nearest neighbours ")
  } else if(metric == "nn-count"){
    lgd_metric <- paste0(avg_par, " number\nof nearest neighbours\nin radius r ")
  } else if(metric == "convhull.volume"){
    lgd_metric <- paste0("Convex hull\nhypervolume ")
  } else if(metric == "sum.variances"){
    lgd_metric <- paste0("Sum of variances")
  } else {
    lgd_metric <- "Diversity"
  }
  
  
  # set up png saving parameters
  if(plot_sr == FALSE){
    png_width <- 420
    png_height <- 300
  } else if(plot_sr == TRUE){
    png_width <- 320
    png_height <- 400
  }
  
  # initialise png saving
  if(save_png == TRUE){
    if(assemb_level == FALSE){
      png(
        here::here(
          "2_Patches", "4_OutputPlots", "3_Spatial_mapping", pam_res, space, pam_type,
          png_filename
        ), 
        width = png_width, height = png_height, units = "mm", pointsize = 24, res = 100,
      )
    } else if(assemb_level == TRUE){
      png(
        here::here(
          "2_Patches", "4_OutputPlots", "3_Spatial_mapping", "3_Assemblage_level_mapping", pam_res, space, pam_type,
          png_filename
        ), 
        width = png_width, height = png_height, units = "mm", pointsize = 24, res = 100,
      )
      }
  }
  
  # set up plotting grid and remove species richness layer of raster, if not plotting
  if(plot_sr == FALSE){
    
    sr_index <- which(layer_names == "species_richness")
    new_layer_names <- layer_names[-sr_index]
    div_rast <- div_rast[[new_layer_names]]
    
    par(mfrow = c(2, 1))
    
  } else if(plot_sr == TRUE){
    
    par(mfrow = c(3, 1))
    
  }
  
  # set font size
  cexval <- 0.7
  
  # get updated names of layers, to apply plotting across
  layers <- terra::names(div_rast)
  
  # just use a for loop, for ease (could also do as an apply function)
  for(layer in layers){
    
    # set up colour breaks, palette and sex-specific legend title
    if(layer == "M"){
      breaks <- div_breaks
      colpal <- divpal
      lgd_title <- paste0(lgd_metric, "(males)")
    } else if(layer == "F"){
      breaks <- div_breaks
      colpal <- divpal
      lgd_title <- paste0(lgd_metric, "(females)")
    } else if(layer == "U"){
      breaks <- div_breaks
      colpal <- divpal
      lgd_title <- paste0(lgd_metric, "(unknown sex)")
    } else if(layer == "species_richness"){
      breaks <- sr_breaks
      colpal <- srpal
      lgd_title <- "species richness"
    }
    
    
    # first plot male diversity raster to get legend only
    terra::plot(
      div_rast, layer, 
      breaks = breaks, 
      col = colpal, 
      #  colNA = "grey",
      legend = TRUE,
      frame = FALSE,
      axes = FALSE, 
      plg = list(title = lgd_title, x = "left", cex = cexval),
      mar = c(1, 2, 1, 2)
    )
    # overlay grey world map so that NA parts of world are grey
    terra::plot(
      world_spatvec, 
      col="grey80", 
      border = NA, 
      add = TRUE,
      axes = "n"
    )
    # overlay male diversity raster again, but without legend
    terra::plot(
      div_rast, layer, 
      breaks = breaks, 
      col = colpal, 
      # colNA = "grey",
      add = TRUE,
      legend = FALSE,
      mar = c(1, 2, 1, 2)
    )
    
  }
  
  # end png saving, if initialised
  if(save_png == TRUE){
    dev.off()
  }
  
  
}

# Functions for calculating diversity by bioregion ----


# get df of chosen metric for each species - for calculating diversity vs global centroid (only works for centroid distance)
metric_vals <- function(pca_data, metric){
  
  # set metric
  metric_get <- set_metric(metric)
  
  # calculate metric
  metric_values <- dispRity::dispRity(pca_data, metric = metric_get)$disparity[[1]][[1]]
  
  # add rownames
  rownames(metric_values) <- rownames(pca_data)
  
  return(metric_values)
  
}


# function to calculate average value of metric for males and females for each ecoregion
# returns [1, 2] matrix that can be inputted to ecoregion data
calc_regions_global_metric <- function(region_index, region_shapes, region_type, metric, div_loss_type,  metric_vals = NULL, pca_data = NULL, null_terrast, avg_par){
  
  # check if combination of parameters is possible
  if(div_loss_type == "global" & metric != "centr-dist"){
    stop(paste("Combination of", div_loss_type, "div_loss_type and", metric, "metric is not possible.
               Please choose new parameters and try again."))
  }
  
  # extract the region
  ecoreg <- region_shapes[region_index, ]
  
  # for debugging purposes (only works if regions == "ecoregions)
  # ecoreg_num <- ecoreg$OBJECTID
  # ecoreg_name <- ecoreg$ECO_NAME
  
  # get number of cells in null raster
  ncells <- length(terra::values(null_terrast))
  
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
    
    # if using global centroid vals
    if(div_loss_type == "global" & metric == "centr-dist"){
      
      # if metric values have already been supplied then just take the average
      if(!is.null(metric_vals)){
        avg_metric <- avg(metric_vals[species_pres, ], avg_type = avg_par, na.rm = TRUE)
      } else
        # if values not supplied, calculate values relative to global centroid
        if(is.null(metric_vals) & !is.null(pca_data)){
          
          # subset PCA data
          region_data <- pca_data[species_pres, ]
          
          # calculate metric values relative to global centroid
          region_metric_vals <- dispRity::dispRity(region_data, metric = set_metric(metric), centroid = rep(0, ncol(region_data)))$disparity[[1]][[1]]
          
          # get average metric value for region
          avg_metric <- avg(region_metric_vals, avg_type = avg_par, na.rm = TRUE)
        }
      
    } else 
      # if using local centroid vals (or other metric)
      if(div_loss_type == "local"){
        
        # subset PCA data
        region_data <- pca_data[species_pres, ]
        
        # calculate metric values relative to global centroid
        region_metric_vals <- dispRity::dispRity(region_data, metric = set_metric(metric))$disparity[[1]][[1]]
        
        # get average metric value for region
        avg_metric <- avg(region_metric_vals, avg_type = avg_par, na.rm = TRUE)
        
      }
    
    
  } else {
    
    # if there are fewer species than threshold, set values = NA
    avg_metric <- NA
    
  }
  
  # set species richness value
  sr <- length(species_pres)
  
  avgs <- c(avg_metric, sr)
  
  if(regions == "ecoregions"){
    ## SET VALUES FOR ECOREGION OBJECTID 207 - "Rock and Ice" - BIOM_NUM 11 - to 
    ## NA - it seems, from looking at Dinerstein et al (2017) and their interactive map
    ## https://ecoregions.appspot.com/ that this should be possibly excluded as it's
    ## basically uninhabited - it has an artificially extremely high species richness
    ## because it spans all of Antarctica, Greenland, parts of the Himalaya, Iceland, 
    ## and Northwestern North America
    if(ecoreg$OBJECTID == 207){
      avgs <- rep(NA, times = 2)
    }
  }
  
  
  return(avgs)
  
}

