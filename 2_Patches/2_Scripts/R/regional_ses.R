
# calculate SES for trimming of a single IUCN level (cumulatively)
calc_iucn_ses <- function(pca_region, region_species, iucn_data, iucn_cat, metric, n_sims, avg_par, parallel_run = FALSE, cluster = NULL){
  
  # set dispRity metric
  # if(metric == "centr-dist"){
  #   metric_get <- get("centroids")
  # } else if (metric == "nn-k"){
  #   metric_get <- get("mean.nn.dist")
  # }else if (metric == "nn-count"){
  #   metric_get <- get("count.neighbours")
  # } else if(metric == "sum.variances"){
  #   metric_get <- list(sum, variances)
  # } else if(metric == "sum.ranges"){
  #   metric_get <- list(sum, ranges)
  # }
  metric_get <- set_metric(metric)
  
  
  # set up vector to store results
  res_vec <- rep(NA, times = 6)
  names(res_vec) <- c("iucn", "species_richness", "raw", "null_mean", "null_sd", "null_se")
  
  
  
  # get list of species within region with threat level of interest trimmed (inclusive - i.e. trims 
  # threat level of interest and any more threatened categories)
  if(iucn_cat == "CR"){
    trimmed_species <- iucn[iucn$species_birdtree %in% region_species & iucn$iucn_cat != "CR", "species_birdtree"]
  } else if(iucn_cat == "EN"){
    trimmed_species <- iucn[iucn$species_birdtree %in% region_species & iucn$iucn_cat != "CR" & iucn$iucn_cat != "EN", "species_birdtree"]
  } else if(iucn_cat == "VU"){
    trimmed_species <- iucn[iucn$species_birdtree %in% region_species & iucn$iucn_cat != "CR" & iucn$iucn_cat != "EN" & iucn$iucn_cat != "VU", "species_birdtree"]
  } else if(iucn_cat == "NT"){
    trimmed_species <- iucn[iucn$species_birdtree %in% region_species & iucn$iucn_cat != "CR" & iucn$iucn_cat != "EN" & iucn$iucn_cat != "VU" & iucn$iucn_cat != "NT", "species_birdtree"]
  }
  
  # check if more than one species in trimmed species list
  # if not, set all diversity metrics to 0 and end function
  if(length(trimmed_species) <= 1){
    
    res_vec["species_richness"] <- length(trimmed_species)
    res_vec["raw"] <- 0
    res_vec["null_mean"] <- NA
    res_vec["null_sd"] <- NA
    res_vec["null_se"] <- NA
    
    return(res_vec)
    
  }
  
  # clip PCA data to only species present in regions (no CR species)
  pca_subset <- pca_region[trimmed_species, ]
  
  # calculate mean of metric of all species in ROI in IUCN cats of interest
  roi_metric_avg <- avg(dispRity(pca_subset, metric = metric_get)$disparity[[1]][[1]], avg_type = avg_par, na.rm = TRUE)
  
  # add species richness and average metric to results vector
  res_vec["species_richness"] <- length(trimmed_species)
  res_vec["raw"] <- roi_metric_avg
  
  # set up matrix to store simulations
  sims <- matrix(NA, nrow=1, ncol = n_sims)
  
  ## For parallelised working
  if(parallel_run == TRUE){
    # set up cluster (note that it's better in terms of RAM to do this within this function, rather than
    # before starting, although it means longer overhead times)
    no_cores <- parallel::detectCores() - 4
    cluster <- parallel::makeCluster(no_cores)
    
    # export necessary objects to the cluster
    parallel::clusterExport(cluster, c("region_species", "trimmed_species", "pca_region", "metric_get", "dispRity", "avg", "avg_par", metric_get), envir = environment())
    
    
    # Calculate SES compared to simulations of random loss of species, using parallelised lapply
    sims[1, ] <- unlist(
      parallel::parLapply(
        cl = cluster,
        1:ncol(sims), 
        function(j) {
          
          coms <- sample(x = region_species, size = length(trimmed_species), replace = FALSE)
          trait_vec <- matrix(pca_region[coms, ], ncol = dim(pca_region)[2])
          
          # Calculate the disparity value for the current simulation (column j)
          disparity_value <- avg(dispRity(trait_vec, metric = metric_get)$disparity[[1]][[1]], avg_type = avg_par, na.rm = TRUE)
          
          return(disparity_value)
        }
      )
    )
    
    # stop cluster
    stopCluster(cluster)
  } else if (parallel_run == FALSE){
    sims[1, ] <- unlist(
      lapply(
        1:ncol(sims), 
        function(j) {
          
          coms <- sample(x = region_species, size = length(trimmed_species), replace = FALSE)
          trait_vec <- matrix(pca_region[coms, ], ncol = dim(pca_region)[2])
          
          # Calculate the disparity value for the current simulation (column j)
          disparity_value <- avg(dispRity(trait_vec, metric = metric_get)$disparity[[1]][[1]], avg_type = avg_par, na.rm = TRUE)
          
          return(disparity_value)
        }
      )
    )
  }
  
  
  # add to results vector
  res_vec["null_mean"] <- mean(sims)
  res_vec["null_sd"] <- sd(sims)
  res_vec["null_se"] <- sd(sims)/sqrt(length(sims))
  
  
  
  # return results
  return(res_vec)
  
}

# calculate SESs for sequential trimming of IUCN levels for a single region
calc_region_ses <- function(region_sf, pca_data, pam, null_raster, iucn, metric, n_sims, regions, append_sf = TRUE, ...){
  
  # show which region is being processed
  if(regions == "biomes"){
    cat("\rProcessing region", region_sf$BIOME_NUM)
  } else if(regions == "ecoregions"){
    cat("\rProcessing region", region_sf$OBJECTID)
  }
  
  
  # set dispRity metric
  if(metric == "centr-dist"){
    metric_get <- get("centroids")
  } else if (metric == "nn-k"){
    metric_get <- get("mean.nn.dist")
  }else if (metric == "nn-count"){
    metric_get <- get("count.neighbours")
  } else if(metric == "sum.variances"){
    metric_get <- list(sum, variances)
  }  else if(metric == "sum.ranges"){
    metric_get <- list(sum, ranges)
  }
  
  # set up results matrix
  results <- matrix(NA, nrow=5, ncol = 7)
  colnames(results) <- c("iucn", "species_richness", "raw", "null_mean", "null_sd", "null_se", "ses")
  
  # create null rast, subsetted to extent of ecoregion
  # first assign values to null_rast to keep track of grid cells
  null_rast_sub <- null_raster
  # get number of grid cells in full raster
  ncells <- length(terra::values(null_raster))
  terra::values(null_rast_sub) <- 1:ncells
  # crop raster to ecoregion extent
  null_rast_sub <- terra::crop(null_rast_sub, region_sf)
  # mask raster to ecoregion extent
  null_rast_sub <- terra::mask(null_rast_sub, region_sf)
  
  # use values of subsetted raster to subset PAM to ROI
  pam_sub <- pam[terra::values(null_rast_sub), ]
  
  # subset pam to only species which are present in the PCA data (this is necessary
  # if we're using all species, not just those with both male and female specimens)
  pam_sub <- subset_pam(pam_sub, get_unique_spp(pca_data))
  
  # get species which are present in subsetted PAM
  if(length(terra::values(null_rast_sub)) < 2){
    # if ROI consists of 1 or fewer grid cells
    species_pres <- names(which(!is.na(pam_sub)))
  } else {
    # if ROI consists of more than 1 grid cell
    species_pres <- names(which(apply(pam_sub, 2, function(x) any(!(is.na(x))))))
  }
  
  # check if there are actually any species present 
  # if no species present, return null results from function
  if(append_sf == TRUE){
    if(length(species_pres) == 0 | is.null(species_pres)){
      results <- rep("nospec_region", times = 5)
      names(results) <- c("ses_all", "ses_CR",  "ses_EN",  "ses_VU",  "ses_NT")
      return(results)
    } else if(length(species_pres) == 1){
      results <- rep("singlespec_region", times = 5)
      names(results) <- c("ses_all", "ses_CR",  "ses_EN",  "ses_VU",  "ses_NT")
      return(results)
    }
  } else if(append_sf == FALSE){
    if(length(species_pres) == 0 | is.null(species_pres)){
      return("No species in region")
    } else if(length(species_pres) == 1){
      return("Only one species in region")
    }
  }
  
  
  # clip PCA data to only species present in region
  pca_region <- pca_data[species_pres, ]
  
  # calculate average of metric of all species in region
  roi_metric_avg <- avg(dispRity(pca_region, metric = metric_get)$disparity[[1]][[1]], avg_type = avg_par, na.rm = TRUE)
  
  # add species richness and average metric to results table
  results[1, "species_richness"] <- length(species_pres)
  results[1, "raw"] <- roi_metric_avg
  
  # Now do the same but for CR species only, comparing to n simulations of random species loss to calculate
  # SES
  iucn_cats <- c("CR", "EN", "VU", "NT")
  
  # # initialise the cluster for parallelisation
  # no_cores <- parallel::detectCores() - 2
  # cl <- parallel::makeCluster(no_cores)
  
  res <- lapply(iucn_cats, calc_iucn_ses,
                pca_region = pca_region, 
                region_species = species_pres, 
                iucn_data = iucn, 
                metric = metric,
                n_sims = n_sims,
                avg_par = avg_par,
                cluster = cl)
  
  # stopCluster(cl)
  
  # populate results table
  results[2:5, 1:6] <- do.call(rbind, res)
  
  # calculate and append SES
  results[2:5, "ses"] <- (results[2:5, "raw"] - results[2:5, "null_mean"])/results[2:5, "null_sd"]
  results[1, "ses"] <- 0
  
  # add iucn category names
  results[, "iucn"] <- c("all", iucn_cats)
  
  # return only the SES column, if appending results to an sf object
  if(append_sf == TRUE){
    
    ses_results <- results[, "ses"]
    ses_results <- as.numeric(ses_results)
    names(ses_results) <- paste("ses", results[, "iucn"], sep = "_")
    
    
    return(ses_results)
    
    # else return the entire results table
  } else{
    
    return(results)
    
  }
  
  
}


# calculate SESs for sequential trimming of IUCN levels for a single taxonomic group
calc_group_ses <- function(group, pca_data, taxonomy, iucn, metric, tax_level, n_sims, ...){
  
  # show which group is being processed
  message(paste0("\rProcessing group ", group))
  
  # set dispRity metric
  metric_get <- set_metric(metric)
  
  # set up results dataframe
  results <- as.data.frame(matrix(NA, nrow=5, ncol = 8))
  colnames(results) <- c("iucn", "n_species", "raw", "null_mean", "null_sd", "null_se", "ses", "n_spec_prop")
  
  # subset taxonomy to only species in pca data
  taxonomy <- taxonomy[taxonomy[, "species"] %in% rownames(pca_data), ]
  
  # get species in group of interest
  species_pres <- taxonomy[taxonomy[, tax_level] == group, "species"]
  
  # subset PCA data to only species present in taxonomic group
  pca_group <- pca_data[species_pres, ]

  # calculate average of metric of all species in group
  group_metric_avg <- avg(dispRity(pca_group, metric = metric_get)$disparity[[1]][[1]], avg_type = avg_par, na.rm = TRUE)
  
  # add species richness and average metric to results table
  results[1, "n_species"] <- length(species_pres)
  results[1, "raw"] <- group_metric_avg
  
  # Now do the same but for CR species only, comparing to n simulations of random species loss to calculate
  # SES
  iucn_cats <- c("CR", "EN", "VU", "NT")
  
  # # initialise the cluster for parallelisation
  # no_cores <- parallel::detectCores() - 2
  # cl <- parallel::makeCluster(no_cores)
  
  res <- lapply(iucn_cats, calc_iucn_ses,
                pca_region = pca_group, 
                region_species = species_pres, 
                iucn_data = iucn, 
                metric = metric,
                n_sims = n_sims,
                avg_par = avg_par,
                cluster = cl)
  
  # stopCluster(cl)
  
  # populate results table
  results[2:5, 1:6] <- do.call(rbind, res)
  
  # calculate and append SES
  results[2:5, "ses"] <- (results[2:5, "raw"] - results[2:5, "null_mean"])/results[2:5, "null_sd"]
  results[1, "ses"] <- 0
  
  # calculate and append proportional species richness
  results[, "n_spec_prop"] <- unlist(lapply(results[, "n_species"], function(n_spec, n_spec_max){
    n_spec_prop <- n_spec / n_spec_max
  }, n_spec_max = max(results[, "n_species"])))
  
  # add iucn category names
  results[, "iucn"] <- c("all", iucn_cats)
  
  # return only the SES column, if appending results to an sf object
  return(results)
    
}
