# Testing assemblage-level convex hull volume calculation


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
  
  # calculate average if output is dimension-level 2
  metric_type <- dispRity::make.metric(metric_get, silent = TRUE)[["type"]]
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

str(pam)
str(sexed_metric_list$M)
str(metric)
str(avg_par)
str(sr_threshold)


div_values <- apply(pam, 1, calc_metric_gridcell, pca_data = sexed_metric_list$M, metric = metric, avg_par = avg_par, min_species = sr_threshold)


# try out hypervolume package
install.packages("hypervolume")
hypervolume::hypervolume(pca_dat)
