# Extracting and comparing measures of diversity from patch colour pattern spaces
# Sensitivity analyses to different ways to clip the data
# 1st November 2024
# Robert MacDonald

# Clear environment
rm(list=ls())

# Load libraries
library(dplyr)
library(dispRity)

## EDITABLE CODE
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Passeriformes"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
# N.B. will need to change the date in the pca_all filename if using usmldbl or usmldblr
# (from 240603 to 240806)
space <- "lab"
# select sex (""allspecimens" or "matchedsex")
# "matchedsex" will subset to only species for which we have at least one
# male and female specimen
sex_match <- "all"
# select metrics c("centr-dist", "nn-k", "nn-count")
# note that nn-k is EXTREMELY slow to run - needs parallelisation (but will probably still
# be too slow to run)
metric <- c("centr-dist", "nn-k", "nn-count")
# choose hyperparameters (number of nearest neighbours to average over, radius to count
# neighbours within, whether or not to return nearest neighbour counts as relative proportion)
nn <- 10
nn_radius <- 100
rad_relative <- TRUE

# Functions ----

# Extract PCA co-ordinates
pca_extract <- function(pca_object){
  
  # extract pca co-ordinates
  pca_vals <- pca_object$x

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

# dispRity-style function to calculate mean distance to given number of nearest neighbours - WITH removal of
# duplicate same-species-different-sex pairs

mean.nn.dist <- function(matrix, rownames = NULL, nn = 5, method = "euclidean", ...) {
  
  ## Set number of neighbours to all if not specified
  if(is.null(nn)){
    nn <- nrow(matrix)
  }
  
  ## Set matrix rownames if specified
  if(!is.null(rownames)){
    rownames(matrix) <- rownames
  }
  
  ## Calculate all pairwise distances
  pair.dists <- as.matrix(dist(matrix, method = method))

  ## Function to calculate nn closest for a single point
  nn_point <- function(one.row, nn, row.index){
    
    # if no species/sex provided just calculate the basic mean nearest neighbour distance
    if(is.null(names(one.row))){
      nn_point_dists <- mean(sort(one.row)[2:(nn + 1)], na.rm = TRUE)
    }
    
    # if species/sex provided, check if more than one instance of same species is in the n nearest neighbours and ignore if so
    else {
      species <- strsplit(names(one.row)[row.index], split = "-")[[1]][1]
      nn_names <- names(sort(one.row)[2:(nn + 1)])
      nn_species <- sapply(strsplit(nn_names, split = "-"), "[", 1)
      # check if same species is nn
      if(any(grepl(species, nn_species))){
        # get position of same species (this will include the specimen of interest)
        pos <- grep(species, names(one.row))
        # get values for these positions (one will be zero - the specimen of interest)
        pos_vals <- one.row[pos]
        names(pos_vals) <- as.character(pos)
        # get position of non-zero same-species specimen (i.e. the opposite sex)
        pos <- as.numeric(names(pos_vals[pos_vals != 0]))
        # set the value for the same species equal to NA
        one.row[pos] <- NA
        # calculate knn mean, removing NAs
        nn_point_dists <- mean(sort(one.row)[2:(nn + 1)], na.rm = TRUE)
        #  nn_point_dists <- mean(c(sort(one.row[2:(pos - 1)]), sort(one.row[(pos + 1):(nn + 2)])))
      }
      else{
        nn_point_dists <- mean(sort(one.row)[2:(nn + 1)], na.rm = TRUE)
      }
    }
    return(nn_point_dists)
  }
  
  # create wrapper function to preserve element names when passing each row to nn_point
  wrapper_function <- function(row.index, matrix, nn) {
    row <- matrix[row.index, ]
    nn_point(row, nn, row.index)
  }
  
  ## Get nn closest for each point
  # mean.nn.dists <- apply(pair.dists, 1, nn_point, nn = nn)
  mean.nn.dists <- sapply(seq_len(nrow(pair.dists)), wrapper_function, matrix = pair.dists, nn = nn)
  
  
  ## Return values
  return(mean.nn.dists)
}

# dispRity-style function to count number of neighbours within a given radius (written by Thomas Guillerme)
count.neighbours <- function(matrix, radius = 1, relative = TRUE, method = "euclidean", ...) {
  ## Check if the matrix is a distance matrix first
  distances <- as.matrix(dist(matrix, method = method))
  ## Set the radius to something if it's a function
  if(is(radius, "function")) {
    radius <- radius(distances)
  }
  ## For each row count how many distances are < radius (minus one is for the diagonal that's = 0)
  counts <- apply(distances, 1, function(one_row, radius) sum(one_row <= radius), radius = radius) - 1
  ## Return the counts
  if(relative) {
    return(unname(counts/ncol(distances)))
  } else {
    return(unname(counts))
  }
}

# Calculate diversity metric for given data
calc_metric <- function(data, disp_metric = c("centr-dist", "nn-k", "nn-count"), ...){
  
  # set dispRity metric
  if(disp_metric == "centr-dist"){
    metric_get <- "centroids"
  } else if (disp_metric == "nn-k"){
    metric_get <- "mean.nn.dist"
  }else if (disp_metric == "nn-count"){
    metric_get <- "count.neighbours"
  }
  
  # coerce data to matrix if not already a matrix
  if(!"matrix" %in% class(data)){
    data <- as.matrix(data)
  }
  
  # coerce to numeric
  class(data) <- "numeric"
  
  # calculate value of metric for each species
  div_values <- dispRity(
    data,
    metric = get(metric_get),
    ...
  )$disparity[[1]][[1]]
  
  # set row and column names
  rownames(div_values) <- rownames(data)
  colnames(div_values) <- disp_metric
  
  return(div_values)
  
}

# Produce text file of metadata for csv file
write_metadata <- function(csv_filename, pca_dat, clade, metric, nn, nn_radius, sex_match, space){
  
  # get species and sex column (to count Nspecies)
  spp_sex <- pca_spp_sex(pca_dat)
  
  # count unique male species
  n_male_spp <- length(unique(spp_sex[spp_sex[, "sex"] == "M", "species"]))
  
  # count unique female species
  n_female_spp <- length(unique(spp_sex[spp_sex[, "sex"] == "F", "species"]))
  
  # count unique unknown sex species
  n_unk_spp <- length(unique(spp_sex[spp_sex[, "sex"] == "U", "species"]))
  
  # define metadata text file filename
  metadata_filename <- sub(".csv", "_metadata.txt", csv_filename)
  
  # write metadata to text file
  # first define subsetting text outside of paste function
  if(sex_match == "matchedsex"){
    subset_text <- "Species list was trimmed to only species for which we have male and female specimens prior to calculation of metrics."
  } else {
    subset_text <- paste("Species subsetting:", sex_match, sep = " ")
  }
  writeLines(
    paste(
      # Line 1 - corresponding csv file
      paste("Metadata to accompany CSV file:", csv_filename, sep = " "),
      # Line 2 - clade/space
      paste("File contains diversity metric information for", clade, "in", space, "space", sep = " "),
      # Line 3 - species subsetting
      subset_text,
      # Line 4 - give number of species of males and females
      paste("Nmales =", n_male_spp, "Nfemales =", n_female_spp, "Nunknownsex =", n_unk_spp, sep = " "),
      # Line 5 - diversity metrics 
      paste("Diversity metrics calculated:", paste(metric, collapse = ", "), sep = " "),
      # Line 6 - hyperparameters
      paste("Hyperparameters: N_nearestneighbours (Knn) =", nn, "; Nearest neighbour radius =", nn_radius, "; Return relative nearest neighbour count =", as.character(rad_relative), sep = " "),
      # Line 7 - timestamp
      paste("Metadata generated on", Sys.time(), sep = " "),
      sep = "\n\n"
      ), 
    con = metadata_filename
  )
  
}



# Load data ----

# load patch data (PCA of whichever colourspace - generated in 02_Patch_Analyse_features.R)
pca_filename <- paste(clade, "patches.231030.PCAcolspaces", space, "240925", "rds", sep = ".")
pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)

# Analysis - metric production ----

# extract PCA co-ordinates
pca_dat <- pca_extract(pca_all)

# add species and sex columns (from rownames)
# only necessary if we're subsetting in some way
if(sex_match != "all"){
  pca_dat <- pca_spp_sex(pca_dat)
  }

# if required, get species for which we have matched male/female data
if(sex_match == "matchedsex"){
  
  # get species to keep
  spp_to_keep <- sex_match_fun(pca_dat)
  
  # subset data to species of interest
  pca_dat <- subset_sex_match(pca_dat, spec_to_keep = spp_to_keep)
  
}


# calculate values of metrics for each species/sex pair
div_metric_vals_list <- list()
div_metric_vals_list <- lapply(metric, function(x)
  
  calc_metric(data = pca_dat, disp_metric = x, nn = nn, rownames = rownames(pca_dat), radius = nn_radius, relative = rad_relative)
  
)

# combine into one matrix
div_metric_vals <- do.call(cbind, div_metric_vals_list)

# write to csv
csv_filename <- here::here(
  "2_Patches", "3_OutputData", "4_Diversity_measures", space,
  paste(clade, "patches", space, sex_match, "diversitymetrics.csv", sep = "_")
)
write.csv(
  div_metric_vals,
  csv_filename
)
# write accompanying metadata to text file
write_metadata(csv_filename, pca_dat, clade, metric, nn, nn_radius, sex_match, space)


# Downstream analysis ----

# choose two metrics of interest
x_metric <- "centr-dist"
y_metric <- "nn-k"

# check if two metrics are correlated
corr_div <- cor.test(div_metric_vals[, x_metric], div_metric_vals[, y_metric])

# inspect results
corr_div

# plot (basic)
plot(div_metric_vals[, x_metric], div_metric_vals[, y_metric])
abline(lm(div_metric_vals[, y_metric] ~ div_metric_vals[, x_metric]))

# fit linear regression model
mod_div <- lm(div_metric_vals[, y_metric] ~ div_metric_vals[, x_metric])

# calculate standardised residuals
resids <- rstandard(mod_div)

# bind to results matrix
div_metric_vals <- cbind(div_metric_vals, std_residuals = resids)

# inspect which birds are either in a particularly sparse area of 
# colour space given their distance to centroid (i.e. species which have
# the highest positive residuals)
head(div_metric_vals[order(-(div_metric_vals[, "std_residuals"])), ], n = 20)
# and those which exist in a particularly dense area of colour space given 
# their distance to centroid (i.e. species with lowest negative residuals)
head(div_metric_vals[order(div_metric_vals[, "std_residuals"]), ], n = 20)

# plot with colour grids
x_axis <- x_metric
y_axis <- y_metric
png(
  here::here(
    "2_Patches", "4_OutputPlots", "2_Diversity_measures", space, 
    paste(clade, "patches", space, sex_match, x_axis, y_axis, "colour_grids.png", sep = "_")
  ),
  width = 4000, height = 2/3 * 4000,
  units = "px",
  pointsize = 50
  
)
# 

xlabel <- x_metric
ylabel <- y_metric
plot(div_metric_vals[, y_axis] ~ div_metric_vals[, x_axis], 
     asp = T, 
     type = "n", 
     xlab = xlabel, ylab = ylabel, las = 1,
#     xlim = 0:max(),  ylim = yrange
)

box(lwd = 2)
rw <- diff(range(div_metric_vals[,1]))/12

for(i in 1:nrow(div_metric_vals)){
  fname <- rownames(div_metric_vals)[i]
  fpng <- png::readPNG(paste0(
    "C:/Users/bop23rxm/Documents/colour_grids_repositioned/", 
    paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
  rasterImage(fpng, 
              xleft = div_metric_vals[i, x_axis] - (rw/15), 
              ybottom = div_metric_vals[i, y_axis]-(rw/9), 
              xright = div_metric_vals[i, x_axis]+(rw/15), 
              ytop = div_metric_vals[i, y_axis]+(rw/9))
  cat(paste0("\r", i, " of ", nrow(div_metric_vals), " processed"))
}
dev.off()
