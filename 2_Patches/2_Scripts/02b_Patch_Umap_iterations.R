# Running UMAP on patch data - multiple iterations and varying parameters
# Robert MacDonald
# 26th March 2024

# clear environment
rm(list=ls())


# Custom functions ----
source(
  here::here(
    "2_Patches", "2_Scripts", "2_BetaVersions", "Patch_Plotting_functions_v1.r"
  )
)

## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Passeriformes"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
# N.B. will need to change the date in the pca_all filename if using usmldbl or usmldblr
# (from 240603 to 240806)
space <- "usml"
# select sex ("M", "F", "All")
sex <- "F"
# select UMAP parameters
nn <- "default"
min_dist <- "default"
## END EDITABLE CODE

# load patch data (PCA of whichever colourspace - generated in 02_Patch_Analyse_features.R)
pca_filename <- paste(clade, "patches.231030.PCAcolspaces", space, "240925", "rds", sep = ".")
pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)

# load taxonomic data
taxo <- read.csv(
  here::here(
    "4_SharedInputData", "BLIOCPhyloMasterTax_2019_10_28.csv"
   )
)


#-------------------------------------------------------------------------------------------------#

# Perform UMAP four times with the default parameters to see if clustering is visually different 
# Plot and save as png
n_umaps <- 4
raw_umaps <- vector("list", length = n_umaps)
date <- gsub("-", "", Sys.Date())

# perform UMAPs
for (i in 1:n_umaps){
  cat("\r", i)
  raw_umaps[[i]] <- umap::umap(pca_all$x, preserve.seed = FALSE)
}

# save
umap_filename <- paste(clade, "patches.pca", space, "UMAPs.iterations", date, "rds", sep = ".")
saveRDS(raw_umaps, 
        file = here::here(
          "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP", 
          umap_filename
          )
        )

# Trim and order taxonomic data to match UMAP data
raw_umaps <- readRDS(file = here::here(
  "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP", "Passeriformes.patches.pca.usml.UMAPs.iterations.20240927.rds"))

# create dataframe from UMAP to manipulate
clean_umaps <- vector("list", length = n_umaps)
for (j in 1:n_umaps){
  clean_umaps[[j]] <- data.frame(raw_umaps[[j]]$layout)
  l <- length(raw_umaps[[j]]$layout[,1])
  for(i in 1:l){
    cat("\r", i)
    split <- strsplit(rownames(clean_umaps[[j]])[i], split = "-")[[1]]
    clean_umaps[[j]]$species[i] <- split[1]
    clean_umaps[[j]]$sex[i] <- split[2]
  }
}

# trim and order taxonomy data to match UMAP data
taxo$specimen <- c(rep(NA, times = length(taxo$TipLabel)))

# get female specimen taxonomic data
fem_taxo <- taxo[paste0(taxo$TipLabel, "-F") %in% rownames(clean_umaps[[1]]), ]
fem_taxo$specimen <- paste0(fem_taxo$TipLabel, "-F")

# get male specimen taxonomic data
male_taxo <- taxo[paste0(taxo$TipLabel, "-M") %in% rownames(clean_umaps[[1]]), ]
male_taxo$specimen <- paste0(male_taxo$TipLabel, "-M")

# get unknown sex specimen taxonomic data
unk_taxo <- taxo[paste0(taxo$TipLabel, "-U") %in% rownames(clean_umaps[[1]]), ]
unk_taxo$specimen <- paste0(unk_taxo$TipLabel, "-U")

# concatenate into one big dataframe with specimen names that match umapDat rownames
taxo_match <- rbind(fem_taxo, male_taxo, unk_taxo)

# remove temp dataframes
rm(fem_taxo, male_taxo, unk_taxo)

# remove species not present in taxonomic data
clean_umaps <- lapply(clean_umaps, function(df){
  df <- df[rownames(df) %in% taxo_match$specimen,]
  return(df)
})

# reorder taxo data to match umap date
taxo_match <- taxo_match[match(rownames(clean_umaps[[1]]), taxo_match$specimen), ]

# Check all specimens match between the two datasets
if(identical(taxo_match$specimen, rownames(clean_umaps[[1]]))){
  print("All specimens in taxonomic and colour data matched")
}

# attach taxon subgroups to UMAP data and rename columns
clean_umaps <- lapply(clean_umaps, function(df){
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$taxon_subgroup <- taxo_match$Taxon_subgroup
  return(df)
})


# plot the four UMAPs together, coloured by taxon subgroup
plot_filename <- paste(clade, "patches.pca", space, "UMAPs.iterations", date, "png", sep = ".")
plot.FourPatch(clean_umaps[[1]],
               clean_umaps[[2]],
               clean_umaps[[3]],
               clean_umaps[[4]],
               taxonCol = TRUE,
               filepath = here::here(
                 "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", "umap_iterations",
                 plot_filename
                 )
               )

# check if identical
identical(clean_umaps[[1]], clean_umaps[[2]])
identical(clean_umaps[[2]], clean_umaps[[3]])
identical(clean_umaps[[3]], clean_umaps[[4]])

#--------------------------------------------------------------------------------#
# Perform UMAP 8 times with different parameters to visually check if clustering is affected
date <- gsub("-", "", Sys.Date())
# set n_neighbours (based on Alam et al., 2024, Nature)
# default is 15
# note that runs with nn > 100 take a long time to run ( ~30 mins for nn = 200)
custom_config <- umap::umap.defaults
nn <- c(25, 35, 45, 55)
nn_umaps <- vector("list", length = length(nn))
for (i in 1:length(nn)){
  cat("\r", i)
  # set custom configuration
  custom_config$n_neighbors <- nn[i]
  nn_umaps[[i]] <- umap::umap(pca_all$x, config = custom_config)
}

# save
umap_filename <- paste(clade, "patches.pca", space, "nnUMAPs", date, "rds", sep = ".")
saveRDS(nn_umaps, 
        file = here::here(
          "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP", 
          umap_filename
        )
      )

# create list of dataframes from UMAPs
nnclean_umaps <- lapply(nn_umaps, function(df){
  df <- data.frame(df$layout)
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$taxon_subgroup <- taxo_match$Taxon_subgroup
  return(df)
})


# plot the four UMAPs together, coloured by taxon subgroup
plot_filename <- paste(clade, "patches.pca", space, "nnUMAPs", date, "png", sep = ".")
plot.FourPatch(nnclean_umaps[[1]],
               nnclean_umaps[[2]],
               nnclean_umaps[[3]],
               nnclean_umaps[[4]],
               taxonCol = T,
               filepath = here::here(
                 "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", "umap_iterations",
                 plot_filename
               )
)


# set minimum distance (again based on Alam et al., 2024, Nature)
# default is 0.1
custom_config <- umap::umap.defaults
min_dist <- c(0.1, 0.25, 0.5, 0.99)
min_dist_umaps <- vector("list", length = length(min_dist))
for (i in 1:length(min_dist)){
  cat("\r", i)
  # set custom configuration
  custom_config$min_dist <- min_dist[i]
  min_dist_umaps[[i]] <- umap::umap(pca_all$x, config = custom_config)
}

umap_filename <- paste(clade, "patches.pca", space, "minDistUMAPs", date, "rds", sep = ".")
saveRDS(min_dist_umaps, 
        file = here::here(
          "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP", 
          plot_filename
        )
      )

# create list of dataframes from UMAPs
min_dist_clean_umaps <- lapply(min_dist_umaps, function(df){
  df <- data.frame(df$layout)
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$taxon_subgroup <- taxo_match$Taxon_subgroup
  return(df)
})

# plot the four UMAPs together, coloured by taxon subgroup
plot_filename <- paste(clade, "patches.pca", space, "minDistUMAPs", date, "png", sep = ".")
plot.FourPatch(min_dist_clean_umaps[[1]],
               min_dist_clean_umaps[[2]],
               min_dist_clean_umaps[[3]],
               min_dist_clean_umaps[[4]],
               taxonCol = T,
               filepath = here::here(
                 "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", "umap_iterations",
                 plot_filename
               )
              )


# function to spot check species in certain ranges

spotCheckUMAP <- function(umap, type, lower1, upper1, lower2, upper2, showValues = FALSE){
  
  # check columns correctly named
  if(!(all(colnames(umap)[1:2] == c("UMAP_1", "UMAP_2")))){
    stop("Please ensure your columns are named 'UMAP_1' and 'UMAP_2'.")
  }
  
  # check if show values
  if(showValues == FALSE){
    
    # check what type of spot check to do
    
    # max of UMAP_1
    if (type == "max1"){
      return(rownames(umap[umap$UMAP_1 == max(umap$UMAP_1),]))
    }
    # max of UMAP_2
    if (type == "max2"){
      return(rownames(umap[umap$UMAP_2 == max(umap$UMAP_2),]))
    }
    # min of UMAP 1
    if (type == "min1"){
      return(rownames(umap[umap$UMAP_1 == min(umap$UMAP_1),]))
    }
    # min of UMAP_2
    if (type == "min2"){
      return(rownames(umap[umap$UMAP_2 == min(umap$UMAP_2),]))
    }
    # within range
    if(type == "range"){
      # check range provided
      if(is.null(lower1) | is.null(upper1) | is.null(lower2) | is.null(upper2)){
        stop("Please provide a range value for all boundaries.")
      }
      # check if any values in range
      if(any(umap$UMAP_1 > lower1 & umap$UMAP_1 < upper1 &
             umap$UMAP_2 > lower2 & umap$UMAP_2 < upper2) == FALSE){
        return("No specimens present in range provided.")
      }
      else {
        return(rownames(umap[umap$UMAP_1 > lower1 & umap$UMAP_1 < upper1 &
                                          umap$UMAP_2 > lower2 & umap$UMAP_2 < upper2,]))
      }
    }
  }
  else if(showValues == TRUE){
    
  }
}


                    
spotCheckUMAP(nnclean_umaps[[1]], type = "range", lower1 = -5, upper1 = 5, lower2 = -5, upper2 = 5)        





cols <- rainbow(length(levels(as.factor(nnclean_umaps[[1]]$taxon_subgroup))))
col_mapping <- cols[match(nnclean_umaps[[1]]$taxon_subgroup, levels(as.factor(nnclean_umaps[[1]]$taxon_subgroup)))]


plot(nnclean_umaps[[1]][,2] ~ nnclean_umaps[[1]][,1],
     asp = T,
     xlab = "UMAP 1", ylab = "UMAP 2", 
     las = 1,
     col = col_mapping); box(lwd=2); title("UMAP of JND xyzlumr colour space. n_neighbour = 5; min dist = 0.1")


#---------------------------------------------------------------------------------------------------------------#
## Perform t-SNE once to obtain basic t-SNE visualisation
## parameters all set to default other than epoch (which controls verbosity, not actual t-SNE parameters)
tsne_base <- tsne::tsne(pca_all$x, 
                        k = 2, 
                        perplexity = 30,
                        max_iter = 1000,
                        epoch = 50)

## Perform t-SNE four times with the default parameters to see if clustering is visually different
## Plot and save as png

# set number of iterations
n_tsne <- 4
raw_tsnes <- vector("list", length = n_tsne)

# perform t_SNEs
for (i in 1:n_tsne){
  cat("\r", i)
  raw_tsnes[[i]] <- tsne::tsne(pca_all$x)
}
                    