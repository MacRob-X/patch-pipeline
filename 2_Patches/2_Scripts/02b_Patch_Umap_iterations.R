# Running UMAP on patch data - multiple iterations and varying parameters
# Robert MacDonald
# 26th March 2024

# clear environment
rm(list=ls())

# load libraries
library(umap)

# Custom functions ----
source(
  here::here(
    "2_Patches", "2_Scripts", "2_BetaVersions", "Patch_Plotting_functions_v1.r"
  )
)

# load patch PCA data
pca_all <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA", 
    "Neoaves.patches.231030.PCAcolspaces.jndxyzlumr.240603.rds"
  )
)

# load taxonomic data
taxo <- readr::read_csv(
  here::here(
    "4_SharedInputData", "BLIOCPhyloMasterTax_2019_10_28.csv"
   )
)


#-------------------------------------------------------------------------------------------------#

# Perform UMAP four times with the default parameters to see if clustering is visually different 
# Plot and save as png
nUMAPs <- 4
rawUMAPs <- vector("list", length = nUMAPs)

# perform UMAPs
for (i in 1:nUMAPs){
  cat("\r", i)
  rawUMAPs[[i]] <- umap(pca_all$x)
}

saveRDS(rawUMAPs, 
        file = here::here(
          "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP", 
          "Neoaves.patches.pca.jndxyzlumr.UMAPs.iterations.240603.rds"
          )
        )

# Trim and order taxonomic data to match UMAP data

# create dataframe from UMAP to manipulate
cleanUMAPs <- vector("list", length = nUMAPs)
for (j in 1:nUMAPs){
  cleanUMAPs[[j]] <- data.frame(rawUMAPs[[j]]$layout)
  l <- length(rawUMAPs[[j]]$layout[,1])
  for(i in 1:l){
    cat("\r", i)
    split <- strsplit(rownames(cleanUMAPs[[j]])[i], split = "-")[[1]]
    cleanUMAPs[[j]]$species[i] <- split[1]
    cleanUMAPs[[j]]$sex[i] <- split[2]
  }
}

# trim and order taxonomy data to match UMAP data
taxo$specimen <- c(rep(NA, times = length(taxo$TipLabel)))

# get female specimen taxonomic data
femTaxo <- taxo[paste0(taxo$TipLabel, "-F") %in% rownames(cleanUMAPs[[1]]), ]
femTaxo$specimen <- paste0(femTaxo$TipLabel, "-F")

# get male specimen taxonomic data
maleTaxo <- taxo[paste0(taxo$TipLabel, "-M") %in% rownames(cleanUMAPs[[1]]), ]
maleTaxo$specimen <- paste0(maleTaxo$TipLabel, "-M")

# get unknown sex specimen taxonomic data
unkTaxo <- taxo[paste0(taxo$TipLabel, "-U") %in% rownames(cleanUMAPs[[1]]), ]
unkTaxo$specimen <- paste0(unkTaxo$TipLabel, "-U")

# concatenate into one big dataframe with specimen names that match umapDat rownames
taxoMatch <- rbind(femTaxo, maleTaxo, unkTaxo)

# remove temp dataframes
rm(femTaxo, maleTaxo, unkTaxo)

# remove species not present in taxonomic data
cleanUMAPs <- lapply(cleanUMAPs, function(df){
  df <- df[rownames(df) %in% taxoMatch$specimen,]
  return(df)
})

# reorder taxo data to match umapDat
taxoMatch <- taxoMatch[match(rownames(cleanUMAPs[[1]]), taxoMatch$specimen), ]

# Check all specimens match between the two datasets
if(identical(taxoMatch$specimen, rownames(cleanUMAPs[[1]]))){
  print("All specimens in taxonomic and colour data matched")
}

# attach taxon subgroups to UMAP data and rename columns
cleanUMAPs <- lapply(cleanUMAPs, function(df){
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$taxon_subgroup <- taxoMatch$Taxon_subgroup
  return(df)
})


# plot the four UMAPs together, coloured by taxon subgroup
plot.FourPatch(cleanUMAPs[[1]],
               cleanUMAPs[[2]],
               cleanUMAPs[[3]],
               cleanUMAPs[[4]],
               taxonCol = TRUE,
               filepath = here::here(
                 "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", "umap_iterations",
                 "Neoaves.patches.pca.jndxyzlumr.UMAPiterations.240603.png"
                 )
               )

# check if identical
identical(cleanUMAPs[[1]], cleanUMAPs[[2]])
identical(cleanUMAPs[[2]], cleanUMAPs[[3]])
identical(cleanUMAPs[[3]], cleanUMAPs[[4]])

#--------------------------------------------------------------------------------#
# Perform UMAP 8 times with different parameters to visually check if clustering is affected

# set n_neighbours (based on Alam et al., 2024, Nature)
# default is 15
# note that runs with nn > 100 take a long time to run ( ~30 mins for nn = 200)
customConfig <- umap.defaults
nn <- c(5, 50, 200, 1000)
nnUMAPs <- vector("list", length = length(nn))
for (i in 1:length(nn)){
  cat("\r", i)
  # set custom configuration
  customConfig$n_neighbors <- nn[i]
  nnUMAPs[[i]] <- umap(pca_all$x, config = customConfig)
}

saveRDS(nnUMAPs, 
        file = here::here(
          "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP", 
          "Neoaves.patches.pca.jndxyzlumr.nnUMAPs.240603.rds"
        )
      )

# create list of dataframes from UMAPs
nncleanUMAPs <- lapply(nnUMAPs, function(df){
  df <- data.frame(df$layout)
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$taxon_subgroup <- taxoMatch$Taxon_subgroup
  return(df)
})


# plot the four UMAPs together, coloured by taxon subgroup
plot.FourPatch(nncleanUMAPs[[1]],
               nncleanUMAPs[[2]],
               nncleanUMAPs[[3]],
               nncleanUMAPs[[4]],
               taxonCol = T,
               filepath = here::here(
                 "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", "umap_iterations",
                 "Neoaves.patches.pca.jndxyzlumr.nnUMAPs.5.50.200.1000.240603.png"
               )
)


# set minimum distance (again based on Alam et al., 2024, Nature)
# default is 0.1
customConfig <- umap.defaults
minDist <- c(0.1, 0.25, 0.5, 0.99)
minDistUMAPs <- vector("list", length = length(minDist))
for (i in 1:length(minDist)){
  cat("\r", i)
  # set custom configuration
  customConfig$min_dist <- minDist[i]
  minDistUMAPs[[i]] <- umap(pca_all$x, config = customConfig)
}

saveRDS(minDistUMAPs, 
        file = here::here(
          "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP", 
          "Neoaves.patches.pca.jndxyzlumr.minDistUMAPs.240604.rds"
        )
      )

# create list of dataframes from UMAPs
mindistcleanUMAPs <- lapply(minDistUMAPs, function(df){
  df <- data.frame(df$layout)
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$taxon_subgroup <- taxoMatch$Taxon_subgroup
  return(df)
})

# plot the four UMAPs together, coloured by taxon subgroup
plot.FourPatch(mindistcleanUMAPs[[1]],
               mindistcleanUMAPs[[2]],
               mindistcleanUMAPs[[3]],
               mindistcleanUMAPs[[4]],
               taxonCol = T,
               filepath = here::here(
                 "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", "umap_iterations",
                 "Neoaves.patches.pca.jndxyzlumr.minDistUMAPs.1.25.50.99.240603.png"
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


                    
spotCheckUMAP(nncleanUMAPs[[1]], type = "range", lower1 = -5, upper1 = 5, lower2 = -5, upper2 = 5)        





cols <- rainbow(length(levels(as.factor(nncleanUMAPs[[1]]$taxon_subgroup))))
col_mapping <- cols[match(nncleanUMAPs[[1]]$taxon_subgroup, levels(as.factor(nncleanUMAPs[[1]]$taxon_subgroup)))]


plot(nncleanUMAPs[[1]][,2] ~ nncleanUMAPs[[1]][,1],
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
                    