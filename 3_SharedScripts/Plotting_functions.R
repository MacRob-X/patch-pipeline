# Functions for visualising PCA- and CNN-generated colour pattern spaces
# 29 February 2024
# Robert MacDonald

library(phytools)

# Basic plotting function ----
# Plots single 2d plot
plot.PixVals <- function(uout, xlabel = "x-axis", ylabel = "y-axis", asp = T) {
  
  # If asp not selected, set x and y axis limits and plot blank canvas
  if(asp != T){
    # Set x and y limits
    xlimit <- c(round(min(uout[, 1]) - 1, digits = 0), round(max(uout[, 1]) + 1, digits = 0))
    ylimit <- c(round(min(uout[, 2]) - 1, digits = 0), round(max(uout[, 2]) + 1, digits = 0))
    # Plot blank canvas
    plot(uout[,2] ~ uout[,1], type="n", 
         xlim = xlimit, ylim = ylimit,
         xlab = xlabel, ylab = ylabel, las = 1); box(lwd=2)
  }
  # if asp selected, plot blank canvas with asp = T (aspect ratio fixed so x and y axis have same scale
  # - this makes the bird pngs the correct shape)
  else if(asp == T){
    plot(uout[,2] ~ uout[,1], type="n", 
         asp = T,
         xlab = xlabel, ylab = ylabel, las = 1); box(lwd=2)
  }
  
  # Add in bird torpedo pngs
  rw <- diff(range(uout[,1]))/12
  for (i in 1:nrow(uout)) {
    cat("\r", i)
    fname <- rownames(uout)[i]
    fname <- paste0(paste(strsplit(fname, split="_")[[1]][1:6], collapse="_"), ".png")
    fpng <- png::readPNG(paste0("./Outputs/pngs/warped/AF4/back/", fname))
    rasterImage(fpng, xleft = uout[i,1]-(rw/2), ybottom = uout[i,2]-(rw/6), xright = uout[i,1]+(rw/2), ytop = uout[i,2]+(rw/6))
  }
  
}

# Plot space with phylogeny mapped on ----
plot.PhyloSpace <- function(uout, tree, xlabel = "x-axis", ylabel = "yaxis", colour_genus = TRUE) {
  
  # Get species names from data to match to tree tips
  species <- rownames(uout)
  species <- strsplit(species, split = "_")
  species <- unique(sapply(species, function(vec) paste(vec[1:2], collapse = "_")))
  
  #Trim tree to match data (remove species that appear in tree but not data)
  # First find how many species to remove from tree (to tell user)
  treedrop <- length(setdiff(tree$tip.label, species))
  tree <- drop.tip(phy, tip = setdiff(tree$tip.label, species))
  
  # Trim data to match tree (remove species that appear in data but not tree)
  # First find how many specimens and species to remove from data (to tell user)
  species <- rownames(uout)
  species <- strsplit(species, split = "_")
  species <- sapply(species, function(vec) paste(vec[1:2], collapse = "_"))
  datdropspecimen <- length(rownames(uout)) - length(rownames(uout[species %in% tree$tip.label, ]))
  datdropspecies <- length(unique(species)) - length(unique(sapply(strsplit(rownames(uout[species %in% tree$tip.label, ]), split = "_"), function(vec) paste(vec[1:2], collapse = "_")))) # This massive line finds the species that are in the data but not the tree, then trims the data to match the tree. It's so complicated because the rownames are the specimen names, not the species names
  # Now trim the dataset to match the tree
  uout <- uout[species %in% tree$tip.label, ]
  # Inform user how many species trimmed from tree, and how many species and specimen images trimmed from data
  print(paste0(treedrop, " species trimmed from tree. ", datdropspecies, " species and ", datdropspecimen, " specimen images trimmed from data."))
  

  
  # Plot blank canvas
  phylomorphospace(tree = tree, X = uout, 
                   xlab = xlabel, ylab = ylabel, 
                   type = "n")
}

# Function to plot 2x2 plot for comparison of different PC axes, pixel values etc.
# Plots four of the basic plots in a grid
plot.FourPix <- function(obj1, obj2, obj3, obj4, filepath, asp = T){
  
  # Initialise png saving
  png(filepath, width = 10, height=15, units = "in", res = 600)
  
  # Set plotting parameters (2x2 grid)
  par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
  
  #Plot graphs
  objects <- list(obj1, obj3, obj2, obj4)
  for(i in 1:4){
    # Set axis labels based on input data (PCA or UMAP - the way this is currently set up, it assumes it's a UMAP
    # unless otherwise specified)
    obj <- objects[[i]]
    if(is.null(colnames(obj)[1])){
      xlabel <- "PCA UMAP axis 1"
      ylabel <- "PCA UMAP axis 2"
    }
    else if(substring(colnames(obj)[1], 1, 2) == "PC"){
      xlabel <- paste0("PC axis", substring(colnames(obj)[1], 3))
      ylabel <- paste0("PC axis", substring(colnames(obj)[2], 3))
    }
    else {
      xlabel <- colnames(obj)[1]
      ylabel <- colnames(obj)[2]
    }
    
    # Create plot
    plot.PixVals(uout = obj, xlabel = xlabel, ylabel = ylabel, asp = asp)
  }
  
  # End png saving
  dev.off()
}
