# Functions for visualising PCA- and CNN-generated colour pattern spaces
# 29 February 2024
# Robert MacDonald


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

#-------------------------------------------------------------------------------------#

# Plot space with phylogeny mapped on ----
plot.PhyloSpace <- function(uout, tree, xlabel = NULL, ylabel = NULL, sex = "all", tiplabs = "off") {
  
  library(phytools)  # for plotting phylomorphospaces
  library(pals)  # for colour palettes
  
  # Trim data to just males or just females, if necessary
  if(sex == "M"){
    specnames <- rownames(uout)
    specnames <- strsplit(specnames, split = "_")
    # convert to dataframe and transpose
    specnames <- data.frame(t(data.frame(specnames)))   # this is ugly - probably a more elegant way to do this
    uout <- uout[which(specnames$X4 == "M"), ]
  }
  else if(sex == "F"){
    specnames <- rownames(uout)
    specnames <- strsplit(specnames, split = "_")
    # convert to dataframe and transpose
    specnames <- data.frame(t(data.frame(specnames)))   # this is ugly - probably a more elegant way to do this
    uout <- uout[which(specnames$X4 == "F"), ]
  }
  else if(sex != "all"){
    stop("Please provide a valid parameter for 'sex' ('M'/'F'/'all')")
    
  }
  
  # Calculate average position for each species
  # Get list of unique species names - note that if M or F is selected, this step becomes redundant, so consider
  # restructuring to improve efficiency
  species <- rownames(uout)
  species <- strsplit(species, split = "_")
  rownames(uout) <- sapply(species, function(vec) paste(vec[1:2], collapse = "_"))
  uniquespecies <- unique(rownames(uout))
  # Match positions of each species name to the specimens in the data
  for(i in 1: length(uout[, 1])){
    # find row numbers of all images of species
    rownums <- which(rownames(uout) == uniquespecies[i])
    # average X and Y positions
    avgX <- sum(uout[, 1][rownums]) / length(uout[, 1][rownums])
    avgY <- sum(uout[, 2][rownums]) / length(uout[, 2][rownums])
    # make values for all instances of species the averages
    uout[, 1][rownums] <- c(rep(avgX, times = length(rownums)))
    uout[, 2][rownums] <- c(rep(avgY, times = length(rownums)))
  }
  # Remove duplicates
  uout <- unique(uout)
  

  #Trim tree to match data (remove species that appear in tree but not data)
  # First find how many species to remove from tree (to tell user)
  treedrop <- length(setdiff(tree$tip.label, uniquespecies))
  tree <- drop.tip(phy, tip = setdiff(tree$tip.label, uniquespecies))
  
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
  
  # Create vector to pass genus name to phylomorphospace() to colour code points
  genera <- sub("_.*", "", rownames(uout))
  names(genera) <- names(uout[,1])
  uniquegenera <- unique(genera)
  
  # Create colour palette
  polychrome <- polychrome()       # from pals package
  tip.cols <- genera
  for(i in 1:length(uniquegenera)){
    tip.cols[tip.cols == uniquegenera[i]] <- polychrome[i]
  }

  
  cols <- c(tip.cols[tree$tip.label], rep("black", tree$Nnode))
  names(cols) <- 1:(length(tree$tip) + tree$Nnode)
  
  # Get x and y axis labels
  if(is.null(xlabel)){
    xlabel <- colnames(uout)[1]
  }
  if(is.null(ylabel)){
    ylabel <- colnames(uout)[2]
  }
  
  
  # Plot phylomorphospace
  phylomorphospace(tree = tree, X = uout, 
                   xlab = xlabel, ylab = ylabel, 
                   control = list(col.node = cols),
                   node.size = c(0, 1.2),
                   label = tiplabs)
  # legend("bottomright", uniquegenera, pt.bg = cols, horiz = TRUE, bty = "n",
  #        pt.cex = 1.5, pch = 21)
  
  # Add convex hulls
  # cols <- setNames(polychrome(n = length(uniquegenera)), sort(unique(getStates(tree))))
  # COLS <- cols[getStates(tree, "tips")]
  # for(i in 1:Ntip(tree)){
  #   label <- tree$tip.label[i]
  #   ii <- which(rownames(genera) == label)
  #   if(length(ii) > 1){
  #     hull <- chull(uout[ii, ])
  #     polygon(uout[ii, ][hull, ], col = COLS[i])
  #   }
  # }
}

#-------------------------------------------------------------------------------------#

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


#-------------------------------------------------------------------------------------#


# Plot four basic patch plots together
plot.FourPatch <- function(obj1, obj2, obj3, obj4, filepath = NULL, asp = T, taxonCol = F, PCobj = NULL){
  
  # If PC object specified, get proportions of variance of each axis
  # Note that for this to work, PCobj must be specified as a prcomp object
  if(!is.null(PCobj)){
    propVar <- round(summary(PCobj)$importance[2, ] * 100, digits = 1)
  }

  
  # Initialise png saving if requested
  if(!is.null(filepath)){
    png(filepath, width = 15, height=15, units = "in", res = 600)
  }
  
  
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
      # get PC axis variance if PCobj specified
      if(!is.null(PCobj)){
        xPropVar <- propVar[as.numeric(substring(colnames(obj)[1], 3))]
        yPropVar <- propVar[as.numeric(substring(colnames(obj)[2], 3))]
        xlabel <- paste0("PC axis ", substring(colnames(obj)[1], 3), " (", xPropVar, "%)")
        ylabel <- paste0("PC axis ", substring(colnames(obj)[2], 3), " (", yPropVar, "%)")
      }
      else {
        xlabel <- paste0("PC axis ", substring(colnames(obj)[1], 3))
        ylabel <- paste0("PC axis ", substring(colnames(obj)[2], 3))
        
      }
    }
    else {
      xlabel <- colnames(obj)[1]
      ylabel <- colnames(obj)[2]
    }
    
    # Set colour mapping for taxon subgroups
    if(taxonCol == TRUE){
      cols <- rainbow(length(levels(as.factor(obj1$taxon_subgroup))))
      col_mapping <- cols[match(obj1$taxon_subgroup, levels(as.factor(obj1$taxon_subgroup)))]
    }
    else {
      col_mapping <- c(rep("black", times = length(obj[,1])))
    }
    
    # Create plot
    # If asp not selected, set x and y axis limits and plot blank canvas
    if(asp != T){
      # Set x and y limits
      xlimit <- c(round(min(obj[, 1]) - 1, digits = 0), round(max(obj[, 1]) + 1, digits = 0))
      ylimit <- c(round(min(obj[, 2]) - 1, digits = 0), round(max(obj[, 2]) + 1, digits = 0))
      # Plot
      plot(obj[,2] ~ obj[,1], 
           xlim = xlimit, ylim = ylimit,
           xlab = xlabel, ylab = ylabel, 
           las = 1,
           col = col_mapping); box(lwd=2)
    }
    # if asp selected, plotwith asp = T (aspect ratio fixed so x and y axis have same scale
    # - this makes the bird pngs the correct shape plus is standard for plotting PCAs etc)
    else if(asp == T){
      plot(obj[,2] ~ obj[,1],
           asp = T,
           xlab = xlabel, ylab = ylabel, 
           las = 1,
           col = col_mapping); box(lwd=2)
    }
  }
  
  # end png saving if required
  if(!is.null(filepath)){
    dev.off()
  }
}


