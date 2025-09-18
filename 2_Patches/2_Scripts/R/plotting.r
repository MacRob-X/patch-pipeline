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
  par(mfcol = c(2, 2), mar = c(4, 4, 1, 1))
  
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
plot.FourPatch <- function(obj1, obj2, obj3, obj4, filepath = NULL, asp = T, taxonCol = NULL){
  
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
      xlabel <- paste0("PC axis ", substring(colnames(obj)[1], 3))
      ylabel <- paste0("PC axis ", substring(colnames(obj)[2], 3))
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


# Plot patches space with colour grids
plot_patch_grids <- function(
    x_axis = NULL, y_axis = NULL,
    data_matrix = NULL,
    umap_obj = NULL, 
    prcomp_obj = NULL,
    colour_grid_path,
    row_names = NULL,
    x_label = NULL, y_label = NULL,
    asp_ratio = c("wrap", "square"),
    save_as = NULL, # currently only "png" and "plot_object" is implemented
    save_path = NULL,
    font_par = NULL,
    thin_number = NULL,
    fourplot = FALSE, fourplot_size = c(10, 10.5)
){
  
  # function to test if axes are atomic, numeric vectors
  is.anv <- function(axis){
    if(is.atomic(axis) & is.numeric(axis) & is.vector(axis)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  # Check only one data type has been provided - stop if not
  if(is.anv(x_axis) & is.anv(y_axis)){
    if(!is.null(data_matrix) | !is.null(umap_obj) | !is.null(prcomp_obj)){
      stop("More than one data type provided.")
    }
  } else if(is.character(x_axis) & is.character(y_axis)){
    if(all(is.null(data_matrix), is.null(umap_obj), is.null(prcomp_obj))){
      stop("No data provided.")
    } else if((!is.null(data_matrix) + !is.null(umap_obj) + !is.null(prcomp_obj)) >= 2){
      stop("More than one data type provided.")
    }
  }
  
  # Check what type of object has been provided and set up plot space
  
  # if axes have been provided as data columns
  if(is.anv(x_axis) & is.anv(y_axis)){
    
    # check if rownames have been provided and match length of data columns 
    # need the rownames to get the png colour grid names
    if(is.null(row_names)){
      stop("No rownames provided.")
    }
    
    plot_space <- cbind(x_axis, y_axis)
    rownames(plot_space) <- row_names
    
    x_axis <- "axis_1"; y_axis <- "axis_2"
    colnames(plot_space) <- c(x_axis, y_axis)
    
    # check if axis labels have been provided
    if(any(is.null(x_label), is.null(y_label))){
      warning("Axis labels not set. Setting to default labels.")
      # set to defaults if not provided
      x_label <- "Axis 1"; y_label <- "Axis 2"
    }
    
  } else if(is.character(x_axis) & is.character(y_axis)){
    # if axes have been provided as characters (e.g. "PC1", "PC2", "UMAP1", "UMAP2 etc)
    
    # Check if data has been provided as a prcomp object, a UMAP object, or as a data matrix
    # and extract data if necessary
    if(!is.null(prcomp_obj)){
      
      plot_space <- prcomp_obj %>% 
        magrittr::extract2("x") %>% 
        as.data.frame()
      
      # get proportions of variance for each axis
      pca_variance <- prcomp_obj %>% 
        summary() %>% 
        magrittr::extract2("importance") %>% 
        as.data.frame()
      
      
    } else if(!is.null(umap_obj)){
      
      plot_space <- umap_obj %>% 
        magrittr::extract2("layout") %>% 
        as.data.frame() %>% 
        rename(
          UMAP1 = V1, UMAP2 = V2
        )
      x_axis <- "UMAP1"; y_axis <- "UMAP2"
 
    } else if(!is.null(data_matrix)){
      
      plot_space <- data_matrix
      
      # check that named axes are in data matrix
      if(!(x_axis %in% colnames(plot_space) & y_axis %in% colnames(plot_space))){
        stop("Chosen axes are not columns in data_matrix.")
      }
      
    }
    
  } else if(any(is.null(x_axis), is.null(y_axis))){
    
    if(!is.null(umap_obj)){
      
      plot_space <- umap_obj %>% 
        magrittr::extract2("layout") %>% 
        as.data.frame() %>% 
        rename(
          UMAP1 = V1, UMAP2 = V2
        )
      x_axis <- "UMAP1"; y_axis <- "UMAP2"
      
    }
    
  }
  
  # if requested, thin to specific number of random species
  if(!is.null(thin_number)){
    
    all_spp <- sapply(strsplit(rownames(plot_space), split = "-"), "[", 1)
    
    unique_spp <- unique(sapply(strsplit(rownames(plot_space), split = "-"), "[", 1))
    
    spp_subset <- sample(unique_spp, thin_number, replace = FALSE)
    
    # subset plot space to only random species
    rows_to_keep <- which(all_spp %in% spp_subset)
    plot_space <- plot_space[rows_to_keep, ]
    
  }
  
  # set axis labels if not yet set
  if(!is.null(umap_obj)){
    x_label <- "UMAP1"
    y_label <- "UMAP2"
  } else if(!is.null(prcomp_obj)) {
    x_label <- paste0(x_axis, " (", round(pca_variance[, x_axis][2], 2) * 100, "% of variance)")
    y_label <- paste0(y_axis, " (", round(pca_variance[, y_axis][2], 2) * 100, "% of variance)")
  } else if(!is.null(data_matrix) & any(is.null(x_label), is.null(y_label))){
    warning("Axis labels not set. Setting to axis names.")
    # set to defaults if not provided
    x_label <- x_axis; y_label <- y_axis
  }
  
  
  # set up saving
  if(!is.null(save_as)){
    if(save_as == "png"){
      
      # set up aspect ratio
      if(asp_ratio == "wrap"){
        # determine axis with greatest range (to determine size of png)
        range_x_axis <- diff(range(plot_space[, x_axis]))
        range_y_axis <- diff(range(plot_space[, y_axis]))
        if(range_x_axis > range_y_axis){
          png_width <- 2500
          png_height <- png_width * (range_y_axis / range_x_axis)
        } else if(range_y_axis > range_x_axis){
          png_height <- 2500
          png_width <- png_height * (range_x_axis / range_y_axis)
        }
      } else if (asp_ratio == "square"){
        png_width <- 2500
        png_height <- 2500
      }
      
      # set up device
      if(is.null(font_par)){
        png(
          save_path,
          width = png_width, height = png_height,
          units = "px",
          pointsize = 40
        )
      } else {
        png(
          save_path,
          width = png_width, height = png_height,
          units = "px",
          pointsize = 40,
          family = font_par
        )
      }
      
    }
  }
  
  # fix aspect ratio as square, if requested
  if(asp_ratio == "square"){
    par(pty = "s")
  }
  
  dev.hold()
  
  plot(plot_space[, y_axis] ~ plot_space[, x_axis], 
       asp = T, 
       type = "n", 
       xlab = "", ylab = "", las = 1,
       # xaxt = "n", yaxt = "n"
       #     xlim = xrange,  ylim = yrange
  )
  title(xlab = x_label, ylab = y_label)
  
  box(lwd = 2)
  # if plotting as a single plot, define colour grid size relative to plot range
  if(fourplot == FALSE){
    rw <- diff(range(plot_space[,1]))/12
    rw_x <- rw / 15
    rw_y <- rw / 9
  } else if(fourplot == TRUE){
    # if plotting as part of a four-panel plot, define grid size as absolute relative to overall plot size
    x_range <- diff(range(plot_space[, x_axis]))
    grid_width_inches <- (fourplot_size[1] / 100) * 1.8
    grid_size_par <- (grid_width_inches / fourplot_size[1]) * x_range
    rw_x <- grid_size_par / 2
    rw_y <- rw_x * (5/3)
  }
  
  
  for(i in 1:nrow(plot_space)){
    fname <- rownames(plot_space)[i]
    fpng <- png::readPNG(paste0(
      colour_grid_path, "/",
      paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
    rasterImage(fpng, 
                xleft = plot_space[i, x_axis] - rw_x, 
                ybottom = plot_space[i, y_axis] - rw_y, 
                xright = plot_space[i, x_axis] + rw_x, 
                ytop = plot_space[i, y_axis] + rw_y)
    cat(paste0("\r", i, " of ", nrow(plot_space), " processed"))
    
    if(isTRUE(all.equal(i/500, as.integer(i/500)))){
      # Force an update to the graphics device after each 500 images
      # This helps prevent memory issues
      if(interactive()) {
        Sys.sleep(0.01)  # Small delay to allow graphics to update
      }
      
      # Free up memory
      gc()
    }
    
  }
  
  dev.flush()
  
  if(!is.null(save_as)){
    dev.off()
  }
  
  
}


# Plot four colour grid plots together (inherits from above colour grid plotting function)
plot_four_cg <- function(
  plot1_data, plot1_x_axis, plot1_y_axis,
  plot2_data, plot2_x_axis, plot2_y_axis,
  plot3_data, plot3_x_axis, plot3_y_axis,
  plot4_data, plot4_x_axis, plot4_y_axis,
  cg_path,
  save_type = "png",  # automatically plot as png - currently the only type supported
  write_folder = NULL,
  file_name,
  png_width = 10, png_height = 10.5, # in inches
  thin_number = NULL
){
  

  if(save_type == "png"){
    
    # if no write path, write to current working directory
    if(is.null(write_folder)){
      write_folder <- getwd()
    }
    
    # set up png plotting
    dpi <- 600
    full_path <- paste(write_folder, file_name, sep = "/")
    png(full_path, width = png_width, height = png_height, units = "in", res = dpi)
    
    # Set plotting parameters (2x2 grid)
    par(mfcol = c(2, 2), mar = c(4, 4, 1, 1))
    
    #Plot graphs
    for(plot_index in 1:4){
      
      # set plot object
      plot_obj <- get(paste0("plot", plot_index, "_data"))
      # set axes
      plot_x_axis <- get(paste0("plot", plot_index, "_x_axis"))
      plot_y_axis <- get(paste0("plot", plot_index, "_y_axis"))
      
      # check plot object class and plot accordingly
      if(class(plot_obj)[1] == "umap"){
        plot_patch_grids(x_axis = plot_x_axis, y_axis = plot_y_axis, umap_obj = plot_obj, colour_grid_path = cg_path, asp_ratio = "square", thin_number = thin_number, fourplot = TRUE, fourplot_size = c(png_width, png_height))
      } else if(class(plot_obj)[1] == "prcomp"){
        plot_patch_grids(x_axis = plot_x_axis, y_axis = plot_y_axis, prcomp_obj = plot_obj, colour_grid_path = cg_path, asp_ratio = "square", thin_number = thin_number, fourplot = TRUE, fourplot_size = c(png_width, png_height))
      } else if(class(plot_obj)[1] == "data.frame" | all(class(plot_obj) == c("matrix", "array"))){
        plot_patch_grids(x_axis = plot_x_axis, y_axis = plot_y_axis, data_matrix = plot_obj, colour_grid_path = cg_path, asp_ratio = "square", thin_number = thin_number, fourplot = TRUE, fourplot_size = c(png_width, png_height))
      }
      

    }
    
    # End png saving
    dev.off()
    
  }
  
  
  
}


# Function to generate heatmaps
# Make a function to plot the PCA loading heatmaps
# 
# Note - requires dplyr, ggplot, tidyterra and raster packages to be loaded

#' Plot PCA loading heatmaps
#' Heatmaps show how pixel/channel combinations are loaded on each PC axis
#' Note - function requires dplyr, ggplot and raster packages to be loaded
#'
#' @param pca_space A prcomp object containing the PCA colour pattern space
#' @param axes PC axes to plot heatmaps for
#' @param pixel_type What do the pixel values represent? One of "srgb", "vurgb", "vrgb", "urgb", "tcsxyz", "lab", "rawusml", "jndxyzlum", "jndxyzlumr", "tcsxyzlum", "tcsxyzlumr", "rawusmldbl", "rawusmldblr"
#' @param write_path Path and filename to write plot to, if saving as a PNG
#' @param rgb_colours Should the plot be plotted in pseudo-RGB colouring? Only possible for pixrl_type == "srgb"
#' @param col_scale_type Should the plot use a unified colour scale for all PC axes, or should each axis have its own colour scale? Individual colour scales for each channel of each axis are also supported
#' @param stack_channels If TRUE, take absolute of each patch loadings and stack channels to give absolute measure of patch important (not compatible with col_scale_type = "axis_channel")
#' @param label_stacked_axes If TRUE and stack_channels == TRUE, labels each stacked channel axis plot with the axis name (e.g. "PC1", "PC2", etc.). If FALSE, stacked channel plot will not have axis labels (useful when combining unstacked and stacked channel plots into a single figure)
#'
#' @return If write_path == NULL, returns a ggplot object of the heatmaps
#' @export
#'
#' @examples
patch_loading_heatmap <- function(pca_space, axes, pixel_type, write_path = NULL, rgb_colours = FALSE, col_scale_type = c("constant", "axis", "axis_channel"), stack_channels = FALSE, label_stacked_axes = TRUE){
  
  # check col_scale_type and stack_channels parameters are compatible
  if(col_scale_type == "axis_channel" & stack_channels == TRUE){
    
    col_scale_type <- "axis"
    
    warning("'axis_channel' colour scale type not compatible with stacked channel heatmap. Setting col_scale_type = 'axis'")
    
  }
  
  # define number of channels based on pixel type
  if(pixel_type == "srgb" | pixel_type == "vurgb" | pixel_type == "vrgb" | pixel_type == "urgb" |
     pixel_type == "tcsxyz" | pixel_type == "lab"){
    n_channels <- 3
  } else if(pixel_type == "rawusml" | pixel_type == "jndxyzlum" | pixel_type == "jndxyzlumr" | 
            pixel_type == "tcsxyzlum" | pixel_type == "tcsxyzlumr") {
    n_channels <- 4
  } else if(pixel_type == "rawusmldbl" | pixel_type == "rawusmldblr"){
    n_channels <- 5
  }
  
  # define channel names
  if(pixel_type == "srgb" | pixel_type == "vurgb" | pixel_type == "vrgb" | pixel_type == "urgb"){
    channel_names <- c("R", "G", "B")
  } else if(pixel_type == "tcsxyz"){
    channel_names <- c("x", "y", "z")
  } else if(pixel_type == "lab"){
    channel_names <- c("L", "a", "b")
  } else if(pixel_type == "rawusml"){
    channel_names <- c("U", "S", "M", "L")
  } else if(pixel_type == "jndxyzlum"){
    channel_names <- c("jndx", "jndy", "jndz", "jndlum")
  } else if(pixel_type == "jndxyzlumr"){
    channel_names <- c("jndx", "jndy", "jndz", "jndlumr")
  } else if(pixel_type == "tcsxyzlum"){
    channel_names <- c("x", "y", "z", "lum")
  } else if(pixel_type == "tcsxyzlumr"){
    channel_names <- c("x", "y", "z", "lumr")
  } else if(pixel_type == "rawusmldbl"){
    channel_names <- c("U", "S", "M", "L", "DBL")
  } else if(pixel_type == "rawusmldblr"){
    channel_names <- c("U", "S", "M", "L", "DBLr")
  }
  
  # create named lists of axes and channels
  names(axes) <- axes
  names(channel_names) <- channel_names
  
  # extract loadings from PCA space
  loadings <- pca_space %>% 
    magrittr::extract2(
      "rotation"
    ) 
  
  # subset to PCs of interest
  loadings <- as.data.frame(loadings[, axes])
  colnames(loadings) <- axes
  
  # set body part names
  body_parts <- c("thr", "cro", "bre", "nap", "bel", "man", "cov", "rum", "fli", "tai")
  # set numbers as names (for plotting later)
  names(body_parts) <- 1:length(body_parts)
  
  # add columns to specify channel and body part from rownames
  loadings <- loadings %>% 
    mutate(
      channel = sapply(strsplit(rownames(.), split = ".", fixed = TRUE), "[", 1),     # need backslash to escape full stop (as automatically recognised as regular expression). Need second backslash to escape first backslash
      body_part = sapply(strsplit(rownames(.), split = ".", fixed = TRUE), "[", 2)
    ) %>% 
    # set levels of body part (to ensure displayed correctly in grid)
    mutate(
      body_part = factor(body_part, levels = body_parts)
    )
  
  # add xy coordinates (used to transform to raster for easier plotting)
  # first get body part matrix in the correct layout - reference matrix to get xy coords
  bp_mat <- matrix(data = body_parts, nrow = 5, ncol = 2, byrow = TRUE)
  # reverse row order in matrix (otherwise grid will be plotted upside down)
  bp_mat <- apply(bp_mat, 2, rev)
  loadings_xy <- loadings %>% 
    mutate(
      x = sapply(body_part, function(val) which(bp_mat == val, arr.ind = TRUE)[2]),
      y = sapply(body_part, function(val) which(bp_mat == val, arr.ind = TRUE)[1])
    )
  # pivot longer
  loadings_xy_long <- loadings_xy %>% 
    tidyr::pivot_longer(
      cols = tidyr::starts_with("PC"),
      names_to = "PC",
      values_to = "loading"
    )
  
  # for non-stacked (individual channel) plot
  if(stack_channels == FALSE){
    
    if(col_scale_type == "constant" | col_scale_type == "axis"){
      ## Plot individual PCs
      ## For use with any pixel type
      
      # create a separate raster for each channel by subsetting
      # first initialise empty list
      n_axes <- length(axes)
      pc_spatrasters <- vector(mode = "list", length = n_axes)
      names(pc_spatrasters) <- axes
      for(axis in axes){
        # get data to use to create raster
        ras_dat <- loadings_xy_long[loadings_xy_long$PC == axis, ]
        # pivot wider to get channels in different layers
        ras_dat <- ras_dat %>% 
          tidyr::pivot_wider(
            names_from = "channel",
            values_from = "loading"
          )
        # order by body part
        ras_dat <- ras_dat[order(ras_dat$body_part), ]
        ras_dat <- ras_dat[, c("x", "y", channel_names)]
        # create raster and convert to spatraster
        pc_raster <- terra::rast(
          raster::rasterFromXYZ(ras_dat)
        )
        # assign to list
        pc_spatrasters[[axis]] <- pc_raster
        
      }
      
      # calculate colour scale
      
      if(col_scale_type == "constant"){
        # use grand min and max for constant colour scale
        #col_lims <- c(min(loadings_xy_long[["loading"]]), max(loadings_xy_long[["loading"]]))
        # use grand min/max - but make uniform about zero
        # find which is larger - absolute of grand min or grand max and set scale based on this
        grand_min <- min(loadings_xy_long[["loading"]])
        grand_max <- max(loadings_xy_long[["loading"]])
        col_max <- max(abs(grand_min), grand_max)
        col_lims <- c(col_max * -1, col_max)
        
      } else if(col_scale_type == "axis"){
        # # use min and max of individual axes for individual colour scales
        # col_lims <- c(min(loadings_xy_long[loadings_xy_long[["PC"]] == axis, "loading"]), max(loadings_xy_long[loadings_xy_long[["PC"]] == axis, "loading"]))
      }
      
      # get limits of x and y coordinates
      x_min <- min(loadings_xy_long$x) - 0.5
      x_max <- max(loadings_xy_long$x) + 0.5
      y_min <- min(loadings_xy_long$y) - 0.5
      y_max <- max(loadings_xy_long$y) + 0.5
      
      # Create plot for each spatraster
      
      # first initialise list to store plots
      pc_plots <- vector(mode = "list", length = n_axes)
      
      # add names of list elements
      names(pc_plots) <- axes
      
      # calculate aspect ratio
      asp <- length(unique(loadings_xy_long$y)) / length(unique(loadings_xy_long$x))
      
      # loop to create each raster
      for(axis in axes){
        
        # create axis plot
        p <- ggplot() + 
          tidyterra::geom_spatraster(data = pc_spatrasters[[axis]]) + 
          xlim(x_min, x_max) + 
          ylim(y_min, y_max) +
          # scale_fill_viridis_c(option = "plasma", limits = col_lims) + 
          facet_grid(cols = vars(lyr)) + 
          theme(aspect.ratio = asp,
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                plot.margin = margin(0, 0, 0, 0))
        
        
        
        
        
        # add colour scale
        if(col_scale_type == "constant"){
          p <- p + 
            scale_fill_viridis_c(option = "plasma", limits = col_lims) + 
            guides(fill = "none") # remove colour legend
        } else {
          p <- p + 
            scale_fill_viridis_c(option = "plasma")
        }
        
        # add to list
        pc_plots[[axis]] <- p
        
      }
      
      if(col_scale_type == "constant"){
        
        # create plot to extract legend
        leg_plot <- ggplot() +
          tidyterra::geom_spatraster(data = pc_spatrasters[[1]]) + 
          facet_wrap(~lyr) + 
          scale_fill_viridis_c(option = "plasma", limits = col_lims) +
          theme(legend.position = "right")
        # extract legend grob to pass to ggarrange
        leg_grob <- ggpubr::get_legend(leg_plot)
        
        # combine individual channel plots into one figure
        loadings_heatmap <- ggpubr::ggarrange(
          plotlist = pc_plots, 
          nrow = n_axes,
          # ncol = n_channels, 
          labels = axes, 
          # label.x = 0.5,
          # label.y = 0.2,
          align = "hv",
          widths = c(1, 1, 1),  # Adjust widths to minimize space
          heights = c(1, 1, 1),
          common.legend = TRUE,
          legend.grob = leg_grob,
          legend = "right"
        )
        
      } else if(col_scale_type == "axis"){
        
        # combine individual channel plots into one figure
        loadings_heatmap <- ggpubr::ggarrange(
          plotlist = pc_plots, 
          nrow = n_axes,
          # ncol = n_channels, 
          labels = axes, 
          # label.x = 0.5,
          # label.y = 0.2,
          align = "hv",
          widths = c(1, 1, 1),  # Adjust widths to minimize space
          heights = c(1, 1, 1),
          common.legend = FALSE,
          legend = "right"
        )
        
      }
      
      return(loadings_heatmap)
    }
    
    # For an individual legend for each channel of each axis
    if(col_scale_type == "axis_channel"){
      
      # Create a template raster
      template_rast <- terra::rast(bp_mat)
      
      # create rasters
      rasters <- lapply(
        axes, function(axis){
          
          lapply(
            channel_names, function(ch){
              
              dat_for_rast <- loadings_xy_long %>% 
                filter(PC == axis, channel == ch) %>% 
                arrange(body_part)
              
              # get values for raster
              vals <- dat_for_rast %>% 
                pull(loading) 
              
              # create raster
              rast <- dat_for_rast %>% 
                select(x, y) %>% 
                as.matrix() %>% 
                terra::rasterize(
                  y = template_rast,
                  values = rep(NA, times = nrow(.))
                )
              # assign values
              terra::values(rast) <- vals
              
              return(rast)
              
            }
          )
          
        }
      )
      
      # combine into a single raster
      rasters <- unlist(rasters, recursive = FALSE)
      heatmap_raster <- terra::rast(rasters)
      
      # plot
      # p <- ggplot() + 
      #   tidyterra::geom_spatraster(data = heatmap_raster) + 
      #   facet_wrap(~lyr, ncol = n_channels) + 
      #   scale_fill_viridis_c(option = "plasma") + 
      #   theme_minimal()
      p <- terra::plot(heatmap_raster, nc = length(channel_names), nr = length(axes), legend = TRUE, col = map.pal("plasma", 100))
      
      return(p)
      
    }
    
  } else if(stack_channels == TRUE){
    
    # for summation of absolute patch loadings across channels
    
    # create a separate raster for each channel by subsetting
    # first initialise empty list
    n_axes <- length(axes)
    pc_spatrasters <- vector(mode = "list", length = n_axes)
    names(pc_spatrasters) <- axes
    for(axis in axes){
      # get data to use to create raster
      ras_dat <- loadings_xy_long[loadings_xy_long$PC == axis, ]
      # pivot wider to get channels in different layers
      ras_dat <- ras_dat %>% 
        tidyr::pivot_wider(
          names_from = "channel",
          values_from = "loading"
        )
      # take the absolute and sum the individual absolute channel loadings
      ras_dat <- ras_dat %>% 
        mutate(
          absolute_loading = rowSums(abs(select(., all_of(unname(channel_names)))))
        ) %>% 
        select(
          - all_of(unname(channel_names))
        )
      # order by body part
      ras_dat <- ras_dat[order(ras_dat$body_part), ]
      ras_dat <- ras_dat[, c("x", "y", "absolute_loading")]
      # create raster and convert to spatraster
      pc_raster <- terra::rast(
        raster::rasterFromXYZ(ras_dat)
      )
      # assign to list
      pc_spatrasters[[axis]] <- pc_raster
      
    }
    
    # calculate colour scale
    
    if(col_scale_type == "constant"){
      all_vals <- lapply(pc_spatrasters, terra::values) %>% unlist()
      # use grand min and max for constant colour scale
      col_lims <- c(
        min(all_vals),
        max(all_vals)
      )
    }
    
    # get limits of x and y coordinates
    x_min <- min(loadings_xy_long$x) - 0.5
    x_max <- max(loadings_xy_long$x) + 0.5
    y_min <- min(loadings_xy_long$y) - 0.5
    y_max <- max(loadings_xy_long$y) + 0.5
    
    # Create plot for each spatraster
    
    # first initialise list to store plots
    pc_plots <- vector(mode = "list", length = n_axes)
    
    # add names of list elements
    names(pc_plots) <- axes
    
    # calculate aspect ratio
    asp <- length(unique(loadings_xy_long$y)) / length(unique(loadings_xy_long$x))
    
    # loop to create each raster
    for(axis in axes){
      
      # create axis plot
      p <- ggplot() + 
        tidyterra::geom_spatraster(data = pc_spatrasters[[axis]]) + 
        xlim(x_min, x_max) + 
        ylim(y_min, y_max) +
        # scale_fill_viridis_c(option = "plasma", limits = col_lims) + 
        facet_grid(cols = vars(lyr)) + 
        theme(aspect.ratio = asp,
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.margin = margin(0, 0, 0, 0))
      
      
      
      
      
      # add colour scale
      if(col_scale_type == "constant"){
        p <- p + 
          scale_fill_distiller(palette = "Reds", direction = 1, limits = col_lims) + 
          guides(fill = "none") # remove colour legend
      } else {
        p <- p + 
          scale_fill_distiller(palette = "Reds", direction = 1)
      }
      
      # add to list
      pc_plots[[axis]] <- p
      
    }
    
    if(col_scale_type == "constant"){
      
      # create plot to extract legend
      leg_plot <- ggplot() +
        tidyterra::geom_spatraster(data = pc_spatrasters[[1]]) + 
        facet_wrap(~lyr) + 
        scale_fill_distiller(palette = "Reds", direction = 1, limits = col_lims) +
        theme(legend.position = "right")
      # extract legend grob to pass to ggarrange
      leg_grob <- ggpubr::get_legend(leg_plot)
      
      # combine individual channel plots into one figure
      # if labelling axes
      if(label_stacked_axes == TRUE){
        loadings_heatmap <- ggpubr::ggarrange(
          plotlist = pc_plots, 
          nrow = n_axes,
          # ncol = n_channels, 
          labels = axes, 
          # label.x = 0.5,
          # label.y = 0.2,
          align = "hv",
          widths = c(1, 1, 1),  # Adjust widths to minimize space
          heights = c(1, 1, 1),
          common.legend = TRUE,
          legend.grob = leg_grob,
          legend = "right"
        )
      } else if(label_stacked_axes == FALSE){ # if not labelling axes
        loadings_heatmap <- ggpubr::ggarrange(
          plotlist = pc_plots, 
          nrow = n_axes,
          # ncol = n_channels, 
          # labels = axes, 
          # label.x = 0.5,
          # label.y = 0.2,
          align = "hv",
          widths = c(1, 1, 1),  # Adjust widths to minimize space
          heights = c(1, 1, 1),
          common.legend = TRUE,
          legend.grob = leg_grob,
          legend = "right"
        )
      }

    } else if(col_scale_type == "axis"){
      
      # combine individual channel plots into one figure
      # if labelling axes
      if(label_stacked_axes == TRUE){
        loadings_heatmap <- ggpubr::ggarrange(
          plotlist = pc_plots, 
          nrow = n_axes,
          # ncol = n_channels, 
          labels = axes, 
          # label.x = 0.5,
          # label.y = 0.2,
          align = "hv",
          widths = c(1, 1, 1),  # Adjust widths to minimize space
          heights = c(1, 1, 1),
          common.legend = FALSE,
          legend = "right"
        )
      } else if(label_stacked_axes == FALSE){
        loadings_heatmap <- ggpubr::ggarrange(
          plotlist = pc_plots, 
          nrow = n_axes,
          # ncol = n_channels, 
          # labels = axes, 
          # label.x = 0.5,
          # label.y = 0.2,
          align = "hv",
          widths = c(1, 1, 1),  # Adjust widths to minimize space
          heights = c(1, 1, 1),
          common.legend = FALSE,
          legend = "right"
        )
      }
      
      
    }
    
    return(loadings_heatmap)
    
  }
  
  
  
}