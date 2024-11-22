# A selection of functions to compute custom metrics of trait space diversity to be implemented 
# using the dispRity package
# Robert MacDonald
# 7th August 2024




#' Calculate mean distance to specified number of nearest neighbours in trait space (knn)
#' 
#' This function calculates the mean distance to a given number of nearest neighbours within a trait space. If species and sex are provided within the rownames (e.g. "Abeillia_abeillei-F") then this function will automatically exclude the distance between the male and female of a species (as this could artificially decrease the mean nearest neighbour distance). If no rownames are provided, the function will include all points in trait space.
#' @param matrix  matrix containing the trait space embeddings/coordinates
#' @param rownames set matrix rownames, if provided
#' @param nn number of nearest neighbours to examine (k). Default is all nearest neighbours (i.e. mean distance to all other points in trait space)
#' @param method type of distance to use - default is Euclidean
#'
#' @return returns a numreic vector (I think) containing the mean distance to k nearest nighbours for each specified point in trait space (i.e. each phenotype)
#' @export
#'
#' @examples implement using dispRity
#' e.g. nn_5 <- as.data.frame(dispRity::dispRity(space_embeddings, metric = mean.nn.dist, nn = 5, rownames = rownames(space_embeddings))$disparity[[1]][[1]])


mean.nn.dist <- function(matrix, rownames = NULL, nn = 5, method = "euclidean") {
  
  ## Set number of neighbours to all if not specified
  if(is.null(nn)){
    nn <- nrow(matrix)
  }
  
  ## Set matrix rownames if specified
  if(!is.null(rownames)){
    rownames(matrix) <- rownames
  }
  
  ## Calculate all pairwise distances
  pair.dists <- as.matrix(vegan::vegdist(matrix, method = method))
  
  ## Function to calculate nn closest for a single point
  nn_point <- function(one.row, nn, row.index){
    
    # if no species/sex provided just calculate the basic mean nearest neighbour distance
    if(is.null(names(one.row))){
      nn_point_dists <- mean(sort(one.row)[2:(nn + 1)])
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
        nn_point_dists <- mean(sort(one.row)[2:(nn + 1)])
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




#' Count number of neighbours within a given radius
#'
#' @param matrix matrix containing the trait space embeddings/coordinates
#' @param radius radius within which to count neighbours (defaults to 1). This can be a number or you can supply a function with which to calculate the radius to use. The function will be applied to a distance matrix generated from the input matrix.
#' @param relative make counts relative to the number of observations? defaults to TRUE
#' @param method type of distance to use - defaults to Euclidean
#'
#' @return
#' @export
#'
#' @examples
count.neighbours <- function(matrix, radius = 1, relative = TRUE, method = "euclidean") {
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

