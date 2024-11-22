# Calculate Participation Ratio of PCA
# 
# Provides a principled way of finding a variance threshold when the ground truth is not known
# i.e. it tells you a reasonable place to cut off your PCA axes without a priori 
# specifying a variance threshold
# 
# Based on https://doi.org/10.1101/214262
# 
# Robert MacDonald
# 7th August 2024

part_ratio <- function(pca_object) {
  
  # Extract eigenvalues from prcomp object
  eigenvals <- pca_object$sdev ^ 2
  
  # Calculate participation ratio
  part_rat <- (sum(eigenvals) ^ 2) / sum(eigenvals ^ 2)
  
  return(part_rat)
  
}
