# Evaluate the best method to estimate the dimensionality of the patch colour pattern space ----
# Use the identified method to estimate dimensionality
# Based on the workflow proposed by Altan et al (2021), Comp Biol
# https://doi.org/10.1371/journal.pcbi.1008591
# Robert MacDonald
# 7th August 2024

# Load libraries ----
library(dplyr)
library(ggplot2)

# Clear environment
rm(list=ls())

# Load custom functions to estimate dimensionality ----
source(
  here::here(
    "3_SharedScripts", "estimate_dimensionality_functions.R"
  )
)

# EDITABLE CODE # ----
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neoaves"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# select sex of interest ("all", "male_female", "male_only", "female_only", "unknown_only")
sex_interest <- "male_female"
# Choose space to work with (usmldbl, usmldblr, xyz, xyzlum, xyzlumr, lab, cie, sRGB, hex, 
# jndxyz, jndxyzlum, jndxyzlumr)
space <- "lab"

# Load data ----

# Load in raw (pre-PCA) colourspaces and extract space of interest
pre_pca_space <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "1_RawColourSpaces", 
    paste(clade, sex_match, "patches.231030.prePCAcolspaces.rds", sep = ".")
  )
) %>% 
  magrittr::extract2(space)

# Load in PCA colourspace
pca_space <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    paste(clade, sex_match,  "patches.231030.PCAcolspaces.rds", sep = ".")
  )
) %>% 
  magrittr::extract2(space)

# Perform analyses ----

# Initialise results table
res <- data.frame(denoised = logical(),
                  method = character(), 
                  method_type = character(), 
                  cutoff = numeric())

## Evaluate linear algorithms ----

### Variance cutoff ----

# Choose variance cutoff (proportion, i.e. 0.9 equals a cutoff of 90%)
varcutoff_prop <- 0.9
var_cutoff <- length(
  summary(pca_space)$importance[3, summary(pca_space)$importance[3, ] < varcutoff_prop]
  ) + 1

# Print result
print(
  paste0(
    "The number of PCA axes to retain, determined by implementing a variance cutoff of ", 
      varcutoff_prop*100, "%, is ", var_cutoff
    )
  )

# Append result to results table
res <- rbind(res, data.frame(denoised = "FALSE",
                             method = "variance_cutoff", 
                             method_type = "linear", 
                             cutoff = var_cutoff))

### Participation ratio ----
partratio_cutoff <- part_ratio(pca_space)
print(
  paste0(
    "The number of PCA axes to retain, determined by participation ratio, is ", 
             partratio_cutoff
    )
  )
# Append result to results table
res <- rbind(res, data.frame(denoised = "FALSE",
                             method = "participation_ratio", 
                             method_type = "linear", 
                             cutoff = partratio_cutoff))

### Parallel analysis ----
# Note that for large sample size, simulated eigenvalues always tend to 1 and
# parallel analysis converges with the classic Kaiser-Guttman rule of retaining 
# the components with eigenvalues > 1 (Guttman, 1954; see Revelle, 2019)
paranal <- psych::fa.parallel(pre_pca_space, fa = "pc", n.iter = 200)
paranal_cutoff <- paranal$ncomp
print(
  paste0(
    "The number of PCA axes to retain, determined by parallel analysis, is ", 
    paranal_cutoff
  )
)
# Append result to results table
res <- rbind(res, data.frame(denoised = "FALSE",
                             method = "parallel_analysis", 
                             method_type = "linear", 
                             cutoff = paranal_cutoff))


## Ok, the above actually isn't necessary - what I need to do is determine if the mapping from
## low-dimensional latent space to high-dimensional full space is linear or non-linear
## I can do this by comparing the variance explained by a linear (PCA 90% variance cutoff) and 
## nonlinear (Joint Autoencoder) denoising model

# Asess how many components to keep using 90% variance cutoff
varcutoff_prop <- 0.9
var_cutoff <- length(
  summary(pca_space)$importance[3, summary(pca_space)$importance[3, ] < varcutoff_prop]
) + 1

# Apply 90% variance cutoff to PCA space
pca_90 <- list()
pca_90$sdev <- pca_space$sdev[1:var_cutoff]
pca_90$rotation <- pca_space$rotation[, 1:var_cutoff]
pca_90$center <- pca_space$center
pca_90$scale <- pca_space$scale
pca_90$x <- pca_space$x[, 1:var_cutoff]

# convert to prcomp class
class(pca_90) <- "prcomp"

# Reconstruct original data from reduced PCA
reconstructed_varcutoff_space <- pca_90$x %*% t(pca_90$rotation)

# add back center and scale, if applicable
if(!is.null(pca_90$center)) {
  reconstructed_varcutoff_space <- sweep(reconstructed_varcutoff_space, 2, pca_90$center, FUN = "+")
}
if(pca_90$scale == TRUE) {
  reconstructed_varcutoff_space <- sweep(reconstructed_varcutoff_space, 2, pca_90$scale, FUN = "*")
}

# inspect reconstructed data
head(reconstructed_varcutoff_space)

# see how far out the reconstruction is by subtracting reconstruction values from original
# values (no need to scale as original data is already scaled to max=1
perf_recon <- (pre_pca_space - reconstructed_varcutoff_space) %>% 
  mutate(
    specimen = rownames(.)
  ) %>% 
  tidyr::pivot_longer(
    cols = c(c(colnames(reconstructed_varcutoff_space))),
    names_to = "body_part",
    values_to = "r"
  )
ggplot(perf_recon, aes(x = r)) + 
  geom_histogram(bins = 40) + 
  facet_wrap(~ body_part)

# calculate proportion of variance in the original data that is captured in the reconstructed 
# data (1 - variance of residuals / variance of original data)
residuals <- pre_pca_space - reconstructed_varcutoff_space
var_orig <- sum(apply(pre_pca_space, 2, var))
var_resid <- sum(apply(residuals, 2, var))
prop_var <- 1 - (var_resid / var_orig)
prop_var

# Proportion of variance explained by a linear model is > 90% - I think this probably implies that 
# colour pattern space for our data is largely linear. This means I can probably use
# a linear dimensionality estimation method. Since Altan et al (2021) find that 
# parallel analysis is the best-performing linear method, I will apply this to the denoised
# (i.e. PCA-reconstructed) dataset to determine dimensionality
# Note that when it comes to looking at the ML embeddings, these are likely to be nonlinear
# so we will have to use a nonlinear dimensionality estimate method
# Note also that in the Altan et al paper they simulate datasets to be linear/nonlinear, and
# they deliberately introduce known amounts of noise. I don't think it is reasonable to 
# assume that noise in our dataset would be of a similar structure to that in their dataset
# (simulated neural firing rates/latent neural signals), so I don't think we should bother 
# with denoising (which in this case would likely just mean the variance cutoff method used here)
# What this means is that we should just use parallel analysis to determine the dimensionality
# of the non-denoised PCA space
# 

paranal_recon <- psych::fa.parallel(reconstructed_varcutoff_space, fa = "pc")
paranal_recon_cutoff <- paranal_recon$ncomp
print(
  paste0(
    "The final number of PCA axes to retain, determined by parallel analysis after denoising, is ", 
    paranal_recon_cutoff
  )
)

# This seems like a very circular method - I'm unsure this is really applicable here
# So from rereading the Altan et al (2021) paper, "Linear dimensionality estimation algorithms may 
# work well for linear datasets, but are likely to overestimate the dimensionality of a manifold 
# arising from a nonlinear mapping between the low-and high-dimensional spaces"
# Since our mapping between the low- and high-dimensional spaces is linear (PCA), there's no
# need to mess about with nonlinear dimensionality estimation algorithms or denoising

# For completeness, try the different linear methods on the denoised space
denoised_pca_space <- prcomp(reconstructed_varcutoff_space)


### Variance cutoff ----

# Choose variance cutoff (proportion, i.e. 0.9 equals a cutoff of 90%)
varcutoff_prop <- 0.9
var_cutoff <- length(
  summary(denoised_pca_space)$importance[3, summary(denoised_pca_space)$importance[3, ] < varcutoff_prop]
) + 1


# Append result to results table
res <- rbind(res, data.frame(denoised = "TRUE",
                             method = "variance_cutoff", 
                             method_type = "linear", 
                             cutoff = var_cutoff))

### Participation ratio ----
partratio_cutoff <- part_ratio(denoised_pca_space)
# Append result to results table
res <- rbind(res, data.frame(denoised = "TRUE",
                             method = "participation_ratio", 
                             method_type = "linear", 
                             cutoff = partratio_cutoff))

### Parallel analysis ----
# Note that for large sample size, simulated eigenvalues always tend to 1
paranal <- psych::fa.parallel(reconstructed_varcutoff_space, fa = "pc")
paranal_cutoff <- paranal$ncomp
# Append result to results table
res <- rbind(res, data.frame(denoised = "TRUE",
                             method = "parallel_analysis", 
                             method_type = "linear", 
                             cutoff = paranal_cutoff))

# Inspect results
res

# Plot results to further inspect
res %>% 
  ggplot(aes(x = method, y = cutoff, colour = method)) + 
  geom_point(size = 2) + 
  facet_wrap(~ denoised, nrow = 2)
