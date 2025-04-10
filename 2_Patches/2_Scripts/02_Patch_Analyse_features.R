
# Patch-wise PCA

rm(list=ls())

library(pavo)
library(ggplot2)
library(tidyverse)

#setwd("~/Dropbox/Projects/Current/Bird_colouration/Colour_pattern_analysis/Patches/")

# Custom functions ----
source("./2_Patches/2_Scripts/R/patch_plotting.R")

# ------------- #

## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neoaves"
# Restrict to only species for which we have male and female data?
mf_restrict <- TRUE

# Input data ----
# (this is the output from 01_Patch_Colourspace_mapping_vX.R)
px_filename <- paste0(clade, ".patches.231030.rawcolspaces.rds")
px <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "1_RawColourSpaces", 
    px_filename
  )
)

# Restrict to only males and females, if requested
if(mf_restrict == TRUE){
  
  # get species for which we have both male AND female data
  males <- unique(px[px[["sex"]] == "M", "species"])
  females <- unique(px[px[["sex"]] == "F", "species"])
  mf_species <- males[!is.na(match(males, females))]
  
  px <- px[px[["species"]] %in% mf_species, ]
  
}


# ------------- #

# Construct colour pattern spaces ----
# the below code doesn't actually change the values of anything, it just changes the format of the data frame to split into
# separate colour spaces and make the values for each dimension for each patch a separate column
# once this is done, the data will be ready to perform PCA to obtain the colour pattern spaces

# set up list of columns to extract from px
spaces_columns <- list(
  usml = c("species","sex","region", "u", "s", "m", "l"),
  usmlraw = c("species","sex","region", "uv", "sw", "mw", "lw"),
  usmldbl = c("species","sex","region","u","s","m", "l", "lum"),
  usmldblr = c("species","sex","region","u","s","m", "l", "lumr"),
  xyz = c("species","sex","region","x","y","z"),
  xyzlum = c("species","sex","region","x","y","z","dbl"),
  lab = c("species","sex","region","L","a","b"),
  ab = c("species","sex","region","a","b"),   # ab (Lab but without the L)
  L = c("species","sex","region","L"),   # L (Lab but only the L)
  cie = c("species","sex","region","cieX","cieY","cieZ"),
  sRGB = c("species","sex","region","sRGB.r","sRGB.g","sRGB.b"),
  hex = c("species","sex","region","hex"),
  jndxyz = c("species","sex","region","x.jnd","y.jnd","z.jnd"),
  jndxyzlum = c("species","sex","region","x.jnd","y.jnd","z.jnd","lum.jnd"),
  jndxyzlumr = c("species","sex","region","x.jndlumr","y.jndlumr","z.jndlumr","lum.jndlumr")
)

# put data into the correct format to perform PCA
# split into separate colour spaces and make the values for each dimension for each patch 
# a separate column
# once this is done, the data will be ready to perform PCA to obtain the colour pattern spaces

# define function to reshape data
reshape_data <- function(space_cols, pixel_data){
  
  dat <- reshape(px[, space_cols], idvar = c("species", "sex"), timevar = c("region"), direction = "wide")
  rownames(dat) <- matrix(apply(dat[,c("species","sex")], 1, paste, collapse='-'), ncol=1)
  return(dat[complete.cases(dat),-c(1:2)])
  
}

# apply function across all space
non_PCA_colourspaces <- lapply(spaces_columns, reshape_data, pixel_data = px)

# Have to do TCS xyzlumr (rescaled lum) separately as it involves rescaling the lum
# dimension
# TCS xyzlumr (rescaled)
vtcs <- data.frame(u=c(1,0,0,0), s=c(0,1,0,0), m=c(0,0,1,0), l=c(0,0,0,1))
d <- max(dist(tcspace(vtcs)[,c("x","y","z")])) # use this value to scale the lum variable (dbl) such that the distance between maximally different lum values (0 and 1) now equals d (= max distance in TCS)
tmp <- px
tmp$dbl <- tmp$dbl * d
dat <- reshape(tmp[,c("species","sex","region","x","y","z","dbl")], idvar = c("species","sex"), timevar = c("region"), direction = "wide")
rownames(dat) <- matrix(apply(dat[,c("species","sex")], 1, paste, collapse='-'), ncol=1)
non_PCA_colourspaces$xyzlumr <- dat[complete.cases(dat),-c(1:2)]




# set filename
if(mf_restrict == TRUE){
  prepca_filename <- paste(clade, "matchedsex", "patches.231030.prePCAcolspaces.rds", sep = ".")
} else{
  prepca_filename <- paste(clade, "allspecimens", "patches.231030.prePCAcolspaces.rds", sep = ".")
}


saveRDS(
  non_PCA_colourspaces,
  here::here(
    "2_Patches", "3_OutputData", "1_RawColourSpaces", 
    prepca_filename
  )
)

# Perform PCA ----

# set filename
if(mf_restrict == TRUE){
  prepca_filename <- paste(clade, "matchedsex", "patches.231030.prePCAcolspaces.rds", sep = ".")
} else{
  prepca_filename <- paste(clade, "allspecimens", "patches.231030.prePCAcolspaces.rds", sep = ".")
}

# reload the raw colour spaces (if necessary)
non_PCA_colourspaces <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "1_RawColourSpaces", 
    prepca_filename
  )
)


# perform PCA for each colour space to get colour pattern spaces (excluding hex space - can't PCA a non-numeric)
spaces <- c("usml", "usmlraw", "usmldbl", "usmldblr", "xyz", "xyzlum", "xyzlumr", "lab", "ab", "L", "cie", "sRGB", "jndxyz", "jndxyzlum", "jndxyzlumr")
pca_spaces <- lapply(spaces, function(space, spaces_list) prcomp(spaces_list[[space]]), spaces_list = non_PCA_colourspaces)
names(pca_spaces) <- spaces

# Save list of PCA colour pattern spaces as RDS file

# set filename
if(mf_restrict == TRUE){
  pca_filename <- paste(clade, "matchedsex", "patches.231030.PCAcolspaces.rds", sep = ".")
} else{
  pca_filename <- paste(clade, "allspecimens", "patches.231030.PCAcolspaces.rds", sep = ".")
}

# save
saveRDS(
  pca_spaces,
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA", 
    pca_filename
  )
)

# The scaled luminance version looks identical to the non-scaled luminance version but actually is ever so slightly different 
# - I've tracked through the whole process and the scaling is working correctly, it just makes essentially no difference to 
# the PCA'd data whether the luminance is scaled or not. Seems to make a very small difference to the UMAP data. Should run everything
# with both sets and check that qualitatively the results are unchanged

# reload PCA (if necessary)
pca_spaces <- readRDS("./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/1_Raw_PCA/Passeriformes.patches.231030.PCAcolspaces.rds")

# plot the first two PCs (just for initial inspection - main visualisation below)
plot(pca_spaces[["lab"]]$x)


# Calculate centroid dists (to see which species is farthest from 'average')

# first initialise dataframe to store results
n_spec <- nrow(pca_spaces[[1]]$x)
cdists <- as.data.frame(matrix(ncol = length(pca_spaces), nrow = n_spec))
colnames(cdists) <- names(pca_spaces)
rownames(cdists) <- rownames(pca_spaces[[1]]$x)

# # calculate centroid distances for species of interest (i.e. each row) in each space
for (i in 1:n_spec) {
  
  cat("\r", i)
  # calculate centroid distances for species of interest in each space
  cdists[i, ] <- unlist(lapply(spaces, function(space, pca_list){
    dist(rbind(rep(0, ncol(pca_list[[space]]$x)), pca_list[[space]]$x[i,]))
  }, pca_list = pca_spaces))
 
}

# get species/sex which are maximally distant in each colour space
rownames(cdists)[cdists$usml == max(cdists$usml)]
rownames(cdists)[cdists$usmldbl == max(cdists$usmldbl)]
rownames(cdists)[cdists$usmldblr == max(cdists$usmldblr)]
rownames(cdists)[cdists$xyz == max(cdists$xyz)]
rownames(cdists)[cdists$xyzlum == max(cdists$xyzlum)]
rownames(cdists)[cdists$xyzlumr == max(cdists$xyzlumr)]
rownames(cdists)[cdists$lab == max(cdists$lab)]
rownames(cdists)[cdists$cie == max(cdists$cie)]
rownames(cdists)[cdists$sRGB == max(cdists$sRGB)]
rownames(cdists)[cdists$jndxyz == max(cdists$jndxyz)]
rownames(cdists)[cdists$jndxyzlum == max(cdists$jndxyzlum)]
rownames(cdists)[cdists$jndxyzlumr == max(cdists$jndxyzlumr)]


# perform UMAP to visualise colour space in two dimensions
umap.xyz <- umap::umap(pca_spaces[["xyz"]]$x)
umap.xyzlum <- umap::umap(pca_spaces[["xyzlum"]]$x)
umap.xyzlumr <- umap::umap(pca_spaces[["xyzlumr"]]$x)
umap.lab <- umap::umap(pca_spaces[["lab"]]$x)
umap.jndxyzlum <- umap::umap(pca_spaces[["jndxyzlum"]]$x)
umap.jndxyzlumr <- umap::umap(pca_spaces[["jndxyzlumr"]]$x)

plot(umap.xyz$layout)
plot(umap.xyzlum$layout)
plot(umap.xyzlumr$layout)   # This looks almost exactly the same as the tcs xyzlum (non-scaled) umap plot but they are very slightly different
plot(umap.lab$layout)
plot(umap.jndxyzlum$layout)
plot(umap.jndxyzlumr$layout)


# Plot PCA axes and UMAP

# choose space of interest
pca_space <- pca_spaces[["lab"]]

plot.FourPatch(obj1 = pca_space$x[, 1:2], 
             obj2 = pca_space$x[, 3:4],
             obj3 = pca_space$x[, 5:6],
             obj4 = umap.lab$layout,
             filepath = "./Plots/featurespaces/patches/patches.lab.pca1-6.umap.back.png",
             asp = T)

# Plot PC axes 1-8 with associated variance proportions
plot.FourPatch(obj1 = pca_space$x[, 1:2], 
               obj2 = pca_space$x[, 3:4],
               obj3 = pca_space$x[, 5:6],
               obj4 = pca_space$x[, 7:8],
               filepath = "./Plots/featurespaces/patches/patches.lab.pca1-8.png",
               PCobj = pca_space,
               asp = T)

# Plot PC axes 9-16 with associated variance proportions
plot.FourPatch(obj1 = pca_space$x[, 9:10], 
               obj2 = pca_space$x[, 11:12],
               obj3 = pca_space$x[, 13:14],
               obj4 = pca_space$x[, 15:16],
               filepath = "./Plots/featurespaces/patches/patches.lab.pca9-16.png",
               PCobj = pca_space,
               asp = T)

# Plot PC axes 17-24 with associated variance proportions
plot.FourPatch(obj1 = pca_space$x[, 17:18], 
               obj2 = pca_space$x[, 19:20],
               obj3 = pca_space$x[, 21:22],
               obj4 = pca_space$x[, 23:24],
               filepath = "./Plots/featurespaces/patches/patches.lab.pca17-24.png",
               PCobj = pca_space,
               asp = T)

# Plot histogram of density for each PC axis (to see spread in each axis)
par(mfrow = c(4, 2), mar = c(4, 4, 1, 1))
for(i in 1:8){    # axis 1-8
  hist(pca.jndxyzlumr$x[, i], breaks = 50,
       xlab = paste0(colnames(pca.jndxyzlumr$x)[i], " (", round(summary(pca.jndxyzlumr)$importance[2, i] * 100, digits = 1), "%)"),
       ylab = "No. of occurrences",
       main = "")
}

# using ggplot2

p <- ggplot(pca.jndxyzlumr$x, aes(x = PC1)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "grey") + 
  theme_grey()
p

# convert PC data to long format tidy data so can plot multiple histograms
# first 8 PCs only
pcLong <- data.frame(pca.jndxyzlumr$x[, 1:8])
# add separate species and sex columns
  for(i in 1:nrow(pcLong)){
    cat("\r", i)
    split <- strsplit(rownames(pcLong)[i], split = "-")[[1]]
    pcLong$species[i] <- split[1]
    pcLong$sex[i] <- split[2]
  }

pcLong <- pivot_longer(pcLong, cols = starts_with("PC"), names_to = "PC_axis", values_to = "value")

# draw histograms

# males, females and unknowns all together
p <- ggplot(pcLong, aes(x = value)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "grey") + 
  theme_grey() + 
  facet_wrap(~ PC_axis, nrow = 4, ncol = 2, scales = "free")
p

# filtered by sex
p <- pcLong %>% 
  filter(sex == "M") %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "grey") + 
  theme_grey() + 
  facet_wrap(~ PC_axis, nrow = 4, ncol = 2, scales = "free")
p

# with both male and female data overlaid
p <- pcLong %>% 
  ggplot(aes(x = value)) +
  geom_histogram(data = pcLong %>% filter(sex == "M"), bins = 50, fill = "skyblue", color = "grey", alpha = 1) + 
  geom_histogram(data = pcLong %>% filter(sex == "F"), bins = 50, fill = "pink", color = "grey", alpha = 0.5) + 
  theme_grey() + 
  facet_wrap(~ PC_axis, nrow = 4, ncol = 2, scales = "free")
p


# ============================= #

# Colour points by family

# create dataframe from UMAP to manipulate
plotData.umap.jndxyzlumr <- data.frame(umap.jndxyzlumr$layout)
l <- length(umap.jndxyzlumr$layout[,1])
for(i in 1:l){
  split <- strsplit(rownames(plotData.umap.jndxyzlumr)[i], split = "-")[[1]]
  plotData.umap.jndxyzlumr$species[i] <- split[1]
  plotData.umap.jndxyzlumr$sex[i] <- split[2]
}


# load taxonomy
taxo <- read.csv("./Data/BLIOCPhyloMasterTax_2016_02_29.csv", strings = F)


# match to species in dataset
taxo <- taxo[taxo$TipLabel %in% px$species, ]

# reorder to match dataset
taxo <- taxo[match(plotData.umap.jndxyzlumr$species, taxo$TipLabel), ]

# create colour mapping
library(pals)
cols <- polychrome((length(unique(taxo$IOCOrder))))
col_mapping <- cols[match(taxo$IOCOrder, unique(taxo$IOCOrder))]

# plot
plot(umap.jndxyzlumr$layout, 
     col = col_mapping,
     xlab = "UMAP 1", ylab = "UMAP 2")
legend("bottomright", 
       legend = levels(as.factor(taxo$IOCOrder)), 
       fill = cols, 
       title = "Order",
       cex = 0.5, pt.cex = 0.5,)
text(umap.jndxyzlumr$layout[,1], 
     umap.jndxyzlumr$layout[,2], 
     labels = taxo$BLFamilyLatin, 
     pos = 1, cex = 0.5)




#-------------------------------------------------------------------------------#
# Use PCA to extrapolate patch values for every value in range of original pixel values
# I'm thinking of this as analogous to creating 'impossible' jaw shapes in morphospace contructions

# Based on https://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com

# get original data
px <- readRDS("./Outputs/features/patches/Passeriformes.patches.231030.processed.240322.rds")

cps.jndxyzlumr <- data.frame(px$x.jndlumr, px$y.jndlumr, px$z.jndlumr, px$lum.jndlumr, 
                             row.names = paste(px$species, px$sex, sep = "-"))

mu = colMeans()
