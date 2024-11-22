
# Patch-wise PCA

rm(list=ls())

library(pavo)
library(ggplot2)
library(tidyverse)

#setwd("~/Dropbox/Projects/Current/Bird_colouration/Colour_pattern_analysis/Patches/")

# Custom functions ----
source("./2_Patches/2_Scripts/2_BetaVersions/Patch_Plotting_functions_v1.r")

# ------------- #

# Input data ----
# (this is the output from 01_Patch_Colourspace_mapping_vX.R)

px <- readRDS("./2_Patches/3_OutputData/1_RawColourspaces/Neoaves.patches.231030.rawcolspaces.processed.240603.rds")


# ------------- #

# Construct colour pattern spaces ----
# the below code doesn't actually change the values of anything, it just changes the format of the data frame to split into
# separate colour spaces and make the values for each dimension for each patch a separate column
# once this is done, the data will be ready to perform PCA to obtain the colour pattern spaces

# TCS xyz
dat <- reshape(px[,c("species","sex","region","x","y","z")], idvar = c("species","sex"), timevar = c("region"), direction = "wide")
rownames(dat) <- matrix(apply(dat[,c("species","sex")], 1, paste, collapse='-'), ncol=1)
cps.xyz <- dat[complete.cases(dat),-c(1:2)]

# TCS xyzlum
dat <- reshape(px[,c("species","sex","region","x","y","z","dbl")], idvar = c("species","sex"), timevar = c("region"), direction = "wide")
rownames(dat) <- matrix(apply(dat[,c("species","sex")], 1, paste, collapse='-'), ncol=1)
cps.xyzlum <- dat[complete.cases(dat),-c(1:2)]

# TCS xyzlumr (rescaled)
vtcs <- data.frame(u=c(1,0,0,0), s=c(0,1,0,0), m=c(0,0,1,0), l=c(0,0,0,1))
d <- max(dist(tcspace(vtcs)[,c("x","y","z")])) # use this value to scale the lum variable (dbl) such that the distance between maximally different lum values (0 and 1) now equals d (= max distance in TCS)
tmp <- px
tmp$dbl <- tmp$dbl * d
dat <- reshape(tmp[,c("species","sex","region","x","y","z","dbl")], idvar = c("species","sex"), timevar = c("region"), direction = "wide")
rownames(dat) <- matrix(apply(dat[,c("species","sex")], 1, paste, collapse='-'), ncol=1)
cps.xyzlumr <- dat[complete.cases(dat),-c(1:2)]

# Lab
dat <- reshape(px[,c("species","sex","region","L","a","b")], idvar = c("species","sex"), timevar = c("region"), direction = "wide")
rownames(dat) <- matrix(apply(dat[,c("species","sex")], 1, paste, collapse='-'), ncol=1)
cps.lab <- dat[complete.cases(dat),-c(1:2)]

# CIE
dat <- reshape(px[,c("species","sex","region","cieX","cieY","cieZ")], idvar = c("species","sex"), timevar = c("region"), direction = "wide")
rownames(dat) <- matrix(apply(dat[,c("species","sex")], 1, paste, collapse='-'), ncol=1)
cps.cie <- dat[complete.cases(dat),-c(1:2)]

# sRGB
dat <- reshape(px[,c("species","sex","region","sRGB.r","sRGB.g","sRGB.b")], idvar = c("species","sex"), timevar = c("region"), direction = "wide")
rownames(dat) <- matrix(apply(dat[,c("species","sex")], 1, paste, collapse='-'), ncol=1)
cps.srgb <- dat[complete.cases(dat),-c(1:2)]

# hex
dat <- reshape(px[,c("species","sex","region","hex")], idvar = c("species","sex"), timevar = c("region"), direction = "wide")
rownames(dat) <- matrix(apply(dat[,c("species","sex")], 1, paste, collapse='-'), ncol=1)
cps.hex <- dat[complete.cases(dat),-c(1:2)]

# jnd xyz
dat <- reshape(px[,c("species","sex","region","x.jnd","y.jnd","z.jnd")], idvar = c("species","sex"), timevar = c("region"), direction = "wide")
rownames(dat) <- matrix(apply(dat[,c("species","sex")], 1, paste, collapse='-'), ncol=1)
cps.jndxyz <- dat[complete.cases(dat),-c(1:2)]

# jnd xyzlum
dat <- reshape(px[,c("species","sex","region","x.jnd","y.jnd","z.jnd","lum.jnd")], idvar = c("species","sex"), timevar = c("region"), direction = "wide")
rownames(dat) <- matrix(apply(dat[,c("species","sex")], 1, paste, collapse='-'), ncol=1)
cps.jndxyzlum <- dat[complete.cases(dat),-c(1:2)]


# jnd xyzlumr
dat <- reshape(px[,c("species","sex","region","x.jndlumr","y.jndlumr","z.jndlumr","lum.jndlumr")], idvar = c("species","sex"), timevar = c("region"), direction = "wide")
rownames(dat) <- matrix(apply(dat[,c("species","sex")], 1, paste, collapse='-'), ncol=1)
cps.jndxyzlumr <- dat[complete.cases(dat),-c(1:2)]

# collate into one list for saving
# these dataframes contain the individual patch colourspaces (e.g. first three columns of first dataframe together make 
# up the belly patch luminance-free tetrahedral colour space)
non_PCA_colourspaces <- list(xyz = cps.xyz,
                             xyzlum = cps.xyzlum,
                             xyzlumr = cps.xyzlumr,
                             lab = cps.lab,
                             cie = cps.cie,
                             sRGB = cps.srgb,
                             hex = cps.hex,
                             jndxyz = cps.jndxyz,
                             jndxyzlum = cps.jndxyzlum,
                             jndxyzlumr = cps.jndxyzlumr)


saveRDS(non_PCA_colourspaces, "./2_Patches/3_OutputData/1_RawColourspaces/Neoaves.patches.231030.prePCAcolspaces.240603.rds")


# Perform PCA ----

# reload the raw colour spaces (if necessary)
non_PCA_colourspaces <- readRDS("./2_Patches/3_OutputData/1_RawColourspaces/Neoaves.patches.231030.prePCAcolspaces.240603.rds")

cps.xyz <- non_PCA_colourspaces$xyz
cps.xyzlum <- non_PCA_colourspaces$xyzlum
cps.xyzlumr <- non_PCA_colourspaces$xyzlumr
cps.lab <- non_PCA_colourspaces$lab
cps.cie <- non_PCA_colourspaces$cie
cps.srgb <- non_PCA_colourspaces$sRGB
cps.hex <- non_PCA_colourspaces$hex
cps.jndxyz <- non_PCA_colourspaces$jndxyz
cps.jndxyzlum <- non_PCA_colourspaces$jndxyzlum
cps.jndxyzlumr <- non_PCA_colourspaces$jndxyzlumr

# perform PCA for each colour space to get colour pattern spaces (excluding hex space - can't PCA a non-numeric)
# TCS xyz
pca.xyz <- prcomp(non_PCA_colourspaces$xyz)
summary(pca.xyz)
saveRDS(pca.xyz, "./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/1_Raw_PCA/Neoaves.patches.231030.PCAcolspaces.xyz.240603.rds")

# xyzlum
pca.xyzlum <- prcomp(non_PCA_colourspaces$xyzlum)
summary(pca.xyzlum)
saveRDS(pca.xyzlum, "./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/1_Raw_PCA/Neoaves.patches.231030.PCAcolspaces.xyzlum.240603.rds")

# xyzlumr
pca.xyzlumr <- prcomp(non_PCA_colourspaces$xyzlumr)
summary(pca.xyzlumr)
saveRDS(pca.xyzlumr, "./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/1_Raw_PCA/Neoaves.patches.231030.PCAcolspaces.xyzlumr.240603.rds")

# lab
pca.lab <- prcomp(non_PCA_colourspaces$lab)
summary(pca.lab)
saveRDS(pca.lab, "./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/1_Raw_PCA/Neoaves.patches.231030.PCAcolspaces.lab.240603.rds")

# cie
pca.cie <- prcomp(non_PCA_colourspaces$cie)
summary(pca.cie)
saveRDS(pca.cie, "./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/1_Raw_PCA/Neoaves.patches.231030.PCAcolspaces.cie.240603.rds")

# srgb
pca.srgb <- prcomp(non_PCA_colourspaces$sRGB)
summary(pca.srgb)
saveRDS(pca.srgb, "./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/1_Raw_PCA/Neoaves.patches.231030.PCAcolspaces.srgb.240603.rds")

# jndxyz
pca.jndxyz <- prcomp(non_PCA_colourspaces$jndxyz)
summary(pca.jndxyz)
saveRDS(pca.jndxyz, "./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/1_Raw_PCA/Neoaves.patches.231030.PCAcolspaces.jndxyz.240603.rds")

# jndxyzlum
pca.jndxyzlum <- prcomp(non_PCA_colourspaces$jndxyzlum)
summary(pca.jndxyzlum)
saveRDS(pca.jndxyzlum, "./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/1_Raw_PCA/Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603.rds")

# jndxyzlumr
pca.jndxyzlumr <- prcomp(non_PCA_colourspaces$jndxyzlumr)
summary(pca.jndxyzlumr)
saveRDS(pca.jndxyzlumr, "./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/1_Raw_PCA/Neoaves.patches.231030.PCAcolspaces.jndxyzlumr.240603.rds")


# The scaled luminance version looks identical to the non-scaled luminance version but actually is ever so slightly different 
# - I've tracked through the whole process and the scaling is working correctly, it just makes essentially no difference to 
# the PCA'd data whether the luminance is scaled or not. Seems to make a very small difference to the UMAP data. Should run everything
# with both sets and check that qualitatively the results are unchanged

# reload PCA (if necessary)
pca.jndxyzlumr <- readRDS("./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/Neoaves.patches.231030.PCAcolspaces.jndxyzlumr.240603.rds")

# plot the first two PCs (just for initial inspection - main visualisation below)
plot(pca.xyz$x)
plot(pca.xyzlum$x)
plot(pca.xyzlumr$x)
plot(pca.lab$x)
plot(pca.cie$x)
plot(pca.srgb$x)
plot(pca.jndxyz$x)
plot(pca.jndxyzlum$x)
plot(pca.jndxyzlumr$x)


# Calculate centroid dists (to see which species is farthest from 'average')

cdists <- data.frame(spec=rownames(pca.jndxyzlumr$x), xyz=NA, xyzlum=NA, xyzlumr=NA, lab=NA)

for (i in 1:nrow(pca.jndxyzlumr$x)) {
  print(i)
  cdists$xyz[i] <- dist(rbind(rep(0, ncol(pca.xyz$x)), pca.xyz$x[i,]))
  cdists$xyzlum[i] <- dist(rbind(rep(0, ncol(pca.xyzlum$x)), pca.xyzlum$x[i,]))
  cdists$xyzlumr[i] <- dist(rbind(rep(0, ncol(pca.xyzlumr$x)), pca.xyzlumr$x[i,]))
  cdists$lab[i] <- dist(rbind(rep(0, ncol(pca.lab$x)), pca.lab$x[i,]))
  cdists$jndxyzlum[i] <- dist(rbind(rep(0, ncol(pca.jndxyzlum$x)), pca.jndxyzlum$x[i,]))
  cdists$jndxyzlumr[i] <- dist(rbind(rep(0, ncol(pca.jndxyzlumr$x)), pca.jndxyzlumr$x[i,]))
}

# get species/sex which are maximally distant in each colour space
cdists$spec[cdists$xyz == max(cdists$xyz)]
cdists$spec[cdists$xyzlum == max(cdists$xyzlum)]
cdists$spec[cdists$xyzlumr == max(cdists$xyzlumr)]
cdists$spec[cdists$lab == max(cdists$lab)]
cdists$spec[cdists$jndxyzlum == max(cdists$jndxyzlum)]
cdists$spec[cdists$jndxyzlumr == max(cdists$jndxyzlumr)]
# and put the entire columns in order
ordered <- cdists$jndxyzlumr
names(ordered) <- cdists$spec
ordered <- order(ordered, decreasing = TRUE)
orderedlum <- cdists[ordered,]

# perform UMAP to visualise colour space in two dimensions
umap.xyz <- umap::umap(pca.xyz$x)
umap.xyzlum <- umap::umap(pca.xyzlum$x)
umap.xyzlumr <- umap::umap(pca.xyzlumr$x)
umap.lab <- umap::umap(pca.lab$x)
umap.jndxyzlum <- umap::umap(pca.jndxyzlum$x)
umap.jndxyzlumr <- umap::umap(pca.jndxyzlumr$x)

plot(umap.xyz$layout)
plot(umap.xyzlum$layout)
plot(umap.xyzlumr$layout)   # This looks almost exactly the same as the tcs xyzlum (non-scaled) umap plot but they are very slightly different
plot(umap.lab$layout)
plot(umap.jndxyzlum$layout)
plot(umap.jndxyzlumr$layout)


# Plot PCA axes and UMAP
plot.FourPatch(obj1 = pca.jndxyzlumr$x[, 1:2], 
             obj2 = pca.jndxyzlumr$x[, 3:4],
             obj3 = pca.jndxyzlumr$x[, 5:6],
             obj4 = umap.jndxyzlumr$layout,
             filepath = "./Plots/featurespaces/patches/patches.jndxyzlumr.pca1-6.umap.back.png",
             asp = T)

# Plot PC axes 1-8 with associated variance proportions
plot.FourPatch(obj1 = pca.jndxyzlumr$x[, 1:2], 
               obj2 = pca.jndxyzlumr$x[, 3:4],
               obj3 = pca.jndxyzlumr$x[, 5:6],
               obj4 = pca.jndxyzlumr$x[, 7:8],
               filepath = "./Plots/featurespaces/patches/patches.jndxyzlumr.pca1-8.png",
               PCobj = pca.jndxyzlumr,
               asp = T)

# Plot PC axes 9-16 with associated variance proportions
plot.FourPatch(obj1 = pca.jndxyzlumr$x[, 9:10], 
               obj2 = pca.jndxyzlumr$x[, 11:12],
               obj3 = pca.jndxyzlumr$x[, 13:14],
               obj4 = pca.jndxyzlumr$x[, 15:16],
               filepath = "./Plots/featurespaces/patches/patches.jndxyzlumr.pca9-16.png",
               PCobj = pca.jndxyzlumr,
               asp = T)

# Plot PC axes 17-24 with associated variance proportions
plot.FourPatch(obj1 = pca.jndxyzlumr$x[, 17:18], 
               obj2 = pca.jndxyzlumr$x[, 19:20],
               obj3 = pca.jndxyzlumr$x[, 21:22],
               obj4 = pca.jndxyzlumr$x[, 23:24],
               filepath = "./Plots/featurespaces/patches/patches.jndxyzlumr.pca17-24.png",
               PCobj = pca.jndxyzlumr,
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
px <- readRDS("./Outputs/features/patches/Neoaves.patches.231030.processed.240322.rds")

cps.jndxyzlumr <- data.frame(px$x.jndlumr, px$y.jndlumr, px$z.jndlumr, px$lum.jndlumr, 
                             row.names = paste(px$species, px$sex, sep = "-"))

mu = colMeans()
