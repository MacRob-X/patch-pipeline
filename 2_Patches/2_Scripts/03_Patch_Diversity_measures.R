# Extracting measures of diversity from patch colour pattern spaces
# 13th March 2024
# Robert MacDonald

# clear environment
rm(list=ls())

# load libraries
library(dispRity)
library(ggfortify)
library(umap)

# create colour palette (colours taken from ggsci palette pal_simpsons(palette = "springfield"))
springfield <- c(HomerYellow = "#FED439", FrinkPink = "#FD8CC1", HomerBlue = "#709AE1", DuffRed = "#C80813", 
                 BartOrange = "#FD7446", BurnsPurple = "#370335", MargeBlue = "#197EC0", LisaOrange = "#F05C3B", 
                 NedGreen = "#46732E", MaggieBlue = "#71D0F5", MargeGreen = "#D5E4A2", BurnsGreen = "#075149", 
                 HomerGrey = "#8A9197", KentRed = "#91331F", BobGreen = "#1A9993", HomerBrown = "#D2AF81")


#----------------------------------------------------------------------------------------------------#
# Functions ----


# Function to trim IUCN and colour data so they match ----
# Requires that the rownames of the colour data are of the form "Genus_species-sex"
# E.g. "Abeillia_abeillei-F"
# Takes as IUCN input the dataframe created by the code in IUCN_API.R

IUCNmatch <- function(colourData, iucnData){
  
  # check if IUCN data contains necessary columns
  if(!("scientific_name" %in% colnames(iucnData))){
    stop("Please provide IUCN data which contains Latin species names in a column titled 'scientific_name'")
  }
  
  # add species name separated by _ to IUCN data
  iucnData$species <- gsub(" ", "_", iucnData$scientific_name)   # note that I'll need to sort out the taxonomy before I do this for real
  iucnData$specimen <- c(rep(NA, times = length(iucnData$species)))
  
  # get female specimen IUCN data
  femIUCN <- iucnData[paste0(iucnData$species, "-F") %in% rownames(colourData), ]
  femIUCN$specimen <- paste0(femIUCN$species, "-F")
  
  # get male specimen IUCN data
  maleIUCN <- iucnData[paste0(iucnData$species, "-M") %in% rownames(colourData), ]
  maleIUCN$specimen <- paste0(maleIUCN$species, "-M")
  
  # get unknown sex specimen IUCN data
  unkIUCN <- iucnData[paste0(iucnData$species, "-U") %in% rownames(colourData), ]
  unkIUCN$specimen <- paste0(unkIUCN$species, "-U")
  
  # concatenate into one big dataframe with specimen names that match colourData rownames
  iucnMatch <- rbind(femIUCN, maleIUCN, unkIUCN)
  
  # trim species categorised as DD, EX, EW, or CR(PE) (after Hughes et al 2022)
  levels(as.factor(iucnMatch$category))
  iucnMatch <- iucnMatch[iucnMatch$category != "DD" & iucnMatch$category != "EX" & iucnMatch$category != "EW" & iucnMatch$category != "CR(PE)", ]
  levels(as.factor(iucnMatch$category))
  
  # remove species not present in IUCN data
  colourData <- colourData[rownames(colourData) %in% iucnMatch$specimen,]
  
  # reorder IUCN data to match colourData
  iucnMatch <- iucnMatch[match(rownames(colourData), iucnMatch$specimen), ]
  
  # Check all specimens match between the two datasets
  if(identical(iucnMatch$specimen, rownames(colourData))){
    print("All specimens in IUCN and colour data matched")
  }
  
  matchedData <- list(iucnMatch, colourData)
  
  return(matchedData)
  
}

#-----------------------------------------------------------------------------------------#

# Load and prepare data ----

# load patch data
pcaAll <- readRDS("./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/Neoaves.patches.231030.PCAcolspaces.jndxyzlumr.240404.rds")

# use first 16 PCs (~97% cumulative percentage of variance)
pcaDat <- pcaAll$x[,1:16]

# load IUCN Red List data
iucn <- readRDS("./4_SharedInputData/IUCN_RedList_data_130324.rds")


#--------------------------------------------------------------#


# Analyses (largely based on code from Hughes et al 2022) ----

# Mean distance to centroid by IUCN category - boxplot

# match iucn and PCA data
matchedData <- IUCNmatch(pcaDat, iucn)
iucnDat <- as.data.frame(matchedData[1])
pcaDat <- as.data.frame(matchedData[2])
rm(matchedData)

# convert data to matrix for analyses
pcaDat <- as.matrix(pcaDat)

# get distances to centroid for each specimen
centroidDistances <- as.data.frame(dispRity(pcaDat, metric = centroids)$disparity[[1]][[1]])
rownames(centroidDistances) <- rownames(pcaDat)
colnames(centroidDistances) <- c("centr_dist")
globDist <- centroidDistances

# attach IUCN categories to distance data
globDist$iucn <- iucnDat$category[match(rownames(globDist), iucnDat$specimen)] # note that since I've already matched the order the match() function is redundant here, but I'm leaving it in as a failsafe

# add duplicate of each specimen with iucn cat "All" (for plotting box plot)
allDist <- globDist
allDist$iucn <- c("All")

all <- rbind(globDist, allDist)

# reorder levels of iucn cat for plotting
all$iucn <- factor(all$iucn, levels = c("CR", "EN", "VU", "NT", "LC", "All"))

# create colour mapping
cols <- springfield[1:length(levels(all$iucn))]

# plot mean distance to centroid of each IUCN category (boxplot)
boxplot(centr_dist ~ iucn, data = all, 
        horizontal = FALSE,
        frame = FALSE,
        notch = TRUE,
        col = cols,
        xlab = "Subset of IUCN categories", ylab = "Distance to centroid")


# ok, so it looks like there might be some kind of difference between the mean distance to centroid (i.e. a measure of colour diversity)
# in different IUCN categories - let's just do a 
# basic ANOVA ----
# to find out

# make a model (will automatically perform an anova with a categorical predictor) with log transformed centroid distances to fulfil assumptions
mod <- lm(log(centr_dist) ~ iucn, data = all)

# check assumptions
autoplot(mod, smooth.colour = NA)

# create the anova table to get F and p values
anova(mod)
# note that I've set the reference level to "CR" here - so this result is telling us there's a significant difference between
# the log of the mean distance to the centroid for critically endangered species and the other categories, but it
# doesn't tell us which categories are driving the difference. Let's do a tukey test to find out

summary(mod)
# looks like LC and ALL species have significantly lower mean distance to centroid than CR species, while EN have significantly
# higher than CR

# convert to aov format
mod_aov <- aov(mod)

tukey_base <- TukeyHSD(mod_aov, ordered = T)
tukey_base

# let's say the pairwise comparisons I'm interested in are All-LC, All-NT, All-VU, All-EN and All-CR
tukey_base$iucn[rownames(tukey_base$iucn) == "All-LC" |
                rownames(tukey_base$iucn) == "NT-All" |
                rownames(tukey_base$iucn) == "VU-All" |
                rownames(tukey_base$iucn) == "EN-All" |
                rownames(tukey_base$iucn) == "CR-All", ]
# well that's interesting - it looks like NT are significantly more diverse in colour than All, as are VU and EN, but CR is not
# (in contrast to what the ANOVA test told us)

# let's check it out for LC-
tukey_base$iucn[rownames(tukey_base$iucn) == "NT-LC" |
                  rownames(tukey_base$iucn) == "CR-LC" |
                  rownames(tukey_base$iucn) == "VU-LC" |
                  rownames(tukey_base$iucn) == "EN-LC", ]

# Just for fun, let's look at distance to centroid vs midpoint of latitudinal range
latData <- read.csv("./Data/data_cooney_etal_latitude.csv", header = TRUE) # from Cooney et al 2022

# get female specimen lat data
femLat <- latData[paste0(latData$Binomial, "-F") %in% rownames(pcaDat), ]
femLat$specimen <- paste0(femLat$Binomial, "-F")
# get male specimen lat data
maleLat <- latData[paste0(latData$Binomial, "-M") %in% rownames(pcaDat), ]
maleLat$specimen <- paste0(maleLat$Binomial, "-M")
# get unknown sex specimen lat data
unkLat <- latData[paste0(latData$Binomial, "-U") %in% rownames(pcaDat), ]
unkLat$specimen <- paste0(unkLat$Binomial, "-U") # looks like there are no species in the lat data that have unknown sex in the patch data
# concatenate into one big dataframe with specimen names that match pcaDat rownames
lats <- rbind(femLat, maleLat)
# remove species not present in lat data
pcaDat <- pcaDat[rownames(pcaDat) %in% lats$specimen,]
# reorder IUCN data to match pcaDat
lats <- lats[match(rownames(pcaDat), lats$specimen), ]
# Check all specimens match between the two datasets
identical(lats$specimen, rownames(pcaDat))
# remove temporary datasets
rm(list = c("femLat", "maleLat", "unkLat"))

# get distances to centroid for each specimen
centroidDistances <- as.data.frame(dispRity(pcaDat, metric = centroids)$disparity[[1]][[1]])
rownames(centroidDistances) <- rownames(pcaDat)
colnames(centroidDistances) <- c("centr_dist")
latDists <- centroidDistances

# attach latitude data to centroid distances data
latDists$midpointLat <- lats$MidpointLat

# plot scatterplot of latitude vs distance to centroid
scatter.smooth(y = latDists$centr_dist, x = latDists$midpointLat, 
     frame = F,
     xlab = "Mean latitude of species range", ylab = "Distance to centroid")

latMod <- lm(centr_dist ~ abs(midpointLat), data = latDists)
# i should really check the assumptions of the model are met here but as I'm just messing around I'll skip it
summary(latMod)
# looks like distance to centroid declines with increasing (absolute) latitude, and the result is highly significant, although
# with an R2 or 0.015, it explains little of the observed variance. So latitude affects uniqueness of colour patternins, but
# only has a very small effect compared to other things



# Colour different IUCN categories on UMAP ----
umapAll <- umap(pcaAll$x)

# create dataframe from UMAP to manipulate
umapMan <- data.frame(umapAll$layout)
l <- length(umapAll$layout[,1])
for(i in 1:l){
  cat("\r", i)
  split <- strsplit(rownames(umapMan)[i], split = "-")[[1]]
  umapMan$species[i] <- split[1]
  umapMan$sex[i] <- split[2]
}

# match IUCN and umap data
matchedData <- IUCNmatch(umapMan, iucn)
iucnDat <- as.data.frame(matchedData[1])
umapDat <- as.data.frame(matchedData[2])
rm(matchedData, umapMan, l, split)

# attach IUCN categories to umap data
umapDat$iucn <- as.factor(iucnDat$category)

# create colour mapping
# create colour mapping
cols <- springfield[1:length(levels(umapDat$iucn))]
col_mapping <- cols[match(umapDat$iucn, unique(umapDat$iucn))]

# plot
plot(umapDat[, 1:2], 
     col = col_mapping,
     xlab = "UMAP 1", ylab = "UMAP 2")
legend("bottomright", 
       legend = levels(umapDat$iucn), 
       fill = cols, 
       title = "IUCN threat level",
       cex = 0.7, pt.cex = 0.7,)

# a little hard to understand - let's try a version where it's all greyed out except the category of interest
# we'll plot one plot for each of the five IUCN categories of interest (plus all together)

# creat transparent grey for the points we're not interested in
transGrey <- gray(level = 0.2, alpha = 0.2)

par(mfrow = c(3, 2))

# Plot one - all together, multicoloured

# create colour mapping
# create colour mapping
cols <- springfield[1:length(levels(umapDat$iucn))]
col_mapping <- cols[match(umapDat$iucn, levels(umapDat$iucn))]

# plot
plot(umapDat[, 1:2], 
     col = col_mapping,
     xlab = "UMAP 1", ylab = "UMAP 2")
legend("bottomright", 
       legend = levels(umapDat$iucn), 
       fill = cols, 
       title = "IUCN threat level",
       cex = 0.7, pt.cex = 0.7,
       box.lty = 0,
       inset = c(0.05, 0.05))

# Plot two - CR
cols <- c("black", rep(transGrey, times = length(levels(umapDat$iucn)) - 1))
col_mapping <- cols[match(umapDat$iucn, levels(umapDat$iucn))]

# plot
plot(umapDat[, 1:2], 
     col = col_mapping,
     xlab = "UMAP 1", ylab = "UMAP 2")
legend("bottomright", 
       legend = levels(umapDat$iucn), 
       fill = cols, 
       title = "IUCN threat level",
       cex = 0.7, pt.cex = 0.7,
       box.lty = 0)

# Plot three - EN
cols <- c(transGrey, "black", rep(transGrey, times = length(levels(umapDat$iucn)) - 2))
col_mapping <- cols[match(umapDat$iucn, levels(umapDat$iucn))]

# plot
plot(umapDat[, 1:2], 
     col = col_mapping,
     xlab = "UMAP 1", ylab = "UMAP 2")
legend("bottomright", 
       legend = levels(umapDat$iucn), 
       fill = cols, 
       title = "IUCN threat level",
       cex = 0.7, pt.cex = 0.7,
       box.lty = 0)

# Plot four - VU
cols <- c(rep(transGrey, times = length(levels(umapDat$iucn)) - 3), 
          "black", 
          rep(transGrey, times = length(levels(umapDat$iucn)) - 3))
col_mapping <- cols[match(umapDat$iucn, levels(umapDat$iucn))]

# plot
plot(umapDat[, 1:2], 
     col = col_mapping,
     xlab = "UMAP 1", ylab = "UMAP 2")
legend("bottomright", 
       legend = levels(umapDat$iucn), 
       fill = cols, 
       title = "IUCN threat level",
       cex = 0.7, pt.cex = 0.7,
       box.lty = 0)

# Plot five - NT
cols <- c(rep(transGrey, times = length(levels(umapDat$iucn)) - 2), 
          "black", 
          rep(transGrey, times = length(levels(umapDat$iucn)) - 4))
col_mapping <- cols[match(umapDat$iucn, levels(umapDat$iucn))]

# plot
plot(umapDat[, 1:2], 
     col = col_mapping,
     xlab = "UMAP 1", ylab = "UMAP 2")
legend("bottomright", 
       legend = levels(umapDat$iucn), 
       fill = cols, 
       title = "IUCN threat level",
       cex = 0.7, pt.cex = 0.7,
       box.lty = 0)

# Plot six - LC
cols <- c(rep(transGrey, times = length(levels(umapDat$iucn)) - 1), 
          "black")
col_mapping <- cols[match(umapDat$iucn, levels(umapDat$iucn))]

# plot
plot(umapDat[, 1:2], 
     col = col_mapping,
     xlab = "UMAP 1", ylab = "UMAP 2")
legend("bottomright", 
       legend = levels(umapDat$iucn), 
       fill = cols, 
       title = "IUCN threat level",
       cex = 0.7, pt.cex = 0.7,
       box.lty = 0)

##### Plot by taxon subgroup 

# taxonomy data
taxo <- read.csv("./4_SharedInputData/BLIOCPhyloMasterTax_2019_10_28.csv")


# create dataframe from UMAP to manipulate
umapDat <- data.frame(umapAll$layout)
l <- length(umapAll$layout[,1])
for(i in 1:l){
  cat("\r", i)
  split <- strsplit(rownames(umapDat)[i], split = "-")[[1]]
  umapDat$species[i] <- split[1]
  umapDat$sex[i] <- split[2]
}

# trim and order taxonomy data to match UMAP data
taxo$specimen <- c(rep(NA, times = length(taxo$TipLabel)))

# get female specimen IUCN data
femTaxo <- taxo[paste0(taxo$TipLabel, "-F") %in% rownames(umapDat), ]
femTaxo$specimen <- paste0(femTaxo$TipLabel, "-F")

# get male specimen IUCN data
maleTaxo <- taxo[paste0(taxo$TipLabel, "-M") %in% rownames(umapDat), ]
maleTaxo$specimen <- paste0(maleTaxo$TipLabel, "-M")

# get unknown sex specimen IUCN data
unkTaxo <- taxo[paste0(taxo$TipLabel, "-U") %in% rownames(umapDat), ]
unkTaxo$specimen <- paste0(unkTaxo$TipLabel, "-U")

# concatenate into one big dataframe with specimen names that match umapDat rownames
taxoMatch <- rbind(femTaxo, maleTaxo, unkTaxo)

# remove temp dataframes
rm(femTaxo, maleTaxo, unkTaxo)

# remove species not present in IUCN data
umapDat <- umapDat[rownames(umapDat) %in% taxoMatch$specimen,]

# reorder taxo data to match umapDat
taxoMatch <- taxoMatch[match(rownames(umapDat), taxoMatch$specimen), ]

# Check all specimens match between the two datasets
if(identical(taxoMatch$specimen, rownames(umapDat))){
  print("All specimens in IUCN and colour data matched")
}

# attach taxon subgroups to UMAP data
umapDat$taxon_subgroup <- taxoMatch$Taxon_subgroup

# create colour mapping
cols <- rainbow(length(levels(as.factor(umapDat$taxon_subgroup))))
col_mapping <- cols[match(umapDat$taxon_subgroup, levels(as.factor(umapDat$taxon_subgroup)))]

# plot
plot(umapDat[, 1:2], 
     col = col_mapping,
     xlab = "UMAP 1", ylab = "UMAP 2",
     asp = T)
legend("top", 
       legend = levels(umapDat$taxon_subgroup), 
       fill = cols, 
       title = "Taxon subgroup",
       cex = 0.5, pt.cex = 0.7,
       box.lty = 0,
       horiz = F)

# Let's try adding convex hulls, for fun
# Function taken from https://chitchatr.wordpress.com/2011/12/30/convex-hull-around-scatter-plot-in-r/
### Plotting function to plot convex hulls
### Filename: Plot_ConvexHull.R
### Notes:
############################################################################

# INPUTS:
# xcoords: x-coordinates of point data
# ycoords: y-coordinates of point data
# lcolor: line color

# OUTPUTS:
# convex hull around data points in a particular color (specified by lcolor)

# FUNCTION:
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}  
# END OF FUNCTION

# function to overlay points for each group on plot
fun_debug <- function(umapDat, cols, pointShapes, legend = T, legPos = "bottomright", sex = NULL, pointLabs = F){
  
  # trim dataset to required sex
  if(sex == "M"){
    umapDat <- umapDat[umapDat$sex == "M",]
  }
  else if(sex == "F"){
    umapDat <- umapDat[umapDat$sex == "F",]
  }
  else if(sex == "U"){
    umapDat <- umapDat[umapDat$sex == "U",]
  }
  
  # make taxon subgroup a factor
  umapDat$taxon_subgroup <- as.factor(umapDat$taxon_subgroup)
  
  # create blank plot
  plot(umapDat[, 1:2], 
       col = col_mapping,
       xlab = "UMAP 1", ylab = "UMAP 2",
       asp = T,
       type = "n")

    # legend positioning variables
  legendX <- par("usr")[2] - 5 # adjust as needed
  legendY <- par("usr")[4] - 0.1 # adjust as needed

    # overlay points, convex hull and legend points for each group
  for(i in 1:length(levels(umapDat$taxon_subgroup))){
    subgroup <- levels(umapDat$taxon_subgroup)[i]
    tmp <- umapDat[umapDat$taxon_subgroup == subgroup, 1:2]
    col <- cols[i]
    points(tmp,
           col = col,
           asp = T,
           pch = pointShapes[i])
    # add point labels if requested
    if(pointLabs == T){
      text(x = tmp$X1, y = tmp$X2, labels = umapDat[umapDat$taxon_subgroup == subgroup,]$species, pos = 3)
    }
    Plot_ConvexHull(xcoord = tmp[, 1],
                    ycoord = tmp[, 2],
                    lcolor = col)
    # add custom legend points
    points(legendX, legendY - i * 0.25,
           pch = pointShapes[i],
           col = col)
    text(legendX + 0.6, legendY - i * 0.25,
         adj = 0, # left align text
         cex = 0.6,
         subgroup)
  }
}

# create list of random point types
set.seed(10) # for reproducibility
pointShapes <- sample(1:8, size = 51, replace = T)   # size = 51 for 51 subtaxa

# plot
fun_debug(umapDat = umapDat, cols = cols, pointShapes = pointShapes, sex = "M", pointLabs = T)


# let's try doing it for passerines only
passDat <- umapDat[taxoMatch$PassNonPass == "PASSERIFORMES", ]
# create colour mapping
cols <- rainbow(length(levels(as.factor(passDat$taxon_subgroup))))
col_mapping <- cols[match(passDat$taxon_subgroup, levels(as.factor(passDat$taxon_subgroup)))]
# create shape mapping
set.seed(11) # for reproducibility
pointShapes <- sample(1:8, size = 17, replace = T)   # size = 17 for 17 subtaxa

fun_debug(umapDat = passDat, cols = cols, pointShapes = pointShapes, sex = "M", pointLabs = F)



# Let's create a bunch of svg files
svg_plot <- function(umapDat, cols, pointShapes, legend = T, legPos = "bottomright", sex = NULL, pointLabs = F){
  
  # trim dataset to required sex
  if(sex == "M"){
    umapDat <- umapDat[umapDat$sex == "M",]
  }
  else if(sex == "F"){
    umapDat <- umapDat[umapDat$sex == "F",]
  }
  else if(sex == "U"){
    umapDat <- umapDat[umapDat$sex == "U",]
  }
  
  # make taxon subgroup a factor
  umapDat$taxon_subgroup <- as.factor(umapDat$taxon_subgroup)
  
  # legend positioning variables
  legendX <- par("usr")[2] + 8 # adjust as needed
  legendY <- par("usr")[4] - 0.1 # adjust as needed
  
  for(i in 1:length(levels(as.factor(umapDat$taxon_subgroup)))){
    
    # get subgroup
    subgroup <- levels(as.factor(umapDat$taxon_subgroup))[i]
    
    # start svg file save
    svg(filename = paste0("./Plots/featurespaces/SubTaxaConvHulls/Passerines", subgroup, ".svg"),
        width = 30, height = 20,
        bg = "transparent")
    
    # create blank plot
    plot(umapDat[, 1:2], 
         col = col_mapping,
         xlab = "UMAP 1", ylab = "UMAP 2",
         asp = T,
         type = "n")
    
    # plot points, hull, legend
    tmp <- umapDat[umapDat$taxon_subgroup == subgroup, 1:2]
    col <- cols[i]
    points(tmp,
           col = col,
           asp = T,
           pch = pointShapes[i])
    
    # add point labels if requested
    if(pointLabs == T){
      text(x = tmp$X1, y = tmp$X2, labels = umapDat[umapDat$taxon_subgroup == subgroup,]$species, pos = 3)
    }
    Plot_ConvexHull(xcoord = tmp[, 1],
                    ycoord = tmp[, 2],
                    lcolor = col)
    # add custom legend points
    points(legendX, legendY - i * 0.25,
           pch = pointShapes[i],
           col = col)
    text(legendX + 0.6, legendY - i * 0.25,
         adj = 0, # left align text
         cex = 0.6,
         subgroup)
    
    dev.off()
  }
}

svg_plot(umapDat = passDat, cols = cols, pointShapes = pointShapesPass, legend = T, sex = "M", pointLabs = T)


#-------------------------------------------------------------------------------------#

## Functional evenness
# The minimum spanning tree distances evenness (from Villeger et al. 2008) - without abundance data
# High FEve means a regular distribution of the traits, while low FEve means clumping or irregular distributions in
# trait space, potentially indicating under-utilization of resources (Schleuter et al., 2010) or the absence of
#  corresponding conditions in the environment.
#  Zheng et al, 2021, Remote Sens. Environ

# working version of PCA data
pcaDat <- pcaAll$x

# match iucn and PCA data
matchedData <- IUCNmatch(pcaDat, iucn)
iucnDat <- as.data.frame(matchedData[1])
pcaDat <- as.data.frame(matchedData[2])
rm(matchedData)

# convert data to matrix for analyses
pcaDat <- as.matrix(pcaDat)

# get overall functionall evenness
minSpanEven <- dispRity::dispRity(pcaDat, metric = dispRity::func.div)

summary(minSpanEven)


## by taxon subgroup

# convert pca data to data frame
pcaDat <- as.data.frame(pcaDat)

# taxonomy data
taxo <- read.csv("./4_SharedInputData/BLIOCPhyloMasterTax_2019_10_28.csv")

# trim and order taxonomy data to match PCA data
taxo$specimen <- c(rep(NA, times = length(taxo$TipLabel)))

# get female specimen IUCN data
femTaxo <- taxo[paste0(taxo$TipLabel, "-F") %in% rownames(pcaDat), ]
femTaxo$specimen <- paste0(femTaxo$TipLabel, "-F")

# get male specimen IUCN data
maleTaxo <- taxo[paste0(taxo$TipLabel, "-M") %in% rownames(pcaDat), ]
maleTaxo$specimen <- paste0(maleTaxo$TipLabel, "-M")

# get unknown sex specimen IUCN data
unkTaxo <- taxo[paste0(taxo$TipLabel, "-U") %in% rownames(pcaDat), ]
unkTaxo$specimen <- paste0(unkTaxo$TipLabel, "-U")

# concatenate into one big dataframe with specimen names that match pcaDat rownames
taxoMatch <- rbind(femTaxo, maleTaxo, unkTaxo)

# remove temp dataframes
rm(femTaxo, maleTaxo, unkTaxo)

# remove species not present in IUCN data
pcaDat <- pcaDat[rownames(pcaDat) %in% taxoMatch$specimen,]

# reorder taxo data to match pcaDat
taxoMatch <- taxoMatch[match(rownames(pcaDat), taxoMatch$specimen), ]

# Check all specimens match between the two datasets
if(identical(taxoMatch$specimen, rownames(pcaDat))){
  print("All specimens in IUCN and colour data matched")
}

# attach taxon subgroups to PCA data
pcaDat$taxon_subgroup <- taxoMatch$Taxon_subgroup


# calculate minimum spanning evenness by taxon subgroup
taxon_subgroups <- unique(pcaDat$taxon_subgroup)
# create placeholder vector
min_span_by_taxon <- vector(
  "numeric",
  length = length(
    taxon_subgroups
  )
)
names(min_span_by_taxon) <- taxon_subgroups

# calculate minimum spanning distance for each subgroup
for(group in taxon_subgroups){
  min_span_by_taxon[group] <- dispRity::dispRity(
    as.matrix(
      pcaDat[pcaDat$taxon_subgroup == group, 1:ncol(pcaDat) - 1]
      ), 
    metric = dispRity::func.eve
    )$disparity[[1]][[1]]
}

saveRDS(
  min_span_by_taxon, 
  file = here::here(
    "2_Patches", "3_OutputData", "4_Diversity_measures",
    "PCA_jndxyzlumr_func_eve_by_taxon_subgroup.RDS"
  )
)

# calculate functional divergence for each group
# # "High and low FDiv values represent broad and narrow niches, respectively; high values indicating differences 
# between individuals or species and low values indicating similarities between individuals or species"
#  - Zheng et al, 2021, Remote Sens. Environ

# create placeholder vector
func_div_by_taxon <- vector(
  "numeric",
  length = length(
    taxon_subgroups
  )
)
names(func_div_by_taxon) <- taxon_subgroups

# calculate minimum spanning distance for each subgroup

for(group in taxon_subgroups){
  func_div_by_taxon[group] <- dispRity::dispRity(
    as.matrix(
      pcaDat[pcaDat$taxon_subgroup == group, 1:ncol(pcaDat) - 1]
    ), 
    metric = dispRity::func.div
  )$disparity[[1]][[1]]
}

saveRDS(
  func_div_by_taxon, 
  file = here::here(
    "2_Patches", "3_OutputData", "4_Diversity_measures",
    "PCA_jndxyzlumr_func_div_by_taxon_subgroup.RDS"
  )
)

# calculate mean nearest neighbour distance for each subgroup 
# note that this will only look for nearest neighbours within the same subgroup

# # create placeholder vector
nn_by_taxon <- vector(
  "numeric",
  length = length(
    taxon_subgroups
  )
)
names(nn_by_taxon) <- taxon_subgroups

# calculate mean nearest neghbour distance for each subgroup

for(group in taxon_subgroups){
  nn_by_taxon[group] <- dispRity::dispRity(
    as.matrix(
      pcaDat[pcaDat$taxon_subgroup == group, 1:ncol(pcaDat) - 1]
    ), 
    metric = dispRity::neighbours
  )$disparity[[1]][[1]] %>% 
    mean()
}


# plot functional divergence and evenness for each subgroup
library(ggplot2)

# count species in each subgroup
n_spec <- pcaDat %>% 
  group_by(taxon_subgroup) %>% 
  count() %>% 
  mutate(taxon_subgroup = janitor::make_clean_names(taxon_subgroup))

# prepare dataframe
div_df <- data.frame(taxon_subgroup = names(janitor::clean_names(min_span_by_taxon)),
                     func_eve = min_span_by_taxon,
                     func_div = func_div_by_taxon,
                     mean_nn = nn_by_taxon) %>% 
  inner_join(n_spec, by = join_by(taxon_subgroup)) %>% 
  rename(
    n_spec = n
  )
rownames(div_df) <- NULL

# pivot longer for plotting
div_df_long <- div_df |> 
  tidyr::pivot_longer(
    cols = c(func_eve, func_div, mean_nn),
    names_to = "metric",
    values_to = "value"
  )

# plot
ggplot() + 
  geom_col(data = div_df_long, 
           aes(x = value, y = taxon_subgroup, colour = taxon_subgroup, fill = taxon_subgroup),
           show.legend = FALSE) + 
  facet_wrap(facets = "metric")

# plot mean_nn against n_spec (to visually assess whether more speciose clades have lower mean_nn and therefore
# more densely pack into colour space)
ggplot() + 
  geom_point(data = div_df, aes(x = n_spec, y = mean_nn, colour = taxon_subgroup, shape = taxon_subgroup)) + 
  scale_shape_manual(values = rep(1:6, times = 100)[1:51])

# same for func_eve against n_spec
ggplot() + 
geom_point(data = div_df, aes(x = n_spec, y = func_eve, colour = taxon_subgroup, shape = taxon_subgroup),
           size = 1.5) + 
  scale_shape_manual(values = rep(1:6, times = 100)[1:51])


# Calculate overall functional evenness separately for males and females
# as might expect males to occupy edges while females occupy centroid

# let's just do it for a small subtaxon first - maybe tyrannidae
tyran_df <- pcaDat %>% 
  filter(taxon_subgroup == "Tyrannidae") %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  )

sexed_tyran_mat <- tyran_df %>% 
  filter(sex == "F") %>% 
  select(-c(taxon_subgroup, species, sex)) %>% 
  as.matrix() 

func_eve_tyran_f <- dispRity::dispRity(
  data = sexed_tyran_mat,
  metric = dispRity::func.eve
  )

func_eve_tyran_m$disparity[[1]][[1]]
func_eve_tyran_f$disparity[[1]][[1]]

func_div_tyran_m$disparity[[1]][[1]]
func_div_tyran_f$disparity[[1]][[1]]


