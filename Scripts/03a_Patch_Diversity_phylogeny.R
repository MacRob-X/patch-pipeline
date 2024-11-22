# Mapping PCA-based diversity measures onto the phylogeny
# Robert MacDonald
# 27th March 2024

# clear environment
rm(list=ls())

# load libraries
library(ape)
library(dispRity)
library(umap)
library(geiger)
library(ggtree)
library(treeio)
library(treeplyr)
library(ggplot2)


# Load data ----

# load patch data (PCA of whichever colourspace - generated in 02_Patch_Analyse_features.R)
pcaAll <- readRDS("./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/Neoaves.patches.231030.PCAcolspaces.jndxyzlumr.240404.rds")

# load first 100 trees of Hackett backbone full trees
phy <- read.tree("./4_SharedInputData/First100_AllBirdsHackett1.tre")


################ SHOULD I TRIM TO MALES/FEMALES ONLY BEFORE I CALCULATE DISTANCE TO CENTROID? ##################
# I'm choosing not to here but maybe speak to Chris about this

#--------------------------------------- SNIPPET 1 ---------------------------------------#


# Calculate distances to centroid for each specimen ----

# convert data to matrix
pcaDat <- as.matrix(pcaAll$x)

# get distances to centroid for each specimen
centroidDistances <- as.data.frame(dispRity(pcaDat[, 1:40], metric = centroids)$disparity[[1]][[1]])
rownames(centroidDistances) <- rownames(pcaDat)
colnames(centroidDistances) <- c("centr_dist")

# prepare data and phylogeny ----

# add species name and sex as columns in centroid data
centroidDistances$species <- c(rep(NA, times = length(centroidDistances[,1])))
centroidDistances$sex <- c(rep(NA, times = length(centroidDistances[,1])))
for(i in 1:length(centroidDistances[,1])){
  cat("\r", i)
  split <- strsplit(rownames(centroidDistances)[i], split = "-")[[1]]
  centroidDistances$species[i] <- split[1]
  centroidDistances$sex[i] <- split[2]
}
rm(split)
# reorder columns
centroidDistances <- centroidDistances[, c("species", "sex", "centr_dist")]

# match centroid data and phylogeny ----

# Choose sex of interest and trim
sexedCentroidDistances <- centroidDistances[centroidDistances$sex == "M", ]

# drop species not in phylogeny
sexedCentroidDistances <- sexedCentroidDistances[sexedCentroidDistances$species %in% phy[[1]]$tip.label, ]

# match phylogeny to species in dataset
phy <- drop.tip.multiPhylo(phy, tip = setdiff(phy[[1]]$tip.label, sexedCentroidDistances$species))


# rename rows as species names and drop species/sex columns
rownames(sexedCentroidDistances) <- sexedCentroidDistances$species
sexedCentroidDistances <- subset(sexedCentroidDistances, select = -c(sex, species))

# reorder centroid distances to match phylogeny
sexedCentroidDistances <- sexedCentroidDistances[match(phy[[1]]$tip.label, rownames(sexedCentroidDistances)), , drop = FALSE]

# check reordering
identical(phy[[1]]$tip.label, rownames(sexedCentroidDistances))

#--------------------------------------- END SNIPPET 1 ---------------------------------------#

#--------------------------------------- SNIPPET 2 ---------------------------------------#


# Calculate mean nearest neightbour distance for each specimen ----

# convert data to matrix
pcaDat <- as.matrix(pcaAll$x)

# get distances to centroid for each specimen
meanNNDistances <- as.data.frame(dispRity(pcaDat[, 1:40], metric = neighbours)$disparity[[1]][[1]])
rownames(meanNNDistances) <- rownames(pcaDat)
colnames(meanNNDistances) <- c("mean_NN_dist")

# prepare data and phylogeny ----

# add species name and sex as columns in mean NN data
meanNNDistances$species <- c(rep(NA, times = length(meanNNDistances[,1])))
meanNNDistances$sex <- c(rep(NA, times = length(meanNNDistances[,1])))
for(i in 1:length(meanNNDistances[,1])){
  cat("\r", i)
  split <- strsplit(rownames(meanNNDistances)[i], split = "-")[[1]]
  meanNNDistances$species[i] <- split[1]
  meanNNDistances$sex[i] <- split[2]
}
rm(split)
# reorder columns
meanNNDistances <- meanNNDistances[, c("species", "sex", "mean_NN_dist")]

# match mean NN data and phylogeny ----

# Choose sex of interest and trim
sexedMeanNNDistances <- meanNNDistances[meanNNDistances$sex == "M", ]

# drop species not in phylogeny
sexedMeanNNDistances <- sexedMeanNNDistances[sexedMeanNNDistances$species %in% phy[[1]]$tip.label, ]

# match phylogeny to species in dataset
phy <- drop.tip.multiPhylo(phy, tip = setdiff(phy[[1]]$tip.label, sexedMeanNNDistances$species))


# rename rows as species names and drop species/sex columns
rownames(sexedMeanNNDistances) <- sexedMeanNNDistances$species
sexedMeanNNDistances <- subset(sexedMeanNNDistances, select = -c(sex, species))

# reorder mean NN distances to match phylogeny
sexedMeanNNDistances <- sexedMeanNNDistances[match(phy[[1]]$tip.label, rownames(sexedMeanNNDistances)), , drop = FALSE]

# check reordering
identical(phy[[1]]$tip.label, rownames(sexedMeanNNDistances))

#--------------------------------------- END SNIPPET 2 ---------------------------------------#

# calculate ancestral state (let's just use the first tree for test purposes) ----
ancStatesBM.M <- fastAnc(phy[[1]], sexedCentroidDistances$centr_dist, CI = TRUE, vars = TRUE)

# produce df to use for plotting
td <- data.frame(node = nodeid(phy[[1]], rownames(sexedCentroidDistances)),
                 centr_dist = sexedCentroidDistances$centr_dist)
nd <- data.frame(node = names(ancStatesBM.M$ace), centr_dist = ancStatesBM.M$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)

# put data and tree together with full_join()
maleCentDistTree <- full_join(phy[[1]], d, by = "node")
rm(td, nd, d)

# get max and min centroid distances for plot colour scale (round to nearest 100)
# allows me to directly compare male and female plots
uplim <- round(max(centroidDistances$centr_dist) / 100) * 100
lowlim <- round(min(centroidDistances$centr_dist) / 100) * 100

# Plot ancestral states on phylogeny with colour according to centroid distance
# Note that I'm not actually interested in the ancestral states but I think it will make the plot look nicer
# and I guess it could show something interesting
p <- ggtree(maleCentDistTree, aes(color = centr_dist),
            ladderize = TRUE, continuous = "colour", size = 0.35) +    # adjust size here to make branches thinner
  scale_color_gradientn(colours = viridisLite::viridis(10), # gradientn allows us to specify how many colours to use in scale
                        limits = c(lowlim, uplim))          # sets limits for colour scale
p <- p + geom_tiplab(offset = 1, size = 0.55) + theme_tree2() + # Adjust the size parameter here to make the tip label text smaller
  theme(panel.grid.major = element_line(color = "black", linewidth = .2),
        panel.grid.minor = element_line(color = "lightgrey", linewidth = .2),
        panel.grid.major.y = element_blank(),
        legend.position.inside = c(0.15, 0.1))
revts(p)  # reverses timescale by setting most recent tip to 0
# save as portrait 100 x 10 inch pdf


#----------------------------------------------------------------------------------------------------------#
# Calculate subgroup mean distance from centroid 

# First run Snippet 1 to get individual sexed centroid distances

# now calculate mean subgroup centroid distances

# taxonomy data
taxo <- read.csv("./4_SharedInputData/BLIOCPhyloMasterTax_2019_10_28.csv")


# trim and order taxonomy data to match centroid data
taxo$specimen <- c(rep(NA, times = length(taxo$TipLabel)))

# get female specimen taxonomic data
femTaxo <- taxo[paste0(taxo$TipLabel, "-F") %in% rownames(centroidDistances), ]
femTaxo$specimen <- paste0(femTaxo$TipLabel, "-F")

# get male specimen taxonomic data
maleTaxo <- taxo[paste0(taxo$TipLabel, "-M") %in% rownames(centroidDistances), ]
maleTaxo$specimen <- paste0(maleTaxo$TipLabel, "-M")

# get unknown sex specimen taxonomic data
unkTaxo <- taxo[paste0(taxo$TipLabel, "-U") %in% rownames(centroidDistances), ]
unkTaxo$specimen <- paste0(unkTaxo$TipLabel, "-U")

# concatenate into one big dataframe with specimen names that match centroidDistances rownames
taxoMatch <- rbind(femTaxo, maleTaxo, unkTaxo)

# if using data for all sexes, remove temp dataframes. Else use individual sexed taxonomic data
#rm(femTaxo, maleTaxo, unkTaxo)

# reorder taxo data to match centroid data
# IF looking at males, use maleTaxo and sexedCentroidDistances
# IF looking at males, use femaleTaxo and sexedCentroidDistances
# IF looking at all specimens, use taxoMatch and centroidDistances
sexedTaxo <- maleTaxo[match(rownames(sexedCentroidDistances), maleTaxo$TipLabel), ]

# Check all specimens match between the two datasets
if(identical(sexedTaxo$TipLabel, rownames(sexedCentroidDistances)) & !any(is.na(sexedTaxo$TipLabel))){
  print("All specimens in taxonomic and centroid data matched")
} else {
  print("Error - specimens in taxonomic and centroid data do not match or there are NAs")
}


# attach sex average centroid distances to taxonomic data
sexedTaxo$centroid_distance <- sexedCentroidDistances$centr_dist

# attach sex average mean NN distances to taxonomic data
sexedTaxo$mean_NN_distance <- sexedMeanNNDistances$mean_NN_dist
# 
# # calculate mean centroid distance of each subgroup and attach
# sexedTaxo$subgroup_centrdist <- rep(NA, times = nrow(sexedTaxo))
# 
# # create dataframe to populate
# subgroups <- data.frame(subgroup = unique(sexedTaxo$Taxon_subgroup), 
#                         centr_dist = rep(NA, times = length(unique(sexedTaxo$Taxon_subgroup))),
#                         centrdist_sd = rep(NA, times = length(unique(sexedTaxo$Taxon_subgroup))))
# # populate dataframe with means of centroid distances for each subgroup
# for (isubgroup in subgroups$subgroup) {
#   # calculate mean subgroup centroid distance and add to dataframe
#   subgroups$centr_dist[subgroups$subgroup == isubgroup] <- mean(sexedTaxo$centroid_distance[sexedTaxo$Taxon_subgroup == isubgroup])
#   # calculate standard deviation of subgroup centroid distance and add to dataframe
#   subgroups$centrdist_sd[subgroups$subgroup == isubgroup] <- sd(sexedTaxo$centroid_distance[sexedTaxo$Taxon_subgroup == isubgroup])
#   # populate main dataframe
#   sexedTaxo$subgroup_centrdist[sexedTaxo$Taxon_subgroup == isubgroup] <- subgroups$centr_dist[subgroups$subgroup == isubgroup]
#   sexedTaxo$subgroup_centrdist_sd[sexedTaxo$Taxon_subgroup == isubgroup] <- subgroups$centrdist_sd[subgroups$subgroup == isubgroup]
# }


# Boxplot of mean subgroup centroid distances ----
# create colour mapping

# Calculate median values for each subgroup
mediansCentroid <- tapply(sexedTaxo$centroid_distance, sexedTaxo$Taxon_subgroup, median)

# Assign colors based on median values using viridis palette
colours <- viridisLite:::viridis(length(mediansCentroid))

# Match colors to subgroups based on their median values
cols <- colours[match(mediansCentroid, sort(mediansCentroid))]

# Set wider margins - bottom, left, top, right
par(mar = c(5, 13, 1, 1), 
    mgp = c(12, 1, 0))

# plot boxplot of mean distance to centroid for each taxon subgroup
boxplot(centroid_distance ~ Taxon_subgroup, data = sexedTaxo,
        horizontal = TRUE,
        frame = FALSE,
        notch = FALSE,
        col = cols,
        xlab = "Taxon subgroup", ylab = "Mean distance to centroid",
        las = 1) # make labels horizontal


# Boxplot of subgroup mean NN distances ----
# Calculate median values for each subgroup
mediansNN <- tapply(sexedTaxo$mean_NN_distance, sexedTaxo$Taxon_subgroup, median)

# Assign colors based on median values using viridis palette
colours <- viridisLite:::viridis(length(mediansNN))

# Match colors to subgroups based on their median values
cols <- colours[match(mediansNN, sort(mediansNN))]

# Set wider margins - bottom, left, top, right
par(mar = c(5, 13, 1, 1), 
    mgp = c(12, 1, 0))

# plot boxplot of mean nearest neighbour distance for each taxon subgroup
boxplot(mean_NN_distance ~ Taxon_subgroup, data = sexedTaxo,
        horizontal = TRUE,
        frame = FALSE,
        notch = FALSE,
        col = cols,
        xlab = "Taxon subgroup", ylab = "Mean nearest neighbour distance",
        las = 1) # make labels horizontal


# Plot scatter plot of nearest neighbour distance and distance to centroid
plot(mean_NN_distance ~ centroid_distance, data = sexedTaxo)

# as ggplot
p <- ggplot(sexedTaxo, aes(x = centroid_distance, y = mean_NN_distance)) + 
  geom_point()
p

# plot medians of each subgroup
summaryData <- sexedTaxo %>% 
  group_by(Taxon_subgroup) %>% 
  summarise(median_NN_distance = median(mean_NN_distance),
            sd_NN_distance = sd(mean_NN_distance),
            se_NN_distance = sd(mean_NN_distance) / sqrt(n()),
            median_centroid_distance = median(centroid_distance),
            sd_centroid_distance = sd(centroid_distance),
            se_centroid_distance = sd(centroid_distance) / sqrt(n()))

# plot summary data
p2 <- ggplot(summaryData, aes(x = median_centroid_distance, y = median_NN_distance,
                              color = Taxon_subgroup, shape = Taxon_subgroup)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = median_NN_distance - sd_NN_distance,
                    ymax = median_NN_distance + sd_NN_distance),
                width = 0.1) + 
  geom_errorbarh(aes(xmin = median_centroid_distance - sd_centroid_distance,
                    xmax = median_centroid_distance + sd_centroid_distance),
                height = 0.1) + 
 geom_text(aes(label = Taxon_subgroup), vjust = -0.5, hjust = -0.05) +   # adds labels to the plot itself
  theme_bw()
p2
