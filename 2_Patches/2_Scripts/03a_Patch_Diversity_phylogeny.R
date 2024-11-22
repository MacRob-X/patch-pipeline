# Mapping PCA-based diversity measures onto the phylogeny
# Robert MacDonald
# 27th March 2024

# clear environment
rm(list=ls())

# load libraries
library(dplyr)
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
pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    "Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603.rds"
  )
)

# load first 100 trees of Hackett backbone full trees
phy <- ape::read.tree(
  here::here(
    "4_SharedInputData", "First100_AllBirdsHackett1.tre"
  )
)

# load taxonomy and rename columns
taxo <- read.csv(
  here::here(
    "4_SharedInputData", "BLIOCPhyloMasterTax_2019_10_28.csv"
  ),
  strings = F
) %>% 
  rename(
    genus = GenusName,
    family = BLFamilyLatin,
    order = IOCOrder,
    taxon_subgroup = Taxon_subgroup
  ) %>% 
  select(
    TipLabel, genus, family, order, taxon_subgroup
  )


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


# ------------------------------------------------------------------------------------------------------#
# Calculate within-group diversity ----
# mean distance to group centroid and mean pairwise distance

# Specify taxonomic level to calculate at ("genus", "family", "order", "taxon_subgroup")
tax_level <- "family"
# Select sex to focus on ("M", "F", "All")
sex_choice <- "M"
## Drop depauperate clades? Specify min. number of species a clade must contain for it to be included
## Minimum is 2, as can't calculate centroid distances or pairwise distances for groups with
## only one species
drop_depaup <- 2
## Plot on log scale?
log_choice <- TRUE



# Convert PCA data to dataframe and append taxonomy
pca_dat <- pca_all %>% 
  magrittr::extract2("x") %>% 
  as.data.frame() %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  ) %>% 
  left_join(taxo, by = join_by("species" == "TipLabel"))

# Filter by sex
if(sex_choice != "All") {
  pca_dat <- pca_dat %>% 
    filter(
      sex == sex_choice
    )
}

# get clades with fewer species than specified minimum
clades_to_keep <- pca_dat %>% 
  count(get(tax_level)) %>% 
  filter(
    n >= drop_depaup
  ) %>% 
  magrittr::extract2(1)

# remove these clades from data
pca_dat <- pca_dat[pca_dat[, tax_level] %in% clades_to_keep, ]


# Get list of unique taxa
taxa <- pca_dat[, tax_level] %>% 
  unique()


# Create empty summary dataframe to populate with diversity metrics (for basic plots)
summary_within_group_diversity <- data.frame(
  taxa = taxa, 
  mean_group_centr_dist = NA, 
  mean_pairwise_dist = NA)
colnames(summary_within_group_diversity) <- c(tax_level, "mean_group_centr_dist", "mean_pairwise_dist")

# Create empty full species-level dataframe to populate with diversity metrics (for phylogenetic tree plots)
within_group_diversity <- data.frame(
  pca_dat[, c("species", "sex", tax_level)],
  group_centr_dist = NA,
  mean_group_centr_dist = NA,
  mean_pairwise_dist = NA
)

# Loop over vector of unique taxa
for (taxon in taxa) {
  
  # Get species data in taxon of interest in matrix form
  dat <- pca_dat[pca_dat[, tax_level] == taxon, ]
  rownames(dat) <- paste(dat$species, dat$sex, sep = "-")
  mat <- as.matrix(dat[, 1:40])
  
  # Calculate distance of each species to within-group centroid and append to main dataframe
  within_group_diversity$group_centr_dist[within_group_diversity[, tax_level] == taxon] <- dispRity::dispRity(
    mat, metric = dispRity::centroids
  )$disparity[[1]][[1]]
  
  # Calculate mean distance to within-group centroid and add to main dataframe
  within_group_diversity$mean_group_centr_dist[within_group_diversity[, tax_level] == taxon] <- mean(
    within_group_diversity$group_centr_dist[within_group_diversity[, tax_level] == taxon]
    )
  
  # Calculate mean distance to within-group centroid and add to summary dataframe
  summary_within_group_diversity$mean_group_centr_dist[summary_within_group_diversity[, tax_level] == taxon] <- mean(
    within_group_diversity$group_centr_dist[within_group_diversity[, tax_level] == taxon]
    )
  
  
  
  # Calculate mean within-group pairwise distance and append to main dataframe
  within_group_diversity$mean_pairwise_dist[within_group_diversity[, tax_level] == taxon] <- mean(
      dispRity::dispRity(
      mat, metric = dispRity::pairwise.dist
    )$disparity[[1]][[1]]
  )
  
  # Add to summary dataframe
  summary_within_group_diversity$mean_pairwise_dist[summary_within_group_diversity[, tax_level] == taxon] <-   within_group_diversity$mean_pairwise_dist[within_group_diversity[, tax_level] == taxon][1]
  
  
}


## Plotting on phylogeny

# Select one tree for plotting purposes
tree <- phy[[1]]

# Get only first row of each taxonomic family/order/tax subgroup
phylo_dat <- within_group_diversity %>% 
  group_by(get(tax_level)) %>% 
  filter(row_number() == 1) %>% 
  as.data.frame() %>% 
  select(
    -ncol(.)
  )
rownames(phylo_dat) <- phylo_dat$species

# Filter tree to one species per taxon
tree <- ape::drop.tip(tree, tip = setdiff(tree$tip.label, rownames(phylo_dat)))

# Reorder diversity data to match phylogeny
phylo_dat <- phylo_dat[match(tree$tip.label, rownames(phylo_dat)), , drop = F]
# Verify identical
identical(rownames(phylo_dat), tree$tip.label)

# Overwrite species column and rownames as taxon
phylo_dat$species <- phylo_dat[, tax_level]
rownames(phylo_dat) <- phylo_dat[, tax_level]

# Overwrite tree tips with taxon only
tree$tip.label <- phylo_dat[, tax_level]

# Get tree trait tip data
td <- data.frame(node = tidytree::nodeid(tree, rownames(phylo_dat)),
                 mean_group_centr_dist = phylo_dat[, "mean_group_centr_dist"],
                 mean_pairwise_dist = phylo_dat[, "mean_pairwise_dist"])
td$node <- as.numeric(td$node)

# put data and tree together
phylo_dat_tree <- full_join(tree, td, by = "node")
rm(td)

## Plot diversity statistics on temrinal branches and tips
# set plotting parameters
if(tax_level == "species"){
  branch_width <- 0.35
  text_size <- 0.55
  line_width <- 0.2
  tip_size <- 1
  folder <- "species_level"
} else if (tax_level == "genus"){
  branch_width <- 0.35
  text_size <- 0.65
  line_width <- 0.3
  tip_size <- 1
  folder <- "genus_level"
} else if (tax_level == "family"){
  branch_width <- 0.45
  text_size <- 1.8
  line_width <- 0.7
  tip_size <- 1
  folder <- "family_level"
} else if (tax_level == "order"){
  branch_width <- 0.35
  text_size <- 2
  line_width <- 1
  tip_size <- 1
  folder <- "order_level"
} else if (tax_level == "taxon_subgroup"){
  branch_width <- 0.7
  text_size <- 2
  line_width <- 2
  tip_size <- 4
  folder <- "taxon_subgroup_level"
}

# Mean distance to group centroid
if(log_choice == TRUE) {
  p <- ggtree::ggtree(phylo_dat_tree, aes(colour = log(mean_group_centr_dist)),
                    layout = "circular",
                    ladderize = TRUE, continuous = "colour", size = branch_width) + 
  geom_tippoint(aes(colour = log(mean_group_centr_dist)), size = tip_size) + 
  geom_tiplab(aes(colour = log(mean_group_centr_dist)), offset = 3, size = text_size) +
  scale_color_viridis_c(option = "D")
} else {
  p <- ggtree::ggtree(phylo_dat_tree, aes(colour = mean_group_centr_dist),
                      layout = "circular",
                      ladderize = TRUE, continuous = "colour", size = branch_width) + 
    geom_tippoint(aes(colour = mean_group_centr_dist), size = tip_size) + 
    geom_tiplab(aes(colour = mean_group_centr_dist), offset = 3, size = text_size) +
    scale_color_viridis_c(option = "D")
  }
p

# Save as SVG
if(log_choice == TRUE) {
  log_filename <- "LOG_"
} else {
  log_filename <- ""
}

svg(
  here::here(
    "2_Patches", "4_OutputPlots", "2_Diversity_measures", "2_Diversity_phylogeny",
    "within_group_diversity", 
    paste0("patch_jndxyzlum_", log_filename, "mean_group_centr_dist", "_", sex_choice, "_", "phylogeny.svg")
  ), 
  width = 6, height = 6
)

p

dev.off()

# Mean within-group pairwise distance
if(log_choice == TRUE) {
  p <- ggtree::ggtree(phylo_dat_tree, aes(colour = log(mean_pairwise_dist)),
                      layout = "circular",
                      ladderize = TRUE, continuous = "colour", size = branch_width) + 
    geom_tippoint(aes(colour = log(mean_pairwise_dist)), size = tip_size) + 
    geom_tiplab(aes(colour = log(mean_pairwise_dist)), offset = 3, size = text_size) +
    scale_color_viridis_c(option = "D")
} else {
  p <- ggtree::ggtree(phylo_dat_tree, aes(colour = mean_pairwise_dist),
                      layout = "circular",
                      ladderize = TRUE, continuous = "colour", size = branch_width) + 
    geom_tippoint(aes(colour = mean_pairwise_dist), size = tip_size) + 
    geom_tiplab(aes(colour = mean_pairwise_dist), offset = 3, size = text_size) +
    scale_color_viridis_c(option = "D")
}
p

# Save as SVG
if(log_choice == TRUE) {
  log_filename <- "LOG_"
} else {
  log_filename <- ""
}

svg(
  here::here(
    "2_Patches", "4_OutputPlots", "2_Diversity_measures", "2_Diversity_phylogeny",
    "within_group_diversity", 
    paste0("patch_jndxyzlum_", log_filename, "mean_pairwise_dist", "_", sex_choice, "_", "phylogeny.svg")
  ), 
  width = 6, height = 6
)

p

dev.off()
