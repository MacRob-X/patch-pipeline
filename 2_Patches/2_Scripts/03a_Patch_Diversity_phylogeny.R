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

## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neognaths"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
# N.B. will need to change the date in the pca_all filename if using usmldbl or usmldblr
# (from 240603 to 240806)
space <- "lab"
# Select number of PC axes to retain ("all" or a number)
axes <- "all"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# select sex of interest ("all", "male_female", "male_only", "female_only", "unknown_only")
sex_interest <- "all"


# Load data ----

# load patch data (PCA of whichever colourspace - generated in 02_Patch_Analyse_features.R)
pca_filename <- paste(clade, sex_match, "patches.250716.PCAcolspaces", "rds", sep = ".")
pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", clade, "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]]

# load taxonomy, rename columns, and convert taxon subgroup to sentence case
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
  ) %>% 
  mutate(
    taxon_subgroup = snakecase::to_sentence_case(taxon_subgroup)
  )


# load first 100 trees of Hackett backbone full trees
phy <- ape::read.tree(
  here::here(
    "4_SharedInputData", "First100_AllBirdsHackett1.tre"
  )
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

# Choose disparity metric to calculate ("centr-dist", "nn-k", "nn-all", "nn-count", "sum.variances", "sum.ranges", "convhull.volume")
metric <- "centr-dist"
# choose number of PC axes to work with ("all" or a number)
axes <- "all"
# Specify taxonomic level to calculate at ("genus", "family", "order", "taxon_subgroup")
tax_level <- "taxon_subgroup"
# Select sex to focus on ("M", "F", "All")
sex_choice <- "All"
# select type of averaging to use ("mean" or "median")
avg_par <- "median"
## Drop depauperate clades? Specify min. number of species a clade must contain for it to be included
## Minimum is 2, as can't calculate centroid distances or pairwise distances for groups with
## only one species
drop_depaup <- 10
## Plot on log scale?
log_choice <- FALSE


# FUNCTIONS

source(
  here::here(
    "2_Patches", "2_Scripts", "R", "mapping.R"
  )
)

library(dispRity)
library(ggplot2)

# convert PCA data to dataframe and append taxonomy
join_taxo <- function(pca_data, taxonomy, spec_taxo = "TipLabel"){
  
  pca_data %>% 
    as.data.frame() %>% 
    mutate(
      species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
      sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
    ) %>% 
    left_join(taxonomy, by = setNames(spec_taxo, "species"))
  
}

# calculate within-group diversity
calc_group_div <- function(sexed_pca_data, tax_level, metric, n_dim, drop_depaup, summary = FALSE, avg_par){
  
  # set dispRity metric
  metric_get <- set_metric(metric)
  
  # get clades with fewer species than specified minimum
  clades_to_keep <- sexed_pca_data %>% 
    count(get(tax_level)) %>% 
    filter(
      n >= drop_depaup
    ) %>% 
    magrittr::extract2(1)
  
  # remove these clades from data
  sexed_pca_data <- sexed_pca_data[sexed_pca_data[, tax_level] %in% clades_to_keep, ]
  
  
  # Get list of unique taxa
  taxa <- sexed_pca_data[, tax_level] %>% 
    unique()
  
  # Create empty full species-level dataframe to populate with diversity metrics (for phylogenetic tree plots)
  within_group_diversity <- data.frame(
    sexed_pca_data[, c("species", "sex", tax_level)],
    group_metric = NA,
    avg_group_metric = NA
  )
  
  if(summary == FALSE){
    
    # Loop over vector of unique taxa
    for (taxon in taxa) {
      
      # Get species data in taxon of interest in matrix form
      dat <- sexed_pca_data[sexed_pca_data[, tax_level] == taxon, ]
      rownames(dat) <- paste(dat$species, dat$sex, sep = "-")
      mat <- as.matrix(dat[, 1:n_dim])
      
      # Calculate distance of each species to within-group centroid and append to main dataframe
      within_group_diversity$group_metric[within_group_diversity[, tax_level] == taxon] <- dispRity::dispRity(
        mat, metric = metric_get
      )$disparity[[1]][[1]]
      
      # Calculate mean within-group metric and add to main dataframe
      within_group_diversity$avg_group_metric[within_group_diversity[, tax_level] == taxon] <- avg(
        within_group_diversity$group_metric[within_group_diversity[, tax_level] == taxon], avg_type = avg_par
      )
      
    }
    
    return(within_group_diversity)
    
  } else if(summary == TRUE){
    
    # Create empty summary dataframe to populate with diversity metrics (for basic plots)
    summary_within_group_diversity <- data.frame(
      taxa = taxa, 
      avg_group_metric = NA
      )
    colnames(summary_within_group_diversity) <- c(tax_level, "avg_group_metric")
    
    # Loop over vector of unique taxa
    for (taxon in taxa) {
      
      # Get species data in taxon of interest in matrix form
      dat <- sexed_pca_data[sexed_pca_data[, tax_level] == taxon, ]
      rownames(dat) <- paste(dat$species, dat$sex, sep = "-")
      mat <- as.matrix(dat[, 1:n_dim])
      
      # Calculate distance of each species to within-group centroid and append to main dataframe
      within_group_diversity$group_metric[within_group_diversity[, tax_level] == taxon] <- dispRity::dispRity(
        mat, metric = metric_get
      )$disparity[[1]][[1]]
      
      # Calculate mean distance to within-group centroid and add to main dataframe
      within_group_diversity$avg_group_metric[within_group_diversity[, tax_level] == taxon] <- avg(
        within_group_diversity$group_metric[within_group_diversity[, tax_level] == taxon], avg_type = avg_par
      )
      
      # Calculate mean distance to within-group centroid and add to summary dataframe
      summary_within_group_diversity$avg_group_metric[summary_within_group_diversity[, tax_level] == taxon] <- avg(
        within_group_diversity$group_metric[within_group_diversity[, tax_level] == taxon], avg_type = avg_par
      )
      
      
    }
    
    return(summary_within_group_diversity)
  }
  
}

# Data preparation and analysis

# restrict PCA to requested number of axes
if(axes != "all"){
  pca_dat <- pca_all[["x"]][, 1:axes]
  n_dim <- axes
} else {
  pca_dat <- pca_all[["x"]]
  n_dim <- ncol(pca_dat)
}

# convert PCA data to dataframe and append taxonomy
pca_dat <- join_taxo(pca_dat, taxo, spec_taxo = "TipLabel")



# get list of sexes (function from mapping code)
# we will use these to iterate over
if(sex_interest == "male_female"){
  sexes <- set_sex_list(sex_interest)
} else{
  # add extra rows for overall (unsexed) diversity
  all_sex_dat <- pca_dat
  all_sex_dat$sex <- paste(all_sex_dat$sex, "All", sep = "_")
  pca_dat <- rbind(pca_dat, all_sex_dat)
  sexes <- c("All", "M", "F")
}


# calculate within-group diversity for each sex separately
within_group_diversity <- lapply(
  sexes, 
  function(sex, pca_data, taxo_level, n_dim, drop_depaup, avg_par, metric){
    
    if(sex != "All"){
      sexed_data <- pca_data[pca_data[, "sex"] == sex, ]
    } else {
      sexed_data <- pca_data[pca_data[, "sex"] %in% c("M_All", "F_All"), ]
      drop_depaup <- drop_depaup * 2 # because the function will think there are twice the number of species when there's a male and female representative
    }
    
    return(calc_group_div(sexed_data, taxo_level, n_dim, drop_depaup, avg_par = avg_par, metric = metric))
    
  }, pca_data = pca_dat, taxo_level = tax_level, n_dim = n_dim, drop_depaup = drop_depaup, avg_par = avg_par, metric = metric
)

# rbind list elements
within_group_diversity <- data.table::rbindlist(within_group_diversity)
# change All sex to All
within_group_diversity[within_group_diversity$sex %in% c("M_All", "F_All"), "sex"] <- "All"


# Boxplots by sex and group

# set aspect ratio
if(clade == "Passeriformes"){
  ar <- 0.5
} else if(clade == "Neoaves"){
  ar <- 0.8
} else if(clade == "Neognaths"){
  ar <- 0.9
}

# Set x-label
if(metric == "centr-dist"){
  x_lab <- "Within-group Centroid Distance"
} else if(metric == "nn-all"){
  x_lab <- "Within-group Pairwise Distance"
} else if(metric == "sum.variances"){
  x_lab <- "Within-group Sum of Variances"
} else if(metric == "displacements"){
  x_lab <- "Relative displacement of group centroid from global centroid"
} else if(metric == "convhull.volume"){
  x_lab <- "Convex hull volume"
} else {
  x_lab <- "Within-group diversity"
}

# convert taxon subgroup to factor with order based on overall metric value
ts_levels <- within_group_diversity %>% 
  filter(
    sex == "All"
  ) %>% 
  arrange(
    desc(avg_group_metric)
  ) %>% 
  select(
    taxon_subgroup
  ) %>% 
  distinct() %>% 
  pull() %>% 
  rev()

within_group_diversity <- within_group_diversity %>% 
  mutate(
    taxon_subgroup = factor(taxon_subgroup, levels = ts_levels),
    sex = factor(sex, levels = sexes)
  )


p <- within_group_diversity %>% 
  # filter(
  #   sex == "F"
  # ) %>% 
  ggplot(aes(y = taxon_subgroup, x = group_metric, fill = avg_group_metric)) + 
  geom_boxplot(notch = FALSE, outliers = FALSE,) + 
  facet_wrap(~ sex, ncol = length(sexes)+1) +
  scale_fill_viridis_c(name = "Median\ncentroid\ndistance") +
  theme_minimal() +
  labs(x = x_lab, y = "Taxon Subgroup")

# save as png

boxplot_filename <- paste0(clade, "_patch_", space, "_", metric, "_", sex_choice, "boxplot.png")
png(
  here::here(
    "2_Patches", "4_OutputPlots", clade, "2_Diversity_measures", "2_Diversity_phylogeny",
    "within_group_diversity", paste(tax_level, "level", sep = "_"), "boxplots",
    boxplot_filename
  ), width = 1200, height = 1200*ar, res = 110
)
print(p)
dev.off()

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
tree$tip.label <- as.character(phylo_dat[, tax_level])

# Get tree trait tip data
td <- data.frame(node = tidytree::nodeid(tree, rownames(phylo_dat)),
                 avg_group_metric = phylo_dat[, "avg_group_metric"])
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

library(ggtree)
# Mean distance to group centroid
if(log_choice == TRUE) {
  p <- ggtree(phylo_dat_tree, aes(colour = log(avg_group_metric)),
                    layout = "rectangular",
                    ladderize = TRUE, continuous = "colour", size = branch_width) + 
  geom_tippoint(aes(colour = log(avg_group_metric)), size = tip_size) + 
  geom_tiplab(aes(colour = log(avg_group_metric)), offset = 1, size = text_size) +
  scale_color_viridis_c(option = "D")
} else {
  p <- ggtree(phylo_dat_tree, aes(colour = avg_group_metric),
                      layout = "rectangular",
                      ladderize = TRUE, continuous = "colour", size = branch_width) + 
    geom_tippoint(aes(colour = avg_group_metric), size = tip_size) + 
    geom_tiplab(aes(colour = avg_group_metric), offset = 1, size = text_size) +
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
    "within_group_diversity", paste(tax_level, "level", sep = "_"), "phylogenies",
    paste0(clade, "_patch_", space, "_", log_filename, "mean_group_centr_dist", "_", sex_choice, "_", "phylogeny.svg")
  ), 
  width = 7, height = 7
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
    "within_group_diversity", paste(tax_level, "level", sep = "_"), "phylogenies",
    paste0(clade, "_patch_", space, "_", log_filename, "mean_pairwise_dist", "_", sex_choice, "_", "phylogeny.svg")
  ), 
  width = 7, height = 7
)

p

dev.off()


# Plot summary statistics for each group on boxplot


# Boxplot of mean subgroup centroid distances 
# create colour mapping

# Calculate median values for each subgroup
medians_centroid <- tapply(within_group_diversity$group_centr_dist, within_group_diversity$taxon_subgroup, median)

# Assign colors based on median values using viridis palette
colours <- viridisLite:::viridis(length(medians_centroid))

# Match colors to subgroups based on their median values
cols <- colours[match(medians_centroid, sort(medians_centroid))]

# Set wider margins - bottom, left, top, right
par(mar = c(5, 13, 1, 1), 
    mgp = c(12, 1, 0))

# plot boxplot of mean distance to centroid for each taxon subgroup
boxplot(group_centr_dist ~ taxon_subgroup, data = within_group_diversity,
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
