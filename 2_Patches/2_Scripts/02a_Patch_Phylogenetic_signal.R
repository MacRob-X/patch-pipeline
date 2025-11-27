# Estimate phylogenetic signal in patch data PCA
# 25th March 2024
# Robert MacDonald

rm(list=ls())

library(motmot)
library(ape)
library(geomorph)
library(ggplot2)
library(dplyr)

## EDITAggplot2## EDITABLE CODE ## ----
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neognaths"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# select sex of interest ("all", "male_female", "male_only", "female_only", "unknown_only")
sex_interest <- "male_female"
# Choose space to work with (usmldbl, usmldblr, xyz, xyzlum, xyzlumr, lab, cie, sRGB, hex, 
# jndxyz, jndxyzlum, jndxyzlumr)
space <- "lab"
# select sex of interest ("M" or "F")
sex <- "F"

# Functions ----

# get list of species (with sex) included in PCA data
get_spp_sex <- function(div_data){
  
  # extract species and sex
  spp_sex <- t(sapply(strsplit(rownames(div_data), split = "-"), "[", 1:2))
  
  # set column names
  colnames(spp_sex) <- c("species", "sex")
  
  return(spp_sex)
  
}


# Load and prepare data ----

# # load phylogeny (full trees 0001 to 1000, Hackett backbone - from birdtree.org)
# hackettTrees <- read.tree("./4_SharedInputData/AllBirdsHackett1.tre")
# 
# # trim to first 100 trees (0001 to 0100)
# phy <- hackettTrees[1:100]
# 
# # save first 100 trees
# write.tree(phy, file = "./4_SharedInputData/First100_AllBirdsHackett1.tre")

# load first 100 trees
phy <- read.tree("./4_SharedInputData/First100_AllBirdsHackett1.tre")

# load PCA patch data 
pca_filename <- paste(clade, sex_match,  "patches.250716.PCAcolspaces.rds", sep = ".")
pcaAll <- readRDS(paste0("./2_Patches/3_OutputData/", clade, "/2_PCA_ColourPattern_spaces/1_Raw_PCA/", pca_filename))[[space]]

# add species name and sex as columns in pca data
pcaDat <- as.data.frame(pcaAll$x)
pcaDat <- cbind(pcaDat, get_spp_sex(pcaDat))

# trim to males/females only
pcaDat <- pcaDat[pcaDat$sex == sex,]

# drop species not in the phylogeny (there shouldn't actually be any of these)
pcaDat <- pcaDat[pcaDat$species %in% phy[[1]]$tip.label, ]

# match phylogeny to species in dataset
phy <- drop.tip(phy, tip = setdiff(phy[[1]]$tip.label, pcaDat$species))

# rename rows as species names and drop species/sex rows
rownames(pcaDat) <- pcaDat$species
pcaDat <- subset(pcaDat, select = -c(sex, species))

#--------------------------------------------------------------------#

# estimate pagel's lambda (phylogenetic signal) for each PC

# dummy variable 1
lambda.col <- vector("list", length = length(phy))
# dummy variable 2
lambda.ml <- rep(list(lambda.col), times = ncol(pcaDat))

for(i in 1:ncol(pcaDat)){
  print(i)
  traitData <- as.matrix(pcaDat[, i])
  rownames(traitData) <- rownames(pcaDat)
  
  # calculate lambda for PC axis i for each tree j in distribution
  ## Note - I might be better off using lapply for this
  for(j in 1:length(phy)){
    cat("\r", j)
    lambda.col[[j]] <-  transformPhylo.ML(phy = phy[[j]], y = traitData, model = "lambda")
  }
  
  lambda.ml[[i]] <- lambda.col
}

saveRDS(lambda.ml, file = "./Outputs/patch.M.jndxyzlumrPCA.lambdaML.rds")

#### Parallelised version

library(parallel)

# Define the function to calculate lambda for a given PC axis across each phylogenetic tree in distribution
calc_lambda <- function(col, data, tree_distrib) {
  traitData <- as.matrix(data[, col])
  rownames(traitData) <- rownames(data)
  parallel::parLapply(cl, 1:length(tree_distrib), function(j) {
    cat("\rProcessing tree", j, "for PC axis", col)
    transformPhylo.ML(phy = tree_distrib[[j]], y = traitData, model = "lambda")
  })
}

# Set up a cluster using the number of cores (4 less than total number of laptop cores)
no_cores <- detectCores() - 4
cl <- makeCluster(no_cores)

# Export the variables and functions needed by the workers
clusterExport(cl, varlist = c("phy", "pcaDat", "transformPhylo.ML", "calc_lambda"))

# Apply the function across all PC axes
lambda.ml <- lapply(
  colnames(pcaDat),
  calc_lambda,
  data = pcaDat, tree_distrib = phy
)
# or do as for loop, which allows displaying which PC it's up to
lambda.ml <- vector("list", length = ncol(pcaDat))
for(i in 1:ncol(pcaDat)) {
  print(paste("Processing PC axis", i))
  lambda.ml[[i]] <- calc_lambda(i, pcaDat, phy)
}

# Stop the cluster after use
stopCluster(cl)

# save
lambda_filepath <- here::here(
  "2_Patches", "3_OutputData", clade, "3_Phylogenetic_signal", space
)
if(!dir.exists(lambda_filepath)){
  dir.create(lambda_filepath, recursive = TRUE)
}
lambda_filename <- paste("pagelsLambda", sex, "hackettTrees100.rds", sep = ".")
saveRDS(lambda.ml, file = paste(lambda_filepath, lambda_filename, sep = "/"))


# Load pagel's lambda for each PC for each sex
lambda.ml <- vector(mode = "list", length = 2)
names(lambda.ml) <- c("M", "F")
for (s in names(lambda.ml)){
  lambda_filename <- paste("pagelsLambda", s, "hackettTrees100.rds", sep = ".")
  lambda.ml[[s]] <- readRDS(paste(lambda_filepath, lambda_filename, sep = "/"))
}
lambda.ml.M <- readRDS(paste(lambda_filepath, lambda_filename, sep = "/"))

# calculate mean lambda and 95% CIs for each PC (i.e. average over the tree distribution)
lambda.ml.means <- vector("list", length = 2)
names(lambda.ml.means) <- c("M", "F")
for(s in names(lambda.ml.means)){
  
  lambda.ml.means.sexed <- vector("list", length = ncol(pcaDat))
  
  for(i in 1:ncol(pcaDat)){
    
    # add mean distribution to list
    lambda.ml.means.sexed[[i]]$means <- c(rep(NA, times = length(phy)))
    for(j in 1:length(phy)){
      lambda.ml.means.sexed[[i]]$means[j] <- lambda.ml[[s]][[i]][[j]]$Lambda[[1]]
    }
    # add grand mean to list
    lambda.ml.means.sexed[[i]]$grandMean <- mean(lambda.ml.means.sexed[[i]]$means)
    # add sd to list
    lambda.ml.means.sexed[[i]]$sd <- sd(lambda.ml.means.sexed[[i]]$means)
    # add 95% CIs to list
    lambda.ml.means.sexed[[i]]$CIs <- c(mean(lambda.ml.means.sexed[[i]]$means + 1.96 * (lambda.ml.means.sexed[[i]]$sd / sqrt(length(phy)))), 
                                  mean(lambda.ml.means.sexed[[i]]$means - 1.96 * (lambda.ml.means.sexed[[i]]$sd / sqrt(length(phy)))))
    
  }
  lambda.ml.means[[s]] <- lambda.ml.means.sexed
}

# convert to dataframe for plotting
lambda.ml.means.df <- data.frame(
  PC = rep(1:ncol(pcaDat), times = 2), 
  sex = c(rep("M", times = ncol(pcaDat)), rep("F", times = ncol(pcaDat))),
  mean = rep(NA, times = ncol(pcaDat) * 2),
  lowerCI = rep(NA, times = ncol(pcaDat) * 2),
  upperCI = rep(NA, times = ncol(pcaDat) * 2)
  )
for(s in c("M", "F")){
  
  for(i in 1:ncol(pcaDat)){
    lambda.ml.means.df[lambda.ml.means.df$sex == s, ]$mean[i] <- lambda.ml.means[[s]][[i]]$grandMean
    lambda.ml.means.df[lambda.ml.means.df$sex == s, ]$lowerCI[i] <- lambda.ml.means[[s]][[i]]$CIs[[1]]
    lambda.ml.means.df[lambda.ml.means.df$sex == s, ]$upperCI[i] <- lambda.ml.means[[s]][[i]]$CIs[[2]]
  }
  
}


# plot lambda plus CIs for each PC
lambda.ml.means.df <- data.frame(PC = 1:ncol(pcaDat))
for(i in 1:ncol(pcaDat)){
  lambda.ml.means.df$mean[i] <- lambda.ml.means[[i]]$grandMean
  lambda.ml.means.df$lowerCI[i] <- lambda.ml.means[[i]]$CIs[[1]]
  lambda.ml.means.df$upperCI[i] <- lambda.ml.means[[i]]$CIs[[2]]
}
plot_filename <- paste("pagelsLambda", sex, "hackettTrees10.png", sep = ".")
png(filename = here::here("2_Patches", "4_OutputPlots", "4_Phylogenetic_signal", space,plot_filename), width = 1000, height = 800)
plot(mean ~ PC, data = lambda.ml.means.df,
     xlab = "PC", ylab = "Mean Lambda across 10 trees",
     ylim = c(0, 1))
arrows(x0 = lambda.ml.means.df$PC,
       y0 = lambda.ml.means.df$lowerCI, 
       x1 =  lambda.ml.means.df$PC,
       y1 = lambda.ml.means.df$upperCI,
       angle = 90, code = 3, length = 0.1)
lines(x = lambda.ml.means.df$PC, y = lambda.ml.means.df$mean)
dev.off()


# ggplot version, with a line for each sex
p <- ggplot(filter(lambda.ml.means.df, sex == "M"), aes(x = PC, y = mean, colour = sex)) + 
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI)) + 
  geom_line(data = filter(lambda.ml.means.df, sex == "F"), aes(x = PC, y = mean, colour = sex)) + 
  geom_point(data = filter(lambda.ml.means.df, sex == "F"), aes(x = PC, y = mean, colour = sex)) + 
  geom_errorbar(data = filter(lambda.ml.means.df, sex == "F"), aes(ymin = lowerCI, ymax = upperCI)) + 
  # facet_wrap(~sex, ncol = 2) +
  labs(x = "Principal Component", y = "Mean \u03bb across 100 trees") + 
  theme_bw() + 
  scale_x_continuous(breaks = 1:nrow(lambda.ml.means.df)) + 
  ylim(0, max(lambda.ml.means.df$mean)) + 
  ggstats::geom_stripped_cols(colour = NA)

plot_filename <- "pagelsLambda_hackettTrees100.png"
ggsave(filename = plot_filename, plot = p, width = 1000, height = 800, units = "px", path = here::here("2_Patches", "4_OutputPlots", clade, "4_Phylogenetic_signal", space), dpi = 150)

screeplot(pcaAll)

# check whether male or female data contains higher phylogenetic signal for each PC axis
lambda.ml.means.M$mean - lambda.ml.means.F$mean
# looks like generally more signal in male data but it's very close (except for PC 26/27)



#--------------------------------------------------------------------#
# Plot overall phylogenetic signal for all PCs
# using generalised Blomberg's K for multidimensional data (Adams, 2014)

# first generate majority-rule consensus tree from all Hackett trees
phy_all <- read.tree("./4_SharedInputData/AllBirdsHackett1.tre")
# MRC tree is the recommended type of consensus tree
# See https://doi.org/10.1093/czoolo/61.6.959

cons_tree <- ape::consensus(phy_all, p = 0.5, rooted = TRUE)

# save
ape::write.tree(cons_tree, file= here::here("4_SharedInputData", "AllBirdsHackett1_majrule_consensustree.txt"))

# toy data to test
smallDat <- pcaDat[sample(1:nrow(pcaDat), 1000, replace = FALSE), 30, drop = FALSE]
smallPhy <- drop.tip(phy[[1]], tip = setdiff(phy[[1]]$tip.label, get_spp_sex(smallDat)[, "species"]))

# Note that the following takes a loooooong time to run (> 24 hours on a high spec laptop)
kMult <- geomorph::physignal(A = as.matrix(smallDat), phy = smallPhy, iter = 999, print.progress = TRUE)

# full dataset
kMult <- geomorph::physignal(A = as.matrix(pcaDat), phy = phy[[1]], iter = 999, print.progress = TRUE)

print(kMult)
plot(kMult)

kmult_filepath <- here::here(
  "2_Patches", "3_OutputData", "3_Phylogenetic_signal", space
)
if(!dir.exists(kmult_filepath)){
  dir.create(kmult_filepath, recursive = TRUE)
}
kmult_filename <- paste("kMult", sex, "hackettTrees1.rds", sep = ".")
saveRDS(kMult, file = paste(kmult_filepath, kmult_filename, sep = "/"))


# Try using multivariate Pagel's Lambda, implemented in motmot
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13343


# Calculate Phylogenetic Autocorrelation Function (ACF) with castor package ----
# https://squidlobster.r-universe.dev/castor/doc/manual.html#get_trait_acf
# castor is specifically designed to work fast with big trees https://doi.org/10.1093/bioinformatics/btx701

# calculate ACF for each PC

# check the tips of the tree match the data
identical(rownames(pcaDat), phy[[1]]$tip.label)

# check the tip order is identical in all trees
identicalValue <- function(x,y) if (identical(x,y)) x else FALSE
check_vals <- lapply(phy, "[[", 4)
Reduce(identicalValue, check_vals) %>% 
  isFALSE() %>% 
  any()
# order of tips is identical in all trees, so only need to match once

# they don't match - reorder the data to match
pcaDat <- pcaDat[match(phy[[1]]$tip.label, rownames(pcaDat)), ]
identical(rownames(pcaDat), phy[[1]]$tip.label)


library(parallel)
library(castor)

# Define the function to calculate ACF for a given PC axis across each phylogenetic tree in distribution
calc_acf <- function(col, data, tree_distrib, parallelise = FALSE) {
  traitData <- data[, col]
  names(traitData) <- rownames(data)
  n_trees <- length(tree_distrib)
  if(parallelise == TRUE){
    parallel::parLapply(cl, 1:n_trees, function(tree_num) {
      cat("\rProcessing tree", tree_num, "for PC axis", col)
      return(get_trait_acf(tree = tree_distrib[[tree_num]], tip_states = traitData))
    })
  } else {
    lapply(1:n_trees, function(tree_num){
      cat("\rProcessing tree", tree_num, "for PC axis", col)
      return(get_trait_acf(tree = tree_distrib[[tree_num]], tip_states = traitData))
    }
    )
  }
  
}


# apply across each tree for each axis in a for loop, which allows displaying which PC it's up to
acf <- vector("list", length = ncol(pcaDat))
for(i in 1:ncol(pcaDat)) {
  print(paste("Processing PC axis", i))
  acf[[i]] <- calc_acf(i, pcaDat, phy, parallelise = FALSE)
}

names(acf) <- colnames(pcaDat)


# transform results from lists to df

# Create all combinations of PC and tree
pc_tree_combs <- expand.grid(
  pc = 1:30,
  tree = 1:100,
  stringsAsFactors = FALSE
)

acf_df <- purrr::pmap_dfr(pc_tree_combs, function(pc, tree) {
  pc_name <- paste0("PC", pc)
  
  # check if successful acf calculation for sub-element
  if (!is.null(acf[[pc_name]][[tree]]) && acf[[pc_name]][[tree]]$success) {
    
    # create dataframe populated with PC/tree combination values if successful
    data.frame(
      phylodist = acf[[pc_name]][[tree]]$phylodistances,
      autocorr = acf[[pc_name]][[tree]]$autocorrelations,
      axis = pc_name,
      tree = tree,
      stringsAsFactors = FALSE
    )
    
  } else {
    
    # return empty dataframe if unsuccessful or missing
    data.frame(
      phylodist = numeric(0),
      autocorr = numeric(0),
      axis = character(0),
      tree = integer(0),
      stringsAsFactors = FALSE
    )
    
  }
})


library(ggplot2)
library(dplyr)

# plot values for a given tree
acf_df %>%
  filter(tree == 1) %>% 
  ggplot(aes(x = phylodist, y = autocorr)) + 
  geom_point() + 
  geom_smooth(method = "lm", formula = y ~ x) + 
  facet_wrap(~ axis)

# plot all values (across trees) for a given axis
acf_df %>%
  filter(axis == "PC2") %>% 
  ggplot(aes(x = phylodist, y = autocorr)) + 
  geom_point() + 
  geom_smooth(method = "lm", formula = y ~ x) + 
  facet_wrap(~ tree)
