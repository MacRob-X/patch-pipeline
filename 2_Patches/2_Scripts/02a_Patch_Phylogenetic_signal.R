# Estimate phylogenetic signal in patch data PCA
# 25th March 2024
# Robert MacDonald

rm(list=ls())

library(motmot)
library(ape)
library(geomorph)

## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Passeriformes"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
space <- "lab"
# select sex of interest ("M" or "F")
sex <- "M"

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

# load phylogeny (full trees 0001 to 1000, Hackett backbone - from birdtree.org)
hackettTrees <- read.tree("./4_SharedInputData/AllBirdsHackett1.tre")

# trim to first 100 trees (0001 to 0100)
phy <- hackettTrees[1:100]

# save first 100 trees
write.tree(phy, file = "./4_SharedInputData/First100_AllBirdsHackett1.tre")

# load first 100 trees
phy <- read.tree("./4_SharedInputData/First100_AllBirdsHackett1.tre")

# load PCA patch data 
pca_filename <- paste(clade, "patches.231030.PCAcolspaces.rds", sep = ".")
pcaAll <- readRDS(paste0("./2_Patches/3_OutputData/2_PCA_ColourPattern_spaces/1_Raw_PCA/", pca_filename))[[space]]

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
lambda.col <- rep(list(transformPhylo.ML(phy = phy[[1]], y = as.matrix(pcaDat[, 1]), model = "lambda")), times = length(phy))
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

# Set up a cluster using the number of cores (4 less than total number of laptop cores)
no_cores <- detectCores() - 4
cl <- makeCluster(no_cores)

# Export the variables and functions needed by the workers
clusterExport(cl, varlist = c("phy", "pcaDat", "transformPhylo.ML"))

# Define the function to calculate lambda for a given PC axis across each phylogenetic tree in distribution
calc_lambda <- function(col, data, tree_distrib) {
  traitData <- as.matrix(data[, col])
  rownames(traitData) <- rownames(data)
  parLapply(cl, 1:length(tree_distrib), function(j) {
    cat("\rProcessing tree", j, "for PC axis", col)
    transformPhylo.ML(phy = tree_distrib[[j]], y = traitData, model = "lambda")
  })
}

# dummy variable for storing results
lambda.ml <- vector("list", length = ncol(pcaDat))

# Apply the function across all PC axes
lambda.ml <- apply(as.matrix(pcaDat), 2, calc_lambda)
for(i in 1:ncol(pcaDat)) {
  print(paste("Processing PC axis", i))
  lambda.ml[[i]] <- calc_lambda(i, pcaDat, phy)
}

# Stop the cluster after use
stopCluster(cl)

# save
saveRDS(lambda.ml, file = paste("./2_Patches/3_OutputData/3_Phylogenetic_signal/patch", sex, space, "lambdaML.rds", sep = "."))


# Load pagel's lambda for each PC
lambda.ml <- readRDS("./2_Patches/3_OutputData/3_Phylogenetic_signal/patch.M.jndxyzlumrPCA.lambdaML.rds")

# calculate mean lambda and 95% CIs for each PC (i.e. average over the tree distribution)
lambda.ml.means <- vector("list", length = ncol(pcaDat))

for(i in 1:ncol(pcaDat)){
  
  # add mean distribution to list
  lambda.ml.means[[i]]$means <- c(rep(NA, times = length(phy)))
  for(j in 1:length(phy)){
    lambda.ml.means[[i]]$means[j] <- lambda.ml[[i]][[j]]$Lambda[[1]]
  }
  # add grand mean to list
  lambda.ml.means[[i]]$grandMean <- mean(lambda.ml.means[[i]]$means)
  # add sd to list
  lambda.ml.means[[i]]$sd <- sd(lambda.ml.means[[i]]$means)
  # add 95% CIs to list
  lambda.ml.means[[i]]$CIs <- c(mean(lambda.ml.means[[i]]$means + 1.96 * (lambda.ml.means[[i]]$sd / sqrt(length(phy)))), 
                                mean(lambda.ml.means[[i]]$means - 1.96 * (lambda.ml.means[[i]]$sd / sqrt(length(phy)))))
  
}





# plot lambda plus CIs for each PC
lambda.ml.means.df <- data.frame(PC = c(1:40))
for(i in 1:length(lambda.ml.means.df$mean)){
  lambda.ml.means.df$mean[i] <- lambda.ml.means[[i]]$grandMean
  lambda.ml.means.df$lowerCI[i] <- lambda.ml.means[[i]]$CIs[[1]]
  lambda.ml.means.df$upperCI[i] <- lambda.ml.means[[i]]$CIs[[2]]
}
plot(mean ~ PC, data = lambda.ml.means.df,
     xlab = "PC", ylab = "Mean Lambda across 100 trees",
     ylim = c(0, 1))
arrows(x0 = lambda.ml.means.df$PC,
       y0 = lambda.ml.means.df$lowerCI, 
       x1 =  lambda.ml.means.df$PC,
       y1 = lambda.ml.means.df$upperCI,
       angle = 90, code = 3, length = 0.1)
lines(x = lambda.ml.means.df$PC, y = lambda.ml.means.df$mean)

screeplot(pcaAll)

# check whether male or female data contains higher phylogenetic signal for each PC axis
lambda.ml.means.M$mean - lambda.ml.means.F$mean
# looks like generally more signal in male data but it's very close (except for PC 17)



#--------------------------------------------------------------------#
# Plot overall phylogenetic signal for all PCs
# using generalised Blomberg's K for multidimensional data (Adams, 2014)

# Matrix version of PCA
pcaDat <- as.matrix(pcaDat)

# Note that the following takes a loooooong time to run (> 24 hours on a high spec laptop)
kMult <- physignal(A = pcaDat, phy = phy[[1]])

saveRDS(kMult, file = "./Outputs/patch.M.jndxyzlumrPCA.kMult.rds")
