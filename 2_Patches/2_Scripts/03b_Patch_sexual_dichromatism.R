# Calculate sexual dichromatism measure from patch data
# Robert MacDonald
# 28th March 2024

# clear environment
rm(list=ls())

# load libraries
library(dispRity)
library(ape)
library(caper)
library(MCMCglmm)
library(phylolm)

# custom functions
source(
  here::here(
    "2_Patches", "2_Scripts", "2_BetaVersions", "Patch_Plotting_functions_v1.r"
  )
)

# Load data ----

# load patch data
pcaAll <- readRDS("./Outputs/features/patches/PCA/Neoaves.patches.pca.jndxyzlumr.rds")



# prepare data and phylogeny ----

# create working version of PCA data
pcaDat <- as.data.frame(pcaAll$x)

# split PCA data into separate dataframes for males and females
# add species name and sex as columns in PCA data
pcaDat$species <- c(rep(NA, times = length(pcaDat[,1])))
pcaDat$sex <- c(rep(NA, times = length(pcaDat[,1])))
for(i in 1:length(pcaDat[,1])){
  cat("\r", i)
  split <- strsplit(rownames(pcaDat)[i], split = "-")[[1]]
  pcaDat$species[i] <- split[1]
  pcaDat$sex[i] <- split[2]
}
rm(split)

# create separate male and female PCA data frames with species as rownames
malePCADat <- pcaDat[pcaDat$sex == "M",]
femalePCADat <- pcaDat[pcaDat$sex == "F",]

# drop species without data for females data
malePCADat <- malePCADat[malePCADat$species %in% femalePCADat$species, ]
# drop species without data for males
femalePCADat <- femalePCADat[femalePCADat$species %in% malePCADat$species, ]

# rename rows as species and drop species and sex columns
rownames(malePCADat) <- malePCADat$species
malePCADat <- subset(malePCADat, select = -c(sex, species))
rownames(femalePCADat) <- femalePCADat$species
femalePCADat <- subset(femalePCADat, select = -c(sex, species))

# reorder female data to match male order
femalePCADat <- femalePCADat[match(rownames(malePCADat), rownames(femalePCADat)), , drop = FALSE]
# check reordering
identical(rownames(malePCADat), rownames(femalePCADat))


# calculate pairwise euclidean distance between males and females of species in PCA colour pattern space
# as measure of sexual dichromatism

# create empty dataframe to populate
eucDistances <- data.frame(distance = rep(NA, times = nrow(malePCADat)), row.names = rownames(malePCADat))

# calculate distance for each male-female pair
for(i in 1:nrow(malePCADat)){
  cat("\r", i)
  # combine each male-female pair into separate dataframe
  df <- rbind(malePCADat[i, ], femalePCADat[i, ])
  # calculate multidimensional euclidean distance between pair and populate distances df
  eucDistances[i,] <- dist(df, method = "euclidean")
}
rm(df)

# plot pairwise euclidean distances as measure of dichromatism on phylogeny ----

# load first 100 trees of Hackett backbone full trees
phy <- read.tree("./Data/First100_AllBirdsHackett1.tre")

# choose one tree for plotting purposes
tree <- phy[[1]]

# trim data to match tree and vice versa (shouldn't really be any species that aren't in the phylogeny anyway
# but best to check)
eucDistances <- eucDistances[rownames(eucDistances) %in% tree$tip.label, , drop = F]
tree <- drop.tip(tree, tip = setdiff(tree$tip.label, rownames(eucDistances)))
# reorder distance data to match phylogeny
eucDistances <- eucDistances[match(tree$tip.label, rownames(eucDistances)), , drop = F]

# calculate ancestral state (for plotting purposes only)
ancDichroBM <- fastAnc(tree, eucDistances$distance, CI = TRUE, vars = TRUE)

# produce df to use for plotting (can log transform eucDistances$distance and ancDichroBM$ace if desired)
td <- data.frame(node = nodeid(tree, rownames(eucDistances)),
                 distance = log(eucDistances$distance))
nd <- data.frame(node = names(ancDichroBM$ace), distance = log(ancDichroBM$ace))
d <- rbind(td, nd)
d$node <- as.numeric(d$node)

# put data and tree together with full_join()
eucDistTree <- full_join(tree, d, by = "node")
rm(td, nd, d)

# Plot ancestral states on phylogeny with colour according to centroid distance
# Note that I'm not actually interested in the ancestral states but I think it will make the plot look nicer
# and I guess it could show something interesting
p <- ggtree(eucDistTree, aes(color = distance),
            ladderize = TRUE, continuous = "colour", size = 0.35) +    # adjust size here to make branches thinner
  scale_color_gradientn(colours = viridisLite::viridis(10)) # gradientn allows us to specify how many colours to use in scale 
p <- p + geom_tiplab(offset = 1, size = 0.55) + theme_tree2() + # Adjust the size parameter here to make the tip label text smaller
  theme(panel.grid.major = element_line(color = "black", linewidth = .2),
        panel.grid.minor = element_line(color = "lightgrey", linewidth = .2),
        panel.grid.major.y = element_blank(),
        legend.position.inside = c(0.15, 0.1))
revts(p)  # reverses timescale by setting most recent tip to 0
# save as portrait 100 x 10 inch pdf




# check if my measure of dichromatism is correlated with a standard TCS xyz measure (from Cooney et al 2019) -----
# load Cooney et al data (Tyrannidae only - individual patch data)
cooneyData <- read.csv("./Data/Cooney_etal_2019/Analysis_Data.csv")

# sum patch dichromatism scores to give total dichromatism score (as in the paper)
# add column for total dichro score

# get species included in study and create dataframe to use
dichroData <- data.frame(species = unique(cooneyData$species), 
                            cooneyDichro = rep(NA, times = length(unique(cooneyData$species))),
                            cooneyMaxDichro = rep(NA, times = length(unique(cooneyData$species))))

# generate total dichro score for each species and also extract max dichro score
for(i in 1:length(dichroData$species)){
  cat("\r", i)
  df <- cooneyData[cooneyData$species == dichroData$species[i],]
  # assign max dichro score
  dichroData$cooneyMaxDichro[i] <- max(df$max.dichro, na.rm = T) # I could equally use min() or even df$cooneyMaxDichro[1] here - there should only be one max dichromatism score for each species
  # calculate and assign total dichro score
  dichroData$cooneyDichro[i] <- sum(df$dichro)
}


# trim my eucDistances data to match and attach
temp <- eucDistances[match(dichroData$species, rownames(eucDistances)), , drop = F]
# check species order matches
identical(dichroData$species, rownames(temp))
# log10 transform and attach eucDistances data (since Cooney dichro data are log10 transformed)
dichroData$jndxyzlumrDichro <- log10(temp$distance)
# remove temp df
rm(temp)

# now check correlation between cooney total dichro and my jndxyzlumr dichro scores
dichroMod <- lm(cooneyDichro ~ jndxyzlumrDichro, data = dichroData)
summary(dichroMod)
# very highly significant correlation (p < 2e-16) but my dichro measure only explains ~ 1/2 variance in cooney dichro (Adj R2 = 0.5377)
dichroMod2 <- lm(jndxyzlumrDichro ~ cooneyDichro, data = dichroData)
summary(dichroMod2)


# test for association between IUCN level and jndxyzlumr dichromatism ----

# load IUCN Red List data
iucn <- readRDS("./Data/IUCN_RedList_data_130324.rds")

# add species name separated by _ to working IUCN dataframe
iucnData <- iucn
iucnData$species <- gsub(" ", "_", iucnData$scientific_name)   # note that I'll need to sort out the taxonomy before I do this for real

# trim species categorised as DD, EX, EW, or CR(PE) (after Hughes et al 2022)
levels(as.factor(iucnData$category))
iucnData <- iucnData[iucnData$category != "DD" & iucnData$category != "EX" & iucnData$category != "EW" & iucnData$category != "CR(PE)", ]
levels(as.factor(iucnData$category))

# remove species not present in IUCN data
jndDichro <- eucDistances[rownames(eucDistances) %in% iucnData$species, , drop = F]

# reorder IUCN data to match dichro Data
iucnData <- iucnData[match(rownames(jndDichro), iucnData$species), ]

# Check all specimens match between the two datasets
if(identical(iucnData$species, rownames(jndDichro))){
  print("All specimens in IUCN and colour data matched")
}

 # attach iucn data to jnd dichro data
jndDichro$iucnCat <- as.factor(iucnData$category)
# manually set levels so that base level is LC
levels(jndDichro$iucnCat) <- c("LC", "NT", "VU", "EN", "CR")

# log10 transform the distance data to improve normality
jndDichro$distance <- log10(jndDichro$distance)

# add species name as a column
jndDichro$species <- rownames(jndDichro)

# let's run a basic single tree MCMCglmm to start with, with iucn category as a multilevel predictor (maybe not appropriate)
# choose one phylogenetic tree
tree <- phy[[1]]

# trim data to match tree and vice versa (shouldn't really be any species that aren't in the phylogeny anyway
# but best to check)
jndDichro <- jndDichro[rownames(jndDichro) %in% tree$tip.label, , drop = F]
tree <- drop.tip(tree, tip = setdiff(tree$tip.label, rownames(jndDichro)))

# Invert covariance matrices (to specify random effects)
animalA <- inverseA(tree)$Ainv

# set a prior (inverse wishart prior for linear regression, residual and phylogenetic variance - nu could also be 0.002)
prior.linear<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))

# run the regression
iucnMod <- MCMCglmm(distance ~ iucnCat,
                    data = jndDichro,
                    random = ~ species,
                    ginverse = list(species = animalA),
                    prior = prior.linear,
                    verbose = T,
                    family = "gaussian",
                    nitt=11000,   # [number of MCMC iterations]
                    thin=10,      # [thinning interval - reduces memory use, and maybe helps to reduce autocorrelation of sample and increase independence]
                    burnin=1000,  # [discard unconverged values]
                    pl=TRUE,      # [save posterior distribution of latent (inferred) variables]
                    pr=TRUE,      # [save posterior distribution of random variables]
                    slice=TRUE)


# try a basic phylolm instead

iucnMod <- phylolm(distance ~ iucnCat, data = jndDichro, phy = tree, model = "BM")
summary(iucnMod)

# let's see if we can make iucn category a pseudocontinuous variable for a bit more clarity
jndDichro$iucnPseudCont <- gsub(pattern = "LC", replacement = 0, x = jndDichro$iucnCat)
jndDichro$iucnPseudCont <- gsub(pattern = "NT", replacement = 1, x = jndDichro$iucnPseudCont)
jndDichro$iucnPseudCont <- gsub(pattern = "VU", replacement = 2, x = jndDichro$iucnPseudCont)
jndDichro$iucnPseudCont <- gsub(pattern = "EN", replacement = 3, x = jndDichro$iucnPseudCont)
jndDichro$iucnPseudCont <- gsub(pattern = "CR", replacement = 4, x = jndDichro$iucnPseudCont)
jndDichro$iucnPseudCont <- as.numeric(jndDichro$iucnPseudCont)

iucnModPC <- phylolm(distance ~ iucnPseudCont, data = jndDichro, phy = tree, model = "BM")
summary(iucnModPC)

# bin iucn into threatened or not threatened binary category
# Create binary IUCN category (LC, NT = 0, EX, CR, EN, VU = 1, DD = NA)
jndDichro$iucnBinCat <- replace(as.character(jndDichro$iucnCat), jndDichro$iucnCat == "LC" | jndDichro$iucnCat == "NT", 0)
jndDichro$iucnBinCat <- replace(as.character(jndDichro$iucnBinCat), jndDichro$iucnBinCat == "EX" | jndDichro$iucnBinCat == "CR" | jndDichro$iucnBinCat == "EN" | jndDichro$iucnBinCat == "VU", 1)

# run logistic regression with binary IUCN category
dichroModBin <- phyloglm(iucnBinCat ~ distance, data = jndDichro, phy = tree)
summary(dichroModBin)

# reverse it
iucnModBin <- phylolm(distance ~ iucnBinCat, data = jndDichro, phy = tree, model = "BM")
summary(iucnModBin)

# reverse the model
dichroModPC <- phylolm(iucnPseudCont ~ distance, data = jndDichro, phy = tree, model = "BM")
summary(dichroModPC)
# looks like there is an extremely small but significant effect (adj R2 = 0.0007313. p = 0.023)


# Check if there's a latitudinal gradient in sexual dichromatism ----
cooney22Data <- read.csv("./Data/Cooney_etal_2022/data_cooney_etal_latitude.csv")

# prepare data

# trim latitude data to match colour data
latData <- cooney22Data[cooney22Data$Binomial %in% rownames(eucDistances), , drop = F]

# trim to males only
latData <- subset(latData, Sex == "M")

# trim colour data to match latitude data
temp <- eucDistances[match(latData$Binomial, rownames(eucDistances)), , drop = F]

# check species order matches
identical(latData$Binomial, rownames(temp))
# attach eucDistances data 
latData$jndxyzlumrDichro <- temp$distance
# remove temp df
rm(temp)

# check correlation between Cooney 2022 qualitative dichromatism score and euclidean distances in jndxyzlumr pattern space
dichroMod22 <- lm(Dichro ~ jndxyzlumrDichro, data = latData)
summary(dichroMod22)
# highly correlated (p < 2e-16, adj R2 = 0.453)

# test for association between latitude and jndxyzlumr dichro (phylogenetically controlled)
# - is latitude a predictor of dichromatism?

# specify prior for a linear regression
prior.linear <- list(G = list(G1 = list(V = 1, nu = 0.02)), R = list(V = 1, nu = 0.02))

# choose one phylogenetic tree
tree <- phy[[1]]

# Invert covariance matrices (to specify random effects)
animalA <- inverseA(tree)$Ainv


# Run the regression using MCMCglmm 
latDichroMod <- MCMCglmm(jndxyzlumrDichro ~ abs(MidpointLat),                    # Specify regression equation
                random = ~ Binomial,                   # Set phylogeny as a random factor
                ginverse = list(Binomial = animalA),  
                prior = prior.linear,                  # Set prior
                verbose = TRUE,                        # Print MH diagnostics to screen
                family = "gaussian",                   # For linear regression
                data = latData,                          # Specify dataset 
                nitt=110000,                            # Number of MCMC iterations
                thin=100,                               # Thinning interval
                burnin=10000,                           # Discard unconverged values
                pl=TRUE,                               # Save posterior distribution of latent (inferred) variables
                pr=TRUE)    

summary(latDichroMod)                
plot(latDichroMod)
