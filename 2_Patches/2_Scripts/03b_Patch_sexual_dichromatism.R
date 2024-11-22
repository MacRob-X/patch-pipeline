# Calculate sexual dichromatism measure from patch data
# Robert MacDonald
# 28th March 2024

# clear environment
rm(list=ls())

# load libraries
library(dplyr)
library(ggplot2)
library(ggtree)

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
pca_all <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    "Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603.rds"
  )
)

# load patch umap
umap <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
    "Neoaves.patches.pca.jndxyzlum.nnUMAPs.25.35.45.55.240603.rds"
  )
) %>% 
  magrittr::extract2(
    1
  )


# prepare data and phylogeny ----

# create working version of PCA data with species and sex as columns
pca_dat <- as.data.frame(pca_all$x) %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
#    sex = janitor::make_clean_names(sapply(strsplit(rownames(.), split = "-"), "[", 2), allow_dupes = T)
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  )


# create separate male and female PCA data frames with species as rownames
male_pca_dat <- pca_dat[pca_dat$sex == "M",]
female_pca_dat <- pca_dat[pca_dat$sex == "F",]

# drop species without data for females data
male_pca_dat <- male_pca_dat[male_pca_dat$species %in% female_pca_dat$species, ]
# drop species without data for males
female_pca_dat <- female_pca_dat[female_pca_dat$species %in% male_pca_dat$species, ]

# rename rows as species and drop species and sex columns
rownames(male_pca_dat) <- male_pca_dat$species
male_pca_dat <- subset(male_pca_dat, select = -c(sex, species))
rownames(female_pca_dat) <- female_pca_dat$species
female_pca_dat <- subset(female_pca_dat, select = -c(sex, species))

# reorder female data to match male order
female_pca_dat <- female_pca_dat[match(rownames(male_pca_dat), rownames(female_pca_dat)), , drop = FALSE]
# check reordering
identical(rownames(male_pca_dat), rownames(female_pca_dat))


# calculate pairwise euclidean distance between males and females of species in PCA colour pattern space
# as measure of sexual dichromatism

# create empty dataframe to populate
dichromatism <- data.frame(distance = rep(NA, times = nrow(male_pca_dat)), row.names = rownames(male_pca_dat))

# calculate distance for each male-female pair
for(i in 1:nrow(male_pca_dat)){
  cat("\r", i)
  # combine each male-female pair into separate dataframe
  df <- rbind(male_pca_dat[i, ], female_pca_dat[i, ])
  # calculate multidimensional euclidean distance between pair and populate distances df
  dichromatism[i,] <- dist(df, method = "euclidean")
}
rm(df)

# plot pairwise euclidean distances as measure of dichromatism on phylogeny ----

# load first 100 trees of Hackett backbone full trees
phy <- read.tree(
  here::here(
    "4_SharedInputData", "First100_AllBirdsHackett1.tre"
  )
)

# choose one tree for plotting purposes
tree <- phy[[1]]

# trim data to match tree and vice versa (shouldn't really be any species that aren't in the phylogeny anyway
# but best to check)
dichromatism <- dichromatism[rownames(dichromatism) %in% tree$tip.label, , drop = F]
tree <- drop.tip(tree, tip = setdiff(tree$tip.label, rownames(dichromatism)))
# reorder distance data to match phylogeny
dichromatism <- dichromatism[match(tree$tip.label, rownames(dichromatism)), , drop = F]

# calculate ancestral state (for plotting purposes only)
dichro_ancest_state_bm <- phytools::fastAnc(tree, dichromatism$distance, CI = TRUE, vars = TRUE)

# produce df to use for plotting (can log transform dichromatism$distance and dichro_ancest_state_bm$ace if desired)
td <- data.frame(node = tidytree::nodeid(tree, rownames(dichromatism)),
                 distance = log(dichromatism$distance))
nd <- data.frame(node = names(dichro_ancest_state_bm$ace), distance = log(dichro_ancest_state_bm$ace))
d <- rbind(td, nd)
d$node <- as.numeric(d$node)

# put data and tree together with full_join()
dichro_tree <- full_join(tree, d, by = "node")
rm(td, nd, d)

# Plot ancestral states on phylogeny with colour according to centroid distance
# Note that I'm not actually interested in the ancestral states but I think it will make the plot look nicer
# and I guess it could show something interesting
p <- ggtree::ggtree(dichro_tree, aes(color = distance),
                    layout = "circular",
            ladderize = TRUE, continuous = "colour", size = 0.35) +    # adjust size here to make branches thinner
  scale_color_gradientn(colours = viridisLite::viridis(10)) # gradientn allows us to specify how many colours to use in scale 
p <- p + ggtree::geom_tiplab(offset = 1, size = 0.55) + theme_tree2() + # Adjust the size parameter here to make the tip label text smaller
  theme(panel.grid.major = element_line(color = "black", linewidth = .2),
        panel.grid.minor = element_line(color = "lightgrey", linewidth = .2),
        panel.grid.major.y = element_blank(),
        legend.position.inside = c(0.15, 0.1))
revts(p)  # reverses timescale by setting most recent tip to 0
# save as portrait 100 x 10 inch pdf


# Plot dichromatism on colour pattern space ----

# add dichromatism score as column in PCA data
dichromatism$species <- rownames(dichromatism)

pca_dichro <- pca_dat %>% 
  filter(
 #   sex == "f",
    sex != "u") %>% 
  right_join(dichromatism, by = "species") %>% 
  rename(dichromatism = distance)
# add rownames to data
rownames(pca_dichro) <- paste(pca_dichro$species, pca_dichro$sex, sep = "-")

pca_dichro %>% 
ggplot(aes(x = PC3, y = PC4, alpha = dichromatism)) + 
  geom_point() + 
  facet_wrap(facets = "sex") + 
  theme_bw()


# Plot colour grids with dichromatism as transparency ----

## Editable Code ##
## Choose sex to plot
sex <- "M"
## Choose axes to plot
pc_axis_1 <- "PC1"
pc_axis_2 <- "PC2"
## End Editable Code ##


# create single-sex dfs
pca_dichro_dat <- pca_dichro[pca_dichro$sex == sex, ]

# get proportions of variance for each axis
pca_variance <- pca_all %>% 
  summary() %>% 
   magrittr::extract2("importance") %>% 
  as.data.frame()


# get x and y ranges
xrange <- c((min(pca_dichro[, pc_axis_1]) - 0.5 ), (max(pca_dichro[, pc_axis_1]) + 0.5 ))
yrange <- c((min(pca_dichro[, pc_axis_2]) - 0.5 ), (max(pca_dichro[, pc_axis_2]) + 0.5 ))

# set axis labels
xlabel <- paste0(pc_axis_1, " (", round(pca_variance[2, pc_axis_1], 2) * 100, "% of variance)")
ylabel <- paste0(pc_axis_2, " (", round(pca_variance[2, pc_axis_2], 2) * 100, "% of variance)")

png(
  here::here(
    "2_Patches", "4_OutputPlots", "4_Dichromatism", "jndxyzlum",
    paste0("patches_jndxyzlum_pca_", pc_axis_1, "_", pc_axis_2, "_dichromatism_", sex, "_colour_grids.png")
  ),
  width = 4020, height = 2568,
  units = "px"
)

plot(pca_dichro_dat[, pc_axis_1] ~ pca_dichro_dat[, pc_axis_2], asp = T, type = "n", xlab = xlabel, ylab = ylabel, las = 1,
     xlim = xrange, ylim = yrange)

# box(lwd = 2)
rw <- diff(range(pca_dichro[,1]))/12

for(i in 1:nrow(pca_dichro_dat)){
  fname <- rownames(pca_dichro_dat)[i]
  fpng <- png::readPNG(paste0(
    "C:/Users/bop23rxm/Documents/colour_grids_repositioned/", 
    paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
  # adjust transparency according to dichromatism
  # set transparency of image as ratio of dichro of interest to max dichro
  alpha <- pca_dichro_dat$dichromatism[i] / max(pca_dichro_dat$dichromatism)
  # add extra channel (alpha) to image for transparency
  fpng <- array(c(fpng, array(alpha, dim = dim(fpng)[1:2])), dim = c(dim(fpng)[1:2], 4))
  rasterImage(fpng, 
              xleft = pca_dichro_dat[, pc_axis_1][i] - (rw/15), 
              ybottom = pca_dichro_dat[, pc_axis_2][i]-(rw/9), 
              xright = pca_dichro_dat[, pc_axis_1][i]+(rw/15), 
              ytop = pca_dichro_dat[, pc_axis_2][i]+(rw/9))
  cat(paste0("\r", i, " of ", nrow(pca_dichro_dat), " processed"))
}

dev.off()



# Same for UMAP

# add dichromatism scores to umap data
umap_dichro <- umap %>% 
  magrittr::extract2(
    "layout"
  ) %>% 
  as.data.frame() %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    #    sex = janitor::make_clean_names(sapply(strsplit(rownames(.), split = "-"), "[", 2), allow_dupes = T)
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  ) %>% 
  filter(
    #   sex == "f",
    sex != "U") %>% 
  right_join(dichromatism, by = "species") %>% 
  rename(
    dichromatism = distance,
    UMAP1 = V1,
    UMAP2 = V2
    )

umap_dichro %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, alpha = dichromatism)) + 
  geom_point() + 
  facet_wrap(facets = "sex") + 
  theme_bw()

# Plot UMAP colour grids with dichromatism as transparency
# add species-sex rownames
rownames(umap_dichro) <- paste(umap_dichro$species, umap_dichro$sex, sep = "-")

# get sex-specific datasets
umap_dichro_m <- umap_dichro[umap_dichro$sex == "M", ]
umap_dichro_f <- umap_dichro[umap_dichro$sex == "F", ]

## Editable Code ##
## Choose sex
sex <- "F"
## End Editable Code ##

# filter dataset based on sex
umap_dat <- umap_dichro[umap_dichro$sex == sex, ]


# get x and y ranges
xrange <- c((min(umap_dichro$UMAP1) - 0.2 ), (max(umap_dichro$UMAP1) + 0.2 ))
yrange <- c((min(umap_dichro$UMAP2) - 0.2 ), (max(umap_dichro$UMAP2) + 0.2 ))

# set axis labels
xlabel <- "UMAP 1"
ylabel <- "UMAP 2"

png(
  here::here(
    "2_Patches", "4_OutputPlots", "4_Dichromatism", "jndxyzlum",
    paste0("patches_jndxyzlum_umap_dichromatism_", sex, "_colour_grids.png")
  ),
  width = 4020, height = 2568,
  units = "px"
)

# change plot background coluor to make white phenotypes easier to see
par(bg = "#F5FFFF")

plot(umap_dat$UMAP1 ~ umap_dat$UMAP2, asp = T, type = "n", xlab = xlabel, ylab = ylabel, las = 1,
     xlim = xrange, ylim = yrange)

# box(lwd = 2)
rw <- diff(range(umap_dichro[,1]))/8

for(i in 1:nrow(umap_dat)){
  fname <- rownames(umap_dat)[i]
  fpng <- png::readPNG(paste0(
    "C:/Users/bop23rxm/Documents/colour_grids_repositioned/", 
    paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
  # adjust transparency according to dichromatism
  # set transparency of image as ratio of dichro of interest to max dichro
  alpha <- umap_dat$dichromatism[i] / max(umap_dat$dichromatism)
  # add extra channel (alpha) to image for transparency
  fpng <- array(c(fpng, array(alpha, dim = dim(fpng)[1:2])), dim = c(dim(fpng)[1:2], 4))
  rasterImage(fpng, 
              xleft = umap_dat$UMAP1[i] - (rw/15), 
              ybottom = umap_dat$UMAP2[i]-(rw/9), 
              xright = umap_dat$UMAP1[i]+(rw/15), 
              ytop = umap_dat$UMAP2[i]+(rw/9))
  cat(paste0("\r", i, " of ", nrow(umap_dat), " processed"))
}

dev.off()


# check if my measure of dichromatism is correlated with a standard TCS xyz measure (from Cooney et al 2019) -----
# load Cooney et al data (Tyrannidae only - individual patch data)
cooney_data <- read.csv(
  here::here(
    "4_SharedInputData", "Cooney_etal_2019", "Analysis_Data.csv"
  )
)

# sum patch dichromatism scores to give total dichromatism score (as in the paper)
# add column for total dichro score

# get species included in study and create dataframe to use
dichro_cooney <- data.frame(species = unique(cooney_data$species), 
                            cooney_dichro = rep(NA, times = length(unique(cooney_data$species))),
                            cooney_max_dichro = rep(NA, times = length(unique(cooney_data$species))))

# generate total dichro score for each species and also extract max dichro score
for(i in 1:length(dichro_cooney$species)){
  cat("\r", i)
  df <- cooney_data[cooney_data$species == dichro_cooney$species[i],]
  # assign max dichro score
  dichro_cooney$cooney_max_dichro[i] <- max(df$max.dichro, na.rm = T) # I could equally use min() or even df$cooney_max_dichro[1] here - there should only be one max dichromatism score for each species
  # calculate and assign total dichro score
  dichro_cooney$cooney_dichro[i] <- sum(df$dichro)
}


# trim my dichromatism data to match and attach
temp <- dichromatism[match(dichro_cooney$species, rownames(dichromatism)), , drop = F]
# check species order matches
identical(dichro_cooney$species, rownames(temp))
# log10 transform and attach dichromatism data (since Cooney dichro data are log10 transformed)
dichro_cooney$colpatspace_dichro <- log10(temp$distance)
# remove temp df
rm(temp)

# now check correlation between cooney total dichro and my jndxyzlumr dichro scores
mod_dichro <- lm(cooney_dichro ~ colpatspace_dichro, data = dichro_cooney)
summary(mod_dichro)
# very highly significant correlation (p < 2e-16) but my dichro measure only explains ~ 1/2 variance in cooney dichro (Adj R2 = 0.5377)
mod_dichro_2 <- lm(colpatspace_dichro ~ cooney_dichro, data = dichro_cooney)
summary(mod_dichro_2)


# test for association between IUCN level and jndxyzlum dichromatism ----

# load IUCN Red List data
iucn <- readr::read_rds(
  here::here(
    "4_SharedInputData", "IUCN_RedList_data_130324.rds"
  )
)

# add species name separated by _ to working IUCN dataframe
iucnData <- iucn
iucnData$species <- gsub(" ", "_", iucnData$scientific_name)   # note that I'll need to sort out the taxonomy before I do this for real

# trim species categorised as DD, EX, EW, or CR(PE) (after Hughes et al 2022)
levels(as.factor(iucnData$category))
iucnData <- iucnData[iucnData$category != "DD" & iucnData$category != "EX" & iucnData$category != "EW" & iucnData$category != "CR(PE)", ]
levels(as.factor(iucnData$category))

# remove species not present in IUCN data
jndDichro <- dichromatism[rownames(dichromatism) %in% iucnData$species, , drop = F]

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

summary(iucnMod)

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
iucnModPC <- phylolm(iucnPseudCont ~ distance, data = jndDichro, phy = tree, model = "BM")
summary(iucnModPC)

# bin iucn into threatened or not threatened binary category
# Create binary IUCN category (LC, NT = 0, EX, CR, EN, VU = 1, DD = NA)
jndDichro$iucnBinCat <- replace(as.character(jndDichro$iucnCat), jndDichro$iucnCat == "LC" | jndDichro$iucnCat == "NT", 0)
jndDichro$iucnBinCat <- replace(as.character(jndDichro$iucnBinCat), jndDichro$iucnBinCat == "EX" | jndDichro$iucnBinCat == "CR" | jndDichro$iucnBinCat == "EN" | jndDichro$iucnBinCat == "VU", 1)

# run logistic regression with binary IUCN category
mod_dichroBin <- phyloglm(iucnBinCat ~ distance, data = jndDichro, phy = tree)
summary(mod_dichroBin)

# reverse it
iucnModBin <- phylolm(distance ~ iucnBinCat, data = jndDichro, phy = tree, model = "BM")
summary(iucnModBin)

# reverse the model
mod_dichroPC <- phylolm(iucnPseudCont ~ distance, data = jndDichro, phy = tree, model = "BM")
summary(mod_dichroPC)
# looks like there is an extremely small but significant effect (adj R2 = 0.0007313. p = 0.023)


# Check if there's a latitudinal gradient in sexual dichromatism ----
cooney22Data <- read.csv(
  here::here(
    "4_SharedInputData", "Cooney_etal_2022", "data_cooney_etal_latitude.csv"
  )
)

# prepare data

# trim latitude data to match colour data
latData <- cooney22Data[cooney22Data$Binomial %in% rownames(dichromatism), , drop = F]

# trim to males only
latData <- subset(latData, Sex == "M")

# trim colour data to match latitude data
temp <- dichromatism[match(latData$Binomial, rownames(dichromatism)), , drop = F]

# check species order matches
identical(latData$Binomial, rownames(temp))
# attach dichromatism data 
latData$colpatspace_dichro <- temp$distance
# remove temp df
rm(temp)

# check correlation between Cooney 2022 qualitative dichromatism score and euclidean distances in jndxyzlumr pattern space
mod_dichro_22 <- lm(Dichro ~ colpatspace_dichro, data = latData)
summary(mod_dichro_22)
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
latmod_dichro <- MCMCglmm(colpatspace_dichro ~ abs(MidpointLat),                    # Specify regression equation
                random = ~ Binomial,                   # Set phylogeny as a random factor
                ginverse = list(Binomial = animalA),  
                prior = prior.linear,                  # Set prior
                verbose = TRUE,                        # Print MH diagnostics to screen
                family = "gaussian",                   # For linear regression
                data = latData,                          # Specify dataset 
                nitt=11000,                            # Number of MCMC iterations
                thin=10,                               # Thinning interval
                burnin=1000,                           # Discard unconverged values
                pl=TRUE,                               # Save posterior distribution of latent (inferred) variables
                pr=TRUE)    

summary(latmod_dichro)                
plot(latmod_dichro)
