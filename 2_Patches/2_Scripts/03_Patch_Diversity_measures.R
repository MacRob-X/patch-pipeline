# Extracting measures of diversity from patch colour pattern spaces
# 13th March 2024
# Robert MacDonald

# clear environment
rm(list=ls())

# load libraries
library(dplyr)
library(ggfortify)

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

match_iucn <- function(colour_data, iucn_data, sex = "all"){
  
  # check if IUCN data contains necessary columns
  if(!("scientific_name" %in% colnames(iucn_data))){
    stop("Please provide IUCN data which contains Latin species names in a column titled 'scientific_name'")
  }
  
  # add species name separated by _ to IUCN data
  iucn_data$species <- gsub(" ", "_", iucn_data$scientific_name)   # note that I'll need to sort out the taxonomy before I do this for real
  iucn_data$specimen <- c(rep(NA, times = length(iucn_data$species)))
  
  if(sex == "F"){
    
    # get female specimen IUCN data
    fem_iucn <- iucn_data[paste0(iucn_data$species, "-F") %in% rownames(colour_data), ]
    fem_iucn$specimen <- paste0(fem_iucn$species, "-F")
    
    # assign to dataframe
    iucn_match <- fem_iucn
    
  } else if (sex == "M") {
    
    # get male specimen IUCN data
    male_iucn <- iucn_data[paste0(iucn_data$species, "-M") %in% rownames(colour_data), ]
    male_iucn$specimen <- paste0(male_iucn$species, "-M")
    
    # assign to dataframe
    iucn_match <- male_iucn
    
  } else if (sex == "all") {
    
    # get female specimen IUCN data
    fem_iucn <- iucn_data[paste0(iucn_data$species, "-F") %in% rownames(colour_data), ]
    fem_iucn$specimen <- paste0(fem_iucn$species, "-F")
    
    # get male specimen IUCN data
    male_iucn <- iucn_data[paste0(iucn_data$species, "-M") %in% rownames(colour_data), ]
    male_iucn$specimen <- paste0(male_iucn$species, "-M")
    
    # get unknown sex specimen IUCN data
    unk_iucn <- iucn_data[paste0(iucn_data$species, "-U") %in% rownames(colour_data), ]
    unk_iucn$specimen <- paste0(unk_iucn$species, "-U")
    
    # concatenate into one big dataframe with specimen names that match colour_data rownames
    iucn_match <- rbind(fem_iucn, male_iucn, unk_iucn)
  }

  # trim species categorised as DD, EX, EW, or CR(PE) (after Hughes et al 2022)
  levels(as.factor(iucn_match$category))
  iucn_match <- iucn_match[iucn_match$category != "DD" & iucn_match$category != "EX" & iucn_match$category != "EW" & iucn_match$category != "CR(PE)", ]
  levels(as.factor(iucn_match$category))
  
  # remove species not present in IUCN data
  colour_data <- colour_data[rownames(colour_data) %in% iucn_match$specimen,]
  
  # reorder IUCN data to match colour_data
  iucn_match <- iucn_match[match(rownames(colour_data), iucn_match$specimen), ]
  
  # Check all specimens match between the two datasets
  if(identical(iucn_match$specimen, rownames(colour_data))){
    print("All specimens in IUCN and colour data matched")
  }
  
  matched_data <- list(iucn_match, colour_data)
  
  return(matched_data)
  
}

# dispRity-style function to calculate mean distance to given number of nearest neighbours
mean.nn.dist <- function(matrix, nn = NULL, method = "euclidean") {
  
  ## Set number of neighbours to all if not specified
  if(is.null(nn)){
    nn <- nrow(matrix)
  }
  
  ## Calculate all pairwise distances
  pair.dists <- as.matrix(vegan::vegdist(matrix, method = method))
  
  ## Get nn closest for each point
  mean.nn.dists <- apply(pair.dists, 1, function(one.row, nn) mean(sort(one.row)[2:nn + 1]), nn = nn)
  
  ## Return values
  return(mean.nn.dists)
}

# dispRity-style function to count number of neighbours within a given radius (written by Thomas Guillerme)
count.neighbours <- function(matrix, radius = 1, relative = TRUE, method = "euclidean") {
  ## Check if the matrix is a distance matrix first
  distances <- as.matrix(dist(matrix, method = method))
  ## Set the radius to something if it's a function
  if(is(radius, "function")) {
    radius <- radius(distances)
  }
  ## For each row count how many distances are < radius (minus one is for the diagonal that's = 0)
  counts <- apply(distances, 1, function(one_row, radius) sum(one_row <= radius), radius = radius) - 1
  ## Return the counts
  if(relative) {
    return(unname(counts/ncol(distances)))
  } else {
    return(unname(counts))
  }
}

#-----------------------------------------------------------------------------------------#

# Load and prepare data ----

# load patch data
pca_all <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    "Neoaves.patches.231030.PCAcolspaces.jndxyzlumr.240603.rds"
  )
)

# load umap patch data (output from 02b_Patch_Umap_iterations.R)
umap_jndxyzlumr <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
    "Neoaves.patches.pca.jndxyzlumr.UMAPs.iterations.240603.rds"
  )
) %>% 
  magrittr::extract2(1) %>% 
  magrittr::extract2("layout") %>% 
  as.data.frame() %>% 
  rename(
    umap_1 = V1, umap_2 = V2
  )

# load IUCN Red List data
iucn <- readr::read_rds(
  here::here(
    "4_SharedInputData", "IUCN_RedList_data_130324.rds"
  )
)


#--------------------------------------------------------------#


# Analyses (largely based on code from Hughes et al 2022) ----

# Mean distance to centroid by IUCN category - boxplot

# get data (could use first 16 PCs (~97% cumulative percentage of variance) if want to simplify)
pca_dat <- pca_all$x

# match iucn and PCA data
matched_data <- match_iucn(pca_dat, iucn, sex = "all")
iucn_dat <- as.data.frame(matched_data[1])
pca_dat <- as.data.frame(matched_data[2])
rm(matched_data)

# convert data to matrix for analyses
pca_dat <- as.matrix(pca_dat)

# get distances to centroid for each species/sex pair
centroid_distances <- as.data.frame(dispRity::dispRity(pca_dat, metric = mean.nn.dist, nn = 5)$disparity[[1]][[1]])
rownames(centroid_distances) <- rownames(pca_dat)
colnames(centroid_distances) <- c("centr_dist")
glob_dist <- centroid_distances

# check density of centroid distances
ggplot() + 
  geom_density(data = centroid_distances, aes(x = centr_dist))
ggplot() + 
  geom_histogram(data = centroid_distances, aes(x = centr_dist), bins = 50)

# attach IUCN categories to distance data
glob_dist <- glob_dist %>% 
  mutate(
    iucn = iucn_dat$category[match(rownames(glob_dist), iucn_dat$specimen)]
  )
# note that since I've already matched the order the match() function is redundant here, but I'm leaving it in as a failsafe

# add duplicate of each specimen with iucn cat "All" (for plotting box plot)
all_dist <- glob_dist
all_dist$iucn <- c("All")

all <- rbind(glob_dist, all_dist)

# add species names and sex as columns
all <- all %>% 
mutate(
  species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
  sex = sapply(strsplit(rownames(.), split = "-"), "[", 2),
)

# reorder levels of iucn cat for plotting
all$iucn <- factor(all$iucn, levels = c("CR", "EN", "VU", "NT", "LC", "All"))

# create colour mapping
cols <- springfield[1:length(levels(all$iucn))]
names(cols) <- levels(all$iucn)

# plot mean distance to centroid of each IUCN category (boxplot)
# ggplot version (can subset by sex)
all %>% 
#  filter(
#  sex == "M",
#  iucn != "All") %>% 
ggplot(aes(x = iucn, y = centr_dist, fill = iucn)) + 
  geom_boxplot(notch = TRUE, show.legend = FALSE) + 
  scale_fill_manual(values = cols) + 
  xlab("Subset of IUCN categories") + ylab("Distance to centroid") +
#  ylim(c(0, 90)) +
  theme_minimal()
  

# ok, so it looks like there might be some kind of difference between the mean distance to centroid (i.e. a measure 
# of colour diversity) in different IUCN categories - let's just do a 
# basic ANOVA ----
# to find out

# reorder levels of iucn cat for modelling (with "LC" as the reference category once 'All' cat is removed)
all$iucn <- factor(all$iucn, levels = c("All", "LC", "CR", "EN", "VU", "NT"))

# make a model (will automatically perform an anova with a categorical predictor) with log transformed centroid distances to fulfil assumptions
mod <- lm(log(centr_dist) ~ iucn, data = all %>% filter(iucn != "All"))

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


# plot IUCN categories on UMAP colourspace
# match iucn and umap data
matched_data <- match_iucn(umap_jndxyzlumr, iucn)
iucn_dat <- as.data.frame(matched_data[1])
umap_dat <- as.data.frame(matched_data[2]) %>% 
  mutate(
    iucn_cat = factor(iucn_dat$category, levels = c("CR", "EN", "VU", "NT", "LC"))
  )
rm(matched_data)

# set colour mapping
cols <- springfield[1:length(levels(umap_dat$iucn_cat))]

# plot
ggplot(umap_dat, aes(x = umap_1, y = umap_2, colour = iucn_cat)) + 
  geom_point() + 
  facet_wrap(~ iucn_cat) + 
  theme_bw()

# and as hexbin
ggplot(umap_dat, aes(x = umap_1, y = umap_2)) + 
  geom_hex() + 
  facet_wrap(~ iucn_cat) + 
  theme_bw()


# plot different IUCN categories in black with others in transparent grey
plots <- vector("list", length =  length(levels(umap_dat$iucn_cat)))

for(i in 1:length(levels(umap_dat$iucn_cat))){
  cat <- levels(umap_dat$iucn_cat)[i]
  plots[i] <- list(
    ggplot(umap_dat, aes(x = umap_1, y = umap_2)) + 
      geom_point(
        data = umap_dat %>% 
          filter(
            iucn_cat != cat
          ),
        aes(x = umap_1, y = umap_2), colour = "skyblue", alpha = 1/8, size = 0.9) + 
      geom_point(
        data = umap_dat %>% 
          filter(
            iucn_cat == cat
          ),
        aes(x = umap_1, y = umap_2), colour = "midnightblue", alpha = 2/3, size = 0.9) + 
      xlab("UMAP axis 1") + ylab("UMAP axis 2") + 
      theme_bw()
  )
}

# plot all on one figure
ggpubr::ggarrange(
  plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], 
  labels = c("A", "B", "C", "D", "E"),
  nrow = 2, ncol = 3
)


## Calculate convex hull volume, sequentially dropping each IUCN category ----
## Note that we will only use the first 8 PCs as otherwise it's too computationally intensive

# attach iucn data to pca data
conv_dat <- pca_dat %>% 
  as.data.frame() %>% 
  mutate(
    iucn = factor(iucn_dat$category[match(rownames(pca_dat), iucn_dat$specimen)], 
                     levels = c("CR", "EN", "VU", "NT", "LC"))
  )

# set up empty list
conv_hull_vols <- data.frame(
  "iucn_subset" = factor(c("All", "-CR", "-EN", "-VU"), levels = c("All", "-CR", "-EN", "-VU")),
  "volume" = c(rep(NA, length = length(levels(conv_dat$iucn)) - 1))
)

# volume of all species
conv_hull_vols$volume[1] <- geometry::convhulln(as.matrix(conv_dat[, 1:8]), 
                                             output.options = TRUE)$vol
# remove CR species
conv_hull_vols$volume[2] <- geometry::convhulln(as.matrix(filter(conv_dat, iucn != "CR")[, 1:8]), 
                                             output.options = TRUE)$vol

# remove CR and EN species
conv_hull_vols$volume[3] <- geometry::convhulln(as.matrix(filter(conv_dat, iucn != "CR" & iucn != "EN")[, 1:8]), 
                                             output.options = TRUE)$vol

# remove CR, EN, and VU species (i.e. leave only NT species)
conv_hull_vols$volume[4] <- geometry::convhulln(as.matrix(filter(conv_dat, iucn != "CR" & iucn != "EN" & iucn != "VU")[, 1:8]), 
                                             output.options = TRUE)$vol

# add column to represent change as proportion of original volume
conv_hull_vols$rel_volume <- conv_hull_vols$volume / max(conv_hull_vols$volume)

# plot results
# create colour mapping
cols <- springfield[1:length(conv_hull_vols)]
names(cols) <- names(conv_hull_vols)

# plot mean distance to centroid of each IUCN category (boxplot)
# ggplot version (can subset by sex)
conv_hull_vols %>% 
  ggplot(aes(x = iucn_subset, y = rel_volume, colour = iucn_subset, size = rel_volume)) + 
  geom_point(notch = TRUE, show.legend = FALSE) + 
  scale_fill_manual(values = cols) + 
  scale_size_continuous(range = c(10, 20)) +
  xlab("Subset of IUCN categories") + ylab("Convex hull hypervolume") +
  #  ylim(c(0, 90)) +
  theme_minimal()





## Just for fun, let's look at distance to centroid vs midpoint of latitudinal range ----
lat_data <- readr::read_csv(
  here::here(
    "4_SharedInputData", "Cooney_etal_2022", "data_cooney_etal_latitude.csv"
  )
)

# get female specimen lat data
fem_lat <- lat_data[paste0(lat_data$Binomial, "-F") %in% rownames(pca_dat), ]
fem_lat$specimen <- paste0(fem_lat$Binomial, "-F")
# get male specimen lat data
male_lat <- lat_data[paste0(lat_data$Binomial, "-M") %in% rownames(pca_dat), ]
male_lat$specimen <- paste0(male_lat$Binomial, "-M")
# get unknown sex specimen lat data
unk_lat <- lat_data[paste0(lat_data$Binomial, "-U") %in% rownames(pca_dat), ]
unk_lat$specimen <- paste0(unk_lat$Binomial, "-U") # looks like there are no species in the lat data that have unknown sex in the patch data
# concatenate into one big dataframe with specimen names that match pca_dat rownames
lats <- rbind(fem_lat, male_lat)
# remove species not present in lat data
pca_dat <- pca_dat[rownames(pca_dat) %in% lats$specimen,]
# reorder IUCN data to match pca_dat
lats <- lats[match(rownames(pca_dat), lats$specimen), ]
# Check all specimens match between the two datasets
identical(lats$specimen, rownames(pca_dat))
# remove temporary datasets
rm(list = c("fem_lat", "male_lat", "unk_lat"))

# get distances to centroid for each specimen
centroid_distances <- as.data.frame(dispRity::dispRity(pca_dat, metric = dispRity::centroids)$disparity[[1]][[1]])
rownames(centroid_distances) <- rownames(pca_dat)
colnames(centroid_distances) <- c("centr_dist")
lat_dists <- centroid_distances

# check centroid distance distribution
ggplot() + 
  geom_histogram(data = centroid_distances, aes(x = centr_dist), bins = 60)
# looks quite skewed, so might want to log-transform
ggplot() + 
  geom_histogram(data = centroid_distances, aes(x = log(centr_dist)), bins = 60)

# attach latitude data to centroid distances data
lat_dists$midpoint_lat <- lats$MidpointLat

# plot scatterplot of latitude vs distance to centroid
scatter.smooth(y = lat_dists$centr_dist, x = lat_dists$midpoint_lat, 
     frame = F,
     xlab = "Mean latitude of species range", ylab = "Distance to centroid")

ggplot(lat_dists, aes(x = midpoint_lat, y = centr_dist)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlab("Mean latitude of species range") + ylab("Distance to centroid")

# with log-transformed centroid distances
ggplot(lat_dists, aes(x = midpoint_lat, y = log(centr_dist))) + 
  geom_point() + 
  geom_smooth(method = "lm")+ 
  xlab("Mean latitude of species range") + ylab("Distance to centroid (log-transformed")

# and plotting absolute value of latitude (to see if species closer to equator have higher distance to centroid)
ggplot(lat_dists, aes(x = abs(midpoint_lat), y = log(centr_dist))) + 
  geom_point() + 
  geom_smooth(method = "lm")+ 
  xlab("Absolute mean latitude of species range") + ylab("Distance to centroid (log-transformed")

latMod <- lm(log(centr_dist) ~ abs(midpoint_lat), data = lat_dists)
# i should really check the assumptions of the model are met here but as I'm just messing around I'll skip it
summary(latMod)
# looks like distance to centroid declines with increasing (absolute) latitude, and the result is highly significant, although
# with an R2 or 0.015, it explains little of the observed variance. So latitude affects uniqueness of colour patterning, but
# only has a very small effect compared to other things



# Colour different IUCN categories on UMAP ----
umapAll <- umap::umap(pca_all$x)

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
matched_data <- match_iucn(umapMan, iucn)
iucn_dat <- as.data.frame(matched_data[1])
umapDat <- as.data.frame(matched_data[2])
rm(matched_data, umapMan, l, split)

# attach IUCN categories to umap data
umapDat$iucn <- as.factor(iucn_dat$category)

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
pca_dat <- pca_all$x

# match iucn and PCA data
matched_data <- match_iucn(pca_dat, iucn)
iucn_dat <- as.data.frame(matched_data[1])
pca_dat <- as.data.frame(matched_data[2])
rm(matched_data)

# convert data to matrix for analyses
pca_dat <- as.matrix(pca_dat)

# get overall functionall evenness
minSpanEven <- dispRity::dispRity(pca_dat, metric = dispRity::func.div)

summary(minSpanEven)


## by taxon subgroup

# convert pca data to data frame
pca_dat <- as.data.frame(pca_dat)

# taxonomy data
taxo <- read.csv("./4_SharedInputData/BLIOCPhyloMasterTax_2019_10_28.csv")

# trim and order taxonomy data to match PCA data
taxo$specimen <- c(rep(NA, times = length(taxo$TipLabel)))

# get female specimen IUCN data
femTaxo <- taxo[paste0(taxo$TipLabel, "-F") %in% rownames(pca_dat), ]
femTaxo$specimen <- paste0(femTaxo$TipLabel, "-F")

# get male specimen IUCN data
maleTaxo <- taxo[paste0(taxo$TipLabel, "-M") %in% rownames(pca_dat), ]
maleTaxo$specimen <- paste0(maleTaxo$TipLabel, "-M")

# get unknown sex specimen IUCN data
unkTaxo <- taxo[paste0(taxo$TipLabel, "-U") %in% rownames(pca_dat), ]
unkTaxo$specimen <- paste0(unkTaxo$TipLabel, "-U")

# concatenate into one big dataframe with specimen names that match pca_dat rownames
taxoMatch <- rbind(femTaxo, maleTaxo, unkTaxo)

# remove temp dataframes
rm(femTaxo, maleTaxo, unkTaxo)

# remove species not present in IUCN data
pca_dat <- pca_dat[rownames(pca_dat) %in% taxoMatch$specimen,]

# reorder taxo data to match pca_dat
taxoMatch <- taxoMatch[match(rownames(pca_dat), taxoMatch$specimen), ]

# Check all specimens match between the two datasets
if(identical(taxoMatch$specimen, rownames(pca_dat))){
  print("All specimens in IUCN and colour data matched")
}

# attach taxon subgroups to PCA data
pca_dat$taxon_subgroup <- taxoMatch$Taxon_subgroup


# calculate minimum spanning evenness by taxon subgroup
taxon_subgroups <- unique(pca_dat$taxon_subgroup)
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
      pca_dat[pca_dat$taxon_subgroup == group, 1:ncol(pca_dat) - 1]
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
      pca_dat[pca_dat$taxon_subgroup == group, 1:ncol(pca_dat) - 1]
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
      pca_dat[pca_dat$taxon_subgroup == group, 1:ncol(pca_dat) - 1]
    ), 
    metric = dispRity::neighbours
  )$disparity[[1]][[1]] %>% 
    mean()
}


# plot functional divergence and evenness for each subgroup
library(ggplot2)

# count species in each subgroup
n_spec <- pca_dat %>% 
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
tyran_df <- pca_dat %>% 
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


