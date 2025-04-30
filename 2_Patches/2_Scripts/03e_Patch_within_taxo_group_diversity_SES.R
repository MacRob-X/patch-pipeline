## Calculate SES for diversity change in each taqxonomic separately
## Robert MacDonald
## 23rd April 2025

# clear environment
rm(list=ls())

# Load libraries
library(parallel)
library(dispRity)
library(dplyr)
library(ggplot2)
library(extrafont)


# Load shared functions ----
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "mapping.R"
  )
)

source(
  here::here(
    "3_SharedScripts", "dispRity_metric_functions.R"
  )
)

# temporary location for the SES calculating functions
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "regional_ses.R"
  )
)


## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Passeriformes"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
# N.B. will need to change the date in the pca_all filename if using usmldbl or usmldblr
# (from 240603 to 240806)
space <- "lab"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# select sex of interest ("all", "male_female", "male_only", "female_only", "unknown_only")
sex_interest <- "male_female"
# select metric ("centr-dist", "nn-k", "nn-count", "nn-all", "sum.variances", "sum.ranges")
metric <- "centr-dist"
# select type of averaging to use ("mean" or "median")
avg_par <- "median"
# select whether to calculate local or global diversity loss (i.e., mean distance to local or global centroid)
div_loss_type <- "local"
# select number of null distributions to generate
n_sims <- 100

# select whther to use liberal, conservative, or nominate IUCN data
# "liberal" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the highest threat level of the multiple species
# "conservative" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the lowest threat level of the multiple species
# "nominate" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the threat level of the BL species that corresponds to the nominate subspecies
iucn_type <- "nominate"

# Specify taxonomic level to calculate at ("genus", "family", "order", "taxon_subgroup")
tax_level <- "taxon_subgroup"
## Drop depauperate clades? Specify min. number of species a clade must contain for it to be included
## Minimum is 2, as can't calculate centroid distances or pairwise distances for groups with
## only one species
drop_depaup <- 10


## Plotting parameters
# Select whether to use binned (based on quantiles) or continuous colour scale
col_scale_type <- "binned"
# If binned, choose the number of quantiles to use
nquants <- 10
## Use viridis_c, viridis turbo palette or custom colour palette?
palette_choice <- "viridis"
## Plot on log scale?
log_choice <- FALSE


# Load data ----

# set PCA filename
pca_filename <- paste(clade, sex_match, "patches.231030.PCAcolspaces", "rds", sep = ".")
# Load PCA (values only)
pca_mat <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]]$x

# load IUCN Red List data
iucn_filename <- paste0("iucn_2024_", iucn_type, ".csv")
iucn <- read.csv(
  here::here(
    "4_SharedInputData", iucn_filename
  )
)[, c("species_birdtree", "iucn_cat")]


# load taxonomy and rename columns
taxo <- read.csv(
  here::here(
    "4_SharedInputData", "BLIOCPhyloMasterTax_2019_10_28.csv"
  ),
  strings = F
) %>% 
  rename(
    species = TipLabel,
    genus = GenusName,
    family = BLFamilyLatin,
    order = IOCOrder,
    taxon_subgroup = Taxon_subgroup
  ) %>% 
  select(
    species, genus, family, order, taxon_subgroup
  )


# Prepare data ----


# get species and sex (in same order as PCA data)
metadata <- as.data.frame(pca_spp_sex(pca_mat, bind_to_data = FALSE))

# add rownames as column (to match to PCA data)
metadata$id <- rownames(pca_mat)

# append taxonomy
metadata <- merge(metadata, taxo, by = "species", all.x = TRUE)

# append IUCN levels
metadata <- merge(metadata, iucn, by.x = "species", by.y = "species_birdtree", all.x = TRUE)

# order to match PCA data
metadata <- metadata[match(rownames(pca_mat), metadata[, "id"]), ]
rownames(metadata) <- rownames(pca_mat)

# generate individual list of sex-specific species values of metric of interest (keep or discard sexes
# depending on value of sex_match)
# first set sex list to iterate over
sexes <- set_sex_list(sex_interest)

# now apply a function to extract named vector lists of values for each sex individually
sexed_pca_list <- lapply(sexes, extract_sex_vals, div_data = pca_mat, assemblage_level = TRUE)
names(sexed_pca_list) <- sexes

# Perform analysis ----

# get list of taxa
groups <- unique(metadata[, tax_level])

# get number of taxa to apply across
n_groups <- length(groups)

# for each taxon, get SES for loss of each IUCN level
# apply function across each sexed dataset
results <- lapply(sexes, function(sex){
  
  pca_sexed <- sexed_pca_list[[sex]]
  
  # apply function across each region (row) in sf object
  lapply(groups, calc_group_ses, 
         pca_data = pca_sexed, 
         taxonomy = taxo, 
         iucn = iucn, 
         metric = metric, 
         tax_level = tax_level, 
         n_sims = n_sims,
         #cluster = cl,
         parallel_run = FALSE
         )
  
})


# name elements of the results list
names(results) <- sexes


# for each sex, process results into dfs to bind to sf data
results_dfs <- lapply(sexes, function(sex) {
  
  # get individual sexed result
  sexed_results <- results[[sex]]
  
  # bind results together into a single dataframe
  sexed_results <- as.data.frame(do.call(rbind, sexed_results))
  
  # add sex column
  sexed_results$sex <- sex
  
  # add group column
  sexed_results$group <- rep(groups, each = nrow(results[[sex]][[1]]))
  
  # convert NaN values to NA (where there are e.g. no CR species to be lost)
  sexed_results <- data.frame(lapply(sexed_results, gsub, pattern = NaN, replacement = NA, fixed = TRUE))
  
})

# bind rows together
results_all <- do.call(rbind, results_dfs)


# Plotting ----

# Make sure numaric arguments are properly converted to numeric
results_all[, c("n_species", "raw", "null_mean", "null_sd", "null_se", "ses", "n_spec_prop")] <- lapply(
  c("n_species", "raw", "null_mean", "null_sd", "null_se", "ses", "n_spec_prop"),
  function(col) as.numeric(results_all[, col])
)
  

# convert IUCN cat lost to factor and adjust levels
results_all[, "iucn"] <- factor(results_all[, "iucn"], levels = c("all", "CR", "EN", "VU", "NT"))
# same for sex (to make plots appear in All-M-F order)
results_all[, "sex"] <- factor(results_all[, "sex"], levels = c("All", "M", "F"))

# get groups with as many or more species than the specified minimum
groups_to_keep <- unique(results_all[results_all[, "iucn"] == "all" & results_all[, "n_species"] >= drop_depaup, "group"])

# convert sex label list and labelling function to change facet labels
sex_vals <- list(
  'All' = "All specimens",
  'M' = "Males",
  'F' = "Females"
)
sex_labeller <- function(variable, value){
  return(sex_vals[value])
}

p <- results_all %>% 
 # filter(sex == "M") %>% 
  # Filter out any NA values in the x variable
  filter(!is.na(n_spec_prop)) %>%
  # filter out any groups with fewer than the requested number of species
  filter(
    group %in% groups_to_keep
    ) %>% 
  ggplot() + 
  geom_point(aes(x = n_spec_prop, y = ses, colour = iucn), size = 3.5, shape = 17) + 
#  xlim(0.6, 1.01) +
  scale_x_reverse() +
  geom_line(aes(x = n_spec_prop, y = ses), alpha = 0.6) +
  scale_color_discrete(type = hcl.colors(n = nlevels(results_all$iucn) + 1, palette = "Terrain")) + 
  xlab("Remaining proportional species richness") + ylab("Standard Effect Size") + 
  labs(colour = "IUCN Status Lost") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = -2, linetype = "dashed") + 
  geom_hline(yintercept = 2, linetype = "dashed") + 
  facet_wrap(~group+sex, ncol = 4) + 
  theme_bw() + 
  theme(text=element_text(size=16,  family="Century Gothic"))


# save plot
png_filename <- paste0(paste(clade, space, tax_level, "SES", metric, avg_par, paste0(n_sims, "sims"), paste0(iucn_type, "iucn"), sex_match, sep = "_"), ".png")

# initialise png saving
png(
  here::here(
    "2_Patches", "4_OutputPlots", "2_Diversity_measures", "2_Diversity_phylogeny", "within_group_diversity", paste0(tax_level, "_level"), "SES",
    png_filename
  ), 
  width = 1450, height = 1550, res = 100, pointsize = 40
)
print(p)
dev.off()
