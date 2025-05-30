# Calculate Phylogenetic Diversity SES reduction with simulated exinctions
# Robert MacDonald
# 30/05/2025

# Note that this uses an outdated R package (PhyloMeasures) that needs to be installed 
# from the CRAN archive as follows
#devtools::install_version("PhyloMeasures")

# Global and regional analyses

# Heavily based on code from Hughes et al 2022 Curr Biol
# https://doi.org/10.1016/j.cub.2022.06.018

# Global SES ----

# clear environment
rm(list=ls())

# Load libraries
library(dplyr)

# Load shared functions ----
source(
  here::here(
    "2_Patches", "2_Scripts", "R", "mapping.R"
  )
)

## EDITABLE CODE ## ----
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Passeriformes"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# select sex of interest ("all", "male_female", "male_only", "female_only", "unknown_only")
sex_interest <- "male_female"
# select type of averaging to use ("mean" or "median")
avg_par <- "median"
# select number of null distributions to generate
n_sims <- 10
# select number of trees to use
n_trees <- 10
# select whther to use liberal, conservative, or nominate IUCN data
# "liberal" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the highest threat level of the multiple species
# "conservative" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the lowest threat level of the multiple species
# "nominate" = Jetz (BirdTree) species that correspond to multiple IUCN (BirdLife) species
# are assigned the threat level of the BL species that corresponds to the nominate subspecies
iucn_type <- "nominate"


# check if tree file exists for specified number of trees - create if not
trimmed_tree_filepath <- here::here(
  "4_SharedInputData", paste0("First", n_trees, "_AllBirdsHackett1.tre")
)
if(!file.exists(trimmed_tree_filepath)){
  full_tree_filepath <- here::here(
    "4_SharedInputData", paste0("AllBirdsHackett1.tre")
    )
  ape::write.tree(
    ape::read.tree(full_tree_filepath)[1:n_trees],
    file = trimmed_tree_filepath
  )
}


# Load data ----

# Load trees
tree_filepath <- here::here("4_SharedInputData", paste0("First", n_trees, "_AllBirdsHackett1.tre"))
trees <- ape::read.tree(tree_filepath)

# Load IUCN Red List data
iucn_filename <- paste0("iucn_2024_", iucn_type, ".csv")
iucn <- read.csv(
  here::here(
    "4_SharedInputData", iucn_filename
  )
)[, c("species_birdtree", "iucn_cat")] %>% 
  rename(
    species = species_birdtree
  )

# Load species contained in our dataset (for specified clade)
# set PCA filename
pca_filename <- paste(clade, sex_match, "patches.231030.PCAcolspaces", "rds", sep = ".")
# Load PCA (species only)
all_species <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[["lab"]] %>% 
  magrittr::extract2("x") %>% 
  get_unique_spp()
  

# Data preparation and analysis ----

# trim IUCN data to only specified species
iucn <- iucn[iucn[["species"]] %in% all_species, ]

# create trimmed lists of bird communities by threat level
iucn_levels <- c("CR", "EN", "VU", "NT", "LC")
trimmed_species_lists <- list()
trimmed_species_lists$all_species <- iucn[["species"]]
trimmed_species_lists$noCR <- iucn[iucn[["iucn_cat"]] != "CR", "species"]
trimmed_species_lists$noEN <- iucn[iucn[["iucn_cat"]] != "CR" & iucn[["iucn_cat"]] != "EN", "species"]
trimmed_species_lists$noVU <- iucn[iucn[["iucn_cat"]] != "CR" & iucn[["iucn_cat"]] != "EN" & iucn[["iucn_cat"]] != "VU", "species"]
trimmed_species_lists$noNT <- iucn[iucn[["iucn_cat"]] != "CR" & iucn[["iucn_cat"]] != "EN" & iucn[["iucn_cat"]] != "VU" & iucn[["iucn_cat"]] != "NT", "species"]

# Create empty community matrix 
tips <- matrix(0, ncol = length(trimmed_species_lists$all_species), nrow = 5)
# Assign species names to column names in community matrix
colnames(tips) <- trimmed_species_lists$all_species
rownames(tips) <- names(trimmed_species_lists)

# define presence/absence in community matrix
for(comm_level in names(trimmed_species_lists)){
  
  tips[comm_level, colnames(tips) %in% trimmed_species_lists[[comm_level]]] <- 1
  
}

# Calculate standardised Phylogenetic Diversity


# Data frame to store ses.pd output (one row for each phylogeny, one column for each community)
ses_out <- data.frame(all_species = rep(NA, n_trees), noCR = NA, noEN = NA, noVU = NA, noNT = NA)
# Data frame to store pd output (one row for each phylogeny, one column for each community)
raw_out <- data.frame(all_species = rep(NA, n_trees), noCR = NA, noEN = NA, noVU = NA, noNT = NA)

# abundance weights for PD calculation (equal abundance)
ab_weights <- rep(1, length(trimmed_species_lists$all_species))
names(ab_weights) <- trimmed_species_lists$all_species

# calculate SES PD and PD
start_time <- Sys.time()
for(tree_number in 1:n_trees){
  
  # trim tree to match full species list
  tree <- ape::drop.tip(trees[[tree_number]], tip = trees[[tree_number]]$tip.label[!trees[[tree_number]]$tip_label %in% trimmed_species_lists$all_species])
  
  # calculate SES PD
  ses_out[, tree_number] <- PhyloMeasures::pd.query(
    tree = tree, 
    matrix = t(tips), 
    null.model = "sequential",
    abundance.weights = ab_weights,
    reps = n_sims,
    standardize = TRUE
    )
  cat("/rCompleted tree ", tree_number, " of ", n_trees)
  
}
end_time <- Sys.time()

end_time - start_time