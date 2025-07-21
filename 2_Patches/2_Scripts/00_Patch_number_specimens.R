# load libraries
library(dplyr)

# clear environment
rm(list=ls())

# load patch pixel values (vRGB; uRGB)
px <- readRDS("./2_Patches/1_InputData/patches.231030.rds")

# load taxonomic data
taxo <- read.csv("./4_SharedInputData/BLIOCPhyloMasterTax_2019_10_28.csv")

# bind and remove any galloanseriformes or palaeognaths
px <- px %>% 
  left_join(
    taxo,
    by = join_by(species == TipLabel)
  ) %>% 
  filter(
    IOCOrder != "STRUTHIONIFORMES", 
    IOCOrder != "RHEIFORMES", 
    IOCOrder != "CASUARIIFORMES", 
    IOCOrder != "APTERYGIFORMES", 
    IOCOrder != "TINAMIFORMES", 
    IOCOrder != "GALLIFORMES", 
    IOCOrder != "ANSERIFORMES"
  )

# count unique species in our data
length(unique(px$species))

# count Neoaves species in Jetz taxonomy
taxo %>% 
  filter(
    IOCOrder != "STRUTHIONIFORMES", 
    IOCOrder != "RHEIFORMES", 
    IOCOrder != "CASUARIIFORMES", 
    IOCOrder != "APTERYGIFORMES", 
    IOCOrder != "TINAMIFORMES", 
    IOCOrder != "GALLIFORMES", 
    IOCOrder != "ANSERIFORMES"
  ) %>% 
  magrittr::extract2("TipLabel") %>% 
  unique() %>% 
  length()

# count Passeriforme species in Jetz taxonomy
taxo %>% 
  filter(
    IOCOrder == "PASSERIFORMES", 
  ) %>% 
  magrittr::extract2("TipLabel") %>% 
  unique() %>% 
  length()

# count species for which we have male and female specimens
# get species for which we have male specimens
male_species <- px %>% 
  filter(
    sex == "M"
  ) %>% 
  select(
    species
  ) %>% 
  unique()
# get species for which we have female specimens
female_species <- px %>% 
  filter(
    sex == "F"
  ) %>% 
  select(
    species
  ) %>% 
  unique()
# check how many species have both male and female data
male_species %>% 
  inner_join(
    female_species
  ) %>% 
  nrow()

# check how many passeriforme species have both male and female data
male_species %>% 
  left_join(
    taxo,
    by = join_by(species == TipLabel)
  ) %>% 
  filter(
    IOCOrder == "PASSERIFORMES"
  ) %>% 
  inner_join(
    female_species
  ) %>% 
  nrow()

# check how percentage coverage of species and genera
n_neoaves_jetz <- taxo %>% 
  filter(
    IOCOrder != "STRUTHIONIFORMES", 
    IOCOrder != "RHEIFORMES", 
    IOCOrder != "CASUARIIFORMES", 
    IOCOrder != "APTERYGIFORMES", 
    IOCOrder != "TINAMIFORMES", 
    IOCOrder != "GALLIFORMES", 
    IOCOrder != "ANSERIFORMES"
  ) %>% 
  magrittr::extract2("TipLabel") %>% 
  unique() %>% 
  length()

n_neoaves_matchedsex <- male_species %>% 
  inner_join(
    female_species
  ) %>% 
  nrow()

# percentage species coverage
(n_neoaves_matchedsex / n_neoaves_jetz) * 100

# check percentage of genera with >= 1 representative in patch data
n_neoaves_genera_jetz <- taxo %>% 
  filter(
    IOCOrder != "STRUTHIONIFORMES", 
    IOCOrder != "RHEIFORMES", 
    IOCOrder != "CASUARIIFORMES", 
    IOCOrder != "APTERYGIFORMES", 
    IOCOrder != "TINAMIFORMES", 
    IOCOrder != "GALLIFORMES", 
    IOCOrder != "ANSERIFORMES"
  ) %>% 
  magrittr::extract2("GenusName") %>% 
  unique() %>% 
  length()

n_neoaves_genera_matchedsex <- px %>% 
  magrittr::extract2("GenusName") %>% 
  unique() %>% 
  length()

(n_neoaves_genera_matchedsex / n_neoaves_genera_jetz) * 100


# Quantify percentage species coverage of IUCN categories ----

## IMPORTANT NOTE 2025-07-18
## THE FOLLOWING DOESN'T WORK, BECAUSE I ONLY HAVE IUCN CATEGORIES FOR THE SPECIES WHICH
## ARE IN THE DATASET
## so the estimate of percentage coverage for each category is an overestimate
## the solution is to redo the taxonomy matching but do it for ALL species present 
## in the taxonomy, rather than restricting it to species for which we have data
## 

# select IUCN type
iucn_type <- "nominate"

# load IUCN data
iucn_filename <- paste0("iucn_2024_", iucn_type, "_patch231030.csv")
iucn <- read.csv(
  here::here(
    "4_SharedInputData", iucn_filename
  )
)[, c("species_birdtree", "iucn_cat")] %>% 
  rename(
    species = species_birdtree
  ) %>% 
  filter(
    !(iucn_cat %in% c("DD", "EW", "EX"))
  )

# TEST
# load old iucn patch data
iucn_filename <- paste0("iucn_2024_", iucn_type, "_patch231030.csv")
iucn_old <- read.csv(
  here::here(
    "4_SharedInputData", iucn_filename
  )
)[, c("species_birdtree", "iucn_cat")] %>% 
  rename(
    species = species_birdtree
  ) %>% 
  filter(
    !(iucn_cat %in% c("DD", "EW", "EX"))
  )

# get species in new but not old iucn (and patch) data
new_species <- iucn %>% 
  filter(
    !(species %in% iucn_old$species)
  )

# get matchedsex species
matchedsex <- male_species %>% 
  inner_join(
    female_species
  )

# check if any of the new species are matched sex
new_species %>% 
  filter(
    species %in% matchedsex$species
  )

# join to matched sex species
matched_iucn_old <- matchedsex %>% 
  left_join(
    iucn_old, by = "species"
  )
# remove 104 species for which we don't have an IUCN category (these are a mixture of species which
# are EX, EW, DD and species which we didn't yet have the data in the old patch data)
matched_iucn_old <- matched_iucn_old %>% 
  filter(!is.na(iucn_cat))

# END TEST

# join to matched sex species
matched_iucn <- male_species %>% 
  inner_join(
    female_species
  ) %>% 
  left_join(
    iucn, by = "species"
  )

# remove 17 species for which we don't have an IUCN category (I've manually checked these on 2025-07-18
# and they're all DD, EX, or EW. N.B. for Zosterops_conspicillatus only the nominate subspecies is extinct)
matched_iucn <- matched_iucn %>% 
  filter(!is.na(iucn_cat))


# get taxonomy with iucn categories (filtered to Neoaves)
taxo_iucn <- taxo %>% 
  filter(
    IOCOrder != "STRUTHIONIFORMES", 
    IOCOrder != "RHEIFORMES", 
    IOCOrder != "CASUARIIFORMES", 
    IOCOrder != "APTERYGIFORMES", 
    IOCOrder != "TINAMIFORMES", 
    IOCOrder != "GALLIFORMES", 
    IOCOrder != "ANSERIFORMES"
  ) %>% 
  left_join(
    iucn, by = join_by("TipLabel" == "species")
  )

# set up contingency table for chisq test
cont_table <- data.frame(
  n_data = rep(NA, 5),
  n_missing = rep(NA, 5)
)
rownames(cont_table) <- c("CR", "EN", "VU", "NT", "LC")

# filter to each IUCN category and print percentage coverage
for(category in c("CR", "EN", "VU", "NT", "LC")){
  
  n_taxo <- taxo_iucn %>% 
    filter(
      iucn_cat == category
    ) %>% 
    nrow()
  
  print(n_taxo)
  
  n_px <- matched_iucn %>% 
    filter(
      iucn_cat == category
    ) %>% 
    nrow()
  
  print(n_px)
  
  pc_cov <- (n_px / n_taxo) * 100
  
  cont_table[category, ] <- c(n_px, (n_taxo - n_px))
  
  print(paste0("Species coverage = ", round(pc_cov, digits = 1), "% for IUCN category ", category))
  
}

# run chi square test to check if % coverage is significantly different among the
# IUCN categories
chisq.test(cont_table)

# pairwise comparisons (with Bonferroni adjustment)
pairwise.prop.test(cont_table[,1], rowSums(cont_table), p.adjust.method = "bonferroni")


# Examine which species are present in the new (2025-07-16) patch dataset (which includes species from
# the Extinct & Endangered cabinet in Tring) but not in the old (2023-10-30) dataset

# load old data, bind  to taxonomy and remove galloanseriformes and palaeognaths
px_old <- readRDS("./2_Patches/1_InputData/patches.231030.rds") %>% 
  left_join(
    taxo,
    by = join_by(species == TipLabel)
  ) %>% 
  filter(
    IOCOrder != "STRUTHIONIFORMES", 
    IOCOrder != "RHEIFORMES", 
    IOCOrder != "CASUARIIFORMES", 
    IOCOrder != "APTERYGIFORMES", 
    IOCOrder != "TINAMIFORMES", 
    IOCOrder != "GALLIFORMES", 
    IOCOrder != "ANSERIFORMES"
  )

# load new data, bind  to taxonomy and remove galloanseriformes and palaeognaths
px_new <- readRDS("./2_Patches/1_InputData/patches.250716.rds") %>% 
  left_join(
    taxo,
    by = join_by(species == TipLabel)
  ) %>% 
  filter(
    IOCOrder != "STRUTHIONIFORMES", 
    IOCOrder != "RHEIFORMES", 
    IOCOrder != "CASUARIIFORMES", 
    IOCOrder != "APTERYGIFORMES", 
    IOCOrder != "TINAMIFORMES", 
    IOCOrder != "GALLIFORMES", 
    IOCOrder != "ANSERIFORMES"
  )

# count species for which we have male and female specimens - OLD patch data
# get species for which we have male specimens
male_species_old <- px_old %>% 
  filter(
    sex == "M"
  ) %>% 
  select(
    species
  ) %>% 
  unique()
# get species for which we have female specimens
female_species_old <- px_old %>% 
  filter(
    sex == "F"
  ) %>% 
  select(
    species
  ) %>% 
  unique()
# check how many species have both male and female data
matched_species_old <- male_species_old %>% 
  inner_join(
    female_species_old
  )

# check how many passeriforme species have both male and female data
matched_passeriformes_old <- male_species_old %>% 
  left_join(
    taxo,
    by = join_by(species == TipLabel)
  ) %>% 
  filter(
    IOCOrder == "PASSERIFORMES"
  ) %>% 
  inner_join(
    female_species_old
  )

# count species for which we have male and female specimens - NEW patch data
# get species for which we have male specimens
male_species_new <- px_new %>% 
  filter(
    sex == "M"
  ) %>% 
  select(
    species
  ) %>% 
  unique()
# get species for which we have female specimens
female_species_new <- px_new %>% 
  filter(
    sex == "F"
  ) %>% 
  select(
    species
  ) %>% 
  unique()
# check how many species have both male and female data
matched_species_new <- male_species_new %>% 
  inner_join(
    female_species_new
  )

# check how many passeriforme species have both male and female data
matched_passeriformes_new <- male_species_new %>% 
  left_join(
    taxo,
    by = join_by(species == TipLabel)
  ) %>% 
  filter(
    IOCOrder == "PASSERIFORMES"
  ) %>% 
  inner_join(
    female_species_new
  )

# isolate the matched sex species which are in the old but not the new data
new_species_neoaves <- matched_species_new %>% 
  filter(
    !(species %in% matched_species_old$species)
  )

# get their IUCN status
new_species_neoaves <- new_species_neoaves %>% 
  left_join(iucn, by = "species")
# WILL HAVE TO REMAKE THE IUCN DATA - it doesn't include any of these new species