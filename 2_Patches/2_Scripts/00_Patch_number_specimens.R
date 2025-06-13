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
    IOCOrder != "GALLIFORMES"
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

# select IUCN type
iucn_type <- "nominate"

# load IUCN data
iucn_filename <- paste0("iucn_2024_", iucn_type, ".csv")
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

# join to matched sex species
matched_iucn <- male_species %>% 
  inner_join(
    female_species
  ) %>% 
  left_join(
    iucn, by = "species"
  )

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

# filter to each IUCN category and print percentage coverage
for(category in c("CR", "EN", "VU", "NT", "LC")){
  
  n_taxo <- taxo_iucn %>% 
    filter(
      iucn_cat == category
    ) %>% 
    nrow()
  
  n_px <- matched_iucn %>% 
    filter(
      iucn_cat == category
    ) %>% 
    nrow()
  
  pc_cov <- (n_px / n_taxo) * 100
  
  print(paste0("Species coverage = ", round(pc_cov, digits = 1), "% for IUCN category ", category))
  
}
