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
# check how many species are different among the male and female species
male_species %>% 
  inner_join(
    female_species
  ) %>% 
  nrow()
