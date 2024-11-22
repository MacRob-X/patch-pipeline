# Fix taxonomy of IUCN data so it matches our Jetz-based taxonomy
# Robert MacDonald
# 30/09/2024

# Code to run each time -----------------------

# load libraries
library(dplyr)

# clear environment
rm(list=ls())

# load IUCN Red List data, subset to class Aves, remove galloanseriformes and palaeognaths
# add column with scientific name separated
# by underscore
iucn_master <- readr::read_rds(
  here::here(
    "4_SharedInputData", "IUCN_RedList_data_130324.rds"
  )
) 

# filter to only Neoaves
iucn_data <- iucn_master %>% 
  filter(
    class_name == "AVES",
    order_name != "STRUTHIONIFORMES", 
 #   IOCOrder != "RHEIFORMES", 
#    IOCOrder != "CASUARIIFORMES", 
#    IOCOrder != "APTERYGIFORMES", 
#    IOCOrder != "TINAMIFORMES", 
    order_name != "GALLIFORMES", 
    order_name != "ANSERIFORMES"
  ) %>% 
  mutate(
    species = sub(" ", "_", scientific_name)
  )

# load raw patch pixel values,
patch_master <- readr::read_rds(
  here::here(
    "2_Patches", "1_InputData", "patches.231030.rds"
  )
)

# load taxonomic data and extract species names
taxo_master <- read.csv(
  here::here(
    "4_SharedInputData", "BLIOCPhyloMasterTax_2019_10_28.csv"
  )
)



# remove gallanseriformes from patch pixel values
patch_data <- patch_master %>% 
  left_join(
    taxo_master,
    by = join_by(species == TipLabel)
  )  %>% 
  filter(
    IOCOrder != "GALLIFORMES"
  ) %>%
  select(
    species, specimen, sex, view, region, coord.x, coord.y, vR, vG, vB, uR, uG, uB, min.r2
  )

# extract patch species names
patch_names <- unique(patch_data$species)


# Using AVONET taxonomy matcher ----
# Downloaded from AVONET supplementary data repository
# https://figshare.com/s/b990722d72a26b5bfead

# load AVONET BirdLife-BirdTree crosswalk
# takes a little while because of the name cleaning
avonet_crosswalk <- read.csv(
  here::here(
    "4_SharedInputData", "avonet", "avonet_v7_birdlife-birdtree-crosswalk.csv"
  )
) %>% 
  rename(
    species_birdlife = Species1,
    species_birdtree = Species3,
    match_type = Match.type,
    match_notes = Match.notes
  ) %>% 
  mutate(
    species_birdlife = sub(" ", "_", species_birdlife),
    species_birdtree = sub(" ", "_", species_birdtree),
    match_type = snakecase::to_snake_case(match_type)
   ) %>% 
  filter( # remove blank rows
    !if_all(everything(), function(x) x == "")
  )

# we can filter out newly described species as these are all species that don't exist in the 
# birdtree taxonomy
# # There are 4 species in the BT taxonomy that are marked as "invalid taxon" and don't have 
# a BL taxonomy name. These are Anthus_longicaudatus, Hypositta_perdita, Lophura_hatinhensis, 
# and Phyllastrephus_leucolepis. None of these species are in the patch data, so we can remove them
# There are also 143 extinct species which have no BT name. Keep these in, as we may get data on 
# extinct species in the future
avonet_crosswalk <- avonet_crosswalk %>% 
  filter(
    match_type != "newly_described_species",
    match_type != "invalid_taxon"
  )

# Check number of species in each mapping type group (1 BL to 1 BT, 1 BL to many BT, many BL to 1 BT,
# extinct)
avonet_crosswalk %>% 
  group_by(match_type) %>% 
  summarise(
    n = n()
  )

## OK - here is where the actual taxonomy matching starts ----

# get full Aves IUCN data
iucn_aves <- iucn_master %>% 
  filter(
    class_name == "AVES",
  ) %>% 
  mutate(
    species = sub(" ", "_", scientific_name),
  ) %>% 
  rename(
    species_birdlife = species
  )

# first need to check that all species for which we have data appear in the AVONET BirdTree taxonomy
# identify any species which don't appear
patch_names[which(!(patch_names %in% avonet_crosswalk$species_birdtree))]
# Campylopterus_curvipennis doesn't appear. This is because it has since been renamed
# "Pampa_curvipennis" (as per IUCN Red List website, accessed 04/09/2024)


# Let's change the BirdTree name from "Pampa_curvipennis" to "Campylopterus_curvipennis"
avonet_crosswalk_checked <- avonet_crosswalk %>% 
  mutate(
    species_birdtree = if_else(species_birdtree == "Pampa_curvipennis", "Campylopterus_curvipennis", species_birdtree)
  )


# Check for species which are in my IUCN data (2023 data) but not in the AVONET crosswalk (
# which uses the HBW-BirdLife Version 5.0 (December 2020) checklist, accessed 18/06/2021)
# this gives the current IUCN name of the species
missing_birdlife_species <- iucn_aves$species_birdlife[which(!(iucn_aves$species_birdlife %in% avonet_crosswalk_checked$species_birdlife))] %>% 
  sort()

# to rectify the missing BirdLife species, I can compare version 5.0 of the HBW-BirdLife digital checklist
# (used by AVONET) with version 8.1 (which I believe is the version used by my IUCN data)
# Just need to replace the IUCN 2023 species names (v8.1) with their v5.0 names, where possible


# use a csv version of the 8.1 checklist and a csv version of the
# 5.0 checklist loaded into R
# load v8.1 and discard rows referring to subspecies
bl_81 <- readr::read_csv(
  here::here(
    "4_SharedInputData", "birdlife_taxonomic_checklists", "HBW-BirdLife_Checklist_v81_Jan24",
    "simplified_rxm_digital_checklist_v81.csv"
  ), 
  skip_empty_rows = TRUE
) %>% 
  filter(
    is.na(sub_spp_id)#,
    # iucn_cat_2023 != "NR"
  ) %>% 
  select(
    -c(subspp_seq, authority, alt_common_names, sub_spp_id)
  ) %>% 
  mutate(
    species_81 = sub(" ", "_", scientific_name)
  )

# load v5.0 and filter out not recognised complexes
bl_5 <- readr::read_csv(
  here::here(
    "4_SharedInputData", "birdlife_taxonomic_checklists", "HBW-BirdLife_Checklist_v5_Dec20",
    "simplified_rxm_digital_checklist_v5.csv"
  ),
  skip_empty_rows = TRUE
) %>% 
  mutate(
    species_5 = sub(" ", "_", scientific_name)
  ) %>% 
  filter(
    iucn_cat_5 != "NR",
    iucn_cat_5 != "UR"
  ) %>% 
  select(
    species_5, sis_rec_id, iucn_cat_5
  )

# check if there are any species which are in the IUCN 2023 data but not the BL v8.1 checklist
iucn_aves$species_birdlife[which(!(iucn_aves$species_birdlife %in% bl_81[bl_81$bl_tax_treat != "NR", ]$species_81))]
# nope, all good - the versions match

# check if there are any species which are in the avonet crosswalk BL list but not the BL v5.0 checklist
avonet_crosswalk_checked$species_birdlife[which(!(avonet_crosswalk_checked$species_birdlife %in% bl_5$species_5))]
# nope, all good - the versions match

# check if there are any species which are in the v5 checklist but only in the v8.1 checklist as
# NR species AND are 1_bl_to_1_bt mapping AND do not appear in the IUCN 2023 data AND are in our 
# patch data - these will not have an 
# IUCN category to assign
# Note that most of these will NOT cause a problem since this list will include species 
# which do have an IUCN category but whose scientific name has changed (e.g. Buettikoferella_bivittata/
# Cincloramphus_bivittatus)
# I could use the 2020 category I guess?
problem_species <- avonet_crosswalk_checked[(avonet_crosswalk_checked$species_birdlife %in% bl_5$species_5[!(bl_5$sis_rec_id %in% bl_81$sis_rec_id[bl_81$bl_tax_treat != "NR"]) & !(bl_5$species_5 %in% iucn_aves$species_birdlife)]) & avonet_crosswalk_checked$match_type == "1_bl_to_1_bt" & avonet_crosswalk_checked$species_birdtree %in% patch_names, ]

# join the birdlife v8.1 checklist to the v5.0 checklist (using SISRecID) and 
# filter to only the species which are missing from the AVONET crosswalk data

bl_checklist_match <- bl_5 %>% 
  full_join(
    bl_81, by = "sis_rec_id"
  ) %>% 
  # filter(
  #   species_81 %in% missing_birdlife_species | species_50 %in% missing_birdlife_species
  # ) %>% 
  select(
    -c(seq, order, family_name, family, subfamily, tribe, common_name, scientific_name, synonyms, taxonomic_sources, spc_rec_id, bl_tax_treat)
  )

# check which 5.0 species don't have an 8.1 name equivalent
bl_checklist_match %>% 
  filter(
    is.na(species_81)
  )

# All 5.0 species have an 8.1 equivalent
 
# identify rows which have duplicated v8.1 species names
dupes_sp81 <- bl_checklist_match[bl_checklist_match$species_81 %in% unique(bl_checklist_match$species_81[duplicated(bl_checklist_match$species_81)]), ]

# remove duplicate bl_81 species which have NA bl_5
duplicated_spp <- unique(bl_checklist_match$species_81[duplicated(bl_checklist_match$species_81)])

for(sp in duplicated_spp){
  dupe_rows <- bl_checklist_match %>% 
    filter(
      species_81 == sp
    )
  # remove row if there's a row with corresponding BL v5 name AND rows with no corresponding BL v5 name
  if(nrow(dupe_rows[!is.na(dupe_rows$species_5), ]) != 0){
    bl_checklist_match <- bl_checklist_match %>% 
      filter(
        !(species_81 == sp & is.na(species_5))
      )
  }
}

# check if any duplicate rows left
dupes_sp81 <- bl_checklist_match[bl_checklist_match$species_81 %in% unique(bl_checklist_match$species_81[duplicated(bl_checklist_match$species_81)]), ]

# can remove the remaining duplicates if their v8.1 IUCN cat is NR and they have no v5 species name
# (as the NR species won't appear in my IUCN data)
duplicated_spp <- unique(bl_checklist_match$species_81[duplicated(bl_checklist_match$species_81)])

for(sp in duplicated_spp){
  dupe_rows <- bl_checklist_match %>% 
    filter(
      species_81 == sp
    )
  # remove rows if v8.1 IUCN cat is NR
  bl_checklist_match <- bl_checklist_match %>% 
    filter(
      !(species_81 == sp & iucn_cat_81 == "NR")
    )
  
}

# check if any duplicate rows left
dupes_sp81 <- bl_checklist_match[bl_checklist_match$species_81 %in% unique(bl_checklist_match$species_81[duplicated(bl_checklist_match$species_81)]), ]
# no duplicated rows left - all good. Let's remove the variables
rm(dupes_sp81, dupe_rows, duplicated_spp, sp)


# now use this data to replace the 2023 IUCN species names with the v5.0 species names
# (don't replace if the v5.0 species name is NA)
iucn_aves_v5 <- iucn_aves %>% 
  left_join(
    bl_checklist_match,
    by = join_by(species_birdlife == species_81)
  ) %>% 
  mutate(
    species_iucn_2023 = species_birdlife
  ) %>% 
  mutate(
    species_birdlife = ifelse(!is.na(species_5), species_5, species_birdlife)
  ) %>% 
  select(
    -species_5
  )
# "species_birdlife" column now has the v5.0 names - original 2023 IUCN names preserved in
# "species_iucn_2023" column



# Join the IUCN data to the avonet crosswalk and add 'final' IUCN category for final category
# which has been manually checked, if necessary. For the one-to-one mapped species, this is the 
# same as the 'initial' IUCN category
avonet_iucn <- iucn_aves_v5 %>% 
  select(
    order_name, family_name, genus_name, species_birdlife, category
  ) %>% 
  full_join(
    avonet_crosswalk_checked,
    by = "species_birdlife"
  ) %>% 
  rename(
    initial_category = category
  ) %>% 
  mutate(
    final_category = NA
  )

# every BirdTree species should now have a corresponding BirdLife species


# filter to only the species for which we have patch data by comparing the patch data species
# with the birdtree species in the updated avonet crosswalk
# When it comes to adding in the data from the second round of image collection (August-October 2024),
# we can just filter this by the new data species so we avoid rechecking the whole thing
avonet_iucn <- avonet_iucn %>% 
  filter(
    species_birdtree %in% patch_names
  )

# Check for rows which are NA for initial category (these are rows for which we potentially
# don't have IUCN data)
na_init_iucn <- avonet_iucn %>% 
  filter(
    is.na(initial_category)
  )

# the important ones of these are the 1_bl_to_1_bt ones, as the others probably have another
# BL species which maps to the BT species
na_init_iucn <- na_init_iucn %>% 
  filter(
    match_type == "1_bl_to_1_bt"
  )

# from manual checking on IUCN Red List and Avibase, these are all species
# which have been lumped by BirdLife between v5.0 Checklist (2020) and v8.1 Checklist (2023)
# now these species are recognised as subspecies of other species
# DECISION: for these species, I can use either the 2020 IUCN category of the species OR
# the 2023 IUCN category of the lumped species
# ---> Assign the 2023 full species category (only differs from the 2020 split species category
## in one of these species in any case)

## Taxonomic sources: Avibase, IUCN Red List
## All of the below 2020 species are now considered subspecies of the 2023 species by IUCN
## Effectively this means that these species should be mapped "1_bl_to_many_bt" (because they're
## split in BT but lumped in BL)
avonet_iucn[avonet_iucn$species_birdtree %in% na_init_iucn$species_birdtree, "match_type"] <- "1_bl_to_many_bt"
## Also need to assign the initial category for these species as the initial category of the
## equivalent 2023 species
## ## Since these species are all in the same genus as their 2023 equivalents, we can also
## assign them the same order, family, and genus
## 
## Cacicus_leucoramphus
## 2020 category: LC
## 2023 category (C. chrysonotus): LC
## MATCH
avonet_iucn[avonet_iucn$species_birdlife == "Cacicus_leucoramphus", "initial_category"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Cacicus_chrysonotus", "initial_category"]
avonet_iucn[avonet_iucn$species_birdlife == "Cacicus_leucoramphus", "order_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Cacicus_chrysonotus", "order_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Cacicus_leucoramphus", "family_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Cacicus_chrysonotus", "family_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Cacicus_leucoramphus", "genus_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Cacicus_chrysonotus", "genus_name"]
## Calendulauda_alopex
## 2020 category: LC
## 2023 category (C. africanoides): LC
## MATCH
avonet_iucn[avonet_iucn$species_birdlife == "Calendulauda_alopex", "initial_category"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Calendulauda_africanoides", "initial_category"]
avonet_iucn[avonet_iucn$species_birdlife == "Calendulauda_alopex", "order_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Calendulauda_africanoides", "order_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Calendulauda_alopex", "family_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Calendulauda_africanoides", "family_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Calendulauda_alopex", "genus_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Calendulauda_africanoides", "genus_name"]
## Caracara_cheriway
## 2020 category: LC
## 2023 category (C. plancus): LC
## MATCH
avonet_iucn[avonet_iucn$species_birdlife == "Caracara_cheriway", "initial_category"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Caracara_plancus", "initial_category"]
avonet_iucn[avonet_iucn$species_birdlife == "Caracara_cheriway", "order_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Caracara_plancus", "order_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Caracara_cheriway", "family_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Caracara_plancus", "family_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Caracara_cheriway", "genus_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Caracara_plancus", "genus_name"]
## Cranioleuca_baroni
## 2020 category: LC
## 2023 category (C. antisiensis):  LC
## MATCH
avonet_iucn[avonet_iucn$species_birdlife == "Cranioleuca_baroni", "initial_category"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Cranioleuca_antisiensis", "initial_category"]
avonet_iucn[avonet_iucn$species_birdlife == "Cranioleuca_baroni", "order_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Cranioleuca_antisiensis", "order_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Cranioleuca_baroni", "family_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Cranioleuca_antisiensis", "family_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Cranioleuca_baroni", "genus_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Cranioleuca_antisiensis", "genus_name"]
## Glaucidium_californicum
## 2020 category: LC
## 2023 category (G. gnoma): LC
## MATCH
avonet_iucn[avonet_iucn$species_birdlife == "Glaucidium_californicum", "initial_category"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Glaucidium_gnoma", "initial_category"]
avonet_iucn[avonet_iucn$species_birdlife == "Glaucidium_californicum", "order_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Glaucidium_gnoma", "order_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Glaucidium_californicum", "family_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Glaucidium_gnoma", "family_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Glaucidium_californicum", "genus_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Glaucidium_gnoma", "genus_name"]
## Mirafra_ashi
## 2020 category: EN
## 2023 category (M. somalica): LC
## NOT MATCH
avonet_iucn[avonet_iucn$species_birdlife == "Mirafra_ashi", "initial_category"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Mirafra_somalica", "initial_category"]
avonet_iucn[avonet_iucn$species_birdlife == "Mirafra_ashi", "order_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Mirafra_somalica", "order_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Mirafra_ashi", "family_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Mirafra_somalica", "family_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Mirafra_ashi", "genus_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Mirafra_somalica", "genus_name"]
## Ramphocelus_costaricensis
## 2020 category: LC
## 2023 category (R. passerinii): LC
## MATCH
avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_costaricensis", "initial_category"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_passerinii", "initial_category"]
avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_costaricensis", "order_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_passerinii", "order_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_costaricensis", "family_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_passerinii", "family_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_costaricensis", "genus_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_passerinii", "genus_name"]
## Ramphocelus_icteronotus
## 2020 category: LC
## 2023 category (R. flammigerus): LC
## MATCH
avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_icteronotus", "initial_category"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_flammigerus", "initial_category"]
avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_icteronotus", "order_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_flammigerus", "order_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_icteronotus", "family_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_flammigerus", "family_name"]
avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_icteronotus", "genus_name"] <-
  avonet_iucn[avonet_iucn$species_birdlife == "Ramphocelus_flammigerus", "genus_name"]



# make final category the same as initial category for the one-to-one mapped species
avonet_iucn <- avonet_iucn %>% 
  mutate(
    final_category = if_else(match_type == "1_bl_to_1_bt", initial_category, final_category)
  )

# we can do the same for species which are mapped as one BirdLife species to many BirdTree species
avonet_iucn <- avonet_iucn %>% 
  mutate(
    final_category = if_else(match_type == "1_bl_to_many_bt", initial_category, final_category)
  )

# we can do the same for extinct species. Note that there are several species which the AVONET
# crosswalk classifies as extinct but are not classified as extinct by IUCN, so we need to assign the 
# initial_category as the final_category, rather than just assigning "EX"
avonet_iucn <- avonet_iucn %>% 
  mutate(
    final_category = if_else(match_type == "extinct", initial_category, final_category)
  )

# check for any duplicate BT species that are not assigned as "many_bl_to_1_bt"
dupes <- avonet_iucn$species_birdtree[duplicated(avonet_iucn$species_birdtree) & avonet_iucn$match_type != "many_bl_to_1_bt"]
avonet_iucn %>% 
  filter(
    species_birdtree %in% dupes
  ) %>% 
  head()

# these are all taxa with imperfect matches where a single BirdTree species is assigned to 
# multiple BirdLife species
# Don't know why these aren't assigned as "many_bl_to_1_bt" - assign as that 
avonet_iucn <- avonet_iucn %>% 
  mutate(
    match_type = if_else(species_birdtree %in% dupes, "many_bl_to_1_bt", match_type)
  )


# filter to just species with many BL species corresponding to 1 BT species
many_bl_1_bt <- avonet_iucn %>% 
  filter(
    match_type == "many_bl_to_1_bt"
  ) 

# for these species, we can automatically assign their final IUCN category as the same as their initial
# IUCN category IF the initial IUCN category is the same for all split BL species mapped to each BT
# species
# loop to do this
bt_spec <- unique(many_bl_1_bt$species_birdtree)
for(spec in bt_spec){
  # check if all BL species have the same IUCN category - assign as final category if so
  iucn_cats <- unique(many_bl_1_bt$initial_category[many_bl_1_bt$species_birdtree == spec])
  if(
    length(iucn_cats) == 1
  ){
    avonet_iucn <- avonet_iucn %>% 
      mutate(
        final_category = if_else(species_birdtree == spec, initial_category, final_category)
      )
  }
}

# check how many "many_bl_to_1_bt" species still have no final IUCN category
avonet_iucn %>% 
  filter(
    match_type == "many_bl_to_1_bt",
    is.na(final_category)
  ) %>% 
  count()
# 636 BL species left - these are species for which the species which are split in the BL taxonomy
# (but are a single species in the BT taxonomy) have different IUCN categories
# These species will need to be assigned a category somehow - conservative (assign lowest cat)?
# liberal (assign highest cat)? or random
# perhaps a good idea to assign the least threatened category to these, as the splitting process
# will probably 'artificially' increase the threat category due to smaller range size

#-------# 
# DEAL WITH THE REMAINING many_bl_to_1_bt" SPECIES

# CONSERVATIVE APPROACH - assign the least threatened category for each of these many_bl_to_1_bt species

# set IUCN cat as factor so we can order it - increasing order of threat level
avonet_iucn <- avonet_iucn %>% 
  mutate(
    initial_category = factor(
      initial_category,
      levels = c("LC", "NT", "VU", "EN", "CR", "EW", "EX", "DD")
    )
  )

# group by BT species, order by initial iucn category, then remove duplicates 
# (keeping first of each group)
to_keep <- avonet_iucn %>% 
  filter(
    match_type == "many_bl_to_1_bt",
    is.na(final_category)
  ) %>% 
  group_by(species_birdtree) %>% 
  arrange(initial_category, .by_group = TRUE) %>% 
  distinct(species_birdtree, .keep_all = TRUE)

avonet_iucn_conservative <- avonet_iucn %>% 
  filter(
    !(match_type == "many_bl_to_1_bt" & is.na(final_category))
  ) %>% 
  bind_rows(to_keep)

# set final category as initial category for these species
avonet_iucn_conservative <- avonet_iucn_conservative %>% 
  mutate(
    final_category = if_else(species_birdtree %in% to_keep$species_birdtree, initial_category, final_category)
  )

# LIBERAL APPROACH - assign the most threatened category for each of these many_bl_to_1_bt species

# set IUCN cat as factor so we can order it - decreasing order of threat level
avonet_iucn <- avonet_iucn %>% 
  mutate(
    initial_category = factor(
      initial_category,
      levels = c("EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD")
    )
  )

# group by BT species, order by initial iucn category, then remove duplicates 
# (keeping first of each group)
to_keep <- avonet_iucn %>% 
  filter(
    match_type == "many_bl_to_1_bt",
    is.na(final_category)
  ) %>% 
  group_by(species_birdtree) %>% 
  arrange(initial_category, .by_group = TRUE) %>% 
  distinct(species_birdtree, .keep_all = TRUE)

avonet_iucn_liberal <- avonet_iucn %>% 
  filter(
    !(match_type == "many_bl_to_1_bt" & is.na(final_category))
  ) %>% 
  bind_rows(to_keep)

# set final category as initial category for these species
avonet_iucn_liberal <- avonet_iucn_liberal %>% 
  mutate(
    final_category = if_else(species_birdtree %in% to_keep$species_birdtree, initial_category, final_category)
  )

# NOMINATE SUBSPECIES APPROACH
# take the threat category of the nominate subspecies as that’s the one we’ve photographed in most cases 

# Filter to only the many_bl_to_1_bt species that haven't already had their final category assigned
# add columns with the BL and BT genus and species names in separate columns
# compare BL species name with BT species name within each row to check if they're the same
# the ones that match are almost certainly the nominate subspecies, so assign their match type as
# "nom_to_nom"
# just to make sure they are indeed the nominate subspecies, flag those rows in which the genus name is different 
# can also assign as "nom_to_nom" those species where the species name is the same but the suffix changes
# from an -us to an -a or vice versa - again, flag these species so they can be manually checked
# 
# flagged rows where the species name is the same but genus name is different (i.e. nom_to_nom, flag) e.g. 
# Ramphiculus_fischeri/Ptilinopus_fischeri are almost certainly synonyms of the same species/subspecies, 
# but flagging them allows me to manually check
# flagged rows where the species name and genus name is different are either
# - rows that can be removed as they are other subspecies of a species for which we have a nom_to_nom match
# - rows that should be retained because there is no easy-to-find nom_to_nom match
# the easiest way to find out if there are many of the latter type is to count the unique birdlife and birdtree
# species names - if they're the same, we already have a nom_to_nom match for all species
find_noms <- avonet_iucn %>% 
  filter(
    match_type == "many_bl_to_1_bt",
    is.na(final_category)
  ) %>% 
  mutate(
    bl_spec = sapply(strsplit(species_birdlife, split = "_"), "[", 2),
    bt_spec = sapply(strsplit(species_birdtree, split = "_"), "[", 2),
    bl_genus = sapply(strsplit(species_birdlife, split = "_"), "[", 1),
    bt_genus = sapply(strsplit(species_birdtree, split = "_"), "[", 1)
  ) %>% 
  mutate(
    match_type = if_else(bl_spec == bt_spec, "nom_to_nom", match_type),
    match_type = if_else(sub("us([^us]*)$", "a\\1", bl_spec) == bt_spec | sub("a([^a]*)$", "us\\1", bl_spec) == bt_spec, "nom_to_nom", match_type)
  ) %>% 
  mutate(
    flag = if_else(bl_genus != bt_genus & match_type == "nom_to_nom", "flag", "not_flag")
  )

# now we just need to check if there are any BT species which don't have a nom_to_nom match type
nom_found_spp <- find_noms %>% 
  filter(
    match_type == "nom_to_nom"
  ) %>% 
  pull(species_birdtree) %>% 
  unique()

find_noms %>% 
  mutate(
    nom_found = if_else(species_birdtree %in% nom_found_spp, "y", "n")
  ) %>%
  filter(
    nom_found == "n"
  ) %>% 
  pull(species_birdtree) %>% 
  unique()

# there are three BirdTree species which we haven't automatically been able to match to 
# a nominate BirdLife species: "Coracina_tenuirostris", "Nectarinia_afra", "Zosterops_palpebrosus"
# I can just manually check these using IUCN Red List website and birdsoftheworld.org (accessed 26/09/2024)
# and then manually assign the nominate subspecies
# Coracina_tenuirostris maps to Edolisoma_tenuirostre
find_noms[find_noms$species_birdlife == "Edolisoma_tenuirostre", "match_type"] <- "nom_to_nom"
# Nectarinia_afra maps to Cinnyris_afer
find_noms[find_noms$species_birdlife == "Cinnyris_afer", "match_type"] <- "nom_to_nom"
# Zosterops_palpebrosus actually already maps to BL species Zosterops_palpebrosus, but
# for some reason the nominate subspecies is marked in the crosswalk as "1_bl_to_many_bt"
# instead of "many_bl_to_1_bt". This means we can just ignore this one, it'll get filtered out in the next step

# now I can cross-reference the avonet_iucn df with the find_noms df (which has all the info regarding
# nominate subspecies)
# first remove all species marked as "many_bl_to_1_bt" which don't already have a final category assigned
# then bind to the find noms data, but filtering out any species which are not the nominate subspecies
avonet_iucn_nominate <- avonet_iucn %>% 
  filter(
    !(match_type == "many_bl_to_1_bt" & is.na(final_category))
  ) %>% 
  bind_rows(
    filter(find_noms, match_type == "nom_to_nom")
  )

# set final category as initial category for these species
avonet_iucn_nominate <- avonet_iucn_nominate %>% 
  mutate(
    final_category = if_else(match_type == "nom_to_nom", initial_category, final_category)
  )


#-------#


# Now we can remove any duplicate BirdTree species, as these should all be "many_bl_to_1_bt" species
# which have had their final category assigned
# just double check that there are no stray duplicates with other match_types
avonet_iucn_conservative$species_birdtree[duplicated(avonet_iucn_conservative$species_birdtree) & avonet_iucn_conservative$match_type != "many_bl_to_1_bt"]
# now remove the duplicate rows
avonet_iucn_conservative <- avonet_iucn_conservative %>% 
  distinct(
    species_birdtree, .keep_all = TRUE
  )
# same for the liberal dataset
avonet_iucn_liberal$species_birdtree[duplicated(avonet_iucn_liberal$species_birdtree) & avonet_iucn_liberal$match_type != "many_bl_to_1_bt"]
# now remove the duplicate rows
avonet_iucn_liberal <- avonet_iucn_liberal %>% 
  distinct(
    species_birdtree, .keep_all = TRUE
  )
# same for nominate dataset
avonet_iucn_nominate$species_birdtree[duplicated(avonet_iucn_nominate$species_birdtree) & avonet_iucn_nominate$match_type != "many_bl_to_1_bt"]
# now remove the duplicate rows
avonet_iucn_nominate <- avonet_iucn_nominate %>% 
  distinct(
    species_birdtree, .keep_all = TRUE
  )

# we should now be DONE with the taxonomy matching ----


# try joining to the patch data, for example
patch_dat <- as.data.frame(patch_names)
patch_dat <- patch_dat %>% 
  left_join(
    avonet_iucn_conservative,
    by = join_by(patch_names == species_birdtree)
  )

# DONE - let's select the relevant columns and save as a csv to be used as IUCN data
iucn_2024_conservative <- avonet_iucn_conservative %>% 
  select(
    order_name, family_name, genus_name, species_birdtree, final_category
  ) %>% 
  rename(
    iucn_cat = final_category
  )

iucn_2024_conservative %>%
  readr::write_csv(
    file = here::here(
    "4_SharedInputData", "iucn_2024_conservative.csv"
  ),
    col_names = TRUE
  )

iucn_2024_liberal <- avonet_iucn_liberal %>% 
  select(
    order_name, family_name, genus_name, species_birdtree, final_category
  ) %>% 
  rename(
    iucn_cat = final_category
  )

iucn_2024_liberal %>%
  readr::write_csv(
    file = here::here(
      "4_SharedInputData", "iucn_2024_liberal.csv"
    ),
    col_names = TRUE
  )

iucn_2024_nominate <- avonet_iucn_nominate %>% 
  select(
    order_name, family_name, genus_name, species_birdtree, final_category
  ) %>% 
  rename(
    iucn_cat = final_category
  )

iucn_2024_nominate %>%
  readr::write_csv(
    file = here::here(
      "4_SharedInputData", "iucn_2024_nominate.csv"
    ),
    col_names = TRUE
  )
