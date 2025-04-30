# Probe PCA space to find best disparity metric
# Largely based on code from T Guillerme's shinyapp
# https://tguillerme.shinyapps.io/moms/


# Clear environment
rm(list=ls())


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


## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Passeriformes"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
# N.B. will need to change the date in the pca_all filename if using usmldbl or usmldblr
# (from 240603 to 240806)
space <- "lab"
# Select number of PC axes to retain ("all" or a number)
axes <- "all"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"
# select sex of interest ("all", "male_female", "male_only", "female_only", "unknown_only")
sex_interest <- "male_female"



# Load data ----

# Load PCA data
pca_filename <- paste(clade, sex_match, "patches.231030.PCAcolspaces", "rds", sep = ".")
pca_dat <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
)[[space]][["x"]]


# Data preparation ----

# restrict PCA to requested number of axes
if(axes != "all"){
  pca_dat <- pca_dat[, 1:axes]
}


# generate individual list of sex-specific species values of metric of interest (keep or discard sexes
# depending on value of sex_match)
# first set sex list to iterate over
sexes <- set_sex_list(sex_interest)
# now apply a function to extract named vector lists of diversity values for each sex individually
sexed_metric_list <- lapply(sexes, extract_sex_vals, div_data = pca_dat, assemblage_level = TRUE)

# name the elements of the list
names(sexed_metric_list) <- sexes

# create function to sample a single PCA space
sample_space <- function(space, sample_size, seed){
  
  # set seed for reproducibility
  set.seed(seed)
  
  # get rownames of random species to retain
  to_retain <- sample(rownames(space), sample_size)
  
  # subset space by rownames and return
  return(space[to_retain, , drop = FALSE])
  
}

# Randomly remove species to limit to 1000 observations
spaces_to_save <- lapply(sexed_metric_list, sample_space, sample_size = 1000, seed = 10)

# save as CSV to upload to moms shinyapp
for(sex in sexes){
  write.csv(spaces_to_save[[sex]], file = paste0(sex, ".", space, ".randomsubset.csv"))
}
lapply(spaces_to_save, function(space) write.csv(space, file = ""))
                         