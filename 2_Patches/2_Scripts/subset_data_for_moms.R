## Create subset of species/sex pairs to upload to moms
## https://tguillerme.shinyapps.io/moms/
## To understand how the different metrics change when species are removed randomly or non-randomly from the space

# load libraries
library(dplyr)


# get PCA space from RDS
pca_jndxyzlum <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA", 
    "Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603.rds"
  )
) %>% 
  magrittr::extract2("x") %>% 
  as.data.frame()


# add species/sex as a column at the start
pca_jndxyzlum <- rownames(pca_jndxyzlum) %>% 
  bind_cols(pca_jndxyzlum)

colnames(pca_jndxyzlum) <- c("specimen", colnames(pca_jndxyzlum[, 2:41]))

# take a random sample of 5000 specimens (to make csv small enough to import to moms)
pca_5000 <- pca_jndxyzlum[-sample(nrow(pca_jndxyzlum), size = 10798), ]

# or actually, a better way to do it might be to convert to int before adding the extra rownames column



# write as csv to import to moms
readr::write_csv(
  pca_5000,
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA", 
    "Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603_5000_species_sample_for_moms.csv"
  ),
  col_names = T
)

