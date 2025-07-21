## Generate colour grids to show patch colours on colour space ----
## Robert MacDonald
## 13th May 2024

## Note that this only needs to be run once to generate the colour grids. The grids can
## then be used to plot using 02_eii_Patch_colour_grids_plot.R

# clear environment
rm(list=ls())

## Load libraries ---- 
library(dplyr)
library(ggplot2)

## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neoaves"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr", etc.)
space <- "lab"
# restrict to only species with male and female data? ("matchedsex" or "allspecimens")
sex_match <- "matchedsex"
# Overwrite previous colour grids in write folder?
overwrite_grids <- FALSE

# load rgb data ----
srgb_dat <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "1_RawColourspaces", 
    paste0(clade, ".", sex_match, ".patches.250716.prePCAcolspaces.rds")
  )
) %>% 
  magrittr::extract2("sRGB")

# load raw (pre-PCA) TCS xyz data
raw_patch_tcs <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "1_RawColourspaces", 
    paste0(clade, ".", sex_match, ".patches.250716.prePCAcolspaces.rds")
  )
) %>% 
  magrittr::extract2("xyz")

# load PCA data (output from 02_Patch_Analyse_features.R)
pca_all <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA", 
    paste0(clade, ".", sex_match, ".patches.250716.PCAcolspaces.rds")
  )
) %>% 
  magrittr::extract2(space)

# load jndxyzlumr umap (output from 02b_Patch_Umap_iterations.R)
umap <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
    paste0(clade, ".patches.", space, ".pca.canonUMAP.rds")
  )
) %>% 
  magrittr::extract2("layout") %>% 
  as.data.frame() %>% 
  rename(
    umap_1 = V1, umap_2 = V2
  )

# or perform UMAP on another colour model PCA dataset
# N.B. best not to do this - better to use the canon UMAP
# umap_xyzlumr <- umap::umap(
#   pca_all$x
# )   %>% 
#   magrittr::extract2("layout") %>% 
#   as.data.frame() %>% 
#   rename(
#     umap_1 = V1, umap_2 = V2
#   )


# get species and sex and tidy up names
srgb_dat <- srgb_dat %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2),
  )  %>% 
  rename(
    bel_r = sRGB.r.bel, bel_g = sRGB.g.bel, bel_b = sRGB.b.bel, 
    bre_r = sRGB.r.bre, bre_g = sRGB.g.bre, bre_b = sRGB.b.bre, 
    cov_r = sRGB.r.cov, cov_g = sRGB.g.cov, cov_b = sRGB.b.cov, 
    cro_r = sRGB.r.cro, cro_g = sRGB.g.cro, cro_b = sRGB.b.cro, 
    fli_r = sRGB.r.fli, fli_g = sRGB.g.fli, fli_b = sRGB.b.fli, 
    man_r = sRGB.r.man, man_g = sRGB.g.man, man_b = sRGB.b.man, 
    nap_r = sRGB.r.nap, nap_g = sRGB.g.nap, nap_b = sRGB.b.nap, 
    rum_r = sRGB.r.rum, rum_g = sRGB.g.rum, rum_b = sRGB.b.rum, 
    tai_r = sRGB.r.tai, tai_g = sRGB.g.tai, tai_b = sRGB.b.tai, 
    thr_r = sRGB.r.thr, thr_g = sRGB.g.thr, thr_b = sRGB.b.thr, 
  )

srgb_dat %>% glimpse()

# convert to tidy data
# get r values
srgb_r <- srgb_dat %>% 
  tidyr::pivot_longer(
    cols = c(bel_r, bre_r, cov_r, cro_r, fli_r, man_r, nap_r, rum_r, tai_r, thr_r),
    names_to = "body_part",
    values_to = "r"
  ) %>% 
  select(species, sex, body_part, r) %>% 
  mutate(
    body_part = stringr::str_replace(body_part, "_r", "")
  )
# get g values
srgb_g <- srgb_dat %>% 
  tidyr::pivot_longer(
    cols = c(bel_g, bre_g, cov_g, cro_g, fli_g, man_g, nap_g, rum_g, tai_g, thr_g),
    names_to = "body_part",
    values_to = "g"
  ) %>% 
  select(species, sex, body_part, g) %>% 
  mutate(
    body_part = stringr::str_replace(body_part, "_g", "")
  )
# get b values
srgb_b <- srgb_dat %>% 
  tidyr::pivot_longer(
    cols = c(bel_b, bre_b, cov_b, cro_b, fli_b, man_b, nap_b, rum_b, tai_b, thr_b),
    names_to = "body_part",
    values_to = "b"
  ) %>% 
  select(species, sex, body_part, b) %>% 
  mutate(
    body_part = stringr::str_replace(body_part, "_b", "")
  )
# combine into one longer dataframe
srgb_long <- srgb_r %>% 
  inner_join(
    srgb_g
  ) %>% 
  inner_join(
    srgb_b
  )

# convert body part to factor so can specify level order (for changing position of body parts in grid)
srgb_long <- srgb_long %>% 
  mutate(
    body_part = factor(body_part, levels = c("thr", "cro", "bre", "nap", "bel", "man", "cov", "rum", "fli", "tai"))
  )
# inspect
srgb_long %>% glimpse()
# remove intermediate dataframes
rm(srgb_r, srgb_g, srgb_b)


# Create colour grids for each specimen
for(i in 1:length(srgb_dat$species)){
  
  # get specimen
  specimen <- paste(srgb_dat[i, "species"], srgb_dat[i, "sex"], sep = "_")
  specimen_rgb <- srgb_long %>% 
    filter(
      species == srgb_dat[i, "species"],
      sex == srgb_dat[i, "sex"]
    ) %>% 
    arrange(body_part)  # this determines the position of each body part in the grid
  
  # if not overwriting previous colour grids, check if there's already a colour grid in the folder
  if(overwrite_grids == FALSE){
    
    # check if grid already exists
    if(file.exists(paste0("C:/Users/bop23rxm/Documents/colour_grids_repositioned/", specimen, ".png"))){
      # skip to next iteration if there's already a grid for this specimen
      next
    }
    
  }
  
  # initialise png saving
  png(
    # here::here(
    #   "2_Patches", "3_OutputData", "5_colour_grids", paste0(specimen, ".png")
    # ),
    paste0("C:/Users/bop23rxm/Documents/colour_grids_repositioned/", specimen, ".png"),   # use local storage as faster
    width = 100, height = 160,
    units = "px"
  )
  
  # Set up the plotting area
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(5, 2)))
  
  # Plot each color in a 5x2 grid without white space
  for (j in 1:10) {
    row <- ((j - 1) %/% 2) + 1
    col <- ((j - 1) %% 2) + 1
    grid::grid.rect(x = (col - 1) / 2, y = (5 - row) / 5,
                    width = 0.5, height = 0.2,
                    just = c("left", "bottom"),
                    gp = grid::gpar(fill = rgb(specimen_rgb$r[j], specimen_rgb$g[j], specimen_rgb$b[j]),
                                    col = NA))
  }
  
  dev.off()
  
  # display where up to 
  cat(paste0("\r", i, " of ", length(srgb_dat$species), " processed"))
}

# Use 02eii_Patch_colour_grids_plot.R to plot the colour grids on the colourspace