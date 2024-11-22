# Reconstruct the theoretical colour pattern phenotype at each point in colour pattern space
# Robert MacDonald
# 7th August 2024

# Load libraries ----
library(dplyr)

# Clear environment
rm(list=ls())

# EDITABLE CODE # ----
# Choose space to work with (usmldbl, usmldblr, xyz, xyzlum, xyzlumr, lab, cie, srgb, hex, 
# jndxyz, jndxyzlum, jndxyzlumr)
space <- "srgb"


# Load in PCA colourspace ----
pca_space <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    paste0("Neoaves.patches.231030.PCAcolspaces.", space, ".240603.rds")
  )
)


# Generate a matrix of every point in colour pattern space ---- 
# Note that this is not feasible for a high resolution with a high number of dimensions
# as the number of data points to generate quickly becomes > number of atoms in known universe
# as a general rule, let's say we can work with max 100 000 generated data points
# which equates to a resolution of 10 for 5 axes
# A resolution of 15 for 4 axes would give a reasonable 50,625 points

# choose number of axes to use

# first get max and min real values of each PC (gamut of actual bird colours)
pc_maxes <- apply(pca_space$x, 2, max)
pc_mins <- apply(pca_space$x, 2, min)

# Let's use the absolute max and absolute min to give our chosen resolution of n points
# across the range
# Let's start from 10% below the minimum range and go up to 10% above
res <- 15
grand_max <- max(pc_maxes[])
grand_min <- min(pc_mins)
new_max <- grand_max + (grand_max/10)
new_min <- grand_min - (grand_min/10)
increment <- (new_max - new_min) / (res-1)

# generate a matrix spanning from the grand min to the grand max, moving in our
# determined increments
n_axes <- 2
npoints <- res ^ n_axes
theor_pca <- matrix(NA, nrow = npoints, ncol = n_axes)


# let's just do it for 2 axes first, as this will simplify the creation of the matrix a lot
n_axes <- 2
# for each point with a given value in PC1, need res points in PC2, so total number of points
# is res ^ n_axes
npoints <- res ^ n_axes
theor_pca <- matrix(NA, nrow = npoints, ncol = n_axes)

# or each point with value new_min in axis 1, need res points in axis 2
j <- 1
for(i in seq(1, npoints, by = res)){
  theor_pca[i:(i+res-1), 1] <- new_min + (j-1)*increment
  theor_pca[i:(i+res-1), 2] <- seq(from = new_min, to = new_max, by = increment)
  j <- j + 1
}

# Need to find a way to generalise the above to n_axes dimensions



# Reconstruct original data from generated PCA data

# First cut off the pca space to match the number of axes selected
# and add our newly generated space data
pca_trunc <- list()
pca_trunc$sdev <- pca_space$sdev[1:n_axes]
pca_trunc$rotation <- pca_space$rotation[, 1:n_axes]
pca_trunc$center <- pca_space$center
pca_trunc$scale <- pca_space$scale
pca_trunc$x <- theor_pca

# convert to prcomp class
class(pca_trunc) <- "prcomp"

# and now back-engineer the space
recon_space <- pca_trunc$x %*% t(pca_trunc$rotation)

# add back center and scale, if applicable
if(!is.null(pca_trunc$center)) {
  recon_space <- sweep(recon_space, 2, pca_trunc$center, FUN = "+")
}
if(pca_trunc$scale == TRUE) {
  recon_space <- sweep(recon_space, 2, pca_trunc$scale, FUN = "*")
}

# inspect reconstructed data
head(recon_space)


# Generate colour pattern grids ----

# First convert USML+DBL, TCS etc to RGB format

# Assuming that data is now in RGB format
srgb_dat <- recon_space

# tidy up names and add rownames of coordinates
rownames(srgb_dat) <- paste(round(theor_pca[, 1], 2), round(theor_pca[, 2], 2), sep = ".")

srgb_dat <- srgb_dat %>% 
  as.data.frame() %>% 
  mutate(
    PC1 = round(theor_pca[, 1], 2),
    PC2 = round(theor_pca[, 2], 2),
  ) %>% 
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


# convert to tidy data
# get r values
srgb_r <- srgb_dat %>% 
  tidyr::pivot_longer(
    cols = c(bel_r, bre_r, cov_r, cro_r, fli_r, man_r, nap_r, rum_r, tai_r, thr_r),
    names_to = "body_part",
    values_to = "r"
  ) %>% 
  select(PC1, PC2, body_part, r) %>% 
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
  select(PC1, PC2, body_part, g) %>% 
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
  select(PC1, PC2, body_part, b) %>% 
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

# Convert any impossible colours (i.e. r, g, or b < 0 or > 1) to hot pink (0.96, 1, 1)
for(i in nrow(srgb_long)){
  if(srgb_long$r[i] < 0 | srgb_long$r[i] > 1 | srgb_long$g[i] < 0 | srgb_long$g[i] > 1 | srgb_long$b[i] < 0 | srgb_long$b[i] > 1) {
    srgb_long$r <- 0.961
    srgb_long$g <- 0.411
    srgb_long$b <- 0.706
  }
}

# Create colour grids for each specimen
for(i in 1:length(srgb_dat$PC1)){
  
  # get specimen
  specimen <- paste(srgb_dat[i, "PC1"], srgb_dat[i, "PC2"], sep = "_")
  specimen_rgb <- srgb_long %>% 
    filter(
      PC1 == srgb_dat[i, "PC1"],
      PC2 == srgb_dat[i, "PC2"]
    ) %>% 
    arrange(body_part)  # this determines the position of each body part in the grid
  
  # initialise png saving
  png(
    here::here(
      "2_Patches", "3_OutputData", "5_colour_grids", "2_Simulated_grids", paste0(specimen, ".png")
    ),
    # paste0("C:/Users/bop23rxm/Documents/colour_grids_repositioned/", specimen, ".png"),   # use local storage as faster
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
  cat(paste0("\r", i, " of ", length(srgb_dat$PC1), " processed"))
}


## Doesn't seem to work for now - it's because literally all the generated patch 
## colours are being converted to hot pink
## I think there must not be enough info in the first two axes to accurately regenerate
## an RGB code that works. Could try to fix this tomorrow
