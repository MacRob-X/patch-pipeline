# Reconstruct the theoretical colour pattern phenotype at each point in colour pattern space
# Robert MacDonald
# 7th August 2024

# Load libraries ----
library(dplyr)

# Clear environment
rm(list=ls())

# EDITABLE CODE # ----
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Passeriformes"
# Choose space to work with (usmldbl, usmldblr, xyz, xyzlum, xyzlumr, lab, cie, srgb, hex, 
# jndxyz, jndxyzlum, jndxyzlumr)
space <- "sRGB"


# Load in PCA colourspace ----
pca_space <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    paste0(clade, ".patches.231030.PCAcolspaces.rds")
  )
)[[space]]


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
res <- 10
grand_max <- max(pc_maxes[])
grand_min <- min(pc_mins)
new_max <- grand_max + (grand_max/10)
new_min <- grand_min - (grand_min/10)
increment <- (new_max - new_min) / (res-1)

# generate a matrix spanning from the grand min to the grand max, moving in our
# determined increments
# for each point with a given value in PC1, need res points in PC2, so total number of points
# is res ^ n_axes
n_axes <- 2
n_points <- res ^ n_axes
theor_pca <- matrix(NA, nrow = n_points, ncol = n_axes)


# # for 2 dimensions
# # for each point with value new_min in axis 1, need res points in axis 2
# j <- 1
# for(i in seq(1, n_points, by = res)){
#   theor_pca[i:(i+res-1), 1] <- new_min + (j-1)*increment
#   theor_pca[i:(i+res-1), 2] <- seq(from = new_min, to = new_max, by = increment)
#   j <- j + 1
# }
# 
# # Need to find a way to generalise the above to n_axes dimensions
# 
# 
# # for 3 dimensions
# j <- 1
# for(i in seq(1, n_points, by = res^2)){
#   theor_pca[i:(i+res^2-1), 3] <- rep(seq(from = new_min, to = new_max, by = increment), times = res)
#   theor_pca[i:(i+res^2-1), 2] <- sort(
#     rep(seq(from = new_min, to = new_max, by = increment), times = res), 
#     decreasing = FALSE
#     )
#   theor_pca[i:(i+res^2-1), 1] <- new_min + (j-1)*increment
#   j <- j + 1
# }
# 
# # for 4 dimensions, need to repeat the 3-dimensional generated space res times - this will
# # provide axes 2-4. Hold axis one constant for each repetition
# j <- 1
# for(i in seq(1, n_points, by = res^3)){
#   theor_pca[i:(i+(res^3)-1), 4] <- rep(seq(from = new_min, to = new_max, by = increment), times = res^2)
#   theor_pca[i:(i+(res^3)-1), 3] <- sort(
#     rep(seq(from = new_min, to = new_max, by = increment), times = res^2), 
#     decreasing = FALSE
#   )
#   theor_pca[i:(i+(res^3)-1), 2] <- new_min + (j-1)*increment
#   j <- j + 1
# }
# 
# # n dimensions
# dim1 <- seq(from = new_min, to = new_max, by = increment)
# dim2 <- cbind(
#   sort(rep(dim1, times = res)),
#   rep(dim1, times = res)
#   
# )
# dim3 <- cbind(
#   sort(rep(dim2[, 1], times = res)),
#   rep(dim2[, 1], times = res),
#   rep(dim2[, 2], times = res)
# )


# generalise the above in a for loop ##################
# same as above, but move left to right across columns (should be simpler)
dim1 <- matrix(
  c(rep(seq(from = new_min, to = new_max, by = increment), times = res^(n_axes-1)),
    rep(NA, times = (res^n_axes)*(n_axes-1))), 
  ncol = n_axes
)
for(i in 2:n_axes){
  
  # make copy of previous dimension
  dim_imin1 <- get(paste0("dim", i-1))
  dim_name <- paste0("dim", i)
  
  # assign new dimension values
  if(i < n_axes){
    pattern <- sort(rep(dim_imin1[1:res^(i-1), i-1], times = res))
    dim_i <-  cbind(
      dim_imin1[, 1:i-1],
      rep(pattern, times = n_points/length(pattern)),
      dim_imin1[, (i+1):n_axes]
    )
  } else if (i == n_axes){
    pattern <- sort(rep(dim_imin1[1:res^(i-1), i-1], times = res))
    dim_i <-  cbind(
      dim_imin1[, 1:i-1],
      rep(pattern, times = n_points/length(pattern))
    )
  }
  
  assign(dim_name, dim_i)
  rm(dim_imin1, dim_i)
}

# assign matrix of generated points to new variable and remove intermediate variables
# N.B. this code will need to be edited each time n_axes is changed
theor_pca <- dim5
rm(dim1, dim2, dim3, dim4, dim5)

#########################################################################################
### ALTERNATIVELY, generate a set number of points spaced randomly between the min and
### max of each individual PC
n_axes <- 10
n_points <- 1000
theor_pca <- matrix(NA, nrow = n_points, ncol = n_axes)

# generate data
for (i in 1:n_axes) {
  # get min and max of axis
  axis_min <- pc_mins[i]
  axis_max <- pc_maxes[i]
  # generate uniform distribution of n points between min and max and assign to 
  # matrix column
 # theor_pca[, i] <- runif(n_points, min = axis_min, max = axis_max)
  
  # OR generate a normal distribution of n points between min and max and assign to
  # matrix column - estimate sd using StdDev = ((max - min)*.997)/6
  # as 6 sigma ~99.7% of population for a normal distribution
  theor_pca[, i] <- rnorm(n_points, sd = ((axis_max - axis_min) * 0.997) / 6 )
}
# assign column names
colnames(theor_pca) <- paste0("PC", 1:n_axes)
# assign 2 random letters and 3 random numbers as rownames (to id each 'specimen')
# code modified from https://stackoverflow.com/questions/42734547/generating-random-strings
row_ids <- function(n_points) {
  a <- do.call(paste0, replicate(1, sample(LETTERS, n_points, TRUE), FALSE))
  paste0(a, sprintf("%03d", sample(9999, n_points, TRUE)), sample(LETTERS, n_points, TRUE))
}
rownames(theor_pca) <- row_ids(n_points)

#########################################################################################

# Reconstruct original data from generated PCA data ----

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
# rownames(srgb_dat) <- paste(round(theor_pca[, 1], 2), round(theor_pca[, 2], 2), sep = ".")

srgb_dat <- srgb_dat %>% 
  as.data.frame() %>% 
  mutate(
    id = rownames(.)
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
  select(id, body_part, r) %>% 
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
  select(id, body_part, g) %>% 
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
  select(id, body_part, b) %>% 
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

# Convert any impossible colours (i.e. r, g, or b < 0 or > 1) to hot pink (0.961, 0.411, 0.706)
pink_r <- 0.961; pink_g <- 0.411; pink_b <- 0.706
# check and convert r < 0
srgb_long[srgb_long$r < 0, "b"] <- pink_b
srgb_long[srgb_long$r < 0, "g"] <- pink_g
srgb_long[srgb_long$r < 0, "r"] <- pink_r
# check and convert g < 0
srgb_long[srgb_long$g < 0, "b"] <- pink_b
srgb_long[srgb_long$g < 0, "r"] <- pink_r
srgb_long[srgb_long$g < 0, "g"] <- pink_g
# check and convert b < 0
srgb_long[srgb_long$b < 0, "r"] <- pink_r
srgb_long[srgb_long$b < 0, "g"] <- pink_g
srgb_long[srgb_long$b < 0, "b"] <- pink_b
# check and convert r > 1
srgb_long[srgb_long$r > 1, "b"] <- pink_b
srgb_long[srgb_long$r > 1, "g"] <- pink_g
srgb_long[srgb_long$r > 1, "r"] <- pink_r
# check and convert g > 1
srgb_long[srgb_long$g > 1, "b"] <- pink_b
srgb_long[srgb_long$g > 1, "r"] <- pink_r
srgb_long[srgb_long$g > 1, "g"] <- pink_g
# check and convert b > 1
srgb_long[srgb_long$b > 1, "r"] <- pink_r
srgb_long[srgb_long$b > 1, "g"] <- pink_g
srgb_long[srgb_long$b > 1, "b"] <- pink_b


# Create colour grids for each specimen
for(i in 1:nrow(srgb_dat)){
  
  # get specimen
  specimen <- srgb_dat[i, "id"]
  specimen_rgb <- srgb_long %>% 
    filter(
      id == srgb_dat[i, "id"]
    ) %>% 
    arrange(body_part)  # this determines the position of each body part in the grid
  
  # initialise png saving (in X drive for space reasons)
  png(
    paste0(
      "X:/cooney_lab/Shared/Rob-MacDonald/simulated_colour_grids/", 
      space, "/", n_axes, "dim_random", "/", specimen, ".png"),
    # here::here(
    #   "2_Patches", "3_OutputData", "5_colour_grids", "2_Simulated_grids", paste0(specimen, ".png")
    # ),
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
  cat(paste0("\r", i, " of ", nrow(srgb_dat), " processed"))
}


## Plot the colour grids on the PC axes ----

# get simulated PCA data and add id to match to colour grid
theor_pca_df <- as.data.frame(theor_pca, row.names = rownames(theor_pca))
colnames(theor_pca_df) <- paste0("PC", 1:n_axes)

# subset to get specific points, e.g. hold 3 PCs steady at min or max while varying other PCs
subset <- theor_pca_df[theor_pca_df$PC3 == max(theor_pca_df$PC3) & 
                         theor_pca_df$PC4 == max(theor_pca_df$PC4) &
                         theor_pca_df$PC5 == max(theor_pca_df$PC5), ]

# get proportions of variance for each axis from original (non-reconstructed) PCA object
pca_variance <- pca_space %>% 
  summary() %>% 
  magrittr::extract2("importance") %>% 
  as.data.frame()


# get x and y ranges
xrange <- c((min(theor_pca_df$PC1) - 0.5 ), (max(theor_pca_df$PC1) + 0.5 ))
yrange <- c((min(theor_pca_df$PC2) - 0.5 ), (max(theor_pca_df$PC2) + 0.5 ))

# set axis labels
xlabel <- paste0("PC1 (", round(pca_variance$PC1[2], 2) * 100, "% of variance)")
ylabel <- paste0("PC2 (", round(pca_variance$PC2[2], 2) * 100, "% of variance)")


png(
  here::here(
    "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", "theoretical_colourspaces",
    space, paste("patches", space, "theoretical_colspace", "PC1", "PC2", "random10dim.png", sep = "_") 
  ),
  width = 4020, height = 2568,
  units = "px"
)


plot(theor_pca_df$PC1 ~ theor_pca_df$PC2, asp = T, type = "n", xlab = xlabel, ylab = ylabel, las = 1,
     xlim = xrange, ylim = yrange)

box(lwd = 2)
rw <- diff(range(theor_pca_df[,1]))/12

# set colour grid location
# gridloc <- here::here(
#   "2_Patches", "3_OutputData", "5_colour_grids", "2_Simulated_grids"
# )
gridloc <- paste0(
  "X:/cooney_lab/Shared/Rob-MacDonald/simulated_colour_grids/", 
  space, "/", n_axes, "dim_random", "/")

for(i in 1:nrow(theor_pca_df)){
  fname <- rownames(theor_pca_df)[i]
  fpng <- png::readPNG(paste0(
    gridloc, "/",
    fname, ".png"))
  rasterImage(fpng, 
              xleft = theor_pca_df$PC1[i] - (rw/(res*(2/3))), 
              ybottom = theor_pca_df$PC2[i]-(rw/(res*(2/5))), 
              xright = theor_pca_df$PC1[i]+(rw/(res*(2/3))), 
              ytop = theor_pca_df$PC2[i]+(rw/(res*(2/5))))
  cat(paste0("\r", i, " of ", nrow(theor_pca_df), " processed"))
}

dev.off()
