## Create visuals showing patch colours on colour space ----
## Robert MacDonald
## 13th May 2024

# clear environment
rm(list=ls())

## Load libraries ---- 
library(dplyr)
library(ggplot2)

## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Passeriformes"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr", etc.)
space <- "lab"

# load rgb data ----
srgb_dat <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "1_RawColourspaces", 
    paste0(clade, ".patches.231030.prePCAcolspaces.rds")
  )
) %>% 
  magrittr::extract2("sRGB")

# load raw (pre-PCA) TCS xyz data
raw_patch_tcs <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "1_RawColourspaces", 
    paste0(clade, ".patches.231030.prePCAcolspaces.rds")
  )
) %>% 
  magrittr::extract2("xyz")

# load PCA data (output from 02_Patch_Analyse_features.R)
pca_all <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA", 
    paste0(clade, ".patches.231030.PCAcolspaces.rds")
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


## Plot UMAP with colour grids as points ----

# get x and y ranges
xrange <- c((min(umap$umap_1) - 0.5 ), (max(umap$umap_1) + 0.5 ))
yrange <- c((min(umap$umap_2) - 0.5 ), (max(umap$umap_2) + 0.5 ))

png(
  here::here(
    "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", space, 
    paste0(clade, "_patches_", space, "_umap_colour_grids.png")
  ),
  width = 4020, height = 2568,
  units = "px"
)

plot(umap[, 1] ~ umap[, 2], asp = T, type="n", xlab="UMAP1", ylab="UMAP2", las=1,
     xlim = xrange, ylim = yrange)

box(lwd = 2)
rw <- diff(range(umap[,1]))/12

for(i in 1:nrow(umap)){
  fname <- rownames(umap)[i]
  fpng <- png::readPNG(paste0(
    "C:/Users/bop23rxm/Documents/colour_grids_repositioned/", 
    paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
  rasterImage(fpng, 
              xleft = umap[i, 1] - (rw/15), 
              ybottom = umap[i, 2]-(rw/9), 
              xright = umap[i, 1]+(rw/15), 
              ytop = umap[i, 2]+(rw/9))
  cat(paste0("\r", i, " of ", nrow(umap), " processed"))
}

dev.off()


## Plot the first two PCA axes with colour grids as points
# get pca data
pca_dat <- pca_all %>% 
  magrittr::extract2("x") %>% 
  as.data.frame()

# get proportions of variance for each axis
pca_variance <- pca_all %>% 
  summary() %>% 
  magrittr::extract2("importance") %>% 
  as.data.frame()


# get x and y ranges
xrange <- c((min(pca_dat$PC1) - 0.5 ), (max(pca_dat$PC1) + 0.5 ))
yrange <- c((min(pca_dat$PC2) - 0.5 ), (max(pca_dat$PC2) + 0.5 ))

# set axis labels
xlabel <- paste0("PC1 (", round(pca_variance$PC1[2], 2) * 100, "% of variance)")
ylabel <- paste0("PC2 (", round(pca_variance$PC2[2], 2) * 100, "% of variance)")

png(
  here::here(
    "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", space, 
    paste0(clade, "_patches_", space, "_pca_pc1_pc2_colour_grids.png")
  ),
  width = 4020, height = 2568,
  units = "px"
)

plot(pca_dat$PC1 ~ pca_dat$PC2, asp = T, type = "n", xlab = xlabel, ylab = ylabel, las = 1,
     xlim = xrange, ylim = yrange)

box(lwd = 2)
rw <- diff(range(pca_dat[,1]))/12

for(i in 1:nrow(pca_dat)){
  fname <- rownames(pca_dat)[i]
  fpng <- png::readPNG(paste0(
    "C:/Users/bop23rxm/Documents/colour_grids_repositioned/", 
    paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
  rasterImage(fpng, 
              xleft = pca_dat$PC1[i] - (rw/15), 
              ybottom = pca_dat$PC2[i]-(rw/9), 
              xright = pca_dat$PC1[i]+(rw/15), 
              ytop = pca_dat$PC2[i]+(rw/9))
  cat(paste0("\r", i, " of ", nrow(pca_dat), " processed"))
}

dev.off()

## Plot raw tetrahedral colour space
# convert to tidy data type
# get species and sex and tidy up names
raw_patch_tcs <- raw_patch_tcs %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2),
  )  %>% 
  rename(
    bel_x = x.bel, bel_y = y.bel, bel_z = z.bel, 
    bre_x = x.bre, bre_y = y.bre, bre_z = z.bre, 
    cov_x = x.cov, cov_y = y.cov, cov_z = z.cov, 
    cro_x = x.cro, cro_y = y.cro, cro_z = z.cro, 
    fli_x = x.fli, fli_y = y.fli, fli_z = z.fli, 
    man_x = x.man, man_y = y.man, man_z = z.man, 
    nap_x = x.nap, nap_y = y.nap, nap_z = z.nap, 
    rum_x = x.rum, rum_y = y.rum, rum_z = z.rum, 
    tai_x = x.tai, tai_y = y.tai, tai_z = z.tai, 
    thr_x = x.thr, thr_y = y.thr, thr_z = z.thr, 
  )

# get x values
tcs_x <- raw_patch_tcs %>% 
  tidyr::pivot_longer(
    cols = c(bel_x, bre_x, cov_x, cro_x, fli_x, man_x, nap_x, rum_x, tai_x, thr_x),
    names_to = "body_part",
    values_to = "x"
  ) %>% 
  select(species, sex, body_part, x) %>% 
  mutate(
    body_part = stringr::str_replace(body_part, "_x", "")
  )
# get y values
tcs_y <- raw_patch_tcs %>% 
  tidyr::pivot_longer(
    cols = c(bel_y, bre_y, cov_y, cro_y, fli_y, man_y, nap_y, rum_y, tai_y, thr_y),
    names_to = "body_part",
    values_to = "y"
  ) %>% 
  select(species, sex, body_part, y) %>% 
  mutate(
    body_part = stringr::str_replace(body_part, "_y", "")
  )
# get z values
tcs_z <- raw_patch_tcs %>% 
  tidyr::pivot_longer(
    cols = c(bel_z, bre_z, cov_z, cro_z, fli_z, man_z, nap_z, rum_z, tai_z, thr_z),
    names_to = "body_part",
    values_to = "z"
  ) %>% 
  select(species, sex, body_part, z) %>% 
  mutate(
    body_part = stringr::str_replace(body_part, "_z", "")
  )
# combine into one longer dataframe
tcs_long <- tcs_x %>% 
  inner_join(
    tcs_y
  ) %>% 
  inner_join(
    tcs_z
  )

rownames(tcs_long) <- paste0(tcs_long$species, tcs_long$sex, )


# get x, y, and z ranges
xrange <- c((min(tcs_long$x) - 0.5 ), (max(tcs_long$x) + 0.5 ))
yrange <- c((min(tcs_long$y) - 0.5 ), (max(tcs_long$y) + 0.5 ))
zrange <- c((min(tcs_long$z) - 0.5 ), (max(tcs_long$z) + 0.5 ))


png(
  here::here(
    "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", "raw_tcs_xyz", "patches_rawtcs_xy_colour_grids.png"
  ),
  width = 4020, height = 2568,
  units = "px"
)

plot(tcs_long$x ~ tcs_long$y, asp = T, type = "n", xlab = "x", ylab = "y", las = 1)

box(lwd = 2)
rw <- diff(range(tcs_long[,"x"]))/12

for(i in 1:nrow(tcs_long)){
  fname <- paste(tcs_long$species[i], tcs_long$sex[i], sep = "-")
  fpng <- png::readPNG(paste0(
    "C:/Users/bop23rxm/Documents/colour_grids_repositioned/", 
    paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
  rasterImage(fpng, 
              xleft = tcs_long$x[i] - (rw/15), 
              ybottom = tcs_long$y[i]-(rw/9), 
              xright = tcs_long$x[i]+(rw/15), 
              ytop = tcs_long$y[i]+(rw/9))
  cat(paste0("\r", i, " of ", nrow(tcs_long), " processed"))
}

dev.off()
