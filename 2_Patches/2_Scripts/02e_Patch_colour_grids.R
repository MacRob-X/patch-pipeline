## Create visuals showing patch colours on colour space ----
## Robert MacDonald
## 13th May 2024

# clear environment
rm(list=ls())

## Load libraries ---- 
library(dplyr)
library(ggplot2)

# load rgb data ----
srgb_dat <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "1_RawColourspaces", "Neoaves.patches.231030.prePCAcolspaces.240404.rds"
  )
) %>% 
  magrittr::extract2("sRGB")

# load umap (output from 02b_Patch_Umap_iterations.R) - specify e.g. jndxyzlum
umap_jndxyzlum <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
    "Neoaves.patches.pca.jndxyzlum.UMAPs.iterations.240603.rds"
  )
) %>% 
  magrittr::extract2(1) %>% 
  magrittr::extract2("layout") %>% 
  as.data.frame() %>% 
  rename(
    umap_1 = V1, umap_2 = V2
  )

# load PCA data (output from 02_Patch_Analyse_features.R) - specify the one you want (e.g. jndxyzlum)
pca_all <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    "Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603.rds"
  )
)

#-------------------------------------------------------------------------------------------------------------------------#
# Creating colour grid pngs (only needs to be run once) ----

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

#-----------------------------------------------------------------------------------------------------------------------------#



## Plot UMAP with colour grids as points ----

# get x and y ranges
xrange <- c((min(umap_jndxyzlum$umap_1) - 0.5 ), (max(umap_jndxyzlum$umap_1) + 0.5 ))
yrange <- c((min(umap_jndxyzlum$umap_2) - 0.5 ), (max(umap_jndxyzlum$umap_2) + 0.5 ))

png(
  here::here(
    "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", "jndxyzlum",
    "patches_jndxyzlum_umap_colour_grids_240603.png"
  ),
  width = 4020, height = 2568,
  units = "px"
)

plot(umap_jndxyzlum[, 1] ~ umap_jndxyzlum[, 2], asp = T, type="n", xlab="UMAP1", ylab="UMAP2", las=1,
     xlim = xrange, ylim = yrange)

box(lwd = 2)
rw <- diff(range(umap_jndxyzlum[,1]))/12

for(i in 1:nrow(umap_jndxyzlum)){
  fname <- rownames(umap_jndxyzlum)[i]
  fpng <- png::readPNG(paste0(
    "C:/Users/bop23rxm/Documents/colour_grids_repositioned/", 
    paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
  rasterImage(fpng, 
              xleft = umap_jndxyzlum[i, 1] - (rw/15), 
              ybottom = umap_jndxyzlum[i, 2]-(rw/9), 
              xright = umap_jndxyzlum[i, 1]+(rw/15), 
              ytop = umap_jndxyzlum[i, 2]+(rw/9))
  cat(paste0("\r", i, " of ", nrow(umap_jndxyzlum), " processed"))
}

dev.off()

## Plot the PCA axes with colour grids as points ----

# get pca data
pca_jndxyzlum <- pca_all %>% 
  magrittr::extract2("x") %>% 
  as.data.frame()

# get proportions of variance for each axis
pca_variance <- pca_all %>% 
  summary() %>% 
  magrittr::extract2("importance") %>% 
  as.data.frame()


# get x and y ranges
xrange <- c((min(pca_jndxyzlum$PC7) - 0.5 ), (max(pca_jndxyzlum$PC7) + 0.5 ))
yrange <- c((min(pca_jndxyzlum$PC8) - 0.5 ), (max(pca_jndxyzlum$PC8) + 0.5 ))

# set axis labels
xlabel <- paste0("PC7 (", round(pca_variance$PC7[2], 2) * 100, "% of variance)")
ylabel <- paste0("PC8 (", round(pca_variance$PC8[2], 2) * 100, "% of variance)")

png(
  here::here(
    "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", "jndxyzlum", 
    "patches_jndxyzlum_pca_PC7_PC8_colour_grids.png"
  ),
  width = 4020, height = 2568,
  units = "px"
)

plot(pca_jndxyzlum$PC7 ~ pca_jndxyzlum$PC8, asp = T, type = "n", xlab = xlabel, ylab = ylabel, las = 1,
     xlim = xrange, ylim = yrange)

box(lwd = 2)
rw <- diff(range(pca_jndxyzlum[,1]))/12

for(i in 1:nrow(pca_jndxyzlum)){
  fname <- rownames(pca_jndxyzlum)[i]
  fpng <- png::readPNG(paste0(
    "C:/Users/bop23rxm/Documents/colour_grids_repositioned/", 
    paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
  rasterImage(fpng, 
              xleft = pca_jndxyzlum$PC7[i] - (rw/15), 
              ybottom = pca_jndxyzlum$PC8[i]-(rw/9), 
              xright = pca_jndxyzlum$PC7[i]+(rw/15), 
              ytop = pca_jndxyzlum$PC8[i]+(rw/9))
  cat(paste0("\r", i, " of ", nrow(pca_jndxyzlum), " processed"))
}

dev.off()


## Plot PCs 1-6 and UMAP 1-2 on the same figure ----

# get pca data
pca_jndxyzlum <- pca_all %>% 
  magrittr::extract2("x") %>% 
  as.data.frame()

# get proportions of variance for each axis
pca_variance <- pca_all %>% 
  summary() %>% 
  magrittr::extract2("importance") %>% 
  as.data.frame()

# Start png plotting

png(
  here::here(
    "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", "jndxyzlum",
    "Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603.pc1-6.umap.colourgrids.png"
  ),
  width = 1000, height = 1000,
  units = "px"
  )

# define rows/columns
par(mfcol=c(2,2),
    mar = c(5.1, 4.1, 2.1, 2.1))


## plot PCs 1-6

# set pcs to use
pcs <- seq(1, 6, by = 2)

for (i in pcs){
  
  # get x and y ranges
  xrange <- c((min(pca_jndxyzlum[, i]) - 0.5 ), (max(pca_jndxyzlum[, i]) + 0.5 ))
  yrange <- c((min(pca_jndxyzlum[, i+1]) - 0.5 ), (max(pca_jndxyzlum[, i+1]) + 0.5 ))
  
  # set axis labels
  xlabel <- paste0("PC", i, " (", round(pca_variance[2, i], 2) * 100, "% of variance)")
  ylabel <- paste0("PC", i+1, " (", round(pca_variance[2, i+1], 2) * 100, "% of variance)")
  
  plot(pca_jndxyzlum[, i] ~ pca_jndxyzlum[, i+1], asp = T, type = "n", xlab = xlabel, ylab = ylabel, las = 1,
       xlim = xrange, ylim = yrange)
  
  rw <- diff(range(pca_jndxyzlum[, i]))/12
  
  for(j in 1:nrow(pca_jndxyzlum)){
    fname <- rownames(pca_jndxyzlum)[j]
    fpng <- png::readPNG(paste0(
      "C:/Users/bop23rxm/Documents/colour_grids_repositioned/", 
      paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
    rasterImage(fpng, 
                xleft = pca_jndxyzlum[j, i] - (rw/10), 
                ybottom = pca_jndxyzlum[j, i+1]-(rw/6), 
                xright = pca_jndxyzlum[j, i]+(rw/10), 
                ytop = pca_jndxyzlum[j, i+1]+(rw/6))
    cat(paste0("\r", j, " of ", nrow(pca_jndxyzlum), " processed"))
  }
  
  
}

## plot UMAP

# get x and y ranges
xrange <- c((min(umap_jndxyzlum$umap_1) - 0.5 ), (max(umap_jndxyzlum$umap_1) + 0.5 ))
yrange <- c((min(umap_jndxyzlum$umap_2) - 0.5 ), (max(umap_jndxyzlum$umap_2) + 0.5 ))

plot(umap_jndxyzlum[, 1] ~ umap_jndxyzlum[, 2], asp = T, type="n", xlab="UMAP 1", ylab="UMAP 2", las=1,
     xlim = xrange, ylim = yrange)

rw <- diff(range(umap_jndxyzlum[,1]))/12

for(i in 1:nrow(umap_jndxyzlum)){
  fname <- rownames(umap_jndxyzlum)[i]
  fpng <- png::readPNG(paste0(
    "C:/Users/bop23rxm/Documents/colour_grids_repositioned/", 
    paste(strsplit(fname, split="-")[[1]][1:2], collapse = "_"), ".png"))
  rasterImage(fpng, 
              xleft = umap_jndxyzlum[i, 1] - (rw/10), 
              ybottom = umap_jndxyzlum[i, 2]-(rw/6), 
              xright = umap_jndxyzlum[i, 1]+(rw/10), 
              ytop = umap_jndxyzlum[i, 2]+(rw/6))
  cat(paste0("\r", i, " of ", nrow(umap_jndxyzlum), " processed"))
}


dev.off()


