# Perform mapping of averaged patch pixel values to cone stimulation values and to various colour spaces
# (raw USML; TCS xyz; TCS xyzlum; TCS xyzlumr; CIE xyz; Lab; sRGB; Hex)
# Main author: Chris Cooney
# Edited by Robert MacDonald
# 04 April 2024

# clear environment
rm(list=ls())

# load libraries ----
library(pavo)
library(rearrr)
library(raster)
library(zoo)
library(rgl)

# Read in cone catch functions -----
# maps pixels in NEF image to cone stimulation values
source("./3_SharedScripts/00_Cone_mapping_functions_v8.R") # pixel vals need to be in range 0-100

# edited pavo function to calculate JND-space coords ----

jnd2xyz.crc <- function (coldistres, center = TRUE, rotate = TRUE, rotcenter = c("mean", "achro"), ref1 = "l", ref2 = "u", axis1 = c(1, 1, 0), axis2 = c(0, 0, 1))  {
  
  pos2 <- function(d12, d13, d23) {
    x3 <- d13
    if (!d23 < d12 + d13) {
      x3 <- x3 * -1
    }
    x3
  }
  pos3 <- function(d12, d13, d23) {
    x3 <- (d13^2 - d23^2 + d12^2)/(2 * d12)
    y3 <- rep(0, 2)
    y3sq <- d13^2 - x3^2
    if (y3sq > 0) {
      y3 <- sqrt(y3sq) * c(1, -1)
    }
    matrix(c(rep(x3, 2), y3), ncol = 2, dimnames = list(NULL, c("x", "y")))
  }
  pos4 <- function(d12, d14, d24, d34) {
    x4 <- (d14^2 - d24^2 + d12^2)/(2 * d12)
    y4 <- ((d14^2 - d34^2 + thirdpointxy["y"]^2 + thirdpointxy["x"]^2)/(2 * thirdpointxy["y"])) - (x4 * (thirdpointxy["x"]/thirdpointxy["y"]))
    z4 <- rep(0, 2)
    z4sq <- d14^2 - x4^2 - y4^2
    if (z4sq > 0) {
      z4 <- sqrt(z4sq) * c(1, -1)
    }
    matrix(c(rep(x4, 2), rep(y4, 2), z4), ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
  }
  ncone <- attr(coldistres, "ncone")
  if (as.numeric(ncone) < 2 || as.numeric(ncone) > 4) {
    stop("only methods for di-, tri- and tetrachromatic models are implemented so far", 
         call. = TRUE)
  }
  references <- attr(coldistres, "resref")
  references <- references[intersect(grep("jnd2xyzrrf", references$patch1, invert = TRUE, fixed = TRUE), grep("jnd2xyzrrf", references$patch2, fixed = TRUE)), ]
  if (all(is.na(coldistres$dL))) {
    coldistres <- coldistres[, !names(coldistres) %in% "dL"]
  }
  combined <- rbind(coldistres, references)
  colmat <- coldist2mat(combined)
  cdmat <- colmat[["dS"]]
  coords <- matrix(NA, nrow = nrow(cdmat), ncol = as.numeric(ncone) - 1, dimnames = list(row.names(cdmat), c("x", "y", "z")[seq(as.numeric(ncone) - 1)]))
  ptnames <- rownames(coords)
  if (ncone == "2") {
    coords[ptnames[1], ] <- 0
    coords[ptnames[2], ] <- cdmat[ptnames[1], ptnames[2]]
    coords[c(ptnames[-(1:2)]), ] <- do.call(rbind, lapply(ptnames[-(1:2)], function(x) { pos2(cdmat[ptnames[1], ptnames[2]], cdmat[ptnames[1], x], cdmat[ptnames[2], x])}))
  }
  if (ncone == "3") {
    coords[ptnames[1], ] <- c(0, 0)
    coords[ptnames[2], ] <- c(cdmat[ptnames[1], ptnames[2]], 0)
    coords[ptnames[3], ] <- pos3(cdmat[ptnames[1], ptnames[2]], cdmat[ptnames[1], ptnames[3]], cdmat[ptnames[2], ptnames[3]])[1, ]
    positions <- lapply(ptnames[-(1:3)], function(x) {
      pos3(cdmat[ptnames[1], ptnames[2]], cdmat[ptnames[1], x], cdmat[ptnames[2], x])
    })
    names(positions) <- ptnames[-(1:3)]
    eucdis <- lapply(positions, function(x) dist(rbind(x, coords[ptnames[3], ]))[c(2, 3)])
    whichdist <- lapply(names(eucdis), function(x) which.min(abs(eucdis[[x]] - cdmat[ptnames[3], x])))
    names(whichdist) <- names(eucdis)
    coords[names(eucdis), ] <- do.call(rbind, lapply(names(eucdis), function(x) positions[[x]][whichdist[[x]], ]))
  }
  if (ncone == "4") {
    coords[ptnames[1], ] <- c(0, 0, 0)
    coords[ptnames[2], ] <- c(cdmat[ptnames[1], ptnames[2]], 0, 0)
    thirdpointxy <- pos3(cdmat[ptnames[1], ptnames[2]], cdmat[ptnames[1], ptnames[3]], cdmat[ptnames[2], ptnames[3]])[1, ]
    coords[ptnames[3], ] <- c(thirdpointxy, 0)
    fourthpointxyz <- pos4(cdmat[ptnames[1], ptnames[2]], cdmat[ptnames[1], ptnames[4]], cdmat[ptnames[2], ptnames[4]], cdmat[ptnames[3], ptnames[4]])[1, ]
    coords[ptnames[4], ] <- fourthpointxyz
    positions <- lapply(ptnames[-(1:4)], function(x) {
      pos4(cdmat[ptnames[1], ptnames[2]], cdmat[ptnames[1], x], cdmat[ptnames[2], x], cdmat[ptnames[3], x])
    })
    names(positions) <- ptnames[-(1:4)]
    eucdis <- lapply(positions, function(x) dist(rbind(x, coords[4, ]))[c(2, 3)])
    whichdist <- lapply(names(eucdis), function(x) which.min(abs(eucdis[[x]] - cdmat[ptnames[4], x])))
    names(whichdist) <- names(eucdis)
    coords[names(eucdis), ] <- do.call(rbind, lapply(names(eucdis), function(x) positions[[x]][whichdist[[x]], ]))
  }
  if ("dL" %in% names(colmat)) {
    ldmat <- colmat[["dL"]]
    coords <- cbind(coords, lum = 0)
    coords[ptnames[1], "lum"] <- 0
    coords[ptnames[2], "lum"] <- ldmat[ptnames[1], ptnames[2]]
    coords[c(ptnames[-(1:2)]), "lum"] <- do.call(rbind, lapply(ptnames[-(1:2)], function(x) {pos2(ldmat[ptnames[1], ptnames[2]], ldmat[ptnames[1], x], ldmat[ptnames[2], x])}))
    # coords[, "lum"] <- coords[, "lum"] - coords["jnd2xyzrrf.achro", "lum"] # edit - unclear why this adjustment is needed
  }
  centroids <- colMeans(coords[grep("jnd2xyzrrf", rownames(coords), invert = TRUE, fixed = TRUE), , drop = FALSE])
  if (center) {
    coords <- sweep(coords, 2, centroids, "-")
  }
  jnd2xyzrrf.ctrd <- colMeans(coords[grep("jnd2xyzrrf", rownames(coords), invert = TRUE, fixed = TRUE), , drop = FALSE])
  chromcoords <- as.data.frame(coords[grep("jnd2xyzrrf", rownames(coords), invert = TRUE, fixed = TRUE), , drop = FALSE])
  refstosave <- as.data.frame(rbind(coords[grep("jnd2xyzrrf", rownames(coords), fixed = TRUE), , drop = FALSE], jnd2xyzrrf.ctrd))
  attr(chromcoords, "class") <- c("colspace", "jnd2xyz", "data.frame")
  attr(chromcoords, "resref") <- refstosave
  if (rotate) {
    rotcenter <- match.arg(rotcenter)
    rotarg <- list(jnd2xyzres = chromcoords, center = rotcenter, ref1 = ref1, ref2 = ref2, axis1 = axis1, axis2 = axis2)
    chromcoords <- do.call(jndrot, rotarg)
  }
  chromcoords
}



#setwd("~/Dropbox/Projects/Current/Bird_colouration/Colour_pattern_analysis/Patches/")

# ------------- #

# Input files ----

# load taxonomy
taxo <- read.csv("./4_SharedInputData/BLIOCPhyloMasterTax_2019_10_28.csv", strings = F)

# load patch pixel values (vRGB; uRGB)
px_master <- readRDS("./2_Patches/1_InputData/patches.231030.rds")

# Choose whether to examine all Neoaves or passerines only
# "neoaves" or "passerines"
group <- "passerines"

if(group == "neoaves"){
  # remove any galloanseriformes or palaeognaths
  px <- px_master |> 
    dplyr::left_join(
      taxo,
      by = dplyr::join_by(species == TipLabel)
    ) |>
    dplyr::filter(
      IOCOrder != "GALLIFORMES"
    ) |>
    dplyr::select(
      species, specimen, sex, view, region, coord.x, coord.y, vR, vG, vB, uR, uG, uB, min.r2
    )
} else if(group == "passerines"){
  # remove any non-passerines
  px <- px_master |> 
    dplyr::left_join(
      taxo,
      by = dplyr::join_by(species == TipLabel)
    ) |>
    dplyr::filter(
      PassNonPass == "PASSERIFORMES"
    ) |>
    dplyr::select(
      species, specimen, sex, view, region, coord.x, coord.y, vR, vG, vB, uR, uG, uB, min.r2
    )
}



# ------------- #

# Aggregate and select colour data ----

# Calculate averages for sexes within species (skip this step if you want data for individual specimens)
px <- aggregate(cbind(vR, vG, vB, uR, uG, uB) ~ species + sex + region, px, mean)

# # Select only species with both Male and Female data (if necessary)
# px <- subset(px, sex %in% c("M","F"))
# px <- subset(px, species%in%names(table(px$species))[table(px$species)==20])

# -------------------------------------------------------------------------------- #

# Map to colourspace ----

# ------------- #
# map - segments
# NOT USE
# px.seg <- data.frame(fun.camA.ideal.segments.ubgr.int2poly1.expanded(vR=px[,"vR"]*100, vG=px[,"vG"]*100, vB=px[,"vB"]*100, uR=px[,"uR"]*100, uB=px[,"uB"]*100)) # map
# px.seg[px.seg <= 0] <- 0.000001
# px.seg[px.seg > 1] <- 1
# tmp <- px.seg[,c("sU","sB","sG","sR")]
# tmp.b <- rowSums(tmp)
# tmp <- tmp / tmp.b # make relative
# colnames(tmp) <- c("S1","S2","S3","S4")
# attr(tmp, "relative") <- TRUE
# attr(tmp, "qcatch") <- "Qi"
# tmp <- colspace(tmp, space = "segment")
# tmp$H[is.na(tmp$H)] <- 0 # correct dark
# tmp$B <- tmp.b
# px.seg <- cbind(px.seg, tmp)

# ------------- #

# map - avian uvs (usml+dbl)
# this produces the raw cone catch values (uv; sw; mw; lw; dbl) and the relative cone catch values (u; s; m; l)
# it does NOT produce a scaled dbl value
### NOT USE - the raw and relative cone catch values are produced in the usml + scaled lum snippet below ###
# px.cone <- data.frame(fun.camA.ideal.uv.int2poly1.expanded(vR=px[,"vR"]*100, vG=px[,"vG"]*100, vB=px[,"vB"]*100, uR=px[,"uR"]*100, uB=px[,"uB"]*100)) # map
# px.cone[px.cone <= 0] <- 0.000001
# px.cone[px.cone > 1] <- 1
# tmp <- data.frame(px.cone[,c("uv","sw","mw","lw")])
# tmp <- tmp/rowSums(tmp) # make relative
# colnames(tmp) <- c("u","s","m","l")
# attr(tmp, "relative") <- TRUE
# attr(tmp, "qcatch") <- "Qi"
# tmp <- colspace(tmp, space = "tcs") # tcsplot(tcs)
# tmp <- tmp[,!colnames(tmp) %in% c("u.r","s.r","m.r","l.r")]
# px.cone <- cbind(px.cone, tmp)
# rm(tmp)

# map - avian uvs (usml + scaled lum)
# this produces the raw cone catch values (uv; sw; mw; lw; dbl) and the relative cone catch values (u; s; m; l)
# plus scaled dbl (luminance) values (lumr) with the distance between theoretical max and min dbl cone stimulations
# the same as the max distance between chromatic cone vertices in TCS xyz space (allowing comparison with the
# relative chromatic cone catch values). It also produces an unscaled luminance values (lum)
# Note that the raw and relative cone catch values are the same as those calculated in the above usml+dbl snippet,

vtcs <- data.frame(u=c(1,0,0,0),
                   s=c(0,1,0,0),
                   m=c(0,0,1,0),
                   l=c(0,0,0,1))
tcs <- tcspace(vtcs)
# max distance in tcs
d <- max(dist(tcs[,c("x","y","z")])) # use this value to scale the lum variable such that the distance between maximally different lum values (0 and 1) now equals d
px.cone.dbl <- data.frame(fun.camA.ideal.uv.int2poly1.expanded(vR=px[,"vR"]*100, vG=px[,"vG"]*100, vB=px[,"vB"]*100, uR=px[,"uR"]*100, uB=px[,"uB"]*100)) # map
px.cone.dbl[px.cone.dbl <= 0] <- 0.000001
px.cone.dbl[px.cone.dbl > 1] <- 1
tmp <- data.frame(px.cone.dbl[,c("uv","sw","mw","lw", "dbl")])
tmp[, 1:4] <- tmp[, 1:4]/rowSums(tmp[, 1:4]) # make relative
colnames(tmp) <- c("u","s","m","l", "lum")
attr(tmp, "relative") <- TRUE
attr(tmp, "qcatch") <- "Qi"
tmp <- colspace(tmp, space = "tcs") # tcsplot(tcs)
tmp <- tmp[,!colnames(tmp) %in% c("u.r","s.r","m.r","l.r")]
# Now scale the lum data so that the max possible size of lum is equal to the max distance between USML TCS vertices
tmp$lumr <- tmp$lum * d
px.cone.dbl <- cbind(px.cone.dbl, tmp)
# note that the above snippet returns usml, xyz, and scaled lum, so no need to have a separate usml or tcs xyz

# map - ciexyz, lab, sRGB, hex
px.cie <- data.frame(fun.camA.ideal.ciexyz.int2poly1.expanded(vR=px[,"vR"]*100, vG=px[,"vG"]*100, vB=px[,"vB"]*100)) # map
px.cie[px.cie <= 0] <- 0.000001
px.cie[px.cie > 1] <- 1
colnames(px.cie) <- paste0("cie", colnames(px.cie))
lab <- convertColor(px.cie, from = "XYZ", to = "Lab", to.ref.white = "E") # "Equal energy"
srgb <- convertColor(lab, from = "Lab", to = "sRGB"); colnames(srgb) <- c("sRGB.r","sRGB.g","sRGB.b")
hex <-  rgb(srgb[,1], srgb[,2], srgb[,3], maxColorValue = 1)
px.cie <- cbind(px.cie, lab, srgb)
px.cie$hex <- hex

# collate 
px <- cbind(px, px.cone.dbl, px.cie)

# map - JND xyzlum
# converts the relative cone catch values to JNDs (allowing for simulataneous analysis and comparison of chromatic and
# luminance aspects of colour)

# efficiently calculate JND-space coords (to be separately for each patch within sex)

px$x.jnd <- NA
px$y.jnd <- NA
px$z.jnd <- NA
px$lum.jnd <- NA


# note: a key issue is that the pavo 'jnd2xyz' function uses the first set of points in the dataset
# as fixed 'anchor points' against which to evaluate the position of all all other points in the datasets.
# Therefore to make points comparable across runs of the function applied to independent datasets or points,
# the same set of anchor points must be used in each case. Also, the function must be prevented from mean-centering
# the resulting jnd scores, which is the only edit to the function (see above). 
# Instead, mean centering can be done on all scores together at the end if required.

anchor_dat <- data.frame(u=c(1,0,0,0), s=c(0,1,0,0), m=c(0,0,1,0), l=c(0,0,0,1), dbl=c(0,1,0,1))
anchor_dat[anchor_dat==0] <- 0.000001
anchor_dat[,c("u","s","m","l")] <- anchor_dat[,c("u","s","m","l")] / rowSums(anchor_dat[,c("u","s","m","l")])

for (i in 1:nrow(px)) {
  
  cat(paste(i, "\r"))
  fpx <- px[i, c("u","s","m","l","dbl"), drop=F]
  fpx <- rbind(anchor_dat, fpx)
  fpx[fpx==0] <- 0.000001
  cd <- suppressMessages({ coldist(fpx, achromatic = T, qcatch = "Qi") })
  cd[cd==0] <- 0.000001
  jnd <- jnd2xyz.crc(cd, center = F, rotate = F)
  px[i,c("x.jnd","y.jnd","z.jnd","lum.jnd")] <- jnd[-(1:nrow(anchor_dat)),]
  
}



# map - JNDxyzlumr 
# converts the relative cone catch values to JNDs (allowing for simulataneous analysis and comparison of chromatic and
# luminance aspects of colour) - this also scales the luminance values so that the max theoretical distance between luminance
# stimulation values is the same as the max theoretical distance between chromatic stimulation values

# coords of tcs vertices
vtcs <- data.frame(u=c(1,0,0,0),
                   s=c(0,1,0,0),
                   m=c(0,0,1,0),
                   l=c(0,0,0,1),
                   lum=c(0,1,0,1))
tcs <- tcspace(vtcs)

# max distance in tcs
d <- max(dist(tcs[,c("x","y","z")]))

# scale lum for comparison
px$dbl.scaled <- px$dbl * d

# efficiently calculate JND-space coords (to be separately for each patch within sex)

px$x.jndlumr <- NA
px$y.jndlumr <- NA
px$z.jndlumr <- NA
px$lum.jndlumr <- NA

# note: a key issue is that the pavo 'jnd2xyz' function uses the first set of points in the dataset
# as fixed 'anchor points' against which to evaluate the position of all all other points in the datasets.
# Therefore to make points comparable across runs of the function applied to independent datasets or points,
# the same set of anchor points must be used in each case. Also, the function must be prevented from mean-centering
# the resulting jnd scores, which is the only edit to the function (see above). 
# Instead, mean centering can be done on all scores together at the end if required.

anchor_dat <- data.frame(u=c(1,0,0,0), s=c(0,1,0,0), m=c(0,0,1,0), l=c(0,0,0,1), dbl.scaled=c(0,1,0,1))
anchor_dat[anchor_dat==0] <- 0.000001
anchor_dat[,c("u","s","m","l")] <- anchor_dat[,c("u","s","m","l")] / rowSums(anchor_dat[,c("u","s","m","l")])

for (i in 1:nrow(px)) {
  
  cat(paste(i, "\r"))
  fpx <- px[i, c("u","s","m","l","dbl.scaled"), drop=F]
  fpx <- rbind(anchor_dat, fpx)
  fpx[fpx==0] <- 0.000001
  cd <- suppressMessages({ coldist(fpx, achromatic = T, qcatch = "Qi") })
  cd[cd==0] <- 0.000001
  jnd <- jnd2xyz.crc(cd, center = F, rotate = F)
  px[i,c("x.jndlumr","y.jndlumr","z.jndlumr","lum.jndlumr")] <- jnd[-(1:nrow(anchor_dat)),]
  
}




# save
saveRDS(px, "./2_Patches/3_OutputData/1_RawColourspaces/Passeriformes.patches.231030.rawcolspaces.rds")



# ============================= #



