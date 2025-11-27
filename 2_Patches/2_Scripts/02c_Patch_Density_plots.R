## Create 2d density plots to visualise density of colourspace fill
## 3rd May 2024
## Robert MacDonald

# Clear environment
rm(list=ls())

## Load libraries ----
library(dplyr)
library(ggplot2)

## Functions ----

# Check collation of UMAP and PCA data is ordered correctly
check_umap_pca <- function(
    df1,
    df2,
    pca
) {
  if(
    all(
      # check umap axis 1 values are correctly ordered (i.e., match the dataset which doesn't
      # include the pca values)
      identical(
        df1 %>% magrittr::extract2("umap_axis_1"),
        df2 %>% magrittr::extract2("umap_axis_1")
      ),
      # check also for umap 2 axis values
      identical(
        df1 %>% magrittr::extract2("umap_axis_2"),
        df2 %>% magrittr::extract2("umap_axis_2")
      ),
      # check that pc1 values match those in the raw pca matrix
      identical(
        df1 %>% magrittr::extract2("PC1"),
        pca %>% as_tibble() %>%  magrittr::extract2("PC1")
      )
    )
  ) {
    TRUE
  }
  else 
    stop("Problem with data collation - UMAP or PCA column values do not match among datasets")
}

# Derive sexual dichromatism score (Euclidean distance between male and female positions of each species
# in PCA colour pattern space)
derive_dichro_score <- function(data) {
  
  # create separate male and female PCA data frames
  males <- data %>% 
    filter(
      sex == "m"
    )
  females<- data %>% 
    filter(
      sex == "f"
    )
  
  # subset to only species which have both male and female data
  males <- males %>% 
    filter(
      species %in% females$species
    )
  females <- females %>% 
    filter(
      species %in% males$species
    )
  
  # check order matches in male and female data sets
  stopifnot(
    identical(males$species, females$species)
  )
  
  
  # create empty df to populate with euclidean distances
  euc_distances <- tibble(
    species = males$species,
    distance = rep(NA, times = nrow(males))
  )
  
  # calculate distance for each male-female pair
  for(i in 1:nrow(males)){
    cat("\r", i)
    # combine each male-female pair into separate dataframe
    tmp <- bind_rows(males[i, ], females[i, ])
    # calculate multidimensional euclidean distance between pair and populate distances df
    euc_distances[i, "distance"] <- as.numeric(dist(tmp, method = "euclidean"))
  }
  
  # attach distances to original data to get matching vector to return
  to_return <- data %>% 
    left_join(
      euc_distances,
      by = "species"
    ) %>% 
    select(distance)
  
  # return distances as vector
  return(to_return$distance)
  
}

# Assert collation of UMAP and PCA data is ordered correctly
assert_umap_pca <- checkmate::makeAssertionFunction(check_umap_pca)

# Load data ----

## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neognaths"
# select type of colour pattern space to use ("jndxyzlum", "usmldbl", "usmldblr")
space <- "lab"
# select whether to use matched sex data (""all" or "matchedsex")
# "matchedsex" will use diversity metrics calculated on a subset of data containing only species
# for which we have both a male and female specimen (and excluding specimens of unknown sex)
sex_match <- "matchedsex"

# PCA data (output from 02_Patch_Analyse_features)
pca_filename <- paste(clade, sex_match, "patches.250716.PCAcolspaces", "rds", sep = ".")
pca_all <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "1_Raw_PCA",
    pca_filename
  )
) %>% 
  magrittr::extract2(space)

# UMAP data (output from 02b_Patch_Umap_iterations.R)
umap_filename <- paste(clade, sex_match, "patches.nn.25.mindist.0.1",  space, "UMAP", "rds", sep = ".")
umap_all <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "2_UMAP",
    umap_filename
  )
) %>% 
  magrittr::extract2("layout")

# taxonomy
taxo_raw <- read.csv(
  here::here(
    "4_SharedInputData", "BLIOCPhyloMasterTax_2019_10_28.csv"
  )
) %>% 
  janitor::clean_names()

# Data munging ----

# create df with species, sex, UMAP values, and taxon subgroups
df <- tibble(
  species = sapply(strsplit(rownames(umap_all), split = "-"), "[", 1),
  sex = sapply(strsplit(rownames(umap_all), split = "-"), "[", 2),
  umap_axis_1 = umap_all[, 1],
  umap_axis_2 = umap_all[, 2]
) %>% 
  left_join(
    taxo_raw,
    by = join_by(species == tip_label)
  ) %>% 
  select(
    species, sex,
    umap_axis_1, umap_axis_2,
    ioc_order, pass_non_pass, taxon_subgroup
  )

# clean character values
df <- df %>% 
  mutate(
    species = snakecase::to_snake_case(species),
    sex = snakecase::to_snake_case(sex),
    ioc_order = snakecase::to_snake_case(ioc_order),
    pass_non_pass = snakecase::to_snake_case(pass_non_pass),
    taxon_subgroup = snakecase::to_snake_case(taxon_subgroup)
  )

df %>% glimpse()

# append raw PCA values
pca_vals <- pca_all %>% 
  magrittr::extract2("x")

df_with_pca <- tibble(
  species = sapply(strsplit(rownames(pca_vals), split = "-"), "[", 1),
  sex = sapply(strsplit(rownames(pca_vals), split = "-"), "[", 2),
) %>% 
  cbind(pca_vals) %>% 
  mutate(
    species = snakecase::to_snake_case(species),
    sex = snakecase::to_snake_case(sex)
  ) %>% 
  inner_join(
    df,
    by = join_by(species, sex)
  ) %>% 
  assert_umap_pca(df, pca_vals)

df_with_pca %>% glimpse()

# derive scores of sexual dichromatism (euclidean distance between male and female of species
# within PCA colour pattern space)
df_with_pca <- df_with_pca %>% 
  mutate(
    dichromatism = derive_dichro_score(.)
  )

# Plotting ----

## Standard scatter plot ----
df %>% 
  ggplot(aes(x = umap_axis_1, y = umap_axis_2, colour = taxon_subgroup)) + 
  geom_point() + 
  theme_bw() 

# males vs females (discarding unknown sex specimens)
df %>% 
  filter(
    sex != "u"
  ) %>% 
  ggplot(aes(x = umap_axis_1, y = umap_axis_2, colour = taxon_subgroup)) + 
  geom_point() + 
  facet_wrap(~ sex) + 
  theme_bw()

# passerines vs non-passerines
df %>% 
  ggplot(aes(x = umap_axis_1, y = umap_axis_2, colour = taxon_subgroup)) + 
  geom_point() + 
  facet_wrap(~ pass_non_pass) + 
  theme_bw()
  

## Plot heatmap-type plots ----

# plot hexbin chart
  ggplot(df, aes(x = umap_axis_1, y = umap_axis_2)) + 
  geom_hex(bins = 80) + 
  scale_fill_viridis_c(option = "plasma") + 
  theme_bw()

# males vs females (discarding unknown sex specimens)
# Also log transform density
df %>%
  filter(
    sex != "u"
  ) %>% 
  ggplot(aes(x = umap_axis_1, y = umap_axis_2)) + 
  geom_hex(bins = 50) + 
  scale_fill_viridis_c(option = "plasma", trans = "log") + 
  facet_wrap(~ sex) + 
  theme_bw()

# passerines vs non-passerines
df %>% 
  ggplot(aes(x = umap_axis_1, y = umap_axis_2)) + 
  geom_hex(bins = 50) + 
  scale_fill_viridis_c(option = "plasma") + 
  facet_wrap(~ pass_non_pass) + 
  theme_bw()


# plot 2d density distribution
df %>% 
  ggplot(aes(x = umap_axis_1, y = umap_axis_2)) + 
  stat_density2d(aes(fill = after_stat(level)), geom = "polygon", n = 50) + 
  theme_bw()

# males vs females (discarding unknown sex specimens)
df %>%
  filter(
    sex != "u"
  ) %>% 
  ggplot(aes(x = umap_axis_1, y = umap_axis_2)) + 
  stat_density2d(aes(fill = after_stat(level)), geom = "polygon", n = 50) + 
  scale_fill_viridis_c(option = "plasma") +
  facet_wrap(~ sex) + 
  theme_bw()

# passerines vs non-passerines
df %>% 
  ggplot(aes(x = umap_axis_1, y = umap_axis_2)) + 
  stat_density2d(aes(fill = after_stat(level)), geom = "polygon", n = 50) + 
  facet_wrap(~ pass_non_pass) + 
  theme_bw()


## Plot dichromatism ----

# first check distribution of dichromatism scores - if highly skewed, might make sense to 
# plot the log-trasformed values
df_with_pca %>% 
  ggplot(aes(x = dichromatism)) + 
  geom_histogram(bins = 100)
# looks very right-skewed, so may be an idea to log-transform prior to plotting or analysis
df_with_pca %>% 
  ggplot(aes(x = log(dichromatism))) + 
  geom_histogram(bins = 100)
# much better

# facet by sex, to see how the colour of male and female specimens of dichromatic species varies
# in terms of position in colourspace
# removing species with dichromatism values below a certain threshold (e.g. 1st quartile) allows
# us to more clearly see how males of highly dichromatic species tend to be at the edge of colour space,
# while females in these species tend to be more central
df_with_pca %>% 
  filter(
 #   pass_non_pass == "nonpasseriformes",
    sex != "u",
    dichromatism >= quantile(dichromatism, na.rm = T)[2],
    !is.na(dichromatism)  ) %>% 
  ggplot(aes(x = umap_axis_1, y = umap_axis_2, alpha = log(dichromatism))) + 
  geom_point(colour = "midnightblue") + 
  facet_wrap(~ sex) + 
  theme_bw()
# we can see that male colouration in species which are dichromatic tends to be towards the top-right
# of the umap colour space
# in contrast, in highly dichromatic species female colouration tends to be more evenly distributed,
# and closer to the centroid
# this suggests that in most dichromatic species, it's the male which has the more unusual colouration

## Plot IUCN category
# load iucn data
iucn <- read.csv(
  here::here(
    "4_SharedInputData", "iucn_2024_nominate.csv"
  )
) %>% 
  mutate(
    species = snakecase::to_snake_case(species_birdtree)
  )

# add to df
df_iucn <- df %>% 
  left_join(iucn, by = "species") %>% 
  mutate(
    iucn_cat = factor(iucn_cat, levels = c("DD", "LC", "NT", "VU", "EN", "CR", "EX"))
  ) %>% 
  left_join(
    df_with_pca,
    by = c("species", "sex")
  )

# plot with transparency according to iucn category
df_iucn %>% 
  filter(
    iucn_cat != "DD",
    sex != "u"
  ) %>% 
  ggplot(aes(x = PC1, y = PC2, alpha = iucn_cat)) + 
  geom_point(colour = "darkred") + 
  facet_wrap(~ sex) + 
  theme_bw()
  
# 

# now plot interactively
# note that we need colour = species here so that the interactive plot will show species name
# when hovered over
p <- df_with_pca %>% 
  filter(
    sex != "u",
    !is.na(dichromatism)
  ) %>% 
    ggplot(aes(x = umap_axis_1, y = umap_axis_2, alpha = dichromatism, colour = species)) + 
    geom_point(show.legend = FALSE) + 
    facet_wrap(~ sex) + 
    theme_bw()
plotly::ggplotly(p)
  
# save interactive plot as html widget
plot_filename <- paste("interactive", clade, "patches", space, "umap_bysex.html", sep = "_")
htmlwidgets::saveWidget(
  plotly::ggplotly(p),
  here::here(
    "2_Patches", "4_OutputPlots", "1_Colourspace_visualisation", 
    plot_filename
  )
)

# Plot individual PC axis distributions ----
library(ggridges)
# first get long version of data
df_pca_long <- df_with_pca %>% 
  tidyr::pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "value") %>% 
  mutate(PC = as.numeric(gsub("PC", "", PC)))

# Plot density plots as ridgeline plot
df_pca_long %>% 
  filter(PC <= 15) %>% 
  ggplot(aes(x = value, y = as.factor(PC), fill = as.factor(PC))) +
  geom_density_ridges()

# PC1 looks pretty non-normal - check the QQ plot
df_with_pca %>% 
  select(PC1) %>% 
  ggplot(aes(sample = PC1)) + 
  stat_qq() + 
  stat_qq_line()

# Same for PCs 1-10
df_pca_long %>% 
  filter(PC <= 10) %>% 
  ggplot(aes(sample = value, colour = factor(PC))) + 
  stat_qq() +
  stat_qq_line() + 
  facet_wrap(~ factor(PC))
  