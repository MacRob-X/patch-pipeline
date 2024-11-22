# 3D UMAP of patch data
# 8th May 2024
# Robert MacDonald

# Clear environment
rm(list=ls())

# Load libraries ----
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

# PCA data (output from 02_Patch_Analyse_features)
pca_jndxyzlumr <- readr::read_rds(
  here::here(
    "2_Patches", "3_OutputData", "2_PCA_ColourPattern_spaces", "Neoaves.patches.231030.PCAcolspaces.jndxyzlumr.240404.rds"
  )
)

# taxonomy
taxo_raw <- read.csv(
  here::here(
    "4_SharedInputData", "BLIOCPhyloMasterTax_2016_02_29.csv"
  )
) %>% 
  janitor::clean_names()

# 

# Perform UMAP ----

# get raw pca values
pca_vals <- pca_jndxyzlumr %>% 
  magrittr::extract2("x")


# perform umap
umap_3d_jndxyzlumr <- umap::umap(pca_vals, n_components = 3)

# Data munging ----

# create df with species, sex, UMAP values, and taxon subgroups
df <- tibble(
  species = sapply(strsplit(rownames(umap_3d_jndxyzlumr$layout), split = "-"), "[", 1),
  sex = sapply(strsplit(rownames(umap_3d_jndxyzlumr$layout), split = "-"), "[", 2),
  umap_axis_1 = umap_3d_jndxyzlumr$layout[, 1],
  umap_axis_2 = umap_3d_jndxyzlumr$layout[, 2],
  umap_axis_3 = umap_3d_jndxyzlumr$layout[, 3],
) %>% 
  left_join(
    taxo_raw,
    by = join_by(species == tip_label)
  ) %>% 
  select(
    species, sex,
    umap_axis_1, umap_axis_2, umap_axis_3,
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


# Visualisation ----

# standard ggplots
ax12 <- df %>% 
  ggplot(aes(x = umap_axis_1, y = umap_axis_2, colour = taxon_subgroup)) + 
  geom_point() + 
  guides(colour = "none") + 
  xlim(-5, 5) + ylim(-5, 5) + 
  theme_bw()
ax23 <- df %>% 
  ggplot(aes(x = umap_axis_2, y = umap_axis_3, colour = taxon_subgroup)) + 
  geom_point() +
  guides(colour = "none") + 
  xlim(-5, 5) + ylim(-5, 5) + 
  theme_bw() 

cowplot::plot_grid(ax12, ax23, labels = "AUTO")

# interactive 3d plot
plotly::plot_ly(
  data = df,
  x = ~ umap_axis_1, y = ~ umap_axis_2, z = ~ umap_axis_3,
  color = ~ species
  
)

# standard plots with dichromatism
ax12_dichro <- df_with_pca %>% 
  filter(
    sex != "u",
    dichromatism >= 10.771,
    !is.na(dichromatism)  ) %>% 
  ggplot(aes(x = umap_axis_1, y = umap_axis_2, alpha = log(dichromatism))) + 
  geom_point(colour = "midnightblue") + 
  facet_wrap(~ sex) + 
  xlim(-5, 5) + ylim(-5, 5) +
  theme_bw()
ax23_dichro <- df_with_pca %>% 
  filter(
    sex != "u",
    dichromatism >= 10.771,
    !is.na(dichromatism)  ) %>% 
  ggplot(aes(x = umap_axis_2, y = umap_axis_3, alpha = log(dichromatism))) + 
  geom_point(colour = "midnightblue") + 
  facet_wrap(~ sex) + 
  xlim(-5, 5) + ylim(-5, 5) +
  theme_bw()
cowplot::plot_grid(ax12_dichro, ax23_dichro, labels = "AUTO", nrow = 2, ncol = 1)

