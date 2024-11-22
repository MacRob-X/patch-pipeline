# Plot and analyse measures of diversity from patch colour pattern spaces
# 13th March 2024
# Robert MacDonald

# clear environment
rm(list=ls())

# load libraries
library(dplyr)
library(ggfortify)


# load patch diversity metric data ----
patch_div <- read.csv(
  here::here(
    "2_Patches", "3_OutputData", "4_Diversity_measures",
    "Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603.diversitymeasures.csv"
  ),
  row.names = 1
) %>% 
  mutate(
    species = sapply(strsplit(rownames(.), split = "-"), "[", 1),
    sex = sapply(strsplit(rownames(.), split = "-"), "[", 2)
  )

patch_div %>% glimpse()

# Plot mean distance to centroid and mean distance to nearest five neighbours by IUCN category ----

# order levels of iucn cat for plotting
patch_div <- patch_div %>% 
  mutate(
    iucn = factor(iucn, levels = c("CR", "EN", "VU", "NT", "LC", "All"))
  )


# create colour palette (colours taken from ggsci palette pal_simpsons(palette = "springfield"))
springfield <- c(HomerYellow = "#FED439", FrinkPink = "#FD8CC1", HomerBlue = "#709AE1", DuffRed = "#C80813", 
                 BartOrange = "#FD7446", BurnsPurple = "#370335", MargeBlue = "#197EC0", LisaOrange = "#F05C3B", 
                 NedGreen = "#46732E", MaggieBlue = "#71D0F5", MargeGreen = "#D5E4A2", BurnsGreen = "#075149", 
                 HomerGrey = "#8A9197", KentRed = "#91331F", BobGreen = "#1A9993", HomerBrown = "#D2AF81")


# create colour mapping
cols <- springfield[1:length(levels(patch_div$iucn))]

# plot mean distance to centroid
centr_dist_plot <- ggplot(data = patch_div, aes(x = iucn, y = centr_dist, fill = iucn)) + 
  geom_boxplot(notch = TRUE, show.legend = FALSE, fill = cols) + 
  xlab("Subset of IUCN categories") + ylab("Mean distance to centroid") +
  theme_bw()

# plot mean distance to nearest five neighbours
nn_5_plot <- ggplot(data = patch_div, aes(x = iucn, y = nn_5_dist, fill = iucn)) + 
  geom_boxplot(notch = TRUE, show.legend = FALSE, fill = cols) + 
  xlab("Subset of IUCN categories") + ylab("Mean distance to five nearest neighbours") +
  theme_bw()

# combine into single plot
diversity_plot_figure <- ggpubr::ggarrange(
  centr_dist_plot, nn_5_plot,
  labels = c("A", "B"),
  ncol = 2, nrow = 1
)
# visualise
diversity_plot_figure

# save to svg
svg(
  here::here(
    "2_Patches", "4_OutputPlots", "2_Diversity_measures",
    "Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603.diversitymeasures.iucn.plot.svg"
  ),
  width = 9, height = 5.25,
)

diversity_plot_figure

dev.off()


# Statistical analysis ----

## Perform ANOVA to check IUCN category and diversity relation ----

### Distance to centroid ----

# subset data to remove duplicate 'All' category rows, reorder levels to make LC the reference category,
# and perform anova
mod_centr_dist <- patch_div %>% 
  filter(
    iucn != "All"
  ) %>% 
  mutate(
    iucn = factor(iucn, levels = c("LC", "CR", "EN", "VU", "NT"))
  ) %>% 
  lm(formula = log(centr_dist) ~ iucn, data = .)

# check assumptions
autoplot(mod_centr_dist, smooth.colour = NA)

# create the anova table to get F and p values
anova(mod_centr_dist)

# inspect coefficients for each category
summary(mod_centr_dist)

# write table to csv
summary(mod_centr_dist) %>% 
  broom::tidy() %>% 
  readr::write_csv(
    file = here::here(
      "2_Patches", "3_OutputData", "4_Diversity_measures",
      "Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603.diversitymeasures.ANOVA.csv"
    )
  )

### Mean distance to five nearest neighbours ----
# subset data to remove duplicate 'All' category rows, reorder levels to make LC the reference category,
# and perform anova
mod_nn5_dist <- patch_div %>% 
  filter(
    iucn != "All"
  ) %>% 
  mutate(
    iucn = factor(iucn, levels = c("LC", "CR", "EN", "VU", "NT"))
  ) %>% 
  lm(formula = log(nn_5_dist) ~ iucn, data = .)

# check assumptions
autoplot(mod_nn5_dist, smooth.colour = NA)

# create the anova table to get F and p values
anova(mod_nn5_dist)

# inspect coefficients for each category
summary(mod_nn5_dist)

# append table to csv
summary(mod_nn5_dist) %>% 
  broom::tidy() %>% 
  readr::write_csv(
    file = here::here(
    "2_Patches", "3_OutputData", "4_Diversity_measures",
    "Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603.diversitymeasures.ANOVA.csv"
  ),
  append = TRUE
)

## Perform ANOVA with IUCN as binary threatened/non-threatened category ----

### Distance to centroid ----

# subset data to remove duplicate 'All' category rows, create binary IUCN category and perform anova
mod_centr_dist_bin <- patch_div %>% 
  filter(
    iucn != "All"
  ) %>% 
  mutate(
    iucn_bin = as.numeric(
      stringr::str_replace_all(iucn, c("LC" = "0", "NT" = "0", "CR" = "1", "EN" = "1", "VU" = "1"))
      )
  ) %>% 
  lm(formula = log(centr_dist) ~ iucn_bin, data = .)

# check assumptions
autoplot(mod_centr_dist_bin, smooth.colour = NA)

# create the anova table to get F and p values
anova(mod_centr_dist_bin)

# inspect coefficients for each category
summary(mod_centr_dist_bin)

# write table to csv
summary(mod_centr_dist_bin) %>% 
  broom::tidy() %>% 
  readr::write_csv(
    file = here::here(
      "2_Patches", "3_OutputData", "4_Diversity_measures",
      "Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603.diversitymeasures.ANOVAbinary.csv"
    )
  )

### Mean distance to five nearest neighbours ----

mod_nn5_dist_bin <- patch_div %>% 
  filter(
    iucn != "All"
  ) %>% 
  mutate(
    iucn_bin = as.numeric(
      stringr::str_replace_all(iucn, c("LC" = "0", "NT" = "0", "CR" = "1", "EN" = "1", "VU" = "1"))
    )
  ) %>% 
  lm(formula = log(nn_5_dist) ~ iucn_bin, data = .)

# check assumptions
autoplot(mod_nn5_dist_bin, smooth.colour = NA)

# create the anova table to get F and p values
anova(mod_nn5_dist_bin)

# inspect coefficients for each category
summary(mod_nn5_dist_bin)

# append table to csv
summary(mod_nn5_dist_bin) %>% 
  broom::tidy() %>% 
  readr::write_csv(
    file = here::here(
      "2_Patches", "3_OutputData", "4_Diversity_measures",
      "Neoaves.patches.231030.PCAcolspaces.jndxyzlum.240603.diversitymeasures.ANOVAbinary.csv"
    ),
    append = TRUE
  )