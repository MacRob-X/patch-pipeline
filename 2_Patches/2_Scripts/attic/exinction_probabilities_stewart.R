# Playing around with extinction probabilities data sent to me by Kerry Stewart
# See paper at https://doi.org/10.1101/2025.01.09.632120
# Robert MacDonald
# 29th April 2025

# clear environment
rm(list=ls())

# Load libraries
library(dplyr)

# function to calculate extinction probabilities of input matrix
# probabilities of extinction by calculating proportion of times each species goes extinct
calc_extinct_prob <- function(threats_matrix){
  
  extinct_prob <- (nrow(threats_matrix) - colSums(threats_matrix)) / nrow(threats_matrix)
  
  return(extinct_prob)
}

# Load in data, convert to matrix
threats_1000 <- read.csv(
  here::here(
    "4_SharedInputData", "all_pamat_all_threats_1000_17122024.csv"
  ), row.names = 1
)

# remove first row (which is just 1 for all species - the "full" scenario)
threats_1000 <- threats_1000[2:nrow(threats_1000), ]

# extract the scenario rownames and unique scenarios
scen_rownames <- sub("_[^_]*$", "", rownames(threats_1000))
scenarios <- unique(scen_rownames)
names(scenarios) <- scenarios

# calculate extinction probabilities under each scenario
extinct_probs <- lapply(
  scenarios, function(scenario, threat_dat, scen_rownames){
    
    # limit the threat data to only the required scenario
    threat_dat <- threat_dat[which(scen_rownames == scenario), ]
    
    # calculate the extinction probability for this scenario
    extinct_prob_scen <- calc_extinct_prob(threat_dat)
    
  }, 
  threat_dat = threats_1000, scen_rownames = scen_rownames
)

# collate into dataframe
extinct_probs_df <- as.data.frame(t(do.call(rbind, extinct_probs)))

# add species rowhttp://127.0.0.1:15011/graphics/plot_zoom_png?width=1078&height=785
extinct_probs_df$species <- rownames(extinct_probs_df)

# Save as CSV for downstream analysis
write.csv(
  extinct_probs_df,
  here::here("4_SharedInputData", "extinction_probabilities_stewart_rxm.csv"), 
  row.names = FALSE
)


# Explore distribution

# load CSV of extinction probabilities
probs <- read.csv(
  here::here("4_SharedInputData", "extinction_probabilities_stewart_rxm.csv")
)

scenarios <- colnames(probs)[1:7]

for(scenario in scenarios){
  hist(log(probs[[scenario]]), breaks = 50, xlab = paste0("Probabilities (", scenario, ")"))
}

# divide into quantiles (for baseline only)
nquants <- 10
prob_quants <- quantile(probs[["none"]][probs[["none"]]>0.01], probs = seq(0, 1, by = 1/nquants), na.rm = TRUE)
