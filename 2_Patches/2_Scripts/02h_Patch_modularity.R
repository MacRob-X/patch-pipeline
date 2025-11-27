# Test for modularity in patch colouration
# 18/11/2025
# Robert MacDonald

# load libraries
library(dplyr)
library(MCMCglmm)


# Functions ----

# set columns of interest based on channel selection
set_resp_cols <- function(pix, channels){
  if(channels != "all"){
    response_cols <- names(pix)[which(grepl(channels, colnames(pix)))]
  } else if(channels == "all"){
    response_cols <- names(pix)
  }
  return(response_cols)
}

## EDITABLE CODE ##
# Select subset of species ("Neoaves" or "Passeriformes")
clade <- "Neognaths"
# Restrict to only species for which we have male and female data?
mf_restrict <- TRUE
# set space
space <- "lab"

# Load data ----

# Load in pre-PCA pixel data
# set filename
if(mf_restrict == TRUE){
  prepca_filename <- paste(clade, "matchedsex", "patches.250716.prePCAcolspaces.rds", sep = ".")
} else{
  prepca_filename <- paste(clade, "allspecimens", "patches.250716.prePCAcolspaces.rds", sep = ".")
}

pix <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", clade, "1_RawColourSpaces", 
    prepca_filename
  )
)[[space]]


# Analysis ----

# change pix names so the patch name is first
names(pix) <- sub("^(\\w+)\\.(\\w+)$", "\\2.\\1", names(pix))

# choose channel to restrict to ("L.", "a.", "b."). Assign "all" if want to cross compare all channels
channel_restrict <- "all"

# restrict to channels of interest
response_cols <- set_resp_cols(pix, channel_restrict)

# get number of response columns
n_rc <- length(response_cols)

# restrict data to response columns
dat <- pix[, response_cols]

# set prior, number of iterations, and burn-in 
prior <- list(
  R = list(V = diag(n_rc), nu = 0.002)  # Residual covariance
)
n_iter <- 130
burn_in <- 30

# specify formula based on response columns (multiple response variables)
model_formula <- as.formula(paste("cbind(", paste(response_cols, collapse = ", "), ") ~ trait"))

# run model
model <- MCMCglmm(model_formula,  # trait-specific means
                  rcov = ~ us(trait):units, # residual covariance matrix unstructured
                  family = rep("gaussian", n_rc),
                  data = dat,
                  prior = prior,
                  verbose = TRUE,
                  nitt = n_iter,
                  burnin = burn_in,
                  thin = 10)

# inspect model
summary(model)

# save model
model_filename <- paste(channel_restrict, n_iter, "it", burn_in, "burnin", "matchedsex_mcmcglmm_corrmat.rds", sep = "_")
saveRDS(
  model,
  here::here(
    "2_Patches", "3_OutputData", clade, "1_RawColourspaces",
    model_filename
  )
)

# reload model if needed
model_filename <- paste(channel_restrict, n_iter, "it", burn_in, "burnin", "matchedsex_mcmcglmm_corrmat.rds", sep = "_")
model <- readRDS(
  here::here(
    "2_Patches", "3_OutputData", clade, "1_RawColourspaces",
    model_filename
  )
)

# Quantify modules ----

# get posterior correlation matrices
post_cor <- posterior.cor(model$VCV)

# can use this to get 95% CIs for each correlation
cor_intervals <- coda::HPDinterval(post_cor, prob = 0.95)

# In a nice matrix format
# Posterior means
cor_mean <- matrix(apply(post_cor, 2, mean), n_rc, n_rc)
rownames(cor_mean) <- colnames(cor_mean) <- response_cols
diag(cor_mean) <- NA


# Lower bound of 95% CI
cor_lower <- matrix(cor_intervals[, "lower"], n_rc, n_rc)
rownames(cor_lower) <- colnames(cor_lower) <- response_cols
diag(cor_lower) <- NA
heatmap(cor_lower)

# Upper bound of 95% CI
cor_upper <- matrix(cor_intervals[, "upper"], n_rc, n_rc)
rownames(cor_upper) <- colnames(cor_upper) <- response_cols
diag(cor_upper) <- NA
heatmap(cor_upper)


# Method 1: clustering dendrogram across all channels and patches
# plot heatmap/dendrogram using posterior means
heatmap(cor_mean)
# save original cor_mean as cor_mean_all (we will overwrite cor_mean later, when verifying
# clusters)
cor_mean_all <- cor_mean

# Method 2: clustering dendrogram for individual channels
cor_mean_a <- cor_mean[grepl(".a", rownames(cor_mean), fixed = T), grepl(".a", colnames(cor_mean), fixed = T)]
cor_mean_b <- cor_mean[grepl(".b", rownames(cor_mean), fixed = T), grepl(".b", colnames(cor_mean), fixed = T)]
cor_mean_L <- cor_mean[grepl(".L", rownames(cor_mean), fixed = T), grepl(".L", colnames(cor_mean), fixed = T)]
heatmap(cor_mean_a)
heatmap(cor_mean_b)
heatmap(cor_mean_L)

# Method 3: average correlations within patches, across channels
block_size <- 3
n_dim <- 30
result_dim <- n_dim / block_size  # result_dim will be 10

# 2. Reshape the correlation matrix into a 4-dimensional array
# Dimensions: (3 rows, 10 blocks_down, 3 columns, 10 blocks_across)
array_blocks <- array(
  data = cor_mean, 
  dim = c(block_size, 
          result_dim, 
          block_size, 
          result_dim)
)

# 3. Permute the dimensions to group the 3x3 blocks correctly
# New dimensions: (3 rows, 3 columns, 10 blocks_down, 10 blocks_across)
# The 3x3 elements are now in dimensions 1 and 2.
permuted_array <- aperm(array_blocks, c(1, 3, 2, 4))

# 4. Apply the mean function to the 3x3 blocks
# We apply the mean over the first two dimensions (1 and 2) 
# for every element in the third and fourth dimensions (the 10x10 blocks).
block_means <- apply(permuted_array, c(3, 4), mean)

# 5. Convert the result back into a 10x10 matrix
averaged_block_matrix <- matrix(
  block_means, 
  nrow = result_dim, 
  ncol = result_dim, 
  byrow = FALSE # This preserves the correct block order
)

# get names
# The original names are in the dimnames of cor_mat
original_names <- rownames(cor_mean)

# Determine the group membership for each original variable.
# E.g., Var_1, Var_2, Var_3 belong to Block 1; Var_4, Var_5, Var_6 belong to Block 2, etc.
group_index <- rep(1:result_dim, each = block_size)

# Group the original names by their block
grouped_names <- split(original_names, group_index)

# Create a concise name for each new block, representing its constituent variables.
# E.g., "Block 1 (Var_1 to Var_3)"
block_names <- sapply(grouped_names, function(names) {
  substr(names, 1, 3)[[1]]
})

# Assign the new names to the resulting matrix
rownames(averaged_block_matrix) <- block_names
colnames(averaged_block_matrix) <- block_names

heatmap(averaged_block_matrix)

# Quantify cluster support ----

# NOTE that this code is ripped straight from Claude - it needs to be checked in the literature
# to ensure this is an appropriate method


#### WARNING - NOT CURRENTLY SUITABLE FOR USE, AS THIS A PRIORI SPECIFIES THE NUMBER OF CLUSTERS EXPECTED

# Convert each posterior sample to a correlation matrix
# CORRECTED: Create array to store all posterior correlation matrices
cor_array <- array(NA, dim = c(n_rc, n_rc, nrow(post_cor)))

for(i in 1:nrow(post_cor)) {
  # Use same logic as your cor_mean calculation
  cor_array[,,i] <- matrix(post_cor[i,], n_rc, n_rc)
  rownames(cor_array[,,i]) <- colnames(cor_array[,,i]) <- response_cols
  # Keep diagonal as 1 for calculations (don't set to NA)
}

# Verify this matches your cor_mean
test_mean <- apply(cor_array, c(1,2), mean)
# Now set diagonal to NA to match your cor_mean format
diag(test_mean) <- NA

# Check match (should be TRUE or very close)
print(all.equal(test_mean, cor_mean, check.attributes = FALSE))
print(paste("Max absolute difference:", 
            max(abs(test_mean - cor_mean), na.rm = TRUE)))

# Visual verification
par(mfrow = c(1,2))
image(cor_mean, main = "Your cor_mean", zlim = c(-1, 1))
image(test_mean, main = "Reconstructed mean", zlim = c(-1, 1))

# Step 3: Perform hierarchical clustering on each posterior sample
cluster_list <- list()

for(i in 1:nrow(post_cor)) {
  # Convert correlation to distance (1 - |r|)
  dist_i <- as.dist(1 - abs(cor_array[,,i]))
  
  # Hierarchical clustering
  cluster_list[[i]] <- hclust(dist_i, method = "average")
}

cat("Created", length(cluster_list), "dendrograms from posterior samples\n")

# Step 4: Calculate co-clustering probabilities
calculate_cocluster_prob <- function(cluster_list, n_clusters, var_names) {
  n_iter <- length(cluster_list)
  n_vars <- length(var_names)
  
  # Matrix to store co-clustering counts
  cocluster_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)
  rownames(cocluster_matrix) <- colnames(cocluster_matrix) <- var_names
  
  for(i in 1:n_iter) {
    # Get cluster assignments for this iteration
    clusters_i <- cutree(cluster_list[[i]], k = n_clusters)
    
    # For each pair, check if they're in the same cluster
    for(j in 1:n_vars) {
      for(k in 1:n_vars) {
        if(clusters_i[j] == clusters_i[k]) {
          cocluster_matrix[j,k] <- cocluster_matrix[j,k] + 1
        }
      }
    }
  }
  
  # Convert counts to probabilities
  cocluster_matrix <- cocluster_matrix / n_iter
  
  return(cocluster_matrix)
}

# Calculate for different numbers of clusters
cocluster_3 <- calculate_cocluster_prob(cluster_list, n_clusters = 3, response_cols)
cocluster_4 <- calculate_cocluster_prob(cluster_list, n_clusters = 4, response_cols)
cocluster_5 <- calculate_cocluster_prob(cluster_list, n_clusters = 5, response_cols)

# Visualize co-clustering probabilities
library(pheatmap)

pheatmap(cocluster_4,
         display_numbers = TRUE,
         number_format = "%.2f",
         main = "Co-clustering Probability (k=4)",
         fontsize_number = 6,
         color = colorRampPalette(c("white", "dodgerblue4"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_method = "average")

# Step 5: Calculate cophenetic correlation distribution
cophenetic_cors <- numeric(nrow(post_cor))

for(i in 1:nrow(post_cor)) {
  dist_i <- as.dist(1 - abs(cor_array[,,i]))
  cophenetic_cors[i] <- cor(dist_i, cophenetic(cluster_list[[i]]))
}

# Summary statistics
cat("\nCophenetic correlation summary:\n")
cat("  Mean:", round(mean(cophenetic_cors), 3), "\n")
cat("  95% CrI:", round(quantile(cophenetic_cors, c(0.025, 0.975)), 3), "\n")

# Step 6: Calculate cluster assignment stability
get_cluster_stability <- function(cluster_list, n_clusters, var_names) {
  n_iter <- length(cluster_list)
  n_vars <- length(var_names)
  
  # Matrix to store cluster assignments
  assignments <- matrix(NA, nrow = n_iter, ncol = n_vars)
  
  for(i in 1:n_iter) {
    assignments[i,] <- cutree(cluster_list[[i]], k = n_clusters)
  }
  
  # For each variable, find modal cluster and its probability
  stability <- data.frame(
    variable = var_names,
    modal_cluster = NA,
    probability = NA
  )
  
  for(j in 1:n_vars) {
    cluster_counts <- table(assignments[,j])
    stability$modal_cluster[j] <- as.numeric(names(which.max(cluster_counts)))
    stability$probability[j] <- max(cluster_counts) / n_iter
  }
  
  # Sort by cluster then by probability
  stability <- stability[order(stability$modal_cluster, -stability$probability),]
  
  return(stability)
}

stability_4 <- get_cluster_stability(cluster_list, n_clusters = 4, response_cols)
print(stability_4)

# Step 7: Create consensus dendrogram
dist_mean <- as.dist(1 - abs(cor_mean))
hc_mean <- hclust(dist_mean, method = "average")

# Plot with cluster boundaries
plot(hc_mean, 
     main = "Consensus Dendrogram\n(Posterior Mean Correlations)",
     xlab = "", 
     sub = "",
     cex = 0.8)
rect.hclust(hc_mean, k = 4, border = "red")


# Same, but for each channel separately

# Extract indices for each channel
L_vars <- grep("\\.L$", response_cols)
a_vars <- grep("\\.a$", response_cols)
b_vars <- grep("\\.b$", response_cols)
# Get variable names for each channel
L_names <- response_cols[L_vars]
a_names <- response_cols[a_vars]
b_names <- response_cols[b_vars]
cor_array_L <- cor_array[L_vars, L_vars, ]
cor_array_a <- cor_array[a_vars, a_vars, ]
cor_array_b <- cor_array[b_vars, b_vars, ]

# Now perform clustering separately for each channel
# L channel
cluster_list_L <- list()
for(i in 1:nrow(post_cor)) {
  dist_i <- as.dist(1 - abs(cor_array_L[,,i]))
  cluster_list_L[[i]] <- hclust(dist_i, method = "average")
}

# a channel
cluster_list_a <- list()
for(i in 1:nrow(post_cor)) {
  dist_i <- as.dist(1 - abs(cor_array_a[,,i]))
  cluster_list_a[[i]] <- hclust(dist_i, method = "average")
}

# b channel
cluster_list_b <- list()
for(i in 1:nrow(post_cor)) {
  dist_i <- as.dist(1 - abs(cor_array_b[,,i]))
  cluster_list_b[[i]] <- hclust(dist_i, method = "average")
}

# Create consensus dendrograms
dist_mean_L <- as.dist(1 - abs(cor_mean_L))
hc_mean_L <- hclust(dist_mean_L, method = "average")

dist_mean_a <- as.dist(1 - abs(cor_mean_a))
hc_mean_a <- hclust(dist_mean_a, method = "average")

dist_mean_b <- as.dist(1 - abs(cor_mean_b))
hc_mean_b <- hclust(dist_mean_b, method = "average")

# Visualize all three
par(mfrow = c(1,3))
plot(hc_mean_L, main = "Lightness (L*) Correlations", xlab = "", sub = "")

plot(hc_mean_a, main = "Red-Green (a*) Correlations", xlab = "", sub = "")
plot(hc_mean_b, main = "Blue-Yellow (b*) Correlations", xlab = "", sub = "")

# Calculate stability for each channel
stability_L <- get_cluster_stability(cluster_list_L, n_clusters = 3, L_names)
stability_a <- get_cluster_stability(cluster_list_a, n_clusters = 3, a_names)
stability_b <- get_cluster_stability(cluster_list_b, n_clusters = 3, b_names)

print("L* channel stability:")
print(stability_L)
print("\na* channel stability:")
print(stability_a)
print("\nb* channel stability:")
print(stability_b)
