# Custom Feature Space Blocking with Stratification for Large Datasets ####

# This function creates training, calibration, and test splits based on
# clustering in feature space. It is designed to handle large datasets
# efficiently by sampling for clustering and then assigning all data points to
# the nearest cluster. It uses the HDBSCAN algorithm for clustering and Gower
# distance for mixed data types.

# The specific use case motivating this function is a dataset of approximately N
# = 800e3, p = 20, with a binary target variable that is highly imbalanced
# (~0.5% positive cases). The intended downstream learner is Random Forest
# Quantile Classifiers (RFQ), potentially with a post-hoc calibration step.

library(readr)
library(dplyr)
library(terra)
library(cluster)
library(dbscan)
library(gower)  # For Gower distance calculations

## Create training dataset ####
rf_global_current <- rast("output/rf_global_pred_regional_current.tif") 
names(rf_global_current) <- "rf_global"
chelsa_nor_250m_cur <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")[[1:19]]
ar50_nor_250m <- rast("output/ar50_250m_EPSG3035.tif")
names(ar50_nor_250m) <- "ar50"
ar50_nor_250m <- as.factor(ar50_nor_250m) # Convert to factor
preds_train <- c(rf_global_current, chelsa_nor_250m_cur, ar50_nor_250m)

presence <- read_csv("output/presence_coords_regional.csv")
absence <- read_csv("output/absence_coords_regional.csv")

# Extract predictor values for presence points
presence_df <- terra::extract(preds_train, presence[c("x", "y")], xy=TRUE, ID=FALSE) |> 
  mutate(response = 1, .before = 1) |>
  drop_na() |> 
  tibble()

# Extract predictor values for absence points
absence_df <- terra::extract(preds_train, absence[c("x", "y")], xy=TRUE, ID=FALSE) |> 
  mutate(response = 0, .before = 1) |>
  drop_na() |> 
  tibble()

# Combine into training dataset
training_data <- bind_rows(presence_df, absence_df) |>
  mutate(response = factor(response), ar50 = factor(ar50)) |> 
  select(-x, -y) # Remove x and y coordinates

# Get final counts
n_final_presence <- sum(training_data$response == "1")
n_final_background <- sum(training_data$response == "0")
n_final_prevalence <- n_final_presence / (n_final_presence + n_final_background)

cat("Final training data - Presences:", n_final_presence, "Background:", n_final_background, "\n")
cat("Prevalence in training data:", round(n_final_prevalence * 100, 3), "%\n")

# Feature weights for feature distances
weights_gini <- read_csv("output/feature_weights_gini.csv")
weights <- pull(weights_gini, weights_gini)
names(weights) <- pull(weights_gini, feature)

## Partition training data into train/calib/test splits ####

create_feature_space_splits_large <- function(data, 
                                              feature_cols,
                                              target_col = "response",
                                              sample_size = 20000,  # Size for clustering
                                              sample_ratio_minority = 0.3,  # Oversample minority in clustering
                                              train_prop = 0.6,
                                              calib_prop = 0.2,
                                              test_prop = 0.2,
                                              min_cluster_size = 10,
                                              cluster_selection_epsilon = 0.0,  # HDBSCAN parameter
                                              minimum_clusters = 10,  # Minimum number of clusters to force if HDBSCAN produces too few
                                              weights = NULL,  # Weights for Gower distance
                                              seed = 42) {
  
  set.seed(seed)
  
  # Step 1: Create stratified sample for clustering
  cat("Creating stratified sample for clustering...\n")
  
  minority_idx <- which(data[[target_col]] == 1)
  majority_idx <- which(data[[target_col]] == 0)
  
  n_minority <- length(minority_idx)
  n_majority <- length(majority_idx)
  
  # Calculate sample sizes
  sample_minority <- min(n_minority, round(sample_size * sample_ratio_minority))
  sample_majority <- min(n_majority, sample_size - sample_minority)
  
  # Sample indices
  sampled_minority_idx <- sample(minority_idx, sample_minority)
  sampled_majority_idx <- sample(majority_idx, sample_majority)
  sampled_idx <- c(sampled_minority_idx, sampled_majority_idx)
  
  # Create sample dataset
  sample_data <- data[sampled_idx, ]
  
  cat(sprintf("Sampled %d observations (%d positive, %d negative) from %d total\n",
              length(sampled_idx), sample_minority, sample_majority, nrow(data)))
  
  # Step 2: Calculate Gower distances on sample
  cat("Calculating Gower distances on sample...\n")
  
  gower_dist_sample <- cluster::daisy(
    x = sample_data[, feature_cols],
    metric = "gower",
    weights = weights
  )
  
  # Step 3: Cluster the sample with HDBSCAN
  cat("Clustering sample with HDBSCAN...\n")
  hdb_result <- hdbscan(gower_dist_sample, 
                        minPts = min_cluster_size,
                        cluster_selection_epsilon = cluster_selection_epsilon)
  
  # Check initial clustering results
  initial_clusters <- length(unique(hdb_result$cluster[hdb_result$cluster > 0]))
  cat(sprintf("HDBSCAN found %d initial clusters\n", initial_clusters))
  
  # If we have too few clusters, use cutree to force more clusters
  if(initial_clusters < minimum_clusters) {
    cat(sprintf("Forcing %d clusters using cutree on HDBSCAN hierarchy...\n", minimum_clusters))
    
    # Use cutree on the hierarchical clustering component
    forced_clusters <- cutree(hdb_result$hc, k = minimum_clusters)
    
    # Create new hdbscan-like result object
    hdb_clusters <- list(
      cluster = forced_clusters,
      minPts = hdb_result$minPts,
      coredist = hdb_result$coredist,
      cluster_scores = hdb_result$cluster_scores,
      membership_prob = hdb_result$membership_prob,
      outlier_scores = hdb_result$outlier_scores,
      hc = hdb_result$hc,
      forced_clustering = TRUE
    )
    
    cat(sprintf("Successfully created %d clusters using cutree\n", length(unique(forced_clusters))))
  } else {
    hdb_clusters <- hdb_result
    hdb_clusters$forced_clustering <- FALSE
  }
  
  # Handle noise points in sample
  if(any(hdb_clusters$cluster == 0)) {
    noise_idx <- which(hdb_clusters$cluster == 0)
    cat(sprintf("Handling %d noise points in sample...\n", length(noise_idx)))
    
    # Extract only needed rows to avoid converting entire matrix
    dist_matrix <- as.matrix(gower_dist_sample)[noise_idx, , drop = FALSE]
    non_noise_idx <- which(hdb_clusters$cluster != 0)
    
    for(j in seq_along(noise_idx)) {
      i <- noise_idx[j]
      if(length(non_noise_idx) > 0) {
        nearest_idx <- non_noise_idx[which.min(dist_matrix[j, non_noise_idx])]
        hdb_clusters$cluster[i] <- hdb_clusters$cluster[nearest_idx]
      }
    }
  }
  
  n_clusters <- length(unique(hdb_clusters$cluster[hdb_clusters$cluster > 0]))
  cat(sprintf("Found %d clusters\n", n_clusters))
  
  # Step 4: Calculate cluster representatives (medoids) for assignment
  cat("Calculating cluster medoids...\n")
  
  cluster_medoids <- list()
  cluster_features <- list()
  
  for(clust in unique(hdb_clusters$cluster)) {
    if(clust == 0) next  # Skip if any noise remains
    
    cluster_mask <- hdb_clusters$cluster == clust
    cluster_indices <- which(cluster_mask)
    
    if(length(cluster_indices) > 1) {
      # Find medoid (most central point) of cluster
      cluster_dist_subset <- as.matrix(gower_dist_sample)[cluster_indices, cluster_indices]
      medoid_idx <- cluster_indices[which.min(rowSums(cluster_dist_subset))]
      cluster_medoids[[as.character(clust)]] <- medoid_idx
      cluster_features[[as.character(clust)]] <- sample_data[medoid_idx, feature_cols]
    } else {
      cluster_medoids[[as.character(clust)]] <- cluster_indices[1]
      cluster_features[[as.character(clust)]] <- sample_data[cluster_indices[1], feature_cols]
    }
  }
  
  # Create reference data for cluster assignment
  medoid_features <- do.call(rbind, cluster_features)
  
  # Step 5: Assign ALL observations to clusters FIRST
  cat("Assigning all observations to clusters...\n")
  
  all_clusters <- integer(nrow(data))
  
  # First, assign the sampled observations
  all_clusters[sampled_idx] <- hdb_clusters$cluster
  
  # Now assign the remaining observations
  unsampled_idx <- setdiff(1:nrow(data), sampled_idx)
  
  if(length(unsampled_idx) > 0) {
    cat(sprintf("Assigning %d unsampled observations to clusters...\n", length(unsampled_idx)))
    
    # Process in batches to manage memory
    batch_size <- 10000
    n_batches <- ceiling(length(unsampled_idx) / batch_size)
    
    # Initialize progress bar
    pb <- txtProgressBar(min = 0, max = n_batches, style = 3)
    
    for(batch in 1:n_batches) {
      # Update progress bar
      setTxtProgressBar(pb, batch)
      
      batch_start <- (batch - 1) * batch_size + 1
      batch_end <- min(batch * batch_size, length(unsampled_idx))
      batch_idx <- unsampled_idx[batch_start:batch_end]
      
      # Calculate distances to medoids
      batch_data <- data[batch_idx, feature_cols]
      
      # Calculate distances using Gower distance
      cluster_assignments_batch <- assign_to_nearest_cluster(
        new_data = batch_data,
        reference_data = medoid_features,
        weights = weights
      )
      
      all_clusters[batch_idx] <- as.integer(names(cluster_features))[cluster_assignments_batch]
    }
    
    # Close progress bar
    close(pb)
  }
  
  # Step 6: Calculate cluster statistics using FULL dataset
  cat("Calculating cluster statistics on full dataset...\n")
  
  cluster_stats <- data.frame(
    cluster = unique(all_clusters[all_clusters > 0]),
    n_total = 0,
    n_presence = 0,
    prevalence = 0
  )
  
  for(i in 1:nrow(cluster_stats)) {
    clust <- cluster_stats$cluster[i]
    cluster_mask <- all_clusters == clust
    cluster_stats$n_total[i] <- sum(cluster_mask)
    cluster_stats$n_presence[i] <- sum(data[cluster_mask, target_col] == 1)
    cluster_stats$prevalence[i] <- cluster_stats$n_presence[i] / cluster_stats$n_total[i]
  }
  
  # Step 7: Assign clusters to train/calib/test splits based on FULL dataset prevalence
  cat("Assigning clusters to train/calib/test splits...\n")
  
  # Sort and assign clusters to splits
  cluster_stats <- cluster_stats[order(cluster_stats$prevalence, decreasing = TRUE), ]
  
  cluster_assignments <- character(nrow(cluster_stats))
  names(cluster_assignments) <- as.character(cluster_stats$cluster)
  
  total_presence <- sum(cluster_stats$n_presence)
  target_train_presence <- round(total_presence * train_prop)
  target_calib_presence <- round(total_presence * calib_prop)
  target_test_presence <- round(total_presence * test_prop)
  
  train_presence <- calib_presence <- test_presence <- 0
  
  for(i in 1:nrow(cluster_stats)) {
    cluster_id <- as.character(cluster_stats$cluster[i])
    cluster_presence <- cluster_stats$n_presence[i]
    
    train_need <- (target_train_presence - train_presence) / max(target_train_presence, 1)
    calib_need <- (target_calib_presence - calib_presence) / max(target_calib_presence, 1)
    test_need <- (target_test_presence - test_presence) / max(target_test_presence, 1)
    
    if(train_need >= calib_need && train_need >= test_need) {
      cluster_assignments[cluster_id] <- "train"
      train_presence <- train_presence + cluster_presence
    } else if(calib_need >= test_need) {
      cluster_assignments[cluster_id] <- "calib"
      calib_presence <- calib_presence + cluster_presence
    } else {
      cluster_assignments[cluster_id] <- "test"
      test_presence <- test_presence + cluster_presence
    }
  }
  
  # Step 8: Create final split assignments
  all_splits <- character(nrow(data))
  
  for(i in 1:nrow(data)) {
    clust <- as.character(all_clusters[i])
    if(clust %in% names(cluster_assignments)) {
      all_splits[i] <- cluster_assignments[clust]
    }
  }
  
  # Step 9: Create results
  train_idx <- which(all_splits == "train")
  calib_idx <- which(all_splits == "calib")
  test_idx <- which(all_splits == "test")
  
  results <- list(
    splits = all_splits,
    clusters = all_clusters,
    cluster_stats = cluster_stats,
    cluster_assignments = cluster_assignments,
    sampled_indices = sampled_idx,
    split_summary = data.frame(
      split = c("train", "calib", "test"),
      n_total = c(length(train_idx), length(calib_idx), length(test_idx)),
      n_presence = c(
        sum(data[train_idx, target_col] == 1),
        sum(data[calib_idx, target_col] == 1),
        sum(data[test_idx, target_col] == 1)
      ),
      prevalence = c(
        sum(data[train_idx, target_col] == 1) / max(length(train_idx), 1),
        sum(data[calib_idx, target_col] == 1) / max(length(calib_idx), 1),
        sum(data[test_idx, target_col] == 1) / max(length(test_idx), 1)
      )
    )
  )
  
  return(results)
}

# Helper function for mixed-type distance calculation using gower distance
assign_to_nearest_cluster <- function(new_data, reference_data, weights = NULL) {
  
  # Calculate Gower distances between new_data and reference_data
  # Using daisy requires combining datasets and extracting relevant distances
  combined_data <- rbind(new_data, reference_data)
  
  distance_obj <- cluster::daisy(
    x = combined_data,
    metric = "gower", 
    weights = weights
  )
  
  # Extract the relevant submatrix
  dist_mat <- as.matrix(distance_obj)
  n_new <- nrow(new_data)
  n_ref <- nrow(reference_data)
  
  # Extract distances from new_data (rows 1:n_new) to reference_data (rows (n_new+1):(n_new+n_ref))
  relevant_distances <- dist_mat[1:n_new, (n_new + 1):(n_new + n_ref), drop = FALSE]
  
  # Find nearest cluster for each observation
  assignments <- apply(relevant_distances, 1, which.min)
  
  return(assignments)
}

## Apply function to create splits ####

# Get feature column names (exclude response variable)
feature_cols <- setdiff(names(training_data), "response")

cat("Feature columns:", paste(feature_cols, collapse = ", "), "\n")

cat("Using feature weights for Gower distance:\n")
print(weights)

set.seed(42)
training_data_N1e5 <- training_data |> 
  sample_n(100000)

# Create feature space splits
splits <- create_feature_space_splits_large(
  data = training_data_N1e5,
  feature_cols = feature_cols,
  target_col = "response",
  sample_size = 10000,  # Cluster on X observations
  sample_ratio_minority = 0.01,  # Aim for X of sample to be positive cases
  min_cluster_size = 20,
  cluster_selection_epsilon = 0,
  minimum_clusters = 15,
  train_prop = 0.6,
  calib_prop = 0.2,
  test_prop = 0.2,
  weights = weights
)

# Display results
print(splits$split_summary)
cat("\nCluster statistics:\n")
print(head(splits$cluster_stats))

## Save splits for downstream modeling ####

# Create split datasets
train_data <- training_data[splits$splits == "train", ]
calib_data <- training_data[splits$splits == "calib", ]
test_data <- training_data[splits$splits == "test", ]

# Save to CSV files
write_csv(train_data, "output/train_data_partitioned.csv")
write_csv(calib_data, "output/calib_data_partitioned.csv")
write_csv(test_data, "output/test_data_partitioned.csv")

# Save split assignments and cluster information
split_info <- data.frame(
  row_id = 1:nrow(training_data),
  split = splits$splits,
  cluster = splits$clusters
)
write_csv(split_info, "output/split_assignments.csv")

cat("\nSaved partitioned datasets to output/ directory\n")