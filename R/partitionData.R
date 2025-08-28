# PAM-based Feature Space Clustering for Data Partitioning ####

# This script implements PAM clustering methodology for data partitioning
# based on the approach from simplePAMResults.R:
# 1) Sample 1e4 rows of training data
# 2) PAM clustering into K clusters with standardized euclidean distances of top X features  
# 3) Assignment of remaining data to nearest cluster medoids

library(readr)
library(dplyr)
library(terra)
library(cluster)
library(RANN)
library(tidyr)
library(ggplot2)
library(sf)

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

# Combine into training dataset (keep x,y coordinates for later use)
training_data_with_coords <- bind_rows(presence_df, absence_df) |>
  mutate(response = factor(response), ar50 = factor(ar50))

training_data <- training_data_with_coords |> 
  select(-x, -y) # Remove x and y coordinates for modeling

# Get final counts
n_final_presence <- sum(training_data$response == "1")
n_final_background <- sum(training_data$response == "0")
n_final_prevalence <- n_final_presence / (n_final_presence + n_final_background)

cat("Final training data - Presences:", n_final_presence, "Background:", n_final_background, "\n")
cat("Prevalence in training data:", round(n_final_prevalence * 100, 3), "%\n")

# Load feature weights
weights_gini <- read_csv("output/feature_weights_gini.csv")
weights <- pull(weights_gini, weights_gini)
names(weights) <- pull(weights_gini, feature)

## PAM-based data partitioning function ####

create_pam_splits <- function(data, 
                              feature_cols,
                              target_col = "response",
                              sample_size = 10000,
                              k_clusters = 10,
                              top_n_features = 5,
                              train_prop = 0.6,
                              calib_prop = 0.2,
                              test_prop = 0.2,
                              feature_weights = NULL,
                              seed = 42) {
  
  set.seed(seed)
  
  # Step 1: Create sample for clustering
  cat("Creating sample for clustering...\n")
  sample_data <- data |> 
    sample_n(sample_size)
  
  sampled_idx <- as.numeric(rownames(sample_data))
  
  cat("Sample prevalence:", 
      round(sum(sample_data[[target_col]] == "1") / nrow(sample_data) * 100, 2), "%\n")
  
  # Step 2: Select top features by weight
  if(!is.null(feature_weights)) {
    top_features <- feature_weights |> 
      sort(decreasing = TRUE) |> 
      head(top_n_features) |>
      names()
    
    cat("Using top", top_n_features, "features:", paste(top_features, collapse = ", "), "\n")
  } else {
    top_features <- feature_cols[1:min(top_n_features, length(feature_cols))]
    cat("No feature weights provided. Using first", length(top_features), "features\n")
  }
  
  # Step 3: Extract feature data for clustering (numeric only)
  clustering_data <- sample_data |>
    select(all_of(top_features)) |>
    select(where(is.numeric))
  
  cat("Clustering data dimensions:", nrow(clustering_data), "x", ncol(clustering_data), "\n")
  
  # Step 4: Run PAM clustering with standardized Euclidean distance
  set.seed(seed)
  pam_result <- pam(
    x = clustering_data,
    k = k_clusters,
    metric = "euclidean", 
    stand = TRUE,  # Standardize features
    pamonce = 5
  )
  
  # Add cluster assignments to training sample
  sample_data$cluster <- factor(pam_result$clustering)
  
  # Show cluster summary
  cluster_summary <- sample_data |>
    group_by(cluster) |>
    summarise(
      n_points = n(),
      n_presence = sum(as.numeric(as.character(.data[[target_col]]))),
      prevalence = n_presence / n_points,
      .groups = 'drop'
    )
  
  cat("\nCluster summary:\n")
  print(cluster_summary)
  
  # Get medoids from PAM result
  medoid_features <- pam_result$medoids |> 
    as_tibble()
  
  # Step 5: Assign ALL observations to clusters using PAM assignment logic
  cat("Assigning all observations to clusters...\n")
  
  all_clusters <- integer(nrow(data))
  
  # First, assign the sampled observations
  all_clusters[sampled_idx] <- pam_result$clustering
  
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
      
      # Extract the same features for batch data
      batch_data_features <- data[batch_idx, ] |>
        select(all_of(top_features)) |>
        select(where(is.numeric))
      
      # Use PAM-style assignment function
      cluster_assignments_batch <- assign_to_pam_clusters(batch_data_features, medoid_features)
      all_clusters[batch_idx] <- cluster_assignments_batch
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

# Fast PAM-based assignment function using RANN for nearest neighbor search
assign_to_pam_clusters <- function(new_data, medoids) {
  
  # Use RANN's nn2() for fast nearest neighbor search
  # Find the 1 nearest neighbor (closest medoid) for each point
  nearest_neighbors <- nn2(
    data = medoids,      # Reference points (medoids)
    query = new_data,    # Query points (new data)
    k = 1,               # Find 1 nearest neighbor
    searchtype = "standard"  # Standard ANN search
  )
  
  # Extract the indices of the nearest medoids
  assignments <- nearest_neighbors$nn.idx[, 1]
  
  return(assignments)
}

## Apply PAM clustering to create splits ####

# Get feature column names (exclude response variable)
feature_cols <- setdiff(names(training_data), "response")

cat("Feature columns:", paste(feature_cols, collapse = ", "), "\n")

cat("Using feature weights for feature selection:\n")
print(weights)

# Create PAM-based feature space splits
splits <- create_pam_splits(
  data = training_data,
  feature_cols = feature_cols,
  target_col = "response",
  sample_size = 10000,  # Sample for clustering
  k_clusters = 10,      # Number of PAM clusters
  top_n_features = 5,   # Top features by weight
  train_prop = 0.6,
  calib_prop = 0.2,
  test_prop = 0.2,
  feature_weights = weights
)

# Display results
print(splits$split_summary)
cat("\nCluster statistics:\n")
print(head(splits$cluster_stats))

## Map visualization ####
presence_idx <- which(training_data_with_coords$response == "1")
background_idx <- sample(which(training_data_with_coords$response == "0"), 100000 - length(presence_idx))
viz_idx <- c(presence_idx, background_idx)

viz_data <- training_data_with_coords[viz_idx, ] |>
  mutate(cluster = factor(splits$clusters[viz_idx]),
         partition = factor(splits$splits[viz_idx])) |>
  arrange(response) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train))

ggplot(viz_data) +
  geom_sf(aes(color = cluster), size = 0.5, alpha = 0.7)

ggplot(viz_data) +
  geom_sf(aes(color = partition), size = 0.5, alpha = 0.7)

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