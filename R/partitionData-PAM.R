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
preds_nor_250m_cur <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")
preds_train <- c(rf_global_current, preds_nor_250m_cur)

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
  mutate(response = factor(response))

training_data <- training_data_with_coords |> 
  select(-x, -y) # Remove x and y coordinates for modeling

# Get final counts
n_final_presence <- sum(training_data$response == "1")
n_final_absence <- sum(training_data$response == "0")
n_final_prevalence <- n_final_presence / (n_final_presence + n_final_absence)

cat("Final training data - Presences:", n_final_presence, "Absence:", n_final_absence, "\n")
cat("Prevalence in training data:", round(n_final_prevalence * 100, 3), "%\n")

# Load feature weights
weights_features <- read_csv("output/weights_feature_data_partitioning.csv")
weights_features <- weights_features |> 
  mutate(feature = forcats::fct_reorder(feature, median_normalized, .fun = first))

# Create plot with error bars
ggplot(weights_features, aes(x = feature, y = median_normalized, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  geom_errorbar(aes(ymin = pmax(0, median_normalized - sd_normalized), 
                    ymax = median_normalized + sd_normalized),
                position = position_dodge(width = 0.9), width = 0.2) +
  coord_flip() +
  labs(title = "Feature Importance from Different Random Forest Methods",
       x = "Feature",
       y = "Normalized Importance",
       fill = "Method") +
  theme_minimal() +
  theme(legend.position = "bottom")

weights <- weights_features |> 
  filter(grepl("Balanced", method)) |> 
  select(feature, median_normalized) |> 
  mutate(
    # Set weights of artype* features to zero to exclude them from clustering
    weight_partitioning = case_when(
    grepl("artype", feature) ~ 0,
    TRUE ~ median_normalized
  )) |> 
  select(feature, weight_partitioning) |>
  arrange(desc(weight_partitioning)) |>
  tibble::deframe()

## Core reusable functions ####

# Select features for clustering based on weights
select_clustering_features <- function(feature_cols, top_n_features, feature_weights = NULL) {
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

  return(top_features)
}

# Run PAM clustering with standardization and return cluster object
run_pam_clustering <- function(sample_data, top_features, k_clusters, target_col, 
                               seed = 42) {

  # Extract feature data for clustering (numeric only)
  clustering_data <- sample_data |>
    select(all_of(top_features)) |>
    select(where(is.numeric))

  cat("Clustering data dimensions:", nrow(clustering_data), "x", ncol(clustering_data), "\n")

  # Run PAM clustering with standardized Euclidean distance
  set.seed(seed)
  pam_result <- pam(
    x = clustering_data,
    k = k_clusters,
    metric = "euclidean",
    stand = TRUE,  # Standardize features
    pamonce = 5
  )

  # Add cluster assignments to sample data
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

  # Calculate standardization parameters for consistent assignment
  standardization_params <- list(
    means = apply(clustering_data, 2, mean),
    sds = apply(clustering_data, 2, sd)
  )

  return(list(
    pam_result = pam_result,
    sample_data = sample_data,
    clustering_data = clustering_data,
    standardization_params = standardization_params
  ))
}

# Fast PAM-based assignment function using standardized distances
assign_to_pam_clusters <- function(new_data, medoids, data_means, data_sds) {

  # Standardize new_data using the same means and SDs as the original clustering
  new_data_std <- scale(new_data, center = data_means, scale = data_sds)
  medoids_std <- scale(medoids, center = data_means, scale = data_sds)

  # Use RANN's nn2() for fast nearest neighbor search on standardized data
  # Find the 1 nearest neighbor (closest medoid) for each point
  nearest_neighbors <- nn2(
    data = medoids_std,      # Reference points (standardized medoids)
    query = new_data_std,    # Query points (standardized new data)
    k = 1,                   # Find 1 nearest neighbor
    searchtype = "standard",  # Standard ANN search
    eps = 0
  )

  # Extract the indices of the nearest medoids
  assignments <- nearest_neighbors$nn.idx[, 1]

  return(assignments)
}

# Assign all observations to clusters using consistent method
assign_all_to_clusters <- function(data, sampled_idx, clustering_result,
                                   top_features, debug = TRUE) {

  cat("Assigning all observations to clusters using consistent nearest-neighbor logic...\n")

  # Get medoids and standardization parameters
  medoid_features <- clustering_result$pam_result$medoids |> as_tibble()
  std_params <- clustering_result$standardization_params

  if(debug) {
    cat("DIAGNOSTIC: Medoids for assignment:\n")
    print(head(medoid_features, 10))
    cat("DIAGNOSTIC: Standardization params - means:", round(std_params$means[1:6], 3), "...\n")
    cat("DIAGNOSTIC: Standardization params - sds:", round(std_params$sds[1:6], 3), "...\n")
    cat("DIAGNOSTIC: Assigning all", nrow(data), "observations (including originally sampled ones)\n")
  }

  # Process ALL observations in batches using assign_to_pam_clusters
  all_clusters <- integer(nrow(data))
  batch_size <- 10000
  n_batches <- ceiling(nrow(data) / batch_size)

  # Track assignments for debugging
  if(debug && n_batches > 0) {
    # Test first small batch to check assignment logic
    test_batch_size <- min(100, nrow(data))
    test_batch_idx <- 1:test_batch_size

    test_data_features <- data[test_batch_idx, ] |>
      select(all_of(top_features)) |>
      select(where(is.numeric))

    cat("DIAGNOSTIC: Test batch features shape:\n", nrow(test_data_features), "x", ncol(test_data_features), "\n")
    cat("DIAGNOSTIC: Test batch feature sample:\n")
    print(head(test_data_features, 10))

    test_assignments <- assign_to_pam_clusters(
      test_data_features,
      medoid_features,
      std_params$means,
      std_params$sds
    )

    cat("DIAGNOSTIC: Test batch cluster assignments:", head(test_assignments, 10), "\n")
    cat("DIAGNOSTIC: Test batch assignment distribution:\n", table(test_assignments), "\n")
  }

  for(batch in 1:n_batches) {

    batch_start <- (batch - 1) * batch_size + 1
    batch_end <- min(batch * batch_size, nrow(data))
    batch_idx <- batch_start:batch_end

    # Extract the same features for batch data
    batch_data_features <- data[batch_idx, ] |>
      select(all_of(top_features)) |>
      select(where(is.numeric))

    # Use PAM assignment for ALL observations
    cluster_assignments_batch <- assign_to_pam_clusters(
      batch_data_features,
      medoid_features,
      std_params$means,
      std_params$sds
    )
    all_clusters[batch_idx] <- cluster_assignments_batch
  }

  if(debug) {
    # Compare with original PAM assignments for verification
    original_assignments <- clustering_result$pam_result$clustering
    assigned_from_nn <- all_clusters[sampled_idx]

    cat("DIAGNOSTIC: Comparing original PAM vs nearest-neighbor for sampled observations:\n")
    cat("DIAGNOSTIC: Original PAM assignments:", head(original_assignments, 10), "\n")
    cat("DIAGNOSTIC: Nearest-neighbor assignments:", head(assigned_from_nn, 10), "\n")
    cat("DIAGNOSTIC: Perfect match:", all(original_assignments == assigned_from_nn), "\n")
    if(!all(original_assignments == assigned_from_nn)) {
      mismatches <- sum(original_assignments != assigned_from_nn)
      cat("DIAGNOSTIC: Number of mismatches:", mismatches, "out of", length(sampled_idx), "\n")
    }
  }

  return(all_clusters)
}

# Calculate cluster statistics on full dataset
calculate_cluster_stats <- function(data, all_clusters, target_col) {

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

  return(cluster_stats)
}

# Assign clusters to splits using dual-constraint optimization
assign_clusters_to_splits <- function(cluster_stats, train_prop, calib_prop, test_prop,
                                      tolerance = 0.05, debug = TRUE) {

  cat("Assigning clusters to train/calib/test splits with dual constraints...\n")

  # Sort clusters by prevalence (highest first) to prioritize presence balance
  cluster_stats <- cluster_stats[order(cluster_stats$prevalence, decreasing = TRUE), ]

  # Initialize tracking variables
  cluster_assignments <- character(nrow(cluster_stats))
  names(cluster_assignments) <- as.character(cluster_stats$cluster)

  # Calculate targets for both presence and total observations
  total_presence <- sum(cluster_stats$n_presence)
  total_observations <- sum(cluster_stats$n_total)

  target_train_presence <- round(total_presence * train_prop)
  target_calib_presence <- round(total_presence * calib_prop)
  target_test_presence <- round(total_presence * test_prop)

  target_train_total <- round(total_observations * train_prop)
  target_calib_total <- round(total_observations * calib_prop)
  target_test_total <- round(total_observations * test_prop)

  # Track current counts
  train_presence <- train_total <- 0
  calib_presence <- calib_total <- 0
  test_presence <- test_total <- 0

  if(debug) {
    cat("DUAL-CONSTRAINT TARGETS:\n")
    cat("Total observations:", total_observations, "| Total presences:", total_presence, "\n")
    cat("Train targets: ", target_train_total, "obs (", train_prop*100, "%), ", target_train_presence, "pres\n")
    cat("Calib targets: ", target_calib_total, "obs (", calib_prop*100, "%), ", target_calib_presence, "pres\n")
    cat("Test targets:  ", target_test_total, "obs (", test_prop*100, "%), ", target_test_presence, "pres\n")
  }

  # Assign clusters using dual-constraint priority scoring
  for(i in 1:nrow(cluster_stats)) {
    cluster_id <- as.character(cluster_stats$cluster[i])
    cluster_presence <- cluster_stats$n_presence[i]
    cluster_total <- cluster_stats$n_total[i]

    # Check capacity constraints (hard stops)
    train_capacity <- (train_total + cluster_total) <= (target_train_total * (1 + tolerance))
    calib_capacity <- (calib_total + cluster_total) <= (target_calib_total * (1 + tolerance))
    test_capacity <- (test_total + cluster_total) <= (target_test_total * (1 + tolerance))

    # Calculate presence need scores (0 = satisfied, 1 = maximum need)
    train_pres_need <- max(0, (target_train_presence - train_presence) / max(target_train_presence, 1))
    calib_pres_need <- max(0, (target_calib_presence - calib_presence) / max(target_calib_presence, 1))
    test_pres_need <- max(0, (target_test_presence - test_presence) / max(target_test_presence, 1))

    # Calculate capacity availability scores (0 = no capacity, 1 = full capacity)
    train_cap_score <- ifelse(train_capacity,
                              1 - (train_total / max(target_train_total, 1)), 0)
    calib_cap_score <- ifelse(calib_capacity,
                              1 - (calib_total / max(target_calib_total, 1)), 0)
    test_cap_score <- ifelse(test_capacity,
                             1 - (test_total / max(target_test_total, 1)), 0)

    # Calculate priority scores (presence need * capacity availability)
    train_priority <- train_pres_need * train_cap_score
    calib_priority <- calib_pres_need * calib_cap_score
    test_priority <- test_pres_need * test_cap_score

    if(debug && i <= 5) {
      cat("Cluster", cluster_id, "(", cluster_total, "obs,", cluster_presence, "pres):\n")
      cat("  Capacities: train=", train_capacity, ", calib=", calib_capacity, ", test=", test_capacity, "\n")
      cat("  Priorities: train=", round(train_priority, 3), ", calib=", round(calib_priority, 3), ", test=", round(test_priority, 3), "\n")
    }

    # Assign to split with highest priority (that has capacity)
    if(train_priority >= calib_priority && train_priority >= test_priority && train_capacity) {
      cluster_assignments[cluster_id] <- "train"
      train_presence <- train_presence + cluster_presence
      train_total <- train_total + cluster_total
      if(debug && i <= 5) cat("  -> Assigned to TRAIN\n")
    } else if(calib_priority >= test_priority && calib_capacity) {
      cluster_assignments[cluster_id] <- "calib"
      calib_presence <- calib_presence + cluster_presence
      calib_total <- calib_total + cluster_total
      if(debug && i <= 5) cat("  -> Assigned to CALIB\n")
    } else if(test_capacity) {
      cluster_assignments[cluster_id] <- "test"
      test_presence <- test_presence + cluster_presence
      test_total <- test_total + cluster_total
      if(debug && i <= 5) cat("  -> Assigned to TEST\n")
    } else {
      # Fallback: assign to split with smallest overshoot
      train_overshoot <- max(0, (train_total + cluster_total) - target_train_total)
      calib_overshoot <- max(0, (calib_total + cluster_total) - target_calib_total)
      test_overshoot <- max(0, (test_total + cluster_total) - target_test_total)

      min_overshoot <- min(train_overshoot, calib_overshoot, test_overshoot)

      if(train_overshoot == min_overshoot) {
        cluster_assignments[cluster_id] <- "train"
        train_presence <- train_presence + cluster_presence
        train_total <- train_total + cluster_total
      } else if(calib_overshoot == min_overshoot) {
        cluster_assignments[cluster_id] <- "calib"
        calib_presence <- calib_presence + cluster_presence
        calib_total <- calib_total + cluster_total
      } else {
        cluster_assignments[cluster_id] <- "test"
        test_presence <- test_presence + cluster_presence
        test_total <- test_total + cluster_total
      }
      if(debug && i <= 5) cat("  -> FALLBACK assignment (all splits near capacity)\n")
    }
  }

  if(debug) {
    cat("\nFINAL ASSIGNMENT RESULTS:\n")
    cat("Train: ", train_total, "/", target_train_total, " obs (", round(train_total/total_observations*100, 1), "%), ",
        train_presence, "/", target_train_presence, " pres\n")
    cat("Calib: ", calib_total, "/", target_calib_total, " obs (", round(calib_total/total_observations*100, 1), "%), ",
        calib_presence, "/", target_calib_presence, " pres\n")
    cat("Test:  ", test_total, "/", target_test_total, " obs (", round(test_total/total_observations*100, 1), "%), ",
        test_presence, "/", target_test_presence, " pres\n")
  }

  return(cluster_assignments)
}

# Assign clusters to CV folds based on prevalence balancing
assign_clusters_to_folds <- function(cluster_stats, n_folds) {

  cat("Assigning clusters to", n_folds, "folds...\n")

  # Sort clusters by prevalence (highest first)
  cluster_stats <- cluster_stats[order(cluster_stats$prevalence, decreasing = TRUE), ]

  fold_assignments <- integer(nrow(cluster_stats))
  names(fold_assignments) <- as.character(cluster_stats$cluster)

  total_presence <- sum(cluster_stats$n_presence)
  fold_presence <- rep(0, n_folds)

  for(i in 1:nrow(cluster_stats)) {
    cluster_id <- as.character(cluster_stats$cluster[i])
    cluster_presence <- cluster_stats$n_presence[i]

    # Find fold with smallest presence count
    best_fold <- which.min(fold_presence)
    fold_assignments[cluster_id] <- best_fold
    fold_presence[best_fold] <- fold_presence[best_fold] + cluster_presence
  }

  return(fold_assignments)
}

## High-level orchestration functions ####

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
                              seed = 42,
                              debug = TRUE) {

  set.seed(seed)

  # Step 1: Create sample for clustering
  cat("Creating sample for clustering...\n")
  sampled_idx <- sample.int(nrow(data), sample_size)
  sample_data <- data |> slice(sampled_idx)

  cat("Sample prevalence:",
      round(sum(sample_data[[target_col]] == "1") / nrow(sample_data) * 100, 2), "%\n")

  if(debug) {
    cat("DIAGNOSTIC: Sample indices range:", min(sampled_idx), "to", max(sampled_idx), "\n")
    presence_in_sample <- sum(sample_data[[target_col]] == "1")
    cat("DIAGNOSTIC: Presences in sample:", presence_in_sample, "out of", nrow(sample_data), "\n")
  }

  # Step 2: Select features and run clustering
  top_features <- select_clustering_features(feature_cols, top_n_features, feature_weights)
  clustering_result <- run_pam_clustering(sample_data, top_features, k_clusters, target_col, seed)

  if(debug) {
    cat("DIAGNOSTIC: Clustering completed. Medoids shape:", nrow(clustering_result$pam_result$medoids), "x", ncol(clustering_result$pam_result$medoids), "\n")
    cat("DIAGNOSTIC: Standardization params - means:", round(clustering_result$standardization_params$means[1:6], 3), "...\n")
    cat("DIAGNOSTIC: Standardization params - sds:", round(clustering_result$standardization_params$sds[1:6], 3), "...\n")
  }

  # Step 3: Assign all observations to clusters
  all_clusters <- assign_all_to_clusters(data, sampled_idx, clustering_result, top_features, debug)

  if(debug) {
    # Check cluster assignments for sampled vs unsampled observations
    sampled_clusters <- all_clusters[sampled_idx]
    unsampled_idx <- setdiff(1:nrow(data), sampled_idx)
    unsampled_clusters <- all_clusters[unsampled_idx]

    cat("DIAGNOSTIC: Sampled observations cluster distribution:\n")
    print(table(sampled_clusters))
    cat("DIAGNOSTIC: Unsampled observations cluster distribution:\n")
    print(table(unsampled_clusters))

    # Check presence distribution across clusters for sampled vs unsampled
    sampled_presence_by_cluster <- tapply(data[sampled_idx, target_col], sampled_clusters, function(x) sum(x == 1))
    unsampled_presence_by_cluster <- tapply(data[unsampled_idx, target_col], unsampled_clusters, function(x) sum(x == 1))

    cat("DIAGNOSTIC: Sampled presences by cluster:\n")
    print(sampled_presence_by_cluster)
    cat("DIAGNOSTIC: Unsampled presences by cluster:\n")
    print(unsampled_presence_by_cluster)
  }

  # Step 4: Calculate cluster statistics and assign to splits
  cluster_stats <- calculate_cluster_stats(data, all_clusters, target_col)
  cluster_assignments <- assign_clusters_to_splits(cluster_stats, train_prop, calib_prop, test_prop)

  # Step 5: Create final split assignments
  all_splits <- character(nrow(data))
  for(i in 1:nrow(data)) {
    clust <- as.character(all_clusters[i])
    if(clust %in% names(cluster_assignments)) {
      all_splits[i] <- cluster_assignments[clust]
    }
  }

  # Step 6: Create results summary
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

create_inner_cv_folds <- function(train_data,
                                  feature_cols,
                                  target_col = "response",
                                  sample_size = 10000,
                                  k_clusters = 6,
                                  n_folds = 3,
                                  top_n_features = 5,
                                  feature_weights = NULL,
                                  seed = 42) {

  set.seed(seed)

  cat("Creating inner CV folds for training data...\n")
  cat("Training data size:", nrow(train_data), "observations\n")

  # Step 1: Create sample for clustering (from training data only)
  cat("Creating sample for CV clustering...\n")
  actual_sample_size <- min(sample_size, nrow(train_data))

  # Sample indices first, then slice
  sampled_idx <- sample.int(nrow(train_data), actual_sample_size)
  sample_data <- train_data |> slice(sampled_idx)

  # Create indexed version for batch processing
  train_data_with_idx <- train_data |>
    mutate(.train_row_id = row_number())

  cat("CV sample size:", actual_sample_size, "\n")
  cat("CV sample prevalence:",
      round(sum(sample_data[[target_col]] == "1") / nrow(sample_data) * 100, 2), "%\n")

  # Step 2: Select features and run clustering
  top_features <- select_clustering_features(feature_cols, top_n_features, feature_weights)
  clustering_result <- run_pam_clustering(sample_data, top_features, k_clusters, target_col, seed)

  # Step 3: Assign all training observations to clusters
  all_cv_clusters <- assign_all_to_clusters(train_data_with_idx, sampled_idx, clustering_result, top_features)

  # Step 4: Calculate cluster statistics and assign to folds
  cv_cluster_stats <- calculate_cluster_stats(train_data, all_cv_clusters, target_col)
  cv_fold_assignments <- assign_clusters_to_folds(cv_cluster_stats, n_folds)

  # Step 5: Create final CV fold assignments
  all_cv_folds <- integer(nrow(train_data))
  for(i in 1:nrow(train_data)) {
    clust <- as.character(all_cv_clusters[i])
    if(clust %in% names(cv_fold_assignments)) {
      all_cv_folds[i] <- cv_fold_assignments[clust]
    }
  }

  # Step 6: Create CV fold summary
  cv_fold_summary <- data.frame(
    cv_fold = 1:n_folds,
    n_total = 0,
    n_presence = 0,
    prevalence = 0
  )

  for(fold in 1:n_folds) {
    fold_idx <- which(all_cv_folds == fold)
    cv_fold_summary$n_total[fold] <- length(fold_idx)
    cv_fold_summary$n_presence[fold] <- sum(train_data[fold_idx, target_col] == 1)
    cv_fold_summary$prevalence[fold] <- cv_fold_summary$n_presence[fold] / max(cv_fold_summary$n_total[fold], 1)
  }

  results <- list(
    cv_folds = all_cv_folds,
    cv_clusters = all_cv_clusters,
    cv_cluster_stats = cv_cluster_stats,
    cv_fold_assignments = cv_fold_assignments,
    cv_fold_summary = cv_fold_summary,
    sampled_indices = sampled_idx
  )

  return(results)
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
  top_n_features = 6,   # Top features by weight
  train_prop = 0.6,
  calib_prop = 0.2,
  test_prop = 0.2,
  feature_weights = weights,
  debug = TRUE          # Enable diagnostic outputs
)

# Display results
print(splits$split_summary)
cat("\nCluster statistics:\n")
print(head(splits$cluster_stats, 10))

## Map visualization ####
training_data_with_coords_outer <- training_data_with_coords |>
  select(response, x, y) |> 
  mutate(cluster = factor(splits$clusters),
         partition = factor(splits$splits)) |> 
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |> 
  vect()
raster_outer <- rasterize(
  training_data_with_coords_outer,
  preds_train, 
  field = c("cluster", "partition"))
plot(raster_outer)

## Apply PAM clustering to create inner CV folds ####

# Get training data indices
train_idx <- which(splits$splits == "train")
train_data <- training_data[train_idx, ]

cat("Creating inner CV folds for training partition...\n")

# Create PAM-based inner CV folds
inner_cv <- create_inner_cv_folds(
  train_data = train_data,
  feature_cols = feature_cols,
  target_col = "response",
  sample_size = 10000,  # Same as outer partitioning
  k_clusters = 6,       # Proportional to 3 CV folds
  n_folds = 3,         # Number of CV folds
  top_n_features = 6,  # Same as outer partitioning
  feature_weights = weights,
  seed = 42
)

# Display results
print(inner_cv$cv_fold_summary)
cat("\nInner CV cluster statistics:\n")
print(head(inner_cv$cv_cluster_stats))

## Map visualization ####
train_data_with_coords_inner <- training_data_with_coords[train_idx, ] |>
  select(response, x, y) |> 
  mutate(cluster = factor(inner_cv$cv_clusters),
         fold = factor(inner_cv$cv_folds)) |> 
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |> 
  vect()
raster_inner <- rasterize(
  train_data_with_coords_inner,
  preds_train, 
  field = c("cluster", "fold"))
plot(raster_inner)

## Save partitioned datasets ####

### Main dataset for modeling ####
# Create main dataset with outer and inner partition columns
# Map inner CV folds back to original training_data indices
inner_folds_full <- rep(NA_integer_, nrow(training_data))
inner_folds_full[train_idx] <- inner_cv$cv_folds

training_data_partitioned <- training_data_with_coords |>
  mutate(
    outer = splits$splits,
    inner = inner_folds_full
  )

write_csv(training_data_partitioned, "output/data_partitioned.csv")

### Partition information in space ####
writeRaster(raster_outer, "output/partition_outer.tif", overwrite=TRUE)
writeRaster(raster_inner, "output/partition_inner.tif", overwrite=TRUE)

# sessionInfo ####

sessioninfo::session_info()
