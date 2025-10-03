
library(readr)
library(dplyr)
library(terra)
library(cluster)
library(tidyr)
library(ggplot2)
library(sf)
library(purrr)
library(twosamples)

## Functions ####

# Quantile cut-based partitioning with probabilistic interleaving/segregation
create_cut_based_partitions <- function(
    data,
    feature_cols,
    target_col,
    n_cuts,
    n_partitions = 5,
    segregation_prob = 0,
    top_n_features = 1,
    feature_weights = NULL,
    seed = 42) {

  cat("Creating", n_partitions, "partitions using", n_cuts, "quantile cuts",
      "(segregation_prob =", segregation_prob, ")...\n")

  # Select top feature(s)
  top_features <- select_clustering_features(
    feature_cols, top_n_features, feature_weights
  )

  # For now, use only the top feature (single variable)
  if (length(top_features) > 1) {
    cat("Note: Using only top feature:", top_features[1], "\n")
  }
  top_feature <- top_features[1]

  # Extract feature values
  feature_values <- data[[top_feature]]

  # Create quantile-based cuts
  quantile_breaks <- quantile(feature_values,
                              probs = seq(0, 1, length.out = n_cuts + 1),
                              na.rm = TRUE)

  # Handle case where quantiles are not unique (can happen with discrete data)
  quantile_breaks <- unique(quantile_breaks)
  actual_cuts <- length(quantile_breaks) - 1

  if (actual_cuts < n_cuts) {
    cat("Warning: Only", actual_cuts, "unique quantile breaks found (requested", n_cuts, ")\n")
  }

  intervals <- cut(feature_values,
                   breaks = quantile_breaks,
                   labels = FALSE,
                   include.lowest = TRUE)

  # Probabilistic interleaving/segregation assignment
  # segregation_prob = 0: Pure interleaving (all intervals distributed round-robin)
  # segregation_prob = 1: Maximum segregation (intervals form contiguous blocks)

  set.seed(seed)
  unique_intervals_ordered <- sort(unique(intervals))
  interval_to_partition <- integer(length(unique_intervals_ordered))

  for (i in seq_along(unique_intervals_ordered)) {
    if (runif(1) < segregation_prob) {
      # SEGREGATE: Join contiguous block (inherit previous partition)
      if (i == 1) {
        interval_to_partition[i] <- sample(1:n_partitions, 1)
      } else {
        interval_to_partition[i] <- interval_to_partition[i - 1]
      }
    } else {
      # INTERLEAVE: Round-robin assignment
      interval_to_partition[i] <- ((i - 1) %% n_partitions) + 1
    }
  }

  # Map intervals to their assigned partitions
  partition_assignments <- interval_to_partition[intervals]

  # Calculate partition statistics
  partition_stats <- data |>
    mutate(partition = factor(partition_assignments)) |>
    group_by(partition) |>
    summarise(
      n_total = n(),
      n_presence = sum(.data[[target_col]] == 1),
      prevalence = n_presence / n_total,
      .groups = "drop"
    )

  cat("\nPartition summary:\n")
  print(partition_stats)

  # Find test partition (closest to overall prevalence)
  overall_prevalence <- sum(partition_stats$n_presence) / sum(partition_stats$n_total)
  prevalence_diffs <- abs(partition_stats$prevalence - overall_prevalence)
  test_partition <- which.min(prevalence_diffs)

  cat("\nTest partition assignment:\n")
  for (p in 1:n_partitions) {
    cat(sprintf(
      "Partition %d: %d obs (%.1f%%), %.3f prevalence%s\n",
      p, partition_stats$n_total[p],
      partition_stats$n_total[p]/sum(partition_stats$n_total)*100,
      partition_stats$prevalence[p],
      ifelse(p == test_partition, " [TEST]", "")
    ))
  }

  # Create outer (test vs train) and inner (CV folds) partitions
  outer_partitions <- ifelse(partition_assignments == test_partition, "test", "train")

  # For inner partitions, renumber non-test partitions as CV folds
  inner_partitions <- rep(NA, length(partition_assignments))
  train_mask <- outer_partitions == "train"
  if (sum(train_mask) > 0) {
    train_partitions <- partition_assignments[train_mask]
    unique_train_partitions <- sort(unique(train_partitions))
    fold_mapping <- setNames(1:length(unique_train_partitions), unique_train_partitions)
    inner_partitions[train_mask] <- fold_mapping[as.character(train_partitions)]
  }

  return(list(
    outer_partitions = outer_partitions,
    inner_partitions = inner_partitions,
    clusters = intervals,  # Return interval assignments for compatibility
    partition_stats = partition_stats,
    test_partition_id = test_partition,
    silhouette = NA  # No silhouette for cut-based approach
  ))
}

# Run CLARA clustering with standardization and return cluster object
run_clara_clustering <- function(data, top_features, k_clusters, target_col,
                                 max_sample_size = 10e3, sample_fraction = 0.25,
                                 samples = 3, seed = 42) {
  # Extract feature data for clustering (numeric only)
  clustering_data <- data |>
    select(all_of(top_features)) |>
    select(where(is.numeric))

  cat("Clustering data dimensions:", nrow(clustering_data), "x", ncol(clustering_data), "\n")

  # Calculate CLARA sample size
  clara_sampsize <- min(max_sample_size, nrow(clustering_data) * sample_fraction)
  cat("CLARA sample size:", clara_sampsize, "\n")

  # Run CLARA clustering with standardized Euclidean distance
  set.seed(seed)
  clara_result <- clara(
    x = clustering_data,
    k = k_clusters,
    metric = "euclidean",
    stand = TRUE,
    rngR = TRUE,
    samples = samples,
    sampsize = clara_sampsize,
    trace = 1
  )

  return(clara_result)
}

# Single clustering approach: consolidate clusters into equal partitions
create_single_clustering_partitions <- function(
    data,
    feature_cols,
    target_col,
    n_clusters,
    n_partitions,
    top_n_features,
    feature_weights = NULL,
    max_sample_size = 10e3,
    sample_fraction = 0.25,
    samples = 3,
    seed = 42) {

  cat("Creating", n_clusters, "CLARA clusters for", n_partitions, "partitions...\n")

  # Select features and run clustering
  top_features <- select_clustering_features(
    feature_cols, top_n_features, feature_weights
  )
  clara_result <- run_clara_clustering(
    data, top_features, n_clusters, target_col, max_sample_size, sample_fraction, samples, seed
  )

  # Calculate cluster statistics
  cluster_stats <- data |>
    mutate(cluster = factor(clara_result$clustering)) |>
    group_by(cluster) |>
    summarise(
      n_total = n(),
      n_presence = sum(.data[[target_col]] == 1),
      prevalence = n_presence / n_total,
      .groups = "drop"
    )

  cat("\nCluster summary:\n")
  print(cluster_stats)

  # Consolidate clusters into equal partitions
  total_obs <- sum(cluster_stats$n_total)
  target_size <- total_obs / n_partitions
  overall_prevalence <- sum(cluster_stats$n_presence) / total_obs

  # Greedy assignment to minimize size deviation
  cluster_assignments <- rep(NA, nrow(cluster_stats))
  partition_sizes <- rep(0, n_partitions)
  partition_presence <- rep(0, n_partitions)

  # Sort clusters by size (largest first) for better packing
  cluster_order <- order(cluster_stats$n_total, decreasing = TRUE)

  for (i in cluster_order) {
    # Find partition with most remaining capacity
    best_partition <- which.min(partition_sizes)
    cluster_assignments[i] <- best_partition
    partition_sizes[best_partition] <- partition_sizes[best_partition] + cluster_stats$n_total[i]
    partition_presence[best_partition] <- partition_presence[best_partition] + cluster_stats$n_presence[i]
  }

  # Calculate partition prevalences and find best test partition
  partition_prevalences <- partition_presence / partition_sizes
  prevalence_diffs <- abs(partition_prevalences - overall_prevalence)
  test_partition <- which.min(prevalence_diffs)

  cat("\nPartition consolidation results:\n")
  for (p in 1:n_partitions) {
    cat(sprintf(
      "Partition %d: %d obs (%.1f%%), %.3f prevalence%s\n",
      p, partition_sizes[p], partition_sizes[p]/total_obs*100,
      partition_prevalences[p],
      ifelse(p == test_partition, " [TEST]", "")
    ))
  }

  # Create final assignments
  all_clusters <- clara_result$clustering
  partition_assignments <- cluster_assignments[all_clusters]

  # Create outer (test vs train) and inner (CV folds) partitions
  outer_partitions <- ifelse(partition_assignments == test_partition, "test", "train")

  # For inner partitions, renumber non-test partitions as CV folds
  inner_partitions <- rep(NA, length(partition_assignments))
  train_mask <- outer_partitions == "train"
  if (sum(train_mask) > 0) {
    train_partitions <- partition_assignments[train_mask]
    unique_train_partitions <- sort(unique(train_partitions))
    fold_mapping <- setNames(1:length(unique_train_partitions), unique_train_partitions)
    inner_partitions[train_mask] <- fold_mapping[as.character(train_partitions)]
  }

  return(list(
    outer_partitions = outer_partitions,
    inner_partitions = inner_partitions,
    clusters = all_clusters,
    cluster_stats = cluster_stats,
    test_partition_id = test_partition,
    silhouette = clara_result$silinfo$avg.width
  ))
}

# Select features for clustering based on weights
select_clustering_features <- function(
    feature_cols, top_n_features, feature_weights = NULL) {
  if (!is.null(feature_weights)) {
    top_features <- feature_weights |>
      sort(decreasing = TRUE) |>
      head(top_n_features) |>
      names()
  } else {
    top_features <- feature_cols[1:min(top_n_features, length(feature_cols))]
  }

  return(top_features)
}

# Evaluate cut-based partitioning for a given n_cuts and segregation_prob
evaluate_cut_partitioning <- function(
    n_cuts,
    segregation_prob = 0,
    data,
    precomputed_distances,
    feature_weights,
    top_n_features = 1,
    featuredist_sample_size = 1e4,
    seed = 42) {

  cat("Evaluating n_cuts =", n_cuts, ", segregation_prob =", segregation_prob, "...\n")

  # Create partitions using cut-based approach
  partition_result <- create_cut_based_partitions(
    data = data,
    feature_cols = setdiff(names(data), "response"),
    target_col = "response",
    n_cuts = n_cuts,
    n_partitions = 5,  # Keep 5 partitions (1 test + 4 CV)
    segregation_prob = segregation_prob,
    top_n_features = top_n_features,
    feature_weights = feature_weights,
    seed = seed
  )

  # Get top feature name (should match precomputed distances)
  top_features <- select_clustering_features(
    setdiff(names(data), "response"),
    top_n_features,
    feature_weights
  )
  top_feature <- top_features[1]

  # Prepare data for CV distance calculation
  # Create combined partition assignments: test partition gets its own ID
  all_partitions <- partition_result$outer_partitions
  combined_folds <- ifelse(all_partitions == "test",
                           0,  # Assign test partition to fold 0
                           partition_result$inner_partitions)

  # Calculate CV distances for this n_cuts across all outer partitions
  cv_dists <- calculate_cv_distances(
    training_data = data,
    variable_name = top_feature,
    cv_folds = combined_folds,
    sample_size = featuredist_sample_size,
    standardize = FALSE,
    seed = seed
  )

  # Debug: Print CV distance statistics
  cat("CV distances - n:", length(cv_dists),
      ", mean:", round(mean(cv_dists, na.rm = TRUE), 4),
      ", range:", round(range(cv_dists, na.rm = TRUE), 4), "\n")

  # Combine all distances and create result structure (unscaled)
  all_dists <- c(precomputed_distances$sample_to_sample,
                 precomputed_distances$prediction_to_sample,
                 cv_dists)

  distance_types <- factor(c(
    rep("sample-to-sample", length(precomputed_distances$sample_to_sample)),
    rep("prediction-to-sample", length(precomputed_distances$prediction_to_sample)),
    rep("CV-distances", length(cv_dists))
  ), levels = c("sample-to-sample", "prediction-to-sample", "CV-distances"))

  # Create feature_dists result structure for compatibility (no scaling)
  feature_dists <- tibble::tibble(
    dist = all_dists,
    what = distance_types,
    dist_type = "feature"
  )

  # Set class and attributes similar to featuredist output
  class(feature_dists) <- c("geodist", class(feature_dists))
  attr(feature_dists, "type") <- "feature"

  # Calculate W_CV statistic
  W_CV <- twosamples::wass_stat(
    feature_dists[feature_dists$what == "CV-distances", "dist"][[1]],
    feature_dists[feature_dists$what == "prediction-to-sample", "dist"][[1]]
  )

  # Print results
  cat("n_cuts =", n_cuts, ", segregation_prob =", segregation_prob,
      "-> W_CV =", round(W_CV, 4), "\n\n")

  return(list(
    n_cuts = n_cuts,
    segregation_prob = segregation_prob,
    partition_result = partition_result,
    W_CV = W_CV,
    silhouette = NA,  # Not applicable for cut-based approach
    top_feature = top_feature,
    feature_dists = feature_dists
  ))
}

# Evaluate clustering for a given k and return W_CV statistic
evaluate_clustering_k <- function(k_clusters, data, precomputed_distances, feature_weights,
                                  top_n_features = 1,
                                  clara_max_sample_size = 10e3,
                                  clara_sample_fraction = 0.25,
                                  clara_samples = 3,
                                  featuredist_sample_size = 1e4,
                                  seed = 42) {

  cat("Evaluating k =", k_clusters, "clusters...\n")

  # Create partitions using existing function
  partition_result <- create_single_clustering_partitions(
    data = data,
    feature_cols = setdiff(names(data), "response"),
    target_col = "response",
    n_clusters = k_clusters,
    n_partitions = 5,  # Keep 5 partitions (1 test + 4 CV)
    top_n_features = top_n_features,
    feature_weights = feature_weights,
    max_sample_size = clara_max_sample_size,
    sample_fraction = clara_sample_fraction,
    samples = clara_samples,
    seed = seed
  )

  # Get top feature name (should match precomputed distances)
  top_features <- select_clustering_features(
    setdiff(names(data), "response"),
    top_n_features,
    feature_weights
  )
  top_feature <- top_features[1]

  # Prepare data for CV distance calculation (all data including test partition)
  # Create combined partition assignments: test partition gets its own ID
  all_partitions <- partition_result$outer_partitions
  combined_folds <- ifelse(all_partitions == "test",
                           0,  # Assign test partition to fold 0
                           partition_result$inner_partitions)

  # Calculate CV distances for this k across all outer partitions
  cv_dists <- calculate_cv_distances(
    training_data = data,
    variable_name = top_feature,
    cv_folds = combined_folds,
    sample_size = featuredist_sample_size,
    standardize = FALSE,
    seed = seed
  )

  # Debug: Print CV distance statistics
  cat("CV distances - n:", length(cv_dists),
      ", mean:", round(mean(cv_dists, na.rm = TRUE), 4),
      ", range:", round(range(cv_dists, na.rm = TRUE), 4), "\n")

  # Combine all distances and create result structure (unscaled)
  all_dists <- c(precomputed_distances$sample_to_sample,
                 precomputed_distances$prediction_to_sample,
                 cv_dists)

  distance_types <- factor(c(
    rep("sample-to-sample", length(precomputed_distances$sample_to_sample)),
    rep("prediction-to-sample", length(precomputed_distances$prediction_to_sample)),
    rep("CV-distances", length(cv_dists))
  ), levels = c("sample-to-sample", "prediction-to-sample", "CV-distances"))

  # Create feature_dists result structure for compatibility (no scaling)
  feature_dists <- tibble::tibble(
    dist = all_dists,
    what = distance_types,
    dist_type = "feature"
  )

  # Set class and attributes similar to featuredist output
  class(feature_dists) <- c("geodist", class(feature_dists))
  attr(feature_dists, "type") <- "feature"

  # Calculate W_CV statistic
  W_CV <- twosamples::wass_stat(
    feature_dists[feature_dists$what == "CV-distances", "dist"][[1]],
    feature_dists[feature_dists$what == "prediction-to-sample", "dist"][[1]]
  )

  # Get silhouette width from partition result (no need to recompute clustering)
  silhouette_avg <- partition_result$silhouette

  # Print results for this k
  cat("k =", k_clusters, "-> W_CV =", round(W_CV, 4), ", silhouette =", round(silhouette_avg, 3), "\n\n")

  return(list(
    k = k_clusters,
    partition_result = partition_result,
    W_CV = W_CV,
    silhouette = silhouette_avg,
    top_feature = top_feature,
    feature_dists = feature_dists
  ))
}

featuredist <- function(training_data, prediction_data, variable_name,
                        cvfold_column = NULL, sample_size = 1e4,
                        seed = NULL, standardize = FALSE, scale_distances = FALSE) {
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Sample training and prediction data
  n_train <- min(nrow(training_data), sample_size)
  n_pred <- min(nrow(prediction_data), sample_size)
  
  train_sampled <- training_data[sample.int(nrow(training_data), n_train), ]
  pred_sampled <- prediction_data[sample.int(nrow(prediction_data), n_pred), ]
  
  # Extract variable values
  train_var <- train_sampled[[variable_name]]
  pred_var <- pred_sampled[[variable_name]]
  
  # Remove NA values
  train_clean <- train_var[!is.na(train_var)]
  pred_clean <- pred_var[!is.na(pred_var)]
  
  # Standardize if requested
  if (standardize) {
    train_mean <- mean(train_clean)
    train_sd <- sd(train_clean)
    if (train_sd > 0) {
      train_clean <- (train_clean - train_mean) / train_sd
      pred_clean <- (pred_clean - train_mean) / train_sd
    }
  }
  
  # Convert to matrices for FNN
  train_matrix <- matrix(train_clean, ncol = 1)
  pred_matrix <- matrix(pred_clean, ncol = 1)
  
  # Calculate sample-to-sample distances
  s2s_dists <- numeric(length(train_clean))
  for (i in seq_along(train_clean)) {
    dists <- FNN::knnx.dist(train_matrix[i, , drop = FALSE], train_matrix, k = 1)
    dists[i] <- NA  # Exclude self-distance
    s2s_dists[i] <- min(dists, na.rm = TRUE)
  }
  
  # Calculate prediction-to-sample distances
  s2p_dists <- numeric(length(pred_clean))
  for (i in seq_along(pred_clean)) {
    dists <- FNN::knnx.dist(pred_matrix[i, , drop = FALSE], train_matrix, k = 1)
    s2p_dists[i] <- min(dists, na.rm = TRUE)
  }
  
  # Combine results
  result <- tibble::tibble(
    dist = c(s2s_dists, s2p_dists),
    what = factor(c(rep("sample-to-sample", length(s2s_dists)),
                    rep("prediction-to-sample", length(s2p_dists))),
                  levels = c("sample-to-sample", "prediction-to-sample",
                             "CV-distances")),
    dist_type = "feature"
  )
  
  # Add CV distances if fold column is provided
  if (!is.null(cvfold_column)) {
    cv_folds <- train_sampled[[cvfold_column]]
    cv_folds <- cv_folds[!is.na(train_var)]  # Match cleaned training data
    
    cv_dists <- numeric(0)
    unique_folds <- unique(cv_folds)
    
    for (fold in unique_folds) {
      test_idx <- which(cv_folds == fold)
      train_idx <- which(cv_folds != fold)
      
      if (length(test_idx) > 0 && length(train_idx) > 0) {
        test_matrix <- matrix(train_clean[test_idx], ncol = 1)
        fold_train_matrix <- matrix(train_clean[train_idx], ncol = 1)
        
        for (i in seq_along(test_idx)) {
          dists <- FNN::knnx.dist(test_matrix[i, , drop = FALSE], fold_train_matrix, k = 1)
          cv_dists <- c(cv_dists, min(dists, na.rm = TRUE))
        }
      }
    }
    
    # Add CV distances to result
    cv_result <- tibble::tibble(
      dist = cv_dists,
      what = factor(rep("CV-distances", length(cv_dists))),
      dist_type = "feature"
    )
    
    result <- dplyr::bind_rows(result, cv_result)
  }
  
  # Scale distances if requested
  if (scale_distances) {
    max_dist <- max(result$dist, na.rm = TRUE)
    if (max_dist > 0) {
      result$dist <- result$dist / max_dist
    }
  }
  
  # Set class and attributes similar to geodist
  class(result) <- c("geodist", class(result))
  attr(result, "type") <- "feature"
  
  # Calculate W statistics
  W_sample <- twosamples::wass_stat(
    result[result$what == "sample-to-sample", "dist"][[1]],
    result[result$what == "prediction-to-sample", "dist"][[1]]
  )
  attr(result, "W_sample") <- W_sample
  
  if (!is.null(cvfold_column)) {
    W_CV <- twosamples::wass_stat(
      result[result$what == "CV-distances", "dist"][[1]],
      result[result$what == "prediction-to-sample", "dist"][[1]]
    )
    attr(result, "W_CV") <- W_CV
  }
  
  return(result)
}

# Calculate prediction-to-sample distances (independent of k)
calculate_prediction_distances <- function(training_data, prediction_data, variable_name,
                                          sample_size = 1e4, standardize = FALSE, seed = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Sample training and prediction data
  n_train <- min(nrow(training_data), sample_size)
  n_pred <- min(nrow(prediction_data), sample_size)

  train_sampled <- training_data[sample.int(nrow(training_data), n_train), ]
  pred_sampled <- prediction_data[sample.int(nrow(prediction_data), n_pred), ]

  # Extract variable values
  train_var <- train_sampled[[variable_name]]
  pred_var <- pred_sampled[[variable_name]]

  # Remove NA values
  train_clean <- train_var[!is.na(train_var)]
  pred_clean <- pred_var[!is.na(pred_var)]

  # Standardize if requested
  if (standardize) {
    train_mean <- mean(train_clean)
    train_sd <- sd(train_clean)
    if (train_sd > 0) {
      train_clean <- (train_clean - train_mean) / train_sd
      pred_clean <- (pred_clean - train_mean) / train_sd
    }
  }

  # Convert to matrices for FNN
  train_matrix <- matrix(train_clean, ncol = 1)
  pred_matrix <- matrix(pred_clean, ncol = 1)

  # Calculate sample-to-sample distances
  s2s_dists <- numeric(length(train_clean))
  for (i in seq_along(train_clean)) {
    dists <- FNN::knnx.dist(train_matrix[i, , drop = FALSE], train_matrix, k = 1)
    dists[i] <- NA  # Exclude self-distance
    s2s_dists[i] <- min(dists, na.rm = TRUE)
  }

  # Calculate prediction-to-sample distances
  s2p_dists <- numeric(length(pred_clean))
  for (i in seq_along(pred_clean)) {
    dists <- FNN::knnx.dist(pred_matrix[i, , drop = FALSE], train_matrix, k = 1)
    s2p_dists[i] <- min(dists, na.rm = TRUE)
  }

  return(list(
    sample_to_sample = s2s_dists,
    prediction_to_sample = s2p_dists,
    train_matrix = train_matrix,
    train_clean = train_clean
  ))
}

# Calculate CV distances only (depends on k via fold assignments)
calculate_cv_distances <- function(training_data, variable_name, cv_folds,
                                  sample_size = 1e4, standardize = FALSE, seed = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Sample training data (consistent with precomputation approach)
  n_train <- min(nrow(training_data), sample_size)
  sample_idx <- sample.int(nrow(training_data), n_train)
  train_sampled <- training_data[sample_idx, ]
  cv_folds_sampled <- cv_folds[sample_idx]

  # Extract variable values
  train_var <- train_sampled[[variable_name]]

  # Remove NA values
  train_clean <- train_var[!is.na(train_var)]
  cv_folds_clean <- cv_folds_sampled[!is.na(train_var)]

  # Standardize if requested (consistent with precomputation)
  if (standardize) {
    train_mean <- mean(train_clean)
    train_sd <- sd(train_clean)
    if (train_sd > 0) {
      train_clean <- (train_clean - train_mean) / train_sd
    }
  }

  cv_dists <- numeric(0)
  unique_folds <- unique(cv_folds_clean)

  # Remove NA folds if any
  unique_folds <- unique_folds[!is.na(unique_folds)]

  for (fold in unique_folds) {
    test_idx <- which(cv_folds_clean == fold)
    train_idx <- which(cv_folds_clean != fold)

    if (length(test_idx) > 0 && length(train_idx) > 0) {
      test_matrix <- matrix(train_clean[test_idx], ncol = 1)
      fold_train_matrix <- matrix(train_clean[train_idx], ncol = 1)

      for (i in seq_along(test_idx)) {
        dists <- FNN::knnx.dist(test_matrix[i, , drop = FALSE], fold_train_matrix, k = 1)
        cv_dists <- c(cv_dists, min(dists, na.rm = TRUE))
      }
    }
  }

  return(cv_dists)
}
