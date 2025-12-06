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
  seed = 42
) {
  cat(
    "Creating",
    n_partitions,
    "partitions using",
    n_cuts,
    "quantile cuts",
    "(segregation_prob =",
    segregation_prob,
    ")...\n"
  )

  # Select top feature(s)
  top_features <- select_clustering_features(
    feature_cols,
    top_n_features,
    feature_weights
  )

  # For now, use only the top feature (single variable)
  if (length(top_features) > 1) {
    cat("Note: Using only top feature:", top_features[1], "\n")
  }
  top_feature <- top_features[1]

  # Extract feature values
  feature_values <- data[[top_feature]]

  # Create quantile-based cuts
  quantile_breaks <- quantile(
    feature_values,
    probs = seq(0, 1, length.out = n_cuts + 1),
    na.rm = TRUE
  )

  # Handle case where quantiles are not unique (can happen with discrete data)
  quantile_breaks <- unique(quantile_breaks)
  actual_cuts <- length(quantile_breaks) - 1

  if (actual_cuts < n_cuts) {
    cat(
      "Warning: Only",
      actual_cuts,
      "unique quantile breaks found (requested",
      n_cuts,
      ")\n"
    )
  }

  intervals <- cut(
    feature_values,
    breaks = quantile_breaks,
    labels = FALSE,
    include.lowest = TRUE
  )

  # Probabilistic interleaving/segregation assignment
  # segregation_prob = 0: Pure interleaving (all intervals distributed round-robin)
  # segregation_prob = 1: Maximum segregation (intervals form contiguous blocks)

  # REVISED: Guarantee all n_partitions are non-empty by assigning first n_partitions intervals
  # one-per-partition, then apply probabilistic logic to remaining intervals

  set.seed(seed)
  unique_intervals_ordered <- sort(unique(intervals))
  n_intervals <- length(unique_intervals_ordered)
  interval_to_partition <- integer(n_intervals)

  # Step 1: Assign first n_partitions intervals (one per partition) to guarantee coverage
  if (n_intervals >= n_partitions) {
    interval_to_partition[1:n_partitions] <- 1:n_partitions

    # Step 2: Apply probabilistic logic to remaining intervals
    if (n_intervals > n_partitions) {
      for (i in (n_partitions + 1):n_intervals) {
        if (runif(1) < segregation_prob) {
          # SEGREGATE: Join contiguous block (inherit previous partition)
          interval_to_partition[i] <- interval_to_partition[i - 1]
        } else {
          # INTERLEAVE: Round-robin assignment
          interval_to_partition[i] <- ((i - 1) %% n_partitions) + 1
        }
      }
    }
  } else {
    # Edge case: fewer intervals than partitions
    # Assign round-robin (some partitions will be empty)
    cat(
      "Warning: Only",
      n_intervals,
      "intervals but",
      n_partitions,
      "partitions requested.\n"
    )
    for (i in seq_along(unique_intervals_ordered)) {
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

  cat("\nPartition assignment:\n")
  for (p in 1:n_partitions) {
    cat(sprintf(
      "Partition %d: %d obs (%.1f%%), %.3f prevalence\n",
      p,
      partition_stats$n_total[p],
      partition_stats$n_total[p] / sum(partition_stats$n_total) * 100,
      partition_stats$prevalence[p]
    ))
  }

  return(list(
    partitions = partition_assignments,
    clusters = intervals, # Return interval assignments for compatibility
    partition_stats = partition_stats,
    silhouette = NA # No silhouette for cut-based approach
  ))
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
  seed = 42
) {
  cat(
    "Evaluating n_cuts =",
    n_cuts,
    ", segregation_prob =",
    segregation_prob,
    "...\n"
  )

  # Create partitions using cut-based approach
  partition_result <- create_cut_based_partitions(
    data = data,
    feature_cols = setdiff(names(data), "response"),
    target_col = "response",
    n_cuts = n_cuts,
    n_partitions = 5, # Keep 5 partitions (1 test + 4 CV)
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

  # Use partition assignments directly as CV folds
  cv_folds <- partition_result$partitions

  # Calculate CV distances across all partitions
  cv_dists <- calculate_cv_distances(
    training_data = data,
    variable_name = top_feature,
    cv_folds = cv_folds,
    sample_size = featuredist_sample_size,
    standardize = FALSE,
    seed = seed
  )

  # Debug: Print CV distance statistics
  cat(
    "CV distances - n:",
    length(cv_dists),
    ", mean:",
    round(mean(cv_dists, na.rm = TRUE), 4),
    ", range:",
    round(range(cv_dists, na.rm = TRUE), 4),
    "\n"
  )

  # Combine all distances and create result structure (unscaled)
  all_dists <- c(
    precomputed_distances$sample_to_sample,
    precomputed_distances$prediction_to_sample,
    cv_dists
  )

  distance_types <- factor(
    c(
      rep("sample-to-sample", length(precomputed_distances$sample_to_sample)),
      rep(
        "prediction-to-sample",
        length(precomputed_distances$prediction_to_sample)
      ),
      rep("CV-distances", length(cv_dists))
    ),
    levels = c("sample-to-sample", "prediction-to-sample", "CV-distances")
  )

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
  cat(
    "n_cuts =",
    n_cuts,
    ", segregation_prob =",
    segregation_prob,
    "-> W_CV =",
    round(W_CV, 4),
    "\n\n"
  )

  return(list(
    n_cuts = n_cuts,
    segregation_prob = segregation_prob,
    partition_result = partition_result,
    W_CV = W_CV,
    silhouette = NA, # Not applicable for cut-based approach
    top_feature = top_feature,
    feature_dists = feature_dists
  ))
}

featuredist <- function(
  training_data,
  prediction_data,
  variable_name,
  cvfold_column = NULL,
  sample_size = 1e4,
  seed = NULL,
  standardize = FALSE,
  scale_distances = FALSE
) {
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
    # Min-max standardization using combined training + prediction data
    combined_values <- c(train_clean, pred_clean)
    min_val <- min(combined_values)
    max_val <- max(combined_values)

    if (max_val > min_val) {
      train_clean <- (train_clean - min_val) / (max_val - min_val)
      pred_clean <- (pred_clean - min_val) / (max_val - min_val)
    } else {
      # All values identical - set to 0.5 (middle of [0,1] range)
      train_clean <- rep(0.5, length(train_clean))
      pred_clean <- rep(0.5, length(pred_clean))
    }
  }

  # Convert to matrices for FNN
  train_matrix <- matrix(train_clean, ncol = 1)
  pred_matrix <- matrix(pred_clean, ncol = 1)

  # Calculate sample-to-sample distances
  s2s_dists <- numeric(length(train_clean))
  for (i in seq_along(train_clean)) {
    dists <- FNN::knnx.dist(
      train_matrix[i, , drop = FALSE],
      train_matrix,
      k = 1
    )
    dists[i] <- NA # Exclude self-distance
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
    what = factor(
      c(
        rep("sample-to-sample", length(s2s_dists)),
        rep("prediction-to-sample", length(s2p_dists))
      ),
      levels = c("sample-to-sample", "prediction-to-sample", "CV-distances")
    ),
    dist_type = "feature"
  )

  # Add CV distances if fold column is provided
  if (!is.null(cvfold_column)) {
    cv_folds <- train_sampled[[cvfold_column]]
    cv_folds <- cv_folds[!is.na(train_var)] # Match cleaned training data

    cv_dists <- numeric(0)
    unique_folds <- unique(cv_folds)

    for (fold in unique_folds) {
      test_idx <- which(cv_folds == fold)
      train_idx <- which(cv_folds != fold)

      if (length(test_idx) > 0 && length(train_idx) > 0) {
        test_matrix <- matrix(train_clean[test_idx], ncol = 1)
        fold_train_matrix <- matrix(train_clean[train_idx], ncol = 1)

        for (i in seq_along(test_idx)) {
          dists <- FNN::knnx.dist(
            test_matrix[i, , drop = FALSE],
            fold_train_matrix,
            k = 1
          )
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
calculate_prediction_distances <- function(
  training_data,
  prediction_data,
  variable_name,
  sample_size = 1e4,
  standardize = FALSE,
  seed = NULL
) {
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
    # Min-max standardization using combined training + prediction data
    combined_values <- c(train_clean, pred_clean)
    min_val <- min(combined_values)
    max_val <- max(combined_values)

    if (max_val > min_val) {
      train_clean <- (train_clean - min_val) / (max_val - min_val)
      pred_clean <- (pred_clean - min_val) / (max_val - min_val)
    } else {
      # All values identical - set to 0.5 (middle of [0,1] range)
      train_clean <- rep(0.5, length(train_clean))
      pred_clean <- rep(0.5, length(pred_clean))
    }
  }

  # Convert to matrices for FNN
  train_matrix <- matrix(train_clean, ncol = 1)
  pred_matrix <- matrix(pred_clean, ncol = 1)

  # Calculate sample-to-sample distances
  s2s_dists <- numeric(length(train_clean))
  for (i in seq_along(train_clean)) {
    dists <- FNN::knnx.dist(
      train_matrix[i, , drop = FALSE],
      train_matrix,
      k = 1
    )
    dists[i] <- NA # Exclude self-distance
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
calculate_cv_distances <- function(
  training_data,
  variable_name,
  cv_folds,
  sample_size = 1e4,
  standardize = FALSE,
  seed = NULL
) {
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
    # Min-max standardization using full training data
    min_val <- min(train_clean)
    max_val <- max(train_clean)

    if (max_val > min_val) {
      train_clean <- (train_clean - min_val) / (max_val - min_val)
    } else {
      # All values identical - set to 0.5 (middle of [0,1] range)
      train_clean <- rep(0.5, length(train_clean))
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
        dists <- FNN::knnx.dist(
          test_matrix[i, , drop = FALSE],
          fold_train_matrix,
          k = 1
        )
        cv_dists <- c(cv_dists, min(dists, na.rm = TRUE))
      }
    }
  }

  return(cv_dists)
}


# Function to calculate G-mean from predictions and true values
# Robust to factor level ordering by explicitly identifying positive/negative classes
calculate_gmean <- function(predicted_class, true_class) {
  # Convert to character to avoid factor level issues
  predicted_class <- as.character(predicted_class)
  true_class <- as.character(true_class)

  # Create confusion matrix
  cm <- table(Actual = true_class, Predicted = predicted_class)

  if (nrow(cm) == 2 && ncol(cm) == 2) {
    # Explicitly identify which row/column corresponds to positive class ("1")
    pos_row <- which(rownames(cm) == "1")
    pos_col <- which(colnames(cm) == "1")
    neg_row <- which(rownames(cm) == "0")
    neg_col <- which(colnames(cm) == "0")

    # Handle case where positive or negative class is missing
    if (
      length(pos_row) == 0 ||
        length(pos_col) == 0 ||
        length(neg_row) == 0 ||
        length(neg_col) == 0
    ) {
      tpr <- 0
      tnr <- 0
      gmean <- 0
    } else {
      # Extract confusion matrix values using explicit indexing
      tp <- cm[pos_row, pos_col] # True positives
      fn <- cm[pos_row, neg_col] # False negatives
      fp <- cm[neg_row, pos_col] # False positives
      tn <- cm[neg_row, neg_col] # True negatives

      # Calculate TPR and TNR
      tpr <- tp / (tp + fn) # Sensitivity
      tnr <- tn / (tn + fp) # Specificity

      # Calculate G-mean
      gmean <- sqrt(tpr * tnr)
    }
  } else {
    # Handle edge cases (only one class present)
    tpr <- 0
    tnr <- 0
    gmean <- 0
  }

  return(list(gmean = gmean, tpr = tpr, tnr = tnr))
}

# Function to calculate comprehensive metrics at a specific threshold
# Robust to factor level ordering and edge cases
get_metrics_at_threshold <- function(predicted_prob, true_class, threshold) {
  # Convert probabilities to class predictions using threshold
  predicted_class <- ifelse(predicted_prob >= threshold, "1", "0")

  # Convert to character to avoid factor level issues
  predicted_class <- as.character(predicted_class)
  true_class <- as.character(true_class)

  # Create confusion matrix
  cm <- table(Actual = true_class, Predicted = predicted_class)

  if (nrow(cm) == 2 && ncol(cm) == 2) {
    # Explicitly identify which row/column corresponds to positive class ("1")
    pos_row <- which(rownames(cm) == "1")
    pos_col <- which(colnames(cm) == "1")
    neg_row <- which(rownames(cm) == "0")
    neg_col <- which(colnames(cm) == "0")

    # Handle case where positive or negative class is missing
    if (
      length(pos_row) == 0 ||
        length(pos_col) == 0 ||
        length(neg_row) == 0 ||
        length(neg_col) == 0
    ) {
      tpr <- 0
      fpr <- 0
      tnr <- 0
      fnr <- 0
      gmean <- 0
    } else {
      # Extract confusion matrix values using explicit indexing
      tp <- cm[pos_row, pos_col] # True positives
      fn <- cm[pos_row, neg_col] # False negatives
      fp <- cm[neg_row, pos_col] # False positives
      tn <- cm[neg_row, neg_col] # True negatives

      # Calculate all metrics
      tpr <- tp / (tp + fn) # Sensitivity = TPR = 1 - FNR
      tnr <- tn / (tn + fp) # Specificity = TNR = 1 - FPR
      fpr <- fp / (fp + tn) # False positive rate = 1 - TNR
      fnr <- fn / (fn + tp) # False negative rate = 1 - TPR

      # Calculate G-mean
      gmean <- sqrt(tpr * tnr)
    }
  } else {
    # Handle edge cases (only one class present in predictions or true labels)
    # Check which class is missing to provide appropriate metrics
    if (nrow(cm) == 1 || ncol(cm) == 1) {
      # Some calculations may still be possible
      actual_classes <- rownames(cm)
      predicted_classes <- colnames(cm)

      if ("1" %in% actual_classes && "1" %in% predicted_classes) {
        # Only positive class present
        tp <- cm["1", "1"]
        tpr <- 1.0 # All positives correctly predicted
        fnr <- 0.0
        fpr <- NA # No negatives to calculate
        tnr <- NA
      } else if ("0" %in% actual_classes && "0" %in% predicted_classes) {
        # Only negative class present
        tn <- cm["0", "0"]
        tnr <- 1.0 # All negatives correctly predicted
        fpr <- 0.0
        tpr <- NA # No positives to calculate
        fnr <- NA
      } else {
        # Complete mismatch
        tpr <- 0
        fpr <- 1
        tnr <- 0
        fnr <- 1
      }
      gmean <- 0 # Cannot calculate meaningful G-mean
    } else {
      # Completely empty or malformed
      tpr <- 0
      fpr <- 0
      tnr <- 0
      fnr <- 0
      gmean <- 0
    }
  }

  return(data.frame(
    TPR = tpr,
    FPR = fpr,
    TNR = tnr,
    FNR = fnr,
    Gmean = gmean
  ))
}

# Function to calculate classification thresholds using ROCR
# Returns three thresholds: default (G-mean optimal or implicit), sens95, spec95
calculate_classification_thresholds <- function(
  predicted_prob,
  true_labels,
  method = c("gmean", "implicit"),
  implicit_class = NULL,
  target_sensitivity = 0.95,
  target_specificity = 0.95
) {
  method <- match.arg(method)

  # Convert to numeric for ROCR
  true_numeric <- as.numeric(as.character(true_labels))

  # Create ROCR prediction object
  pred_obj <- ROCR::prediction(predicted_prob, true_numeric)

  # Extract performance metrics
  # Note: ROCR orders thresholds from high to low (Inf -> -Inf)
  # TPR increases as threshold decreases (low -> high index)
  # TNR decreases as threshold decreases (low -> high index)
  tpr_values <- ROCR::performance(pred_obj, "tpr")@y.values[[1]]
  tnr_values <- ROCR::performance(pred_obj, "tnr")@y.values[[1]]
  all_thresholds <- ROCR::performance(pred_obj, "tpr")@x.values[[1]]

  # 1. Calculate default threshold based on method
  if (method == "gmean") {
    # G-mean optimal threshold (for BRF)
    gmean_values <- sqrt(tpr_values * tnr_values)
    optimal_idx <- which.max(gmean_values)
    threshold_default <- all_thresholds[optimal_idx]
  } else if (method == "implicit") {
    # Implicit threshold from class predictions (for RFQ)
    if (is.null(implicit_class)) {
      stop("implicit_class must be provided when method = 'implicit'")
    }
    predicted_positives <- implicit_class == "1"
    if (sum(predicted_positives) > 0) {
      threshold_default <- min(predicted_prob[predicted_positives])
    } else {
      threshold_default <- 0.5 # Fallback
    }
  }

  # 2. High sensitivity threshold: TPR >= target_sensitivity
  # Get FIRST threshold where TPR >= target (highest threshold = max specificity)
  sens_candidates <- which(tpr_values >= target_sensitivity)
  if (length(sens_candidates) > 0) {
    threshold_sens <- all_thresholds[sens_candidates[1]]
  } else {
    # If target sensitivity not achievable, use highest sensitivity
    threshold_sens <- all_thresholds[which.max(tpr_values)]
    cat(
      "Warning:",
      target_sensitivity * 100,
      "% sensitivity not achievable, using max TPR =",
      round(max(tpr_values), 3),
      "\n"
    )
  }

  # 3. High specificity threshold: TNR >= target_specificity
  # Get LAST threshold where TNR >= target (lowest threshold = max sensitivity)
  spec_candidates <- which(tnr_values >= target_specificity)
  if (length(spec_candidates) > 0) {
    threshold_spec <- all_thresholds[spec_candidates[length(spec_candidates)]]
  } else {
    # If target specificity not achievable, use highest specificity
    threshold_spec <- all_thresholds[which.max(tnr_values)]
    cat(
      "Warning:",
      target_specificity * 100,
      "% specificity not achievable, using max TNR =",
      round(max(tnr_values), 3),
      "\n"
    )
  }

  return(list(
    threshold_default = threshold_default,
    threshold_sens = threshold_sens,
    threshold_spec = threshold_spec
  ))
}

## Model training functions ####

# Function to train a BRF model with threshold optimization
train_brf_model <- function(
  cv_fold,
  mtry,
  nodesize,
  job_id,
  train_data,
  test_partition = 1
) {
  cat(
    "Job",
    job_id,
    "- BRF Fold",
    cv_fold,
    "mtry:",
    mtry,
    "nodesize:",
    nodesize,
    "...\n"
  )

  # Split data by partition (exclude test_partition from all tuning)
  # train_fold: all partitions except cv_fold and test_partition
  # val_fold: only cv_fold partition
  train_fold <- train_data |>
    filter(partition != test_partition, partition != cv_fold) |>
    select(-partition) |>
    as.data.frame()

  val_fold <- train_data |>
    filter(partition == cv_fold) |>
    select(-partition) |>
    as.data.frame()

  # Convert response to factor with explicit levels (0, 1) for consistency
  train_fold$response <- factor(
    as.character(train_fold$response),
    levels = c("0", "1")
  )
  val_fold$response <- factor(
    as.character(val_fold$response),
    levels = c("0", "1")
  )

  # Ensure ar50 is properly handled as factor if it exists
  if ("ar50" %in% names(train_fold)) {
    train_fold$ar50 <- as.factor(as.character(train_fold$ar50))
    val_fold$ar50 <- as.factor(as.character(val_fold$ar50))
  }

  # Calculate class sizes for balanced sampling
  n_presence <- sum(train_fold$response == "1")
  n_absence <- sum(train_fold$response == "0")
  balanced_size <- min(n_presence, n_absence)

  # Train BRF model with balanced sampling
  model <- randomForest::randomForest(
    formula = response ~ .,
    data = train_fold,
    ntree = 1000,
    mtry = mtry,
    nodesize = nodesize,
    sampsize = c(balanced_size, balanced_size), # Balanced sampling
    replace = TRUE,
    importance = FALSE,
    do.trace = FALSE
  )

  # Get predictions on training fold for threshold optimization
  train_pred_prob <- predict(model, newdata = train_fold, type = "prob")[, "1"]
  train_true <- as.numeric(as.character(train_fold$response))

  # Optimize threshold using ROCR on training folds
  pred_obj <- ROCR::prediction(train_pred_prob, train_true)

  # Calculate G-mean for all possible thresholds
  tpr_values <- ROCR::performance(pred_obj, "tpr")@y.values[[1]]
  tnr_values <- ROCR::performance(pred_obj, "tnr")@y.values[[1]]
  gmean_values <- sqrt(tpr_values * tnr_values)

  # Find optimal threshold
  optimal_idx <- which.max(gmean_values)
  optimal_threshold <- ROCR::performance(pred_obj, "tpr")@x.values[[1]][
    optimal_idx
  ]

  # Get predictions on validation fold
  val_pred_prob <- predict(model, newdata = val_fold, type = "prob")[, "1"]
  val_pred_class <- ifelse(val_pred_prob >= optimal_threshold, "1", "0")
  val_true <- as.character(val_fold$response)

  # Calculate G-mean on validation set with optimal threshold
  metrics <- calculate_gmean(val_pred_class, val_true)
  gmean <- metrics$gmean
  tpr <- metrics$tpr
  tnr <- metrics$tnr

  # Create result row
  result <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    model_type = "BRF",
    cv_fold = cv_fold,
    mtry = mtry,
    nodesize = nodesize,
    gmean = gmean,
    tpr = tpr,
    tnr = tnr,
    stringsAsFactors = FALSE
  )

  # Append to CSV immediately
  write_csv(result, progress_file, append = TRUE)
  cat("  Job", job_id, "completed - G-mean:", round(gmean, 4), "\n")

  return(result)
}

# Function to train a RFQ model (existing approach)
train_rfq_model <- function(
  cv_fold,
  mtry,
  nodesize,
  job_id,
  train_data,
  test_partition = 1
) {
  cat(
    "Job",
    job_id,
    "- RFQ Fold",
    cv_fold,
    "mtry:",
    mtry,
    "nodesize:",
    nodesize,
    "...\n"
  )

  # Split data by partition (exclude test_partition from all tuning)
  # train_fold: all partitions except cv_fold and test_partition
  # val_fold: only cv_fold partition
  train_fold <- train_data |>
    filter(partition != test_partition, partition != cv_fold) |>
    select(-partition) |>
    as.data.frame() # Convert tibble to data.frame for imbalanced()

  val_fold <- train_data |>
    filter(partition == cv_fold) |>
    select(-partition) |>
    as.data.frame() # Convert tibble to data.frame for imbalanced()

  # Convert response to factor with explicit levels (0, 1) for consistency
  train_fold$response <- factor(
    as.character(train_fold$response),
    levels = c("0", "1")
  )
  val_fold$response <- factor(
    as.character(val_fold$response),
    levels = c("0", "1")
  )

  # Ensure ar50 is properly handled as factor if it exists
  if ("ar50" %in% names(train_fold)) {
    train_fold$ar50 <- as.factor(as.character(train_fold$ar50))
    val_fold$ar50 <- as.factor(as.character(val_fold$ar50))
  }

  # Train model with current parameters
  model <- randomForestSRC::imbalanced(
    formula = response ~ .,
    data = train_fold,
    ntree = 3000,
    mtry = mtry,
    nodesize = nodesize,
    importance = FALSE, # Skip importance for speed during tuning
    do.trace = FALSE
  )

  # Predict on validation set
  pred <- predict(model, newdata = val_fold)

  # Calculate G-mean on validation set
  pred_class <- pred$class
  true_class <- val_fold$response

  # Calculate G-mean using common function
  metrics <- calculate_gmean(pred_class, true_class)
  gmean <- metrics$gmean
  tpr <- metrics$tpr
  tnr <- metrics$tnr

  # Create result row
  result <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    model_type = "RFQ",
    cv_fold = cv_fold,
    mtry = mtry,
    nodesize = nodesize,
    gmean = gmean,
    tpr = tpr,
    tnr = tnr,
    stringsAsFactors = FALSE
  )

  # Append to CSV immediately
  write_csv(result, progress_file, append = TRUE)
  cat("  Job", job_id, "completed - G-mean:", round(gmean, 4), "\n")

  return(result)
}

# Unified function to train either model type
train_single_model <- function(
  cv_fold,
  model_type,
  mtry,
  nodesize,
  job_id,
  train_data,
  test_partition = 1
) {
  if (model_type == "BRF") {
    return(train_brf_model(
      cv_fold,
      mtry,
      nodesize,
      job_id,
      train_data,
      test_partition
    ))
  } else if (model_type == "RFQ") {
    return(train_rfq_model(
      cv_fold,
      mtry,
      nodesize,
      job_id,
      train_data,
      test_partition
    ))
  } else {
    stop("Unknown model_type: ", model_type)
  }
}
