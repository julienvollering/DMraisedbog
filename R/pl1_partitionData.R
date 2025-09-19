# CLARA-based Feature Space Clustering for Data Partitioning ####

# This script implements CLARA clustering methodology for data partitioning
# into K clusters with standardized euclidean distances of top X features

library(readr)
library(dplyr)
library(terra)
library(cluster)
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
presence_df <- terra::extract(preds_train, presence[c("x", "y")],
  xy = TRUE,
  ID = FALSE
) |>
  mutate(response = 1, .before = 1) |>
  drop_na() |>
  tibble()

# Extract predictor values for absence points
absence_df <- terra::extract(
  preds_train, absence[c("x", "y")],
  xy = TRUE,
  ID = FALSE
) |>
  mutate(response = 0, .before = 1) |>
  drop_na() |>
  tibble()

# Combine into training dataset (keep x,y coordinates for later use)
data_with_coords <- bind_rows(presence_df, absence_df) |>
  mutate(response = factor(response))

data <- data_with_coords |>
  select(-x, -y) # Remove x and y coordinates for modeling

# Get final counts
n_final_presence <- sum(data$response == "1")
n_final_absence <- sum(data$response == "0")
n_final_prevalence <- n_final_presence / (n_final_presence + n_final_absence)

cat("Final training data - Presences:", n_final_presence, "Absence:", n_final_absence, "\n")
cat("Prevalence in training data:", round(n_final_prevalence * 100, 3), "%\n")

# Load feature weights
weights_features <- read_csv("output/pl1/weights_feature_data_partitioning.csv")
weights_features <- weights_features |>
  mutate(feature = forcats::fct_reorder(
    feature, median_normalized,
    .fun = first
  ))

weights <- weights_features |>
  filter(grepl("Balanced", method)) |>
  select(feature, median_normalized) |>
  mutate(
    # Set weights of artype* features to zero to exclude them from clustering
    weight_partitioning = case_when(
      grepl("artype", feature) ~ 0,
      TRUE ~ median_normalized
    )
  ) |>
  select(feature, weight_partitioning) |>
  arrange(desc(weight_partitioning)) |>
  tibble::deframe()

## Functions ####

### Core functions ####

# Run CLARA clustering with standardization and return cluster object
run_clara_clustering <- function(data, top_features, k_clusters, target_col,
                                 seed = 42) {
  # Extract feature data for clustering (numeric only)
  clustering_data <- data |>
    select(all_of(top_features)) |>
    select(where(is.numeric))

  cat("Clustering data dimensions:", nrow(clustering_data), "x", ncol(clustering_data), "\n")

  # Run CLARA clustering with standardized Euclidean distance
  set.seed(seed)
  clara_result <- clara(
    x = clustering_data,
    k = k_clusters,
    metric = "euclidean",
    # Standardize features:
    stand = TRUE,
    rngR = TRUE,
    # Number of samples to be drawn from dataset -- docs recommend 10-100:
    samples = 10,
    # Observations in each sample -- docs recommended > 40 + 2*k:
    # 10e3 is about 1% of data, SAC means nearby points are similar.
    # Need to capture full spatial variation.
    sampsize = min(10e3, nrow(clustering_data) * 0.05), # Should be 10e3!
    trace = 1
  )

  # Add cluster assignments to data
  data$cluster <- factor(clara_result$clustering)

  # Show cluster summary
  cluster_summary <- data |>
    group_by(cluster) |>
    summarise(
      n_points = n(),
      n_presence = sum(as.numeric(as.character(.data[[target_col]]))),
      prevalence = n_presence / n_points,
      .groups = "drop"
    )

  cat("\nCluster summary:\n")
  print(cluster_summary)

  return(list(
    clara_result = clara_result,
    data_with_clusters = data,
    clustering_data = clustering_data
  ))
}

# Assign clusters to partitions using dual-constraint optimization
assign_clusters_to_partitions <- function(
    cluster_stats, partition_props, partition_names,
    tolerance = 0.05) {
  n_partitions <- length(partition_props)

  cat("Assigning clusters to", n_partitions, "partitions with dual constraints...\n")

  # Sort clusters by prevalence (highest first) to prioritize presence balance
  cluster_stats <- cluster_stats[order(cluster_stats$prevalence, decreasing = TRUE), ]

  # Initialize tracking variables
  cluster_assignments <- character(nrow(cluster_stats))
  names(cluster_assignments) <- as.character(cluster_stats$cluster)

  # Calculate targets for both presence and total observations
  total_presence <- sum(cluster_stats$n_presence)
  total_observations <- sum(cluster_stats$n_total)

  target_presence <- round(total_presence * partition_props)
  target_total <- round(total_observations * partition_props)
  names(target_presence) <- partition_names
  names(target_total) <- partition_names

  # Track current counts
  current_presence <- rep(0, n_partitions)
  current_total <- rep(0, n_partitions)
  names(current_presence) <- partition_names
  names(current_total) <- partition_names

  cat("DUAL-CONSTRAINT TARGETS:\n")
  cat("Total observations:", total_observations, "| Total presences:", total_presence, "\n")
  for (i in 1:n_partitions) {
    cat(sprintf(
      "%s targets: %d obs (%.1f%%), %d pres\n",
      partition_names[i], target_total[i], partition_props[i] * 100, target_presence[i]
    ))
  }

  # Assign clusters using dual-constraint priority scoring
  for (i in 1:nrow(cluster_stats)) {
    cluster_id <- as.character(cluster_stats$cluster[i])
    cluster_presence <- cluster_stats$n_presence[i]
    cluster_total <- cluster_stats$n_total[i]

    # Check capacity constraints (hard stops)
    capacity <- (current_total + cluster_total) <= (target_total * (1 + tolerance))

    # Calculate presence need scores (0 = satisfied, 1 = maximum need)
    pres_need <- pmax(0, (target_presence - current_presence) / pmax(target_presence, 1))

    # Calculate capacity availability scores (0 = no capacity, 1 = full capacity)
    cap_score <- ifelse(capacity, 1 - (current_total / pmax(target_total, 1)), 0)

    # Calculate priority scores (presence need * capacity availability)
    priority <- pres_need * cap_score

    # Assign to partition with highest priority (that has capacity)
    best_partition <- which.max(priority)
    if (capacity[best_partition]) {
      cluster_assignments[cluster_id] <- partition_names[best_partition]
      current_presence[best_partition] <- current_presence[best_partition] + cluster_presence
      current_total[best_partition] <- current_total[best_partition] + cluster_total
    } else {
      # Fallback: assign to partition with smallest overshoot
      overshoot <- pmax(0, (current_total + cluster_total) - target_total)
      best_partition <- which.min(overshoot)

      cluster_assignments[cluster_id] <- partition_names[best_partition]
      current_presence[best_partition] <- current_presence[best_partition] + cluster_presence
      current_total[best_partition] <- current_total[best_partition] + cluster_total
    }
  }

  cat("\nFINAL ASSIGNMENT RESULTS:\n")
  for (i in 1:n_partitions) {
    cat(sprintf(
      "%s: %d/%d obs (%.1f%%), %d/%d pres\n",
      partition_names[i], current_total[i], target_total[i],
      round(current_total[i] / total_observations * 100, 1),
      current_presence[i], target_presence[i]
    ))
  }

  return(cluster_assignments)
}

# Coordinate clustering, assignment, and results
create_clara_partitions <- function(
    data,
    feature_cols,
    target_col,
    k_clusters,
    top_n_features,
    partition_props,
    partition_names,
    feature_weights = NULL,
    tolerance = 0.05,
    seed = 42) {
  cat("Creating CLARA clusters for data partitioning...\n")

  # Select features and run clustering
  top_features <- select_clustering_features(
    feature_cols, top_n_features,
    feature_weights
  )
  clustering_result <- run_clara_clustering(
    data, top_features, k_clusters,
    target_col, seed
  )

  all_clusters <- clustering_result$clara_result$clustering
  cluster_stats <- calculate_group_stats(
    data, all_clusters, target_col, "cluster",
    verbose = TRUE
  )
  cluster_assignments <- assign_clusters_to_partitions(
    cluster_stats, partition_props, partition_names, tolerance
  )
  all_partitions <- cluster_assignments[as.character(all_clusters)]
  partition_stats <- calculate_group_stats(
    data, all_partitions, target_col, "partition"
  )

  results <- list(
    partitions = all_partitions,
    partition_stats = partition_stats,
    clusters = all_clusters,
    cluster_stats = cluster_stats,
    cluster_assignments = cluster_assignments
  )

  return(results)
}

### Helper functions ####

# Select features for clustering based on weights
select_clustering_features <- function(
    feature_cols, top_n_features, feature_weights = NULL) {
  if (!is.null(feature_weights)) {
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


# Concise function to calculate group statistics using tidyverse
calculate_group_stats <- function(
    data, groups, target_col, group_col_name, verbose = FALSE) {
  if (verbose) cat("Calculating group statistics...\n")

  data |>
    mutate(group_var = groups) |>
    group_by(group_var) |>
    summarise(
      n_total = n(),
      n_presence = sum(.data[[target_col]] == 1),
      prevalence = n_presence / n_total,
      .groups = "drop"
    ) |>
    rename(!!group_col_name := group_var)
}

### High-level interface functions ####

create_outer <- function(
    data,
    feature_cols,
    target_col = "response",
    k_clusters = 10,
    top_n_features = 5,
    train_prop = 0.6,
    calib_prop = 0.2,
    test_prop = 0.2,
    feature_weights = NULL,
    seed = 42) {
  result <- create_clara_partitions(
    data = data,
    feature_cols = feature_cols,
    target_col = target_col,
    k_clusters = k_clusters,
    top_n_features = top_n_features,
    partition_props = c(train_prop, calib_prop, test_prop),
    partition_names = c("train", "calib", "test"),
    feature_weights = feature_weights,
    seed = seed
  )

  return(result)
}

create_inner <- function(
    train_data,
    feature_cols,
    target_col = "response",
    k_clusters = 6,
    n_folds = 3,
    top_n_features = 5,
    feature_weights = NULL,
    seed = 42) {
  cat("Creating inner CV folds for training data using CLARA...\n")
  cat("Training data size:", nrow(train_data), "observations\n")

  # Create equal proportions for CV folds
  fold_props <- rep(1 / n_folds, n_folds)
  fold_names <- as.character(1:n_folds)

  result <- create_clara_partitions(
    data = train_data,
    feature_cols = feature_cols,
    target_col = target_col,
    k_clusters = k_clusters,
    top_n_features = top_n_features,
    partition_props = fold_props,
    partition_names = fold_names,
    feature_weights = feature_weights,
    seed = seed
  )

  return(result)
}

## Apply CLARA clustering to create outer splits ####

# Get feature column names (exclude response variable)
feature_cols <- setdiff(names(data), "response")

cat("Feature columns:", paste(feature_cols, collapse = ", "), "\n")

cat("Using feature weights for feature selection:\n")
print(weights)

# Create CLARA-based feature space splits
outer <- create_outer(
  data = data,
  feature_cols = feature_cols,
  target_col = "response",
  k_clusters = 10, # Number of CLARA clusters
  top_n_features = 6, # Top features by weight
  train_prop = 0.6,
  calib_prop = 0.2,
  test_prop = 0.2,
  feature_weights = weights
)

# Display results
print(outer$partition_stats)
print(outer$cluster_stats)

## Map visualization ####
data_with_coords_outer <- data_with_coords |>
  select(response, x, y) |>
  mutate(
    cluster = factor(outer$clusters),
    partition = factor(outer$partitions)
  ) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
  vect()
raster_outer <- rasterize(
  data_with_coords_outer,
  preds_train,
  field = c("cluster", "partition")
)
plot(raster_outer)

## Apply CLARA clustering to create inner CV folds ####

# Get training data indices
train_idx <- which(outer$partitions == "train")
train_data <- data[train_idx, ]

cat("Creating inner CV folds for training partition...\n")

# Create CLARA-based inner CV folds
inner <- create_inner(
  train_data = train_data,
  feature_cols = feature_cols,
  target_col = "response",
  k_clusters = 6, # Proportional to 3 CV folds
  n_folds = 3, # Number of CV folds
  top_n_features = 6, # Same as outer partitioning
  feature_weights = weights,
  seed = 42
)

# Display results
print(inner$partition_stats)
print(inner$cluster_stats)

## Map visualization ####
data_with_coords_inner <- data_with_coords[train_idx, ] |>
  select(response, x, y) |>
  mutate(
    cluster = factor(inner$clusters),
    fold = factor(inner$partitions)
  ) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
  vect()
raster_inner <- rasterize(
  data_with_coords_inner,
  preds_train,
  field = c("cluster", "fold")
)
plot(raster_inner)

## Save partitioned datasets ####

### Main dataset for modeling ####
# Create main dataset with outer and inner partition columns
# Map inner CV folds back to original data indices
inner_folds_full <- rep(NA_integer_, nrow(data))
inner_folds_full[train_idx] <- inner$partitions

data_partitioned <- data_with_coords |>
  mutate(
    outer = outer$partitions,
    inner = inner_folds_full
  )

write_csv(data_partitioned, "output/pl1/data_partitioned.csv")

### Partition information in space ####
writeRaster(raster_outer, "output/pl1/partition_outer.tif", overwrite = TRUE)
writeRaster(raster_inner, "output/pl1/partition_inner.tif", overwrite = TRUE)

# sessionInfo ####

sessioninfo::session_info()
