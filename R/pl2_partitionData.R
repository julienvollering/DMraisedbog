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

## Read modeling frame ####
mf <- read_csv("output/pl2/modeling_frame_regional.csv")

raster_current <- rast("output/pl2/scenario_current.tif")

mf_current_with_coords <- mf |> 
  filter(scenario == "current") |> 
  select(-scenario)
mf_current <- mf_current_with_coords |>
  select(-x, -y) # Remove x and y coordinates for modeling

n_presence <- sum(mf_current$response == "1")
n_absence <- sum(mf_current$response == "0")
prevalence <- n_presence / (n_presence + n_absence)

cat("Final training data - Presences:", n_presence, "Absence:", n_absence, "\n")
cat("Prevalence in training data:", round(prevalence * 100, 2), "%\n")

## Load feature weights ####
weights_features <- read_csv("output/pl2/weights_feature_data_partitioning.csv")
weights_features <- weights_features |>
  mutate(feature = forcats::fct_reorder(
    feature, median,
    .fun = first
  ))

weights <- weights_features |>
  select(feature, median) |>
  arrange(desc(median)) |>
  tibble::deframe()

## Functions ####

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
    stand = TRUE,
    rngR = TRUE,
    samples = 10,
    sampsize = min(10e3, nrow(clustering_data) * 0.25),
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
    seed = 42) {

  cat("Creating", n_clusters, "CLARA clusters for", n_partitions, "partitions...\n")

  # Select features and run clustering
  top_features <- select_clustering_features(
    feature_cols, top_n_features, feature_weights
  )
  clara_result <- run_clara_clustering(
    data, top_features, n_clusters, target_col, seed
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
    test_partition_id = test_partition
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

    cat("Using top", top_n_features, "features:", paste(top_features, collapse = ", "), "\n")
  } else {
    top_features <- feature_cols[1:min(top_n_features, length(feature_cols))]
    cat("No feature weights provided. Using first", length(top_features), "features\n")
  }

  return(top_features)
}

## Apply single CLARA clustering for partitioning ####

# Get feature column names (exclude response variable)
feature_cols <- setdiff(names(mf_current), "response")

cat("Using feature weights for feature selection:\n")
print(weights)

# Create single CLARA clustering with consolidation
partitioning_result <- create_single_clustering_partitions(
  data = mf_current,
  feature_cols = feature_cols,
  target_col = "response",
  n_clusters = 5,        # Number of CLARA clusters
  n_partitions = 5,      # Number of partitions: 1 outer test, remaining inner CV
  top_n_features = 1,    # Top features by weight
  feature_weights = weights
)

# Display results
cat("\nOuter partition summary (test vs train):\n")
outer_summary <- table(partitioning_result$outer_partitions)
print(outer_summary)

cat("\nInner partition summary (CV folds):\n")
inner_summary <- table(partitioning_result$inner_partitions, useNA = "always")
print(inner_summary)

## Map visualization ####
mf_current_with_coords_partitioned <- mf_current_with_coords |>
  select(response, x, y) |>
  mutate(
    cluster = factor(partitioning_result$clusters),
    outer = factor(partitioning_result$outer_partitions),
    inner = partitioning_result$inner_partitions
  ) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(raster_current)) |>
  vect()

raster <- rasterize(
  mf_current_with_coords_partitioned,
  raster_current,
  field = c("cluster", "outer", "inner")
)
plot(raster["outer"])
plot(raster["inner"])

## Save partitioned datasets ####

mf_partitioned <- mf_current_with_coords |>
  mutate(
    scenario = "current",
    outer = partitioning_result$outer_partitions,
    inner = partitioning_result$inner_partitions,
    .before = 1
  ) |>
  bind_rows(filter(mf, scenario != "current")) |>
  select(scenario, outer, inner, response, everything())

write_csv(mf_partitioned, "output/pl2/modeling_frame_regional_partitioned.csv",
          append = FALSE)

writeRaster(raster, "output/pl2/partition.tif", overwrite = TRUE)

# sessionInfo ####

sessioninfo::session_info()
