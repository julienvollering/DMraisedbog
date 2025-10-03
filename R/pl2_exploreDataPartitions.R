# Explore Data Partitions ####

# This script explores the data partitions created by pl2_partitionDataByFeatures.R
# by calculating nearest neighbor distance (NND) ECDF plots and W statistics for
# the top weighted features

library(readr)
library(dplyr)
library(terra)
library(tidyr)
library(ggplot2)
library(sf)
library(purrr)
library(furrr)
library(twosamples)

## Configuration ####

# Parallelization
n_cores <- 4 # Physical cores (not logical processors)
plan(multisession, workers = n_cores)

# Feature distance calculation subsampling
# Higher values = more accurate W_CV estimates but slower computation
# Runtime scales quadratically with sample size
featuredist_sample_size <- 1e4 # Max samples for W_CV calculation

# Number of top features to explore
n_top_features <- 10

## Read modeling frame and partitions ####

mf_partitioned <- read_csv("output/pl2/modeling_frame_regional_partitioned.csv")

# Extract current scenario with partitions
mf_current_with_partitions <- mf_partitioned |>
  filter(scenario == "current") |>
  select(-scenario)

# Create combined fold column for featuredist
combined_folds <- ifelse(mf_current_with_partitions$outer == "test",
  0, # Assign test partition to fold 0
  mf_current_with_partitions$inner
)

mf_current <- mf_current_with_partitions |>
  select(-x, -y, -outer, -inner) |> # Remove coordinates and partition columns
  mutate(cv_fold = combined_folds) # Add combined fold column

# Extract future scenario for prediction distances
toJoin <- mf_partitioned |>
  filter(scenario == "future") |>
  select(-scenario, -outer, -inner)

mf_future_presence <- mf_partitioned |>
  filter(response == 1) |>
  select(x, y) |>
  left_join(toJoin, by = c("x", "y")) |>
  select(-x, -y)

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

cat("\nTop", n_top_features, "weighted features for partitioning:\n")
print(head(weights, n_top_features))

## Functions ####

source("R/functions.R")

## Calculate distances for top features ####

# Get top N features
top_features <- names(head(weights, n_top_features))

cat("\nCalculating distances for top", n_top_features, "features...\n")

# Calculate distances for each feature using a single featuredist call
feature_distance_results <- future_map(top_features, function(feature_name) {
  cat("\nProcessing feature:", feature_name, "\n")

  # Calculate all distances using featuredist with cv_fold column
  dists <- featuredist(
    training_data = mf_current,
    prediction_data = mf_future_presence,
    variable_name = feature_name,
    cvfold_column = "cv_fold",
    sample_size = featuredist_sample_size,
    seed = 42,
    standardize = FALSE,
    scale_distances = TRUE
  )

  # Add feature name to result
  result <- dists |>
    mutate(feature = feature_name)

  # Extract W statistics from attributes
  W_sample <- attr(dists, "W_sample")
  W_CV <- attr(dists, "W_CV")

  cat("  W_sample =", round(W_sample, 4), "\n")
  cat("  W_CV =", round(W_CV, 4), "\n")

  # Preserve attributes
  attr(result, "W_sample") <- W_sample
  attr(result, "W_CV") <- W_CV

  return(result)
}, .options = furrr_options(seed = TRUE))

# Combine all results
all_feature_distances <- bind_rows(feature_distance_results)

## Print W statistics summary ####

cat("\n=== W Statistics Summary ===\n")
w_stats_summary <- map_dfr(feature_distance_results, function(result) {
  tibble(
    feature = unique(result$feature),
    W_sample = attr(result, "W_sample"),
    W_CV = attr(result, "W_CV")
  )
}) |>
  arrange(W_CV)

print(w_stats_summary)

## Visualization: ECDF plots for all features ####

# Create individual ECDF plots for each feature
for (feature_name in top_features) {
  p <- all_feature_distances |>
    filter(feature == feature_name) |>
    ggplot(aes(x = dist, color = what)) +
    stat_ecdf(linewidth = 1) +
    coord_cartesian(xlim = c(0, NA)) +
    labs(
      title = paste("Nearest Neighbor Distance ECDF:", feature_name),
      x = "Distance (scaled)",
      y = "ECDF",
      color = "Distance type"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  print(p)
}

# sessionInfo ####

sessioninfo::session_info()
