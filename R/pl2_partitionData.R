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
library(purrr)
library(twosamples)

## Configuration: Subsampling for Performance Control ####

# Run tests?
run_tests <- FALSE

# CLARA clustering subsampling
# Higher values = more accurate clustering but slower computation
# Memory usage scales with sample size
clara_max_sample_size <- 1e4      # Max samples for CLARA algorithm
clara_sample_fraction <- 0.05     # Fraction of data to sample (if smaller than max)
clara_samples <- 1                # Number of sample sets to draw and evaluate

# Feature distance calculation subsampling
# Higher values = more accurate W_CV estimates but slower computation
# Runtime scales quadratically with sample size
featuredist_sample_size <- 3e4     # Max samples for featuredist W_CV calculation

# Evaluation range
# More k values = more thorough search but longer runtime
# Runtime scales linearly with number of k values
k_range <- seq(5, 50, by = 20)     # Range of k_clusters to evaluate

# Performance guidance:
# - For exploration: use current defaults
# - For production: consider increasing sample sizes
# - For speed: reduce sample sizes or k_range
# - Memory issues: reduce clara_max_sample_size

## Read modeling frame ####
mf <- read_csv("output/pl2/modeling_frame_regional.csv")

raster_current <- rast("output/pl2/scenario_current.tif")

mf_current_with_coords <- mf |> 
  filter(scenario == "current") |> 
  select(-scenario)
mf_current <- mf_current_with_coords |>
  select(-x, -y) # Remove x and y coordinates for modeling

toJoin <- mf |> 
  filter(scenario == "future") |> 
  select(-scenario)
mf_future_presence <- mf |> 
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

## Functions ####

source("R/functions.R")

# Unit test: Verify prediction-to-sample distances are consistent
if (run_tests) {
  
  cat("Running unit test: prediction-to-sample distance consistency...\n")
  test_distances_1 <- calculate_prediction_distances(
    training_data = mf_current,
    prediction_data = mf_future_presence,
    variable_name = top_feature,
    sample_size = nrow(mf_current),  # Use full dataset
    standardize = FALSE,
    seed = 42
  )
  
  test_distances_2 <- calculate_prediction_distances(
    training_data = mf_current,
    prediction_data = mf_future_presence,
    variable_name = top_feature,
    sample_size = nrow(mf_current),  # Use full dataset
    standardize = FALSE,
    seed = 42
  )
  
  # Check if prediction-to-sample distances are identical
  pred_distances_match <- identical(test_distances_1$prediction_to_sample,
                                    test_distances_2$prediction_to_sample)
  sample_distances_match <- identical(test_distances_1$sample_to_sample,
                                      test_distances_2$sample_to_sample)
  
  cat("Prediction-to-sample distances identical:", pred_distances_match, "\n")
  cat("Sample-to-sample distances identical:", sample_distances_match, "\n")
  
  if (pred_distances_match && sample_distances_match) {
    cat("✓ Unit test PASSED: Distance calculations are deterministic\n\n")
  } else {
    cat("✗ Unit test FAILED: Distance calculations are not consistent\n")
    cat("This suggests sampling or randomization issues in distance calculation\n\n")
  }
  
}

## Systematic evaluation of k_clusters using W_CV minimization ####

# Get feature column names (exclude response variable)
feature_cols <- setdiff(names(mf_current), "response")

cat("Using feature weights for feature selection:\n")
print(weights)

cat("Evaluating k_clusters range:", paste(k_range, collapse = ", "), "\n")

# Pre-compute prediction distances (independent of k)
cat("Pre-computing prediction distances...\n")
top_features <- select_clustering_features(
  feature_cols, 1, weights
)

precomputed_distances <- calculate_prediction_distances(
  training_data = mf_current,
  prediction_data = mf_future_presence,
  variable_name = top_feature,
  sample_size = featuredist_sample_size,
  standardize = FALSE,
  seed = 42
)
cat("Pre-computation complete.\n\n")

# Systematic evaluation using tibble and purrr
cat("Starting systematic evaluation of k_clusters...\n")
clustering_evaluation <- tibble(k = k_range) |>
  mutate(
    results = map(k, ~evaluate_clustering_k(
      k_clusters = .x,
      data = mf_current,
      precomputed_distances = precomputed_distances,
      feature_weights = weights,
      top_n_features = 1,
      clara_max_sample_size = clara_max_sample_size,
      clara_sample_fraction = clara_sample_fraction,
      clara_samples = clara_samples,
      featuredist_sample_size = featuredist_sample_size,
      seed = 42
    )),
    W_CV = map_dbl(results, ~.x$W_CV),
    silhouette = map_dbl(results, ~.x$silhouette),
    top_feature = map_chr(results, ~.x$top_feature)
  )

# Display evaluation results
cat("\nK-clusters evaluation results:\n")
evaluation_summary <- clustering_evaluation |>
  select(k, W_CV, silhouette, top_feature) |>
  arrange(W_CV)
print(evaluation_summary)

# Visualization of evaluation results
p1 <- ggplot(clustering_evaluation, aes(x = k, y = W_CV)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 2) +
  geom_vline(xintercept = clustering_evaluation$k[which.min(clustering_evaluation$W_CV)],
             linetype = "dashed", color = "red", alpha = 0.7) +
  labs(
    title = "W_CV vs Number of Clusters",
    subtitle = paste("Optimal k =", clustering_evaluation$k[which.min(clustering_evaluation$W_CV)]),
    x = "Number of clusters (k)",
    y = "W_CV (scaled feature distance)"
  ) +
  theme_minimal()

print(p1)

plot(clustering_evaluation$results[[1]]$feature_dists, stat = "ecdf")
clustering_evaluation$results[[1]]$feature_dists |> 
  group_by(what) |> 
  summarise(dist = median(dist))
clustering_evaluation$results[[1]]$feature_dists |> 
  ggplot(aes(x = dist, color = what)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(0, 100))
plot(clustering_evaluation$results[[2]]$feature_dists, stat = "ecdf")
clustering_evaluation$results[[2]]$feature_dists |> 
  group_by(what) |> 
  summarise(dist = median(dist))
clustering_evaluation$results[[2]]$feature_dists |> 
  ggplot(aes(x = dist, color = what)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(0, 100))
plot(clustering_evaluation$results[[3]]$feature_dists, stat = "ecdf")
clustering_evaluation$results[[3]]$feature_dists |> 
  group_by(what) |> 
  summarise(dist = median(dist))
clustering_evaluation$results[[3]]$feature_dists |> 
  ggplot(aes(x = dist, color = what)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(0, 100))


# Find optimal k (minimum W_CV)
optimal_k <- clustering_evaluation$k[which.min(clustering_evaluation$W_CV)]
cat("\nOptimal k_clusters (minimum W_CV):", optimal_k, "\n")

# Save evaluation results
write_csv(clustering_evaluation |> select(-results),
          "output/pl2/k_clusters_evaluation.csv", append = FALSE)

# Extract optimal partitioning result
optimal_result <- clustering_evaluation |>
  filter(k == optimal_k) |>
  pull(results) |>
  pluck(1)

plot(optimal_result$feature_dists, stat = "ecdf")

partitioning_result <- optimal_result$partition_result

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
