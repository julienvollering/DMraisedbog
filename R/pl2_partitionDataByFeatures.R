# Quantile Cut-based Feature Space Partitioning ####

# This script implements quantile cut-based partitioning for data splitting
# into 5 partitions using interleaved assignment across feature value ranges
# Optimizes n_cuts parameter to minimize W_CV (Wasserstein distance between
# CV-distances and prediction-to-sample distances)

library(readr)
library(dplyr)
library(terra)
library(tidyr)
library(ggplot2)
library(sf)
library(purrr)
library(furrr)
library(twosamples)
library(CAST)

## Configuration ####

# Parallelization
n_cores <- 4  # Physical cores (not logical processors)
plan(multisession, workers = n_cores)

# Run tests?
run_tests <- FALSE

# Feature distance calculation subsampling
# Higher values = more accurate W_CV estimates but slower computation
# Runtime scales quadratically with sample size
featuredist_sample_size <- 3e4     # Max samples for W_CV calculation

# Evaluation ranges for 2D grid search
# n_cuts: Controls granularity of feature space division
# segregation_prob: Controls balance between interleaving (0) and segregation (1)

n_cuts_range <- seq(500, 2000, by = 250)  # Coarser grid for 2D search
segregation_prob_range <- seq(0.25, 0.9, by = 0.05)  # 0 = pure interleaving, 0.9 = high segregation

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

## Systematic evaluation of n_cuts using W_CV minimization ####

# Get feature column names (exclude response variable)
feature_cols <- setdiff(names(mf_current), "response")

cat("Using feature weights for feature selection:\n")
print(weights)

cat("Evaluating n_cuts range:", paste(n_cuts_range, collapse = ", "), "\n")
cat("Evaluating segregation_prob range:", paste(segregation_prob_range, collapse = ", "), "\n")

# Pre-compute prediction distances (independent of n_cuts and segregation_prob)
cat("Pre-computing prediction distances...\n")
top_feature <- select_clustering_features(
  feature_cols,
  1,
  weights
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

# 2D grid search over n_cuts and segregation_prob (parallelized)
cat("Starting parallel 2D grid search over n_cuts × segregation_prob...\n")
cat("Using", n_cores, "cores\n")
partitioning_evaluation <- expand_grid(
  n_cuts = n_cuts_range,
  segregation_prob = segregation_prob_range
) |>
  mutate(
    results = future_map2(n_cuts, segregation_prob, ~evaluate_cut_partitioning(
      n_cuts = .x,
      segregation_prob = .y,
      data = mf_current,
      precomputed_distances = precomputed_distances,
      feature_weights = weights,
      top_n_features = 1,
      featuredist_sample_size = featuredist_sample_size,
      seed = 42
    ), .options = furrr_options(seed = TRUE)),
    W_CV = map_dbl(results, ~.x$W_CV),
    top_feature = map_chr(results, ~.x$top_feature)
  )

# Display evaluation results
cat("\nTop 10 parameter combinations by W_CV:\n")
evaluation_summary <- partitioning_evaluation |>
  select(n_cuts, segregation_prob, W_CV, top_feature) |>
  arrange(W_CV) |>
  head(10)
print(evaluation_summary)

# Visualization: 2D heatmap
p1 <- ggplot(partitioning_evaluation, aes(x = n_cuts, y = segregation_prob, fill = W_CV)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  geom_point(data = partitioning_evaluation |> slice_min(W_CV, n = 1),
             color = "red", size = 4, shape = 21, stroke = 2) +
  labs(
    title = "W_CV across n_cuts × segregation_prob",
    subtitle = "Red point = optimal combination",
    x = "Number of quantile cuts (n_cuts)",
    y = "Segregation probability",
    fill = "W_CV"
  ) +
  theme_minimal()

print(p1)

# Line plot: W_CV vs segregation_prob for each n_cuts
p2 <- ggplot(partitioning_evaluation, aes(x = segregation_prob, y = W_CV, color = factor(n_cuts))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "W_CV vs Segregation Probability",
    x = "Segregation probability",
    y = "W_CV (Wasserstein distance)",
    color = "n_cuts"
  ) +
  theme_minimal()

print(p2)

# Find optimal parameters (minimum W_CV)
optimal_row <- partitioning_evaluation |> slice_min(W_CV, n = 1)
optimal_n_cuts <- optimal_row$n_cuts
optimal_segregation_prob <- optimal_row$segregation_prob
cat("\nOptimal parameters (minimum W_CV):\n")
cat("  n_cuts =", optimal_n_cuts, "\n")
cat("  segregation_prob =", optimal_segregation_prob, "\n")
cat("  W_CV =", round(optimal_row$W_CV, 4), "\n")

dists_optimal <- partitioning_evaluation$results[[which.min(partitioning_evaluation$W_CV)]]$feature_dists
plot(dists_optimal, stat = "ecdf")
dists_optimal |>
  group_by(what) |>
  summarise(dist = median(dist))
dists_optimal |>
  ggplot(aes(x = dist, color = what)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(0, 10))

# Save evaluation results
write_csv(partitioning_evaluation |> select(-results),
          "output/pl2/n_cuts_evaluation.csv", append = FALSE)

# Extract optimal partitioning result
optimal_result <- optimal_row$results[[1]]

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
    interval = factor(partitioning_result$clusters),
    outer = factor(partitioning_result$outer_partitions),
    inner = partitioning_result$inner_partitions
  ) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(raster_current)) |>
  vect()

raster <- rasterize(
  mf_current_with_coords_partitioned,
  raster_current,
  field = c("interval", "outer", "inner")
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
