# Spatial Blocking for Data Partitioning ####

# This script implements spatial blocking methodology for data partitioning
# using blockCV::cv_spatial with systematic evaluation of block sizes

library(readr)
library(dplyr)
library(terra)
library(tidyr)
library(ggplot2)
library(sf)
library(purrr)
library(twosamples)
library(blockCV)
library(CAST)

## Configuration: Subsampling for Performance Control ####

# Run tests?
run_tests <- FALSE

# Feature distance calculation subsampling
# Higher values = more accurate W_CV estimates but slower computation
# Runtime scales quadratically with sample size
featuredist_sample_size <- 3e4     # Max samples for featuredist W_CV calculation

# Evaluation range
# More block sizes = more thorough search but longer runtime
# Runtime scales linearly with number of block sizes
block_sizes_km <- seq(10, 150, by = 20)  # Range of block sizes (in km) to evaluate

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

# Get top feature for distance calculations
top_feature <- names(weights)[1]

## Functions ####

source("R/functions.R")

# Create spatial blocking partitions using blockCV
create_spatial_blocking_partitions <- function(
    data_with_coords,
    block_size_m,
    target_col = "response",
    n_partitions = 5,
    seed = 42) {

  cat("Creating spatial blocks with size =", block_size_m, "m...\n")

  # Convert data to sf object
  data_sf <- data_with_coords |>
    st_as_sf(coords = c("x", "y"), crs = st_crs(raster_current))

  # Create spatial blocks using blockCV
  set.seed(seed)
  spatial_cv <- cv_spatial(
    x = data_sf,
    column = target_col,
    size = block_size_m,
    k = n_partitions,
    selection = "random",
    iteration = 100,
    progress = FALSE
  )

  # Extract fold assignments
  fold_assignments <- spatial_cv$folds_ids

  # Calculate fold statistics
  fold_stats <- data_with_coords |>
    mutate(fold = factor(fold_assignments)) |>
    group_by(fold) |>
    summarise(
      n_total = n(),
      n_presence = sum(.data[[target_col]] == 1),
      prevalence = n_presence / n_total,
      .groups = "drop"
    )

  cat("\nFold summary:\n")
  print(fold_stats)

  # Find fold with prevalence closest to overall prevalence
  overall_prevalence <- sum(fold_stats$n_presence) / sum(fold_stats$n_total)
  prevalence_diffs <- abs(fold_stats$prevalence - overall_prevalence)
  test_fold <- as.numeric(as.character(fold_stats$fold[which.min(prevalence_diffs)]))

  cat("\nFold assignments:\n")
  for (f in sort(unique(fold_assignments))) {
    fold_data <- fold_stats[fold_stats$fold == f, ]
    cat(sprintf(
      "Fold %d: %d obs (%.1f%%), %.3f prevalence%s\n",
      f, fold_data$n_total, fold_data$n_total/sum(fold_stats$n_total)*100,
      fold_data$prevalence,
      ifelse(f == test_fold, " [TEST]", "")
    ))
  }

  # Create outer (test vs train) and inner (CV folds) partitions
  outer_partitions <- ifelse(fold_assignments == test_fold, "test", "train")

  # For inner partitions, renumber non-test folds as CV folds
  inner_partitions <- rep(NA, length(fold_assignments))
  train_mask <- outer_partitions == "train"
  if (sum(train_mask) > 0) {
    train_folds <- fold_assignments[train_mask]
    unique_train_folds <- sort(unique(train_folds))
    fold_mapping <- setNames(1:length(unique_train_folds), unique_train_folds)
    inner_partitions[train_mask] <- fold_mapping[as.character(train_folds)]
  }

  # Calculate prevalence balance metric (analogous to silhouette for clustering)
  prevalence_balance <- 1 - sd(fold_stats$prevalence) / mean(fold_stats$prevalence)

  return(list(
    outer_partitions = outer_partitions,
    inner_partitions = inner_partitions,
    folds = fold_assignments,
    fold_stats = fold_stats,
    test_fold_id = test_fold,
    prevalence_balance = prevalence_balance,
    spatial_cv = spatial_cv
  ))
}

# Evaluate spatial blocking for a given block size and return W_CV statistic
evaluate_blocksize <- function(block_size_km, data, data_with_coords,
                               precomputed_distances, top_feature,
                               featuredist_sample_size = 1e4,
                               seed = 42) {

  cat("Evaluating block size =", block_size_km, "km...\n")

  # Convert km to meters
  block_size_m <- block_size_km * 1000

  # Create spatial partitions
  partition_result <- create_spatial_blocking_partitions(
    data_with_coords = data_with_coords,
    block_size_m = block_size_m,
    target_col = "response",
    n_partitions = 5,  # Keep 5 partitions (1 test + 4 CV)
    seed = seed
  )

  # Prepare data for CV distance calculation (all data including test partition)
  # Create combined partition assignments: test partition gets its own ID
  all_partitions <- partition_result$outer_partitions
  combined_folds <- ifelse(all_partitions == "test",
                          0,  # Assign test partition to fold 0
                          partition_result$inner_partitions)

  # Calculate CV distances for this block size across all outer partitions
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

  # Get prevalence balance from partition result
  prevalence_balance <- partition_result$prevalence_balance

  # Print results for this block size
  cat("block_size =", block_size_km, "km -> W_CV =", round(W_CV, 4),
      ", prevalence_balance =", round(prevalence_balance, 3), "\n\n")

  return(list(
    block_size_km = block_size_km,
    partition_result = partition_result,
    W_CV = W_CV,
    prevalence_balance = prevalence_balance,
    top_feature = top_feature,
    feature_dists = feature_dists
  ))
}

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

## Systematic evaluation of block sizes using W_CV minimization ####

cat("Using top feature for distance calculation:", top_feature, "\n")
cat("Evaluating block sizes (km):", paste(block_sizes_km, collapse = ", "), "\n")

# Pre-compute prediction distances (independent of block size)
cat("Pre-computing prediction distances...\n")
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
cat("Starting systematic evaluation of block sizes...\n")
blocking_evaluation <- tibble(block_size_km = block_sizes_km) |>
  mutate(
    results = map(block_size_km, ~evaluate_blocksize(
      block_size_km = .x,
      data = mf_current,
      data_with_coords = mf_current_with_coords,
      precomputed_distances = precomputed_distances,
      top_feature = top_feature,
      featuredist_sample_size = featuredist_sample_size,
      seed = 42
    )),
    W_CV = map_dbl(results, ~.x$W_CV),
    prevalence_balance = map_dbl(results, ~.x$prevalence_balance),
    top_feature = map_chr(results, ~.x$top_feature)
  )

# Display evaluation results
cat("\nBlock size evaluation results:\n")
evaluation_summary <- blocking_evaluation |>
  select(block_size_km, W_CV, prevalence_balance, top_feature) |>
  arrange(W_CV)
print(evaluation_summary)

# Visualization of evaluation results
p1 <- ggplot(blocking_evaluation, aes(x = block_size_km, y = W_CV)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 2) +
  geom_vline(xintercept = blocking_evaluation$block_size_km[which.min(blocking_evaluation$W_CV)],
             linetype = "dashed", color = "red", alpha = 0.7) +
  labs(
    title = "W_CV vs Spatial Block Size",
    subtitle = paste("Optimal block size =", blocking_evaluation$block_size_km[which.min(blocking_evaluation$W_CV)], "km"),
    x = "Spatial block size (km)",
    y = "W_CV (scaled feature distance)"
  ) +
  theme_minimal()

print(p1)

# Plot distance distributions for first few block sizes
for (i in 1:min(3, length(blocking_evaluation$results))) {
  cat("\nBlock size:", blocking_evaluation$block_size_km[i], "km\n")

  p <- blocking_evaluation$results[[i]]$feature_dists |>
    ggplot(aes(x = dist, color = what)) +
    stat_ecdf() +
    coord_cartesian(xlim = c(0, 20)) +
    labs(
      title = paste("Distance distributions - Block size", blocking_evaluation$block_size_km[i], "km"),
      x = "Feature distance",
      y = "ECDF"
    ) +
    theme_minimal()

  print(p)

  median_dists <- blocking_evaluation$results[[i]]$feature_dists |>
    group_by(what) |>
    summarise(median_dist = median(dist))
  print(median_dists)
}

# Find optimal block size (minimum W_CV)
optimal_block_size <- blocking_evaluation$block_size_km[which.min(blocking_evaluation$W_CV)]
cat("\nOptimal block size (minimum W_CV):", optimal_block_size, "km\n")

dists_optimal <- blocking_evaluation$results[[which.min(blocking_evaluation$W_CV)]]$feature_dists
plot(dists_optimal, stat = "ecdf")
dists_optimal |> 
  group_by(what) |> 
  summarise(dist = median(dist))
dists_optimal |> 
  ggplot(aes(x = dist, color = what)) +
  stat_ecdf() + 
  coord_cartesian(xlim = c(0, 10)) 

# Save evaluation results
write_csv(blocking_evaluation |> select(-results),
          "output/pl2/block_size_evaluation.csv", append = FALSE)

# Extract optimal partitioning result
optimal_result <- blocking_evaluation |>
  filter(block_size_km == optimal_block_size) |>
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
    fold = factor(partitioning_result$folds),
    outer = factor(partitioning_result$outer_partitions),
    inner = partitioning_result$inner_partitions
  ) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(raster_current)) |>
  vect()

raster <- rasterize(
  mf_current_with_coords_partitioned,
  raster_current,
  field = c("fold", "outer", "inner")
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

write_csv(mf_partitioned, "output/pl2/modeling_frame_regional_partitioned_spatial.csv",
          append = FALSE)

writeRaster(raster, "output/pl2/partition_spatial.tif", overwrite = TRUE)

# sessionInfo ####

sessioninfo::session_info()