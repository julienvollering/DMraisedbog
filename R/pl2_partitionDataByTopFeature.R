# Top-Feature Partitioning ####

# This script implements quantile cut-based partitioning using the most important
# feature. Optimizes n_cuts × segregation_prob to minimize W_CV.

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
n_cores <- 4
plan(multisession, workers = n_cores)

# Feature distance calculation subsampling
featuredist_sample_size <- 3e4

# Evaluation ranges for 2D grid search (configurable)
n_cuts_range <- seq(1000, 2000, by = 500)
segregation_prob_range <- seq(0.5, 0.75, by = 0.05)

## Read modeling frame ####
mf <- read_csv("output/pl2/modeling_frame_regional.csv")

raster_current <- rast("output/pl2/scenario_current.tif")

mf_current_with_coords <- mf |>
  filter(scenario == "current") |>
  select(-scenario)
mf_current <- mf_current_with_coords |>
  select(-x, -y)

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

cat("Training data - Presences:", n_presence, "Absences:", n_absence, "\n")
cat("Prevalence:", round(prevalence * 100, 2), "%\n")

## Load feature weights and identify most important feature ####
weights_features <- read_csv("output/pl2/weights_feature_data_partitioning.csv")
weighting_method <- "Balanced Random Forest"

weights <- weights_features |>
  filter(method == weighting_method) |>
  select(feature, median) |>
  arrange(desc(median)) |>
  tibble::deframe()

most_important_feature <- names(weights)[1]

cat("\nMost important feature:", most_important_feature, "\n")

## Functions ####
source("R/functions.R")

## Pre-compute prediction distances on most important feature ####

cat("\nPre-computing prediction distances on", most_important_feature, "...\n")

precomputed_distances <- calculate_prediction_distances(
  training_data = mf_current,
  prediction_data = mf_future_presence,
  variable_name = most_important_feature,
  sample_size = featuredist_sample_size,
  standardize = FALSE,
  seed = 42
)

cat("Pre-computation complete.\n\n")

## 2D grid search over n_cuts and segregation_prob ####

cat("Starting parallel 2D grid search over n_cuts × segregation_prob...\n")
cat("Using", n_cores, "cores\n")

# Create single-feature weights vector for compatibility
feature_weights <- setNames(
  weights[most_important_feature],
  most_important_feature
)

# Check for existing evaluation results
eval_file <- "output/pl2/n_cuts_evaluation_topfeature.csv"
if (file.exists(eval_file)) {
  cat("Loading existing evaluation results from:", eval_file, "\n")
  existing_results <- read_csv(eval_file, show_col_types = FALSE)
  cat("Found", nrow(existing_results), "previous evaluations\n")
} else {
  existing_results <- tibble(
    n_cuts = numeric(),
    segregation_prob = numeric(),
    W_CV = numeric()
  )
  cat("No existing results found, starting fresh\n")
}

# Create grid of all combinations to evaluate
all_combinations <- expand_grid(
  n_cuts = n_cuts_range,
  segregation_prob = segregation_prob_range
)

# Identify which combinations still need evaluation
combinations_to_run <- all_combinations |>
  anti_join(existing_results, by = c("n_cuts", "segregation_prob"))

cat("Parameter grid:", nrow(all_combinations), "total combinations\n")
cat("Already evaluated:", nrow(existing_results), "combinations\n")
cat("To evaluate:", nrow(combinations_to_run), "combinations\n\n")

if (nrow(combinations_to_run) > 0) {
  # Run only unevaluated combinations
  new_evaluations <- combinations_to_run |>
    mutate(
      results = future_map2(
        n_cuts,
        segregation_prob,
        ~ evaluate_cut_partitioning(
          n_cuts = .x,
          segregation_prob = .y,
          data = mf_current,
          precomputed_distances = precomputed_distances,
          feature_weights = feature_weights,
          top_n_features = 1,
          featuredist_sample_size = featuredist_sample_size,
          seed = 42
        ),
        .options = furrr_options(seed = TRUE)
      ),
      W_CV = map_dbl(results, ~ .x$W_CV)
    )

  # Save new results immediately (without results column)
  new_results_summary <- new_evaluations |>
    select(n_cuts, segregation_prob, W_CV)

  write_csv(new_results_summary, eval_file, append = file.exists(eval_file))
  cat("Saved", nrow(new_results_summary), "new evaluations to", eval_file, "\n")

  # Combine all results for analysis
  partitioning_evaluation <- new_evaluations
  all_results_summary <- bind_rows(existing_results, new_results_summary)
} else {
  cat("All combinations already evaluated, skipping grid search\n")
  partitioning_evaluation <- tibble(
    n_cuts = numeric(),
    segregation_prob = numeric(),
    results = list(),
    W_CV = numeric()
  )
  all_results_summary <- existing_results
}

## Display and visualize results ####

cat("\nTop 10 parameter combinations by W_CV (across all runs):\n")
evaluation_summary <- all_results_summary |>
  arrange(W_CV) |>
  head(10)
print(evaluation_summary)

# Heatmap
p1 <- ggplot(
  all_results_summary,
  aes(x = n_cuts, y = segregation_prob, fill = W_CV)
) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  geom_point(
    data = all_results_summary |> slice_min(W_CV, n = 1),
    color = "red",
    size = 4,
    shape = 21,
    stroke = 2
  ) +
  labs(
    title = paste(
      "W_CV across n_cuts × segregation_prob (Most Important:",
      most_important_feature,
      ")"
    ),
    subtitle = "Red point = optimal combination",
    x = "Number of quantile cuts (n_cuts)",
    y = "Segregation probability",
    fill = "W_CV"
  ) +
  theme_minimal()

print(p1)

# Line plot
p2 <- ggplot(
  all_results_summary,
  aes(x = segregation_prob, y = W_CV, color = factor(n_cuts))
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = paste(
      "W_CV vs Segregation Probability (Most Important:",
      most_important_feature,
      ")"
    ),
    x = "Segregation probability",
    y = "W_CV",
    color = "n_cuts"
  ) +
  theme_minimal()

print(p2)

## Extract optimal results ####

optimal_params <- all_results_summary |> slice_min(W_CV, n = 1)
optimal_n_cuts <- optimal_params$n_cuts
optimal_segregation_prob <- optimal_params$segregation_prob
optimal_w_cv <- optimal_params$W_CV

cat("\nOptimal parameters (minimum W_CV across all runs):\n")
cat("  n_cuts =", optimal_n_cuts, "\n")
cat("  segregation_prob =", optimal_segregation_prob, "\n")
cat("  W_CV =", round(optimal_w_cv, 4), "\n")

# Re-run optimal combination to get full results
cat("\nRe-running optimal combination to generate partitions...\n")
optimal_result <- evaluate_cut_partitioning(
  n_cuts = optimal_n_cuts,
  segregation_prob = optimal_segregation_prob,
  data = mf_current,
  precomputed_distances = precomputed_distances,
  feature_weights = feature_weights,
  top_n_features = 1,
  featuredist_sample_size = featuredist_sample_size,
  seed = 42
)

# ECDF plot
plot(optimal_result$feature_dists, stat = "ecdf")

## Save results ####

partitioning_result <- optimal_result$partition_result

cat("\nPartition summary:\n")
print(table(partitioning_result$partitions))

## Map visualization ####
mf_current_with_coords_partitioned <- mf_current_with_coords |>
  select(response, x, y) |>
  mutate(
    interval = factor(partitioning_result$clusters),
    partition = factor(partitioning_result$partitions)
  ) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(raster_current)) |>
  vect()

raster <- rasterize(
  mf_current_with_coords_partitioned,
  raster_current,
  field = c("interval", "partition")
)

plot(raster["partition"])

## Save partitioned datasets ####

mf_partitioned <- mf_current_with_coords |>
  mutate(
    scenario = "current",
    partition = partitioning_result$partitions,
    .before = 1
  ) |>
  bind_rows(filter(mf, scenario != "current")) |>
  select(scenario, partition, response, everything())

write_csv(
  mf_partitioned,
  "output/pl2/modeling_frame_regional_partitioned_topfeature.csv",
  append = FALSE
)

writeRaster(raster, "output/pl2/partition_topfeature.tif", overwrite = TRUE)

# sessionInfo ####

sessioninfo::session_info()
