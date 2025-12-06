# Explore Feature Space Distances Between Current and Future Conditions ####

# This script calculates W_sample for each feature to quantify distributional
# shifts between current and future conditions.
# W_sample compares sample-to-sample distances vs prediction-to-sample distances.

library(readr)
library(dplyr)
library(purrr)
library(furrr)
library(ggplot2)
library(twosamples)

source("R/functions.R")

## Configuration ####

# Parallelization
n_cores <- 4
plan(multisession, workers = n_cores)

# Feature distance calculation subsampling
featuredist_sample_size <- 1e4

# Random seed for reproducibility
seed <- 42

## Read modeling frame ####
mf <- read_csv("output/pl2/modeling_frame_regional.csv")

mf_current <- mf |>
  filter(scenario == "current") |>
  select(-scenario, -x, -y)

# Get prediction domain: future conditions at presence locations
# These are the locations we care most about for our research questions
toJoin <- mf |>
  filter(scenario == "future") |>
  select(-scenario)

mf_future_presence <- mf |>
  filter(response == 1) |>
  select(x, y) |>
  left_join(toJoin, by = c("x", "y")) |>
  select(-x, -y)

cat("Training data:", nrow(mf_current), "observations\n")
cat(
  "Prediction locations (future @ presences):",
  nrow(mf_future_presence),
  "\n"
)

## Load and select feature weights ####
weights_features <- read_csv("output/pl2/weights_feature_data_partitioning.csv")

weighting_method <- "Balanced Random Forest"

weights <- weights_features |>
  filter(method == weighting_method) |>
  select(feature, median) |>
  arrange(desc(median)) |>
  tibble::deframe()

## Calculate W_sample for all features ####

feature_names <- names(weights)

cat(
  "Calculating W_sample for",
  length(feature_names),
  "features using",
  n_cores,
  "cores...\n"
)

# Calculate feature distances and W_sample for each feature (parallelized)
w_sample_results <- future_map_dfr(
  feature_names,
  function(feat) {
    # Calculate distances using featuredist function
    dists <- featuredist(
      training_data = mf_current,
      prediction_data = mf_future_presence,
      variable_name = feat,
      cvfold_column = NULL,
      sample_size = featuredist_sample_size,
      seed = seed,
      standardize = TRUE,
      scale_distances = FALSE
    )

    # Extract W_sample and calculate summary statistics
    s2s_dists <- dists[dists$what == "sample-to-sample", "dist"][[1]]
    p2s_dists <- dists[dists$what == "prediction-to-sample", "dist"][[1]]

    tibble(
      feature = feat,
      W_sample = attr(dists, "W_sample"),
      s2s_median = median(s2s_dists, na.rm = TRUE),
      p2s_median = median(p2s_dists, na.rm = TRUE)
    )
  },
  .options = furrr_options(seed = TRUE)
)

# Add feature importance weights to results
w_sample_results <- w_sample_results |>
  left_join(
    tibble(feature = names(weights), importance = unname(weights)),
    by = "feature"
  ) |>
  arrange(desc(W_sample))

cat("\nW_sample results:\n")
print(w_sample_results, n = Inf)

## Save results ####
write_csv(
  w_sample_results,
  "output/pl2/Wsample_by_feature.csv",
  append = FALSE
)

## Visualizations ####

# W_sample vs Importance scatter plot
p <- ggplot(w_sample_results, aes(x = importance, y = W_sample)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(
    aes(label = feature),
    hjust = -0.1,
    vjust = 0.5,
    size = 3,
    check_overlap = TRUE
  ) +
  labs(
    title = "Feature Importance vs Distribution Shift",
    x = "Feature Importance",
    y = "W_sample (standardized)"
  ) +
  theme_minimal()

print(p)

# ECDF plots for *most shifted* and *most important* features
most_shifted_feature <- w_sample_results |>
  slice_max(W_sample, n = 1) |>
  pull(feature)

top_feature <- w_sample_results |>
  slice_max(importance, n = 1) |>
  pull(feature)

selected_features <- unique(c(most_shifted_feature, top_feature))

distance_plots <- map(selected_features, function(feat) {
  dists <- featuredist(
    training_data = mf_current,
    prediction_data = mf_future_presence,
    variable_name = feat,
    cvfold_column = NULL,
    sample_size = featuredist_sample_size,
    seed = seed,
    standardize = FALSE,
    scale_distances = FALSE
  )

  w_val <- round(attr(dists, "W_sample"), 4)
  imp_val <- round(weights[feat], 4)

  ggplot(dists, aes(x = dist, color = what)) +
    stat_ecdf(linewidth = 1) +
    coord_cartesian(xlim = c(0, quantile(dists$dist, 0.99, na.rm = TRUE))) +
    labs(
      title = feat,
      x = "Distance (unstandardized)",
      y = "Cumulative Probability",
      color = "Distance Type"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
})

walk(distance_plots, print)

# sessionInfo ####

sessioninfo::session_info()
