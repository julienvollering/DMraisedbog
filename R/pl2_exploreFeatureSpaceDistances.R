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

## Configuration ####

# Parallelization
n_cores <- 4
plan(multisession, workers = n_cores)

# Feature distance calculation subsampling
featuredist_sample_size <- 3e4

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
cat("Prediction locations (future @ presences):", nrow(mf_future_presence), "\n")

## Load feature weights ####
weights_features <- read_csv("output/pl2/weights_feature_data_partitioning.csv")

weights <- weights_features |>
  select(feature, median) |>
  arrange(desc(median)) |>
  tibble::deframe()

## Load functions ####
source("R/functions.R")

## Calculate W_sample for all features ####

feature_names <- names(weights)

cat("Calculating W_sample for", length(feature_names), "features using", n_cores, "cores...\n")

# Calculate feature distances and W_sample for each feature (parallelized)
w_sample_results <- future_map_dfr(feature_names, function(feat) {
  # Calculate distances using featuredist function
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

  # Extract W_sample and calculate summary statistics
  s2s_dists <- dists[dists$what == "sample-to-sample", "dist"][[1]]
  p2s_dists <- dists[dists$what == "prediction-to-sample", "dist"][[1]]

  tibble(
    feature = feat,
    W_sample = attr(dists, "W_sample"),
    s2s_median = median(s2s_dists, na.rm = TRUE),
    p2s_median = median(p2s_dists, na.rm = TRUE)
  )
}, .options = furrr_options(seed = TRUE))

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
write_csv(w_sample_results, "output/pl2/w_sample_by_feature.csv", append = FALSE)

## Visualizations ####

# W_sample vs Importance scatter plot
p <- ggplot(w_sample_results, aes(x = importance, y = W_sample)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(aes(label = feature), hjust = -0.1, vjust = 0.5, size = 3, check_overlap = TRUE) +
  labs(
    title = "Feature Importance vs Distribution Shift",
    x = "Feature Importance",
    y = "W_sample"
  ) +
  theme_minimal()

print(p)

# ECDF plots for most shifted and top 5 most important features
most_shifted_feature <- w_sample_results |>
  slice_max(W_sample, n = 1) |>
  pull(feature)

top_5_features <- w_sample_results |>
  slice_max(importance, n = 5) |>
  pull(feature)

selected_features <- unique(c(most_shifted_feature, top_5_features))

cat("\nGenerating ECDF plots for:", paste(selected_features, collapse = ", "), "\n")

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
      subtitle = paste0("W_sample = ", w_val, ", Importance = ", imp_val),
      x = "Distance",
      y = "Cumulative Probability",
      color = "Distance Type"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
})

walk(distance_plots, print)

## PCA analysis on top 5 features ####

cat("\n=== PCA Analysis ===\n")
cat("Performing PCA on top 5 most important features...\n")

# Fit PCA on current training data
pca_data_current <- mf_current |>
  select(all_of(top_5_features))

pca_result <- prcomp(pca_data_current, center = TRUE, scale. = TRUE)

cat("\nPCA summary:\n")
print(summary(pca_result))

# Add PC1 and PC2 to current data
mf_current_pca <- mf_current |>
  mutate(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2]
  )

# Project future data into PC space
pca_data_future <- mf_future_presence |>
  select(all_of(top_5_features))

future_pcs <- predict(pca_result, newdata = pca_data_future)[, 1:2]
mf_future_presence_pca <- mf_future_presence |>
  mutate(
    PC1 = future_pcs[, 1],
    PC2 = future_pcs[, 2]
  )

# Visualize PCA
pca_df <- as.data.frame(pca_result$x[, 1:2]) |>
  mutate(response = as.factor(mf_current$response))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = response)) +
  geom_point(alpha = 0.5) +
  labs(
    title = "PCA of Top 5 Features (Current)",
    x = "PC1",
    y = "PC2",
    color = "Response"
  ) +
  theme_minimal()

print(p_pca)

# Combined current + future
pca_df_future <- mf_future_presence_pca |>
  select(PC1, PC2, response) |>
  mutate(response = as.factor(response))

pca_combined <- bind_rows(
  current = pca_df,
  future = pca_df_future,
  .id = "scenario"
)

p_pca_combined <- ggplot(pca_combined, aes(x = PC1, y = PC2, color = response, shape = scenario)) +
  geom_point(alpha = 0.5) +
  labs(
    title = "PCA of Top 5 Features (Current + Future)",
    x = "PC1",
    y = "PC2",
    color = "Response",
    shape = "Scenario"
  ) +
  theme_minimal()

print(p_pca_combined)

# Calculate feature distances for PC1
cat("\nCalculating feature distances for PC1...\n")
dists_pc1 <- featuredist(
  training_data = mf_current_pca,
  prediction_data = mf_future_presence_pca,
  variable_name = "PC1",
  cvfold_column = NULL,
  sample_size = featuredist_sample_size,
  seed = seed,
  standardize = FALSE,
  scale_distances = FALSE
)

w_sample_pc1 <- attr(dists_pc1, "W_sample")

cat("\nPC1 W_sample:", round(w_sample_pc1, 4), "\n")

# ECDF plot for PC1
p_pc1 <- ggplot(dists_pc1, aes(x = dist, color = what)) +
  stat_ecdf(linewidth = 1) +
  coord_cartesian(xlim = c(0, quantile(dists_pc1$dist, 0.99, na.rm = TRUE))) +
  labs(
    title = "PC1 Distance Distributions",
    subtitle = paste0("W_sample = ", round(w_sample_pc1, 4)),
    x = "Distance",
    y = "Cumulative Probability",
    color = "Distance Type"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_pc1)

# Save PCA results for use in partitioning
saveRDS(list(
  pca_result = pca_result,
  top_5_features = top_5_features,
  w_sample_pc1 = w_sample_pc1
), "output/pl2/pca_results.rds")

# sessionInfo ####

sessioninfo::session_info()
