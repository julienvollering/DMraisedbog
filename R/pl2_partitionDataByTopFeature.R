# Top-Feature Partitioning ####

# This script implements presence-based partitioning using the most important
# feature. Sorts presences along the feature and divides into k groups.
# Absences are assigned based on overlap with partition boundaries.
# Absences in gaps between partitions are dropped.

library(readr)
library(dplyr)
library(terra)
library(ggplot2)
library(sf)

source("R/functions.R")

## Configuration ####

# Number of partitions
k_partitions <- 5

# Minimum presences per partition
min_presences <- 50

# Random seed for reproducibility
partition_seed <- 42

## Read modeling frame ####
mf <- read_csv("output/pl2/modeling_frame_regional.csv")

raster_current <- rast("output/pl2/scenario_current.tif")

mf_current_with_coords <- mf |>
  filter(scenario == "current") |>
  select(-scenario)

mf_current <- mf_current_with_coords |>
  select(-x, -y)

n_presence <- sum(mf_current$response == 1)
n_absence <- sum(mf_current$response == 0)
prevalence <- n_presence / (n_presence + n_absence)

cat("\n=== DATA SUMMARY ===\n")
cat("Training data - Presences:", n_presence, "Absences:", n_absence, "\n")
cat("Prevalence:", round(prevalence * 100, 2), "%\n")
cat("Requested partitions:", k_partitions, "\n")
cat("Minimum presences per partition:", min_presences, "\n")
cat("Total presences needed:", k_partitions * min_presences, "\n")

## Load feature weights and identify most important feature ####
weights_features <- read_csv("output/pl2/weights_feature_data_partitioning.csv")
weighting_method <- "Balanced Random Forest"

weights <- weights_features |>
  filter(method == weighting_method) |>
  select(feature, median) |>
  arrange(desc(median)) |>
  tibble::deframe()

most_important_feature <- names(weights)[1]

cat("\n=== TOP FEATURE ===\n")
cat("Most important feature:", most_important_feature, "\n")
cat("Feature importance:", round(weights[most_important_feature], 4), "\n")

## Apply presence-based partitioning ####

cat("\n=== PARTITIONING ===\n")
cat("Running presence-based partitioning algorithm...\n")

partitioning_result <- partition_by_presence_sorting(
  data = mf_current,
  feature_name = most_important_feature,
  k = k_partitions,
  min_pres = min_presences,
  seed = partition_seed
)

# Display results
cat("\nPartition summary:\n")
cat("  Partitions created:", partitioning_result$k, "\n")
cat("  Feature used:", partitioning_result$feature_name, "\n")
cat("  Observations dropped (in gaps):", partitioning_result$n_dropped, "\n\n")

cat("Presences per partition:\n")
print(partitioning_result$n_presences)

cat("\nAbsences per partition:\n")
print(partitioning_result$n_absences)

cat("\nTotal observations per partition:\n")
partition_totals <- partitioning_result$n_presences +
  partitioning_result$n_absences
names(partition_totals) <- paste0("P", 1:k_partitions)
print(partition_totals)

## Feature distribution visualization ####

cat("\n=== VISUALIZATION ===\n")

# Histogram of feature with partition boundaries
partition_boundaries <- data.frame(
  partition = 1:k_partitions,
  x_lower = sapply(partitioning_result$partition_list, function(p) p$x_lower),
  x_upper = sapply(partitioning_result$partition_list, function(p) p$x_upper)
)

cat("Partition boundaries on feature", most_important_feature, ":\n")
print(partition_boundaries)

# Plot feature distribution by partition
mf_current_plot <- mf_current |>
  mutate(partition = factor(partitioning_result$partitions))

p_dist <- ggplot(
  mf_current_plot |> filter(!is.na(partition)),
  aes(x = .data[[most_important_feature]], fill = partition)
) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  facet_wrap(~partition, ncol = 1) +
  labs(
    title = paste("Distribution of", most_important_feature, "by Partition"),
    x = most_important_feature,
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_dist)

# Boxplot by partition and response
p_box <- ggplot(
  mf_current_plot |> filter(!is.na(partition)),
  aes(
    x = partition,
    y = .data[[most_important_feature]],
    fill = factor(response)
  )
) +
  geom_boxplot() +
  labs(
    title = paste(most_important_feature, "by Partition and Response"),
    x = "Partition",
    y = most_important_feature,
    fill = "Response"
  ) +
  theme_minimal()

print(p_box)

## Map visualization ####

cat("\nCreating spatial visualization...\n")

mf_current_with_coords_partitioned <- mf_current_with_coords |>
  select(response, x, y) |>
  mutate(
    partition = factor(partitioning_result$partitions)
  ) |>
  filter(!is.na(partition)) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(raster_current)) |>
  vect()

raster <- rasterize(
  mf_current_with_coords_partitioned,
  raster_current,
  field = "partition"
)

plot(
  raster["partition"],
  main = paste(
    "Spatial Distribution of Partitions\n(Feature:",
    most_important_feature,
    ")"
  ),
  col = rainbow(k_partitions)
)

## Save partitioned datasets ####

cat("\n=== SAVING RESULTS ===\n")

mf_partitioned <- mf_current_with_coords |>
  mutate(
    scenario = "current",
    partition = partitioning_result$partitions,
    .before = 1
  ) |>
  bind_rows(filter(mf, scenario != "current")) |>
  select(scenario, partition, response, everything())

output_csv <- "output/pl2/modeling_frame_regional_partitioned_topfeature.csv"
output_tif <- "output/pl2/partition_topfeature.tif"

write_csv(mf_partitioned, output_csv, append = FALSE)
cat("Saved partitioned modeling frame to:", output_csv, "\n")

writeRaster(raster, output_tif, overwrite = TRUE)
cat("Saved partition raster to:", output_tif, "\n")

cat("\nPartitioning complete!\n")
cat("  Total observations:", nrow(mf_current), "\n")
cat(
  "  Observations in partitions:",
  sum(!is.na(partitioning_result$partitions)),
  "\n"
)
cat("  Observations dropped:", partitioning_result$n_dropped, "\n")
cat(
  "  Drop rate:",
  round(100 * partitioning_result$n_dropped / nrow(mf_current), 2),
  "%\n"
)

# sessionInfo ####

sessioninfo::session_info()
