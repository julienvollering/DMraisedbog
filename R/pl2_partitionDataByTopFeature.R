# Top-Feature Partitioning ####

# PURPOSE: Cuts the frame into k CV partitions by sorting bog presences along the top materially-shifting feature (bio10), assigning other classes by range overlap.

# This script implements presence-based partitioning using the most important
# feature. Sorts presences along the feature and divides into k groups.
# Absences are assigned by a complete tiling of the feature axis: interior cuts at the
# midpoint between adjacent presence groups, terminal intervals open to +/- Inf. Nothing
# is dropped. Each row also carries `envelope_side` -- its position relative to the
# presence envelope (`inside` / `below` / `above`) -- so skill can be reported both
# inside-envelope and tail-inclusive without refitting anything.

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

# The class that anchors the partitions. With the 3-class response, only bog rows sort
# the feature axis; nonpeat and otherpeat rows are assigned by range overlap.
presence_level <- "bog"

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

n_presence <- sum(mf_current$response == presence_level)
prevalence <- n_presence / nrow(mf_current)

cat("\n=== DATA SUMMARY ===\n")
cat("Training rows:", nrow(mf_current), "\n")
print(table(mf_current$dataset, mf_current$response))
cat("Prevalence of", presence_level, ":", round(prevalence * 100, 2), "%\n")
cat("Requested partitions:", k_partitions, "\n")
cat("Minimum presences per partition:", min_presences, "\n")
cat("Total presences needed:", k_partitions * min_presences, "\n")

## Load feature weights and identify most important feature ####
weights_features <- read_csv("output/pl2/weights_feature_data_partitioning.csv")
weighting_method <- "Balanced Random Forest"

# Restricted to predictors that MOVE MATERIALLY under the scenario (`material`, written
# by pl2_weightFeaturesDataPartitioning.R: projection-dynamic AND a mean shift of at least
# 0.5 SD of the predictor's own spread). The partitions exist to make the pairwise CV an
# extrapolation test along the axis the projection actually travels, and two filters are
# needed to get there:
#
#  - the unrestricted ruler's top feature is `slope`, static terrain that is identical in
#    the current and future frames, so folds cut on it differ in terrain and not climate;
#  - the top merely-dynamic feature is `bio02`, which changes in nearly every cell but
#    shifts only 0.26 SD, so folds cut on it separate the data along a gradient the
#    projection barely traverses.
#
# The bar (17 of 32 dynamic features clear it) leaves `bio10` as the top-VI candidate,
# which is also the cut axis of the section 1.3 sweep -- so pl2's CV and pl3's sweep sit
# on the same gradient rather than on two unrelated ones.
weights <- weights_features |>
  filter(method == weighting_method, material) |>
  select(feature, median) |>
  arrange(desc(median)) |>
  tibble::deframe()

most_important_feature <- names(weights)[1]

cat("\n=== TOP FEATURE ===\n")
cat("Most important materially-shifting feature:", most_important_feature, "\n")
cat("Feature importance:", round(weights[most_important_feature], 4), "\n")
cat("Candidates considered:", length(weights), "of", nrow(filter(weights_features, method == weighting_method)), "\n")

## Apply presence-based partitioning ####

cat("\n=== PARTITIONING ===\n")
cat("Running presence-based partitioning algorithm...\n")

partitioning_result <- partition_by_presence_sorting(
  data = mf_current,
  feature_name = most_important_feature,
  k = k_partitions,
  min_pres = min_presences,
  seed = partition_seed,
  presence_level = presence_level
)

# Display results
cat("\nPartition summary:\n")
cat("  Partitions created:", partitioning_result$k, "\n")
cat("  Feature used:", partitioning_result$feature_name, "\n\n")

# The tiling is complete, so every row lands in exactly one partition unless the feature
# itself is NA. A non-zero count here is therefore a real signal, and stops the run.
stopifnot(partitioning_result$n_dropped == 0)

cat("Rows per partition and class:\n")
print(partitioning_result$n_by_class)

# Every partition has to carry all three classes, or a fold trained on it silently
# becomes a two-class problem and the contrast the partition exists to test is absent.
if (any(partitioning_result$n_by_class == 0)) {
  cat("\nWARNING: at least one partition x class cell is empty (see table above)\n")
}

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

# Boxplot by partition and class
p_box <- ggplot(
  mf_current_plot |> filter(!is.na(partition)),
  aes(
    x = partition,
    y = .data[[most_important_feature]],
    fill = response
  )
) +
  geom_boxplot() +
  labs(
    title = paste(most_important_feature, "by Partition and Response"),
    x = "Partition",
    y = most_important_feature,
    fill = "Class"
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
    envelope_side = partitioning_result$envelope_side,
    .before = 1
  ) |>
  bind_rows(filter(mf, scenario != "current")) |>
  select(scenario, partition, envelope_side, response, everything())

output_csv <- "output/pl2/modeling_frame_regional_partitioned_topfeature.csv"

write_csv(mf_partitioned, output_csv, append = FALSE)
cat("Saved partitioned modeling frame to:", output_csv, "\n")
cat("  Rows partitioned:", sum(!is.na(partitioning_result$partitions)),
    "of", nrow(mf_current), "\n")

# Row accounting: the partitioned frame is the modelling frame with two columns added --
# no row gained, none lost -- and every current row carries a partition.
stopifnot(
  nrow(mf_partitioned) == nrow(mf),
  sum(!is.na(partitioning_result$partitions)) == nrow(mf_current)
)

# sessionInfo ####

sessioninfo::session_info()
