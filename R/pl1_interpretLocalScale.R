# Interpret Local Scale Model: Gain/Loss of Suitability in Lyngstad Polygons ####

# This script analyzes how Lyngstad polygons gain or lose binary suitability
# under future climate scenarios using the local scale model predictions.
# Uses the 95% sensitivity threshold from evaluateLocalScale.R to convert
# probability predictions to binary suitability classifications.

# For each Lyngstad polygon, calculates:
# - Complete loss: all cells within polygon become unsuitable in future
# - Partial loss: at least one cell within polygon becomes unsuitable in future
# - Complete gain: all cells within polygon become suitable in future
# - Partial gain: at least one cell within polygon becomes suitable in future

# Results are summarized by polygon count and area to assess landscape-scale
# impacts of climate change on raised bog suitability.

library(readr)
library(dplyr)
library(tidyr)
library(sf)
library(terra)

## Step 1: Load data and models ####

# Load existing local model predictions
current_predictions <- rast("output/pl1/rf_local_pred_current_all_layers.tif")
future_predictions <- rast("output/pl1/rf_local_pred_future_all_layers.tif")

crs_common <- st_crs(future_predictions)

# Load Lyngstad polygons
lyngstad_polygons <- st_read("data/DMraisedbog.gpkg", layer = "lyngstad-MTYPE_A") |>
  st_transform(crs_common)

# Note: Calibration already applied in the prediction files
# Band 1: p_raw (raw probabilities)
# Band 2: p_cal (calibrated probabilities)
# Band 3: class_raw (raw binary classification)
# Band 4: class_cal (calibrated binary classification)

cat("Data loaded successfully.\n")
cat("Lyngstad polygons:", nrow(lyngstad_polygons), "\n")
cat("Current predictions:", nlyr(current_predictions), "layers\n")
cat("Future predictions:", nlyr(future_predictions), "layers\n")

## Step 2: Extract 95% sensitivity threshold ####

cat("\nExtracting 95% sensitivity threshold from evaluation results...\n")

# Read performance comparison to get the threshold
performance_data <- read_csv("output/pl1/test_performance_comparison.csv")

# Extract 95% sensitivity threshold
threshold_95sens <- performance_data |>
  filter(threshold_type == "95_sensitivity") |>
  pull(threshold) |>
  unique()

if(length(threshold_95sens) != 1) {
  stop("Could not find unique 95% sensitivity threshold in performance data")
}

cat("95% sensitivity threshold:", threshold_95sens, "\n")

## Step 3: Extract calibrated probabilities by polygon using terra::extract ####

cat("\nExtracting calibrated probabilities by polygon...\n")

# Extract calibrated probabilities (Band 2) for current and future scenarios
cat("Extracting current scenario probabilities...\n")
current_extract <- terra::extract(
  current_predictions[[2]], lyngstad_polygons,
  cells = TRUE, xy = TRUE) |> 
  drop_na(p_cal)

cat("Extracting future scenario probabilities...\n")
future_extract <- terra::extract(
  future_predictions[[2]], lyngstad_polygons,
  cells = TRUE, xy = TRUE) |> 
  drop_na(p_cal)

cat("Current extractions:", nrow(current_extract), "values\n")
cat("Future extractions:", nrow(future_extract), "values\n")

## Step 4: Calculate suitability changes by polygon ####

cat("\nCalculating suitability changes by polygon...\n")

# Function to analyze polygon changes using extracted values 

# current_data <- current_extract; future_data <- future_extract; polygons <- lyngstad_polygons; threshold <- threshold_95sens
analyze_polygon_changes <- function(current_data, future_data, polygons, threshold) {

  # Apply threshold to create binary suitability
  current_data$suitable <- current_data[["p_cal"]] >= threshold  # Band 2 values
  future_data$suitable <- future_data[["p_cal"]] >= threshold

  # Combine current and future data with scenario identifier
  combined_data <- bind_rows(
    current_data |> filter(!is.na(p_cal)) |> mutate(scenario = "current"),
    future_data |> filter(!is.na(p_cal)) |> mutate(scenario = "future")
  )

  # Calculate summary statistics by polygon
  results <- combined_data |>
    group_by(ID, scenario) |>
    summarise(
      n_cells = n(),
      n_suitable = sum(suitable, na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = scenario,
      values_from = c(n_cells, n_suitable)
    ) |>
    mutate(
      # Use current n_cells (should be same for both scenarios)
      n_cells = coalesce(n_cells_current, n_cells_future),
      # Calculate change categories
      unchanged = n_suitable_current == n_suitable_future,
      complete_loss = n_suitable_current == n_cells & n_suitable_future == 0,
      complete_gain = n_suitable_current == 0 & n_suitable_future == n_cells,
      partial_loss = !complete_loss & n_suitable_current > n_suitable_future,
      partial_gain = !complete_gain & n_suitable_current < n_suitable_future,
    ) |>
    select(
      id_polygon = ID,
      n_cells,
      n_suitable_current,
      n_suitable_future,
      unchanged,
      complete_loss,
      partial_loss,
      complete_gain,
      partial_gain
    ) |>
    # Check exactly one change category is TRUE per row
    rowwise() |>
    mutate(n_categories = sum(c(unchanged, complete_loss, partial_loss, complete_gain, partial_gain))) |>
    ungroup()

  if (any(results$n_categories != 1)) {
    stop("Invalid change categories detected")
  }

  results <- results |> select(-n_categories)

  return(results)
}

# Analyze polygon changes using extracted data
polygon_changes <- analyze_polygon_changes(
  current_extract,
  future_extract,
  lyngstad_polygons,
  threshold_95sens
)

cat("Analysis completed for", nrow(polygon_changes), "polygons with prediction data\n")

## Step 5: Generate summary statistics ####

cat("\nGenerating summary statistics...\n")

polygon_summary <- polygon_changes |> 
  count(unchanged, complete_loss, partial_loss, complete_gain, partial_gain) |> 
  mutate(change_type = case_when(
    unchanged ~ "Unchanged",
    complete_loss ~ "Complete Loss",
    partial_loss ~ "Partial Loss",
    complete_gain ~ "Complete Gain",
    partial_gain ~ "Partial Gain",
    TRUE ~ "Other"
  )) |>
  select(change_type, n) |>
  mutate(percent = round(n / sum(n) * 100, 2))
  
polygon_changes |> 
  group_by(unchanged, complete_loss, partial_loss, complete_gain, partial_gain) |>
  count()

## Step 6: Print results ####

cat("\n=== SUMMARY: LYNGSTAD POLYGON SUITABILITY CHANGES ===\n")
cat("Threshold used: 95% sensitivity =", threshold_95sens, "\n\n")

cat("POLYGON COUNT SUMMARY:\n")
print(polygon_summary)

## Step 7: Parallel analysis with AOA masking ####

cat("\n=== PARALLEL ANALYSIS: AOA-MASKED SUITABILITY CHANGES ===\n")

# Load Area of Applicability (AOA) mask
aoa_mask <- rast("output/pl1/aoa_future_local_scale.tif")

cat("AOA mask loaded successfully.\n")
cat("AOA mask dimensions:", dim(aoa_mask), "\n")

# Lyngstad polygons within AOA
lyngstad_polygons_aoa <- terra::extract(aoa_mask, lyngstad_polygons) |> 
  group_by(ID) |> 
  summarize(within_aoa = all(AOA == 1))
lyngstad_polygons_aoa <- lyngstad_polygons |>
  bind_cols(lyngstad_polygons_aoa)

lyngstad_polygons_aoa |> 
  st_drop_geometry() |> 
  group_by(within_aoa) |>
  count()

# Mask the prediction rasters to AOA (not strictly necessary since lyngstad_polygons_aoa)
current_predictions_aoa <- mask(current_predictions, aoa_mask, maskvalues = 0)
future_predictions_aoa <- mask(future_predictions, aoa_mask, maskvalues = 0)

# Extract calibrated probabilities using AOA-masked rasters
cat("Extracting current scenario probabilities (AOA-masked)...\n")
current_extract_aoa <- terra::extract(
  current_predictions_aoa[[2]], 
  filter(lyngstad_polygons_aoa, within_aoa == TRUE),
  cells = TRUE, xy = TRUE) |> 
  drop_na(p_cal) # should be zero rows if within_aoa == TRUE

cat("Extracting future scenario probabilities (AOA-masked)...\n")
future_extract_aoa <- terra::extract(
  future_predictions_aoa[[2]],
  filter(lyngstad_polygons_aoa, within_aoa == TRUE),
  cells = TRUE, xy = TRUE) |> 
  drop_na(p_cal) # should be zero rows if within_aoa == TRUE

cat("Current AOA extractions:", nrow(current_extract_aoa), "values\n")
cat("Future AOA extractions:", nrow(future_extract_aoa), "values\n")

# Analyze polygon changes using AOA-masked data
cat("Calculating AOA-masked suitability changes by polygon...\n")
polygon_changes_aoa <- analyze_polygon_changes(
  current_extract_aoa,
  future_extract_aoa,
  lyngstad_polygons,
  threshold_95sens
)

cat("AOA analysis completed for", nrow(polygon_changes_aoa), "polygons inside AOA\n")

# Generate AOA summary statistics
cat("Generating AOA summary statistics...\n")

polygon_summary_aoa <- polygon_changes_aoa |>
  count(unchanged, complete_loss, partial_loss, complete_gain, partial_gain) |>
  mutate(change_type = case_when(
    unchanged ~ "Unchanged",
    complete_loss ~ "Complete Loss",
    partial_loss ~ "Partial Loss",
    complete_gain ~ "Complete Gain",
    partial_gain ~ "Partial Gain",
    TRUE ~ "Other"
  )) |>
  select(change_type, n) |>
  mutate(percent = round(n / sum(n) * 100, 2))

# Print AOA results
cat("\n=== AOA-MASKED RESULTS ===\n")
cat("Threshold used: 95% sensitivity =", threshold_95sens, "\n\n")

cat("AOA-MASKED POLYGON COUNT SUMMARY:\n")
print(polygon_summary_aoa)

## Step 8: Compare full vs AOA-masked results ####

cat("\n=== COMPARISON: FULL vs AOA-MASKED ANALYSIS ===\n")

# Create comparison summary
comparison_summary <- bind_rows(
  polygon_summary |> mutate(analysis = "Full"),
  polygon_summary_aoa |> mutate(analysis = "AOAmasked")
) |>
  select(analysis, change_type, n, percent) |>
  pivot_wider(
    names_from = analysis,
    values_from = c(n, percent),
    names_sep = "_"
  )

cat("COMPARISON TABLE (Full vs AOA-masked):\n")
print(comparison_summary)

# sessionInfo ####

sessioninfo::session_info()
