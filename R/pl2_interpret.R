# Spatial Interpretation of Current and Future Predictions ####

# This script reads prediction rasters for current and future scenarios and
# summarizes them spatially across Lyngstad survey polygons (raised bogs).
# Analyzes threshold crossing patterns using training data prevalence threshold.

library(readr)
library(dplyr)
library(tidyr)
library(sf)
library(terra)

## Configuration ####

# Select which partitioning schemes to process
partitioning_schemes <- c("importance", "pca", "maxshift")

## Load spatial data ####

# Read Lyngstad raised bog polygons
cat("Loading Lyngstad raised bog polygons...\n")
lyngstad <- st_read("data/DMraisedbog.gpkg", layer = "lyngstad-MTYPE_A", quiet = TRUE)
cat("Loaded", nrow(lyngstad), "polygons\n")

# Transform to EPSG:3035 to match modeling CRS
lyngstad_proj <- st_transform(lyngstad, crs = 3035)

# Convert to SpatVector for terra operations
lyngstad_vect <- vect(lyngstad_proj)

## Load prevalence threshold ####

# Read modeling frame to calculate prevalence (same as pl2_predict.R)
mf <- read_csv("output/pl2/modeling_frame_regional.csv", show_col_types = FALSE)

train_data <- mf |>
  filter(scenario == "current") |>
  select(response)

prevalence <- mean(train_data$response == 1)
cat("Training data prevalence (threshold):", round(prevalence, 4), "\n\n")

## Main loop over partitioning schemes ####

all_data <- list()

# scheme <- partitioning_schemes[1]  # For testing

for (scheme in partitioning_schemes) {
  cat("Processing scheme:", toupper(scheme), "\n")

  ### Load predictions ####

  pred_file <- paste0("output/pl2/rf_local_pred_", scheme, ".tif")
  pred_rasters <- rast(pred_file)

  ### Extract zonal statistics for both scenarios ####

  zonal <- terra::extract(
    pred_rasters,
    lyngstad_vect,
    fun = mean) |>
    as_tibble()

  ### Analyze change ####
  
  change <- zonal |>
    mutate(
      change = future - current,
      positive_current = current >= prevalence,
      positive_future = future >= prevalence,
      transition = case_when(
        positive_current & positive_future ~ "remain_positive",
        !positive_current & !positive_future ~ "remain_negative",
        !positive_current & positive_future ~ "become_positive",
        positive_current & !positive_future ~ "become_negative"
      )
    ) |> 
    select(ID, current, future, change, transition)

  # Store data for combined CSV
  all_data[[scheme]] <- change

  ### Create output polygon layer ####

  # Add to polygon layer
  lyngstad_predicted <- lyngstad_proj
  lyngstad_predicted <- bind_cols(lyngstad_predicted, change) |> 
    mutate(scheme = scheme, .before = current)

  # Save polygon layer
  st_write(
    lyngstad_predicted,
    dsn = "output/pl2/lyngstad_predictions.gpkg",
    layer = scheme,
    delete_layer = TRUE,
    quiet = TRUE
  )

  ### Print summary statistics ####

  n_polygons <- nrow(change)
  transition_counts <- change |>
    count(transition) |>
    mutate(prop = n / n_polygons)

  print(transition_counts)
}

# sessionInfo ####

sessioninfo::session_info()
