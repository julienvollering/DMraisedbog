# Spatial Interpretation of Current and Future Predictions ####

# This script reads prediction rasters for current and future scenarios and
# summarizes them spatially across Lyngstad survey polygons (raised bogs).
# Analyzes threshold crossing patterns using training data prevalence threshold.

library(readr)
library(dplyr)
library(tidyr)
library(sf)
library(terra)

## Load spatial data ####

cat("Loading Lyngstad raised bog polygons...\n")
lyngstad <- st_read(
  "data/DMraisedbog.gpkg",
  layer = "lyngstad-MTYPE_A",
  quiet = TRUE
)
cat("Loaded", nrow(lyngstad), "polygons\n")

lyngstad_proj <- st_transform(lyngstad, crs = 3035)
lyngstad_vect <- vect(lyngstad_proj)

## Load prevalence threshold ####

mf <- read_csv("output/pl2/modeling_frame_regional.csv", show_col_types = FALSE)

train_data <- mf |>
  filter(scenario == "current") |>
  select(response)

prevalence <- mean(train_data$response == 1)
cat("Training data prevalence (threshold):", round(prevalence, 4), "\n\n")

## Main loop over models ####

models <- c("brf", "rfq")
all_data <- list()

for (model in models) {
  cat("Processing model:", toupper(model), "\n")

  ### Load predictions ####

  pred_file <- paste0("output/pl2/rf_local_pred_", model, ".tif")
  pred_rasters <- rast(pred_file)

  ### Extract zonal statistics for both scenarios ####

  zonal <- terra::extract(
    pred_rasters,
    lyngstad_vect,
    fun = mean
  ) |>
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

  all_data[[model]] <- change

  ### Create output polygon layer ####

  lyngstad_predicted <- lyngstad_proj
  lyngstad_predicted <- bind_cols(lyngstad_predicted, change) |>
    mutate(model = model, .before = current)

  st_write(
    lyngstad_predicted,
    dsn = "output/pl2/lyngstad_predictions.gpkg",
    layer = model,
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
