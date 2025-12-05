library(terra)
library(sf)
library(rnaturalearth)
library(tidyverse)

# Clean up old terra temporary files from previous runs
terra::tmpFiles(remove = TRUE)

# Function to get memory usage of all objects
memory_usage <- function() {
  obj_names <- ls(envir = .GlobalEnv)
  obj_sizes <- sapply(obj_names, function(x) {
    object.size(get(x, envir = .GlobalEnv))
  })
  # Convert to MB and sort
  obj_sizes_mb <- round(obj_sizes / 1024^2, 2)
  sort(obj_sizes_mb, decreasing = TRUE)
}

# Global model ####

# Plan for Predictor Stack Preparation for Global (EU) Model
#
# Objective: Prepare CHELSA bioclim+ and Trace21k variables for nested SDM global model
# - Input: CHELSA variables at 30 arcsec (~1km) resolution in WGS84
# - Output: Single GeoTIFF with bands for all predictors at 5000m resolution in EPSG:3035, masked to EU+Norway

## Load CHELSA data ####

### Past climate ####
chelsa_past_files <- list.files(
  "data/CHELSA/1981-2010",
  pattern = "CHELSA_.*\\.tif$",
  full.names = TRUE
) %>%
  sort()

chelsa_past_stack <- rast(chelsa_past_files)
names(chelsa_past_stack) <- stringr::str_extract(
  names(chelsa_past_stack),
  "(?<=CHELSA_).+(?=_1981)"
)

### Future climate ####
chelsa_future_files <- list.files(
  "data/CHELSA/2071-2100",
  pattern = "CHELSA_.*\\.tif$",
  full.names = TRUE
) %>%
  sort()

chelsa_future_stack <- rast(chelsa_future_files)
names(chelsa_future_stack) <- stringr::str_extract(
  names(chelsa_future_stack),
  "(?<=CHELSA_gfdl-esm4_ssp585_).+(?=_2071)"
)

### Paleo-derived predictors ####

# Load paleo-derived predictor files (excluding n_consecutive_icefree which is correlated)
paleo_files <- list.files(
  "data/CHELSA/paleo_derived",
  pattern = "^paleo_.*\\.tif$",
  full.names = TRUE
) %>%
  # Exclude the correlated variable
  setdiff(., grep("n_consecutive", ., value = TRUE)) %>%
  sort()

cat("Loading", length(paleo_files), "paleo-derived predictors:\n")
cat(paste(basename(paleo_files), collapse = "\n"), "\n")

paleo_stack <- rast(paleo_files)

# Extract clean layer names (remove _EUextent_EPSG4326.tif suffix)
paleo_names <- basename(paleo_files) %>%
  str_remove("_EUextent_EPSG4326\\.tif$")
names(paleo_stack) <- paleo_names

### Check scaling across variables ####

# Compare scaling between past and future climate stacks
past_stats <- chelsa_past_stack[1e7 + 1:1e6] |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  group_by(variable) |>
  summarise(median = median(value, na.rm = TRUE), .groups = "drop")
future_stats <- chelsa_future_stack[1e7 + 1:1e6] |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  group_by(variable) |>
  summarise(median = median(value, na.rm = TRUE), .groups = "drop")
checkvars <- left_join(
  past_stats,
  future_stats,
  by = "variable",
  suffix = c("_past", "_future")
) |>
  mutate(median_ratio = abs(median_future / median_past)) |>
  filter(median_ratio > 5 | median_ratio < 0.2) |> # Identify variables with possible scaling differences (factor >5)
  pull(variable)

# Approximate Norway bounding box in WGS84
norway_extent <- ext(c(
  xmin = 4,
  xmax = 32,
  ymin = 57,
  ymax = 72
))

plot(chelsa_past_stack[[checkvars]], ext = norway_extent)
plot(chelsa_future_stack[[checkvars]], ext = norway_extent)
# No variables appear to have different scaling

## Create EU + Norway mask ####

# Get European countries
europe_countries <- ne_countries(
  continent = "europe",
  scale = 10,
  returnclass = "sf"
)

# Get EU member states (approximate list - may need updating)
eu_countries <- c(
  "Austria",
  "Belgium",
  "Bulgaria",
  "Croatia",
  "Cyprus",
  "Czechia",
  "Denmark",
  "Estonia",
  "Finland",
  "France",
  "Germany",
  "Greece",
  "Hungary",
  "Ireland",
  "Italy",
  "Latvia",
  "Lithuania",
  "Luxembourg",
  "Malta",
  "Netherlands",
  "Poland",
  "Portugal",
  "Romania",
  "Slovakia",
  "Slovenia",
  "Spain",
  "Sweden"
)

# Filter to EU countries + Norway
target_countries <- europe_countries %>%
  filter(name %in% c(eu_countries, "Norway"))

# Transform to EPSG:3035 (European Equal Area)
target_countries <- st_transform(target_countries, crs = "EPSG:3035")

# Create combined polygon for masking
mask_polygon <- target_countries %>%
  st_union() %>%
  st_sf()

## Create template grid and project to EPSG:3035 ####

# Define European extent in EPSG:3035 coordinates
# Approximate bounding box covering EU + Norway in EPSG:3035
europe_ext_3035 <- ext(c(
  xmin = 2000000,
  xmax = 7000000,
  ymin = 1000000,
  ymax = 5500000
))

# Create 5000m resolution template grid in EPSG:3035
template_5km <- rast(europe_ext_3035, resolution = 5000, crs = "EPSG:3035")

# Define European extent in WGS84 for cropping before reprojection
# Approximate bounding box covering EU + Norway in WGS84
europe_ext_wgs84 <- ext(c(xmin = -10, xmax = 35, ymin = 35, ymax = 72))

# Crop to European extent first (in WGS84) to reduce memory requirements
chelsa_cropped_wgs84 <- crop(chelsa_past_stack, europe_ext_wgs84)

# Project CHELSA data from WGS84 to EPSG:3035 at 5000m resolution
chelsa_3035 <- project(
  x = chelsa_cropped_wgs84,
  y = template_5km,
  method = "bilinear"
)

## Mask to EU + Norway ####

# Mask to country boundaries (removes sea and non-target countries)
chelsa_masked <- mask(chelsa_3035, mask_polygon)

## Process paleo-derived predictors ####

# Crop to European extent (in WGS84) to reduce memory requirements
paleo_cropped_wgs84 <- crop(paleo_stack, europe_ext_wgs84)

# Project to EPSG:3035 at 5000m resolution
paleo_3035 <- project(
  x = paleo_cropped_wgs84,
  y = template_5km,
  method = "bilinear"
)

# Mask to country boundaries
paleo_masked <- mask(paleo_3035, mask_polygon)

# Combine CHELSA and paleo-derived predictors
chelsa_masked <- c(chelsa_masked, paleo_masked)

cat("Combined predictor stack has", nlyr(chelsa_masked), "layers\n")

## Fix threshold-based variables with truncated ranges ####

# Pre-correction coverage check
cat("\nChecking coverage before correction (global model):\n")
all_counts <- global(chelsa_masked, "notNA")
reference_count <- max(all_counts$notNA)
cat("  Reference coverage (maximum across layers):", reference_count, "cells\n")

reduced_coverage <- all_counts$notNA < reference_count
if (any(reduced_coverage)) {
  cat("  Variables with reduced coverage compared to maximum:\n")
  for (i in which(reduced_coverage)) {
    missing <- reference_count - all_counts$notNA[i]
    cat(sprintf(
      "    - %s: %d cells (%d cells missing)\n",
      names(chelsa_masked)[i],
      all_counts$notNA[i],
      missing
    ))
  }
}

# Variables like gdd5, gdd10, swe, gst are truncated to positive values,
# creating NAs in cold regions where they should logically be 0.
# Convert NA to 0 only in cells where all other layers have valid data.

threshold_vars <- c("gdd10", "gst", "swe")
vars_present <- names(chelsa_masked)[names(chelsa_masked) %in% threshold_vars]

if (length(vars_present) > 0) {
  # Create reference mask: cells where ALL non-threshold layers have data
  non_threshold_idx <- which(!names(chelsa_masked) %in% threshold_vars)
  reference_layers <- chelsa_masked[[non_threshold_idx]]
  all_valid_mask <- all(!is.na(reference_layers))

  # Convert NA to 0 for threshold variables where reference layers have data
  cat(
    "\nConverting NA to 0 for threshold variables where reference layers have data...\n"
  )
  for (var in vars_present) {
    var_layer <- chelsa_masked[[var]]
    na_mask <- is.na(var_layer)
    cells_to_fix <- na_mask & all_valid_mask
    n_changed <- global(cells_to_fix, "sum", na.rm = TRUE)$sum

    if (n_changed > 0) {
      chelsa_masked[[var]] <- ifel(cells_to_fix, 0, var_layer)
      cat(sprintf("  - %s: %d cells converted from NA to 0\n", var, n_changed))
    } else {
      cat(sprintf("  - %s: no cells needed conversion\n", var))
    }
  }
  cat("\n")
}

## Write output ####
writeRaster(
  chelsa_masked,
  filename = "output/predictors_global_5km_EUNorway_EPSG3035.tif",
  overwrite = TRUE,
  names = names(chelsa_masked)
)

# Quick visualization check
if (interactive()) {
  plot(
    chelsa_masked[[c(
      which(names(chelsa_masked) == "bio01"),
      which(names(chelsa_masked) == "paleo_years_icefreeland")
    )]]
  )
}

# Regional model ####

# Plan for Predictor Stack Preparation for Regional (Norway) Model
#
# Objective: Prepare predictors for nested SDM regional model at 250m resolution
# - CHELSA bioclim+ variables (current + future scenarios)
# - Paleo-derived predictors (Trace21k-derived variables, same for current + future)
# - Terrain variables (elevation, slope) from DTM50
# - AR50 land cover classes (artype: 10,20,30,50,60,70,81)
# - Output: EPSG:3035, 250m resolution, mainland Norway extent

## Create template grid in EPSG:3035 ####

# List all DTM tiles to define Norway extent
dtm_files <- list.files(
  "data/DTM50_UTM33_20250613",
  pattern = "\\.tif$",
  full.names = TRUE
)

# Load and mosaic DTM tiles in original UTM33 CRS
dtm_tiles <- map(dtm_files, rast)
dtm_mosaic_utm33 <- do.call(mosaic, dtm_tiles)

# Aggregate DTM from 50m to 250m in original UTM33 CRS
# This avoids distortion from resampling before derivative calculation
dtm_250m_utm33 <- aggregate(dtm_mosaic_utm33, fact = 5, fun = "mean")

# Define terrain-derived variables in extensible structure
# Each entry: name, terra::terrain parameter (v), unit (if applicable)
# To add new variables: append to this list (e.g., aspect, TRI, TPI)
terrain_variables <- list(
  list(name = "elevation", type = "base", source = dtm_250m_utm33),
  list(name = "slope", type = "terrain", v = "slope", unit = "degrees")
)

# Calculate all terrain variables at 250m resolution in UTM33
terrain_layers_utm33 <- list()
for (var_def in terrain_variables) {
  if (var_def$type == "base") {
    # Base DTM layer (elevation)
    layer <- var_def$source
    names(layer) <- var_def$name
  } else if (var_def$type == "terrain") {
    # Derived using terra::terrain()
    layer <- terrain(dtm_250m_utm33, v = var_def$v, unit = var_def$unit)
    names(layer) <- var_def$name
  }
  terrain_layers_utm33[[var_def$name]] <- layer
  cat("Calculated", var_def$name, "at 250m in UTM33\n")
}

# Create terrain stack in UTM33
terrain_stack_utm33 <- rast(terrain_layers_utm33)

# Get Norway extent in UTM33, then transform to EPSG:3035
norway_ext_utm33 <- ext(dtm_mosaic_utm33)
norway_ext_3035 <- project(
  norway_ext_utm33,
  from = crs(dtm_mosaic_utm33),
  to = "EPSG:3035"
)

# Create 250m template using same origin as global model (template_5km)
# Get origin from global template (5km grid)
global_origin_x <- xmin(template_5km)
global_origin_y <- ymin(template_5km)

# Calculate grid-aligned extent using the same origin
# Align Norway extent to 250m grid with same origin as 5km grid
xmin_aligned <- global_origin_x +
  floor((norway_ext_3035$xmin - global_origin_x) / 250) * 250
xmax_aligned <- global_origin_x +
  ceiling((norway_ext_3035$xmax - global_origin_x) / 250) * 250
ymin_aligned <- global_origin_y +
  floor((norway_ext_3035$ymin - global_origin_y) / 250) * 250
ymax_aligned <- global_origin_y +
  ceiling((norway_ext_3035$ymax - global_origin_y) / 250) * 250

norway_ext_aligned <- ext(c(
  xmin = xmin_aligned,
  xmax = xmax_aligned,
  ymin = ymin_aligned,
  ymax = ymax_aligned
))

# Create 250m resolution template grid with aligned extent and shared origin
template_250m <- rast(norway_ext_aligned, resolution = 250, crs = "EPSG:3035")

## Process AR50 land cover ####

# Load artype 60 (raised bogs) raster created in QGIS (50m resolution, 0/1 values)
ar50_50m <- rast("output/ar50_artype_layers/ar50_50m_EPSG3035_artype60.tif")

# Aggregate to 250m using mean (5x5 cells -> continuous [0,1] values)
ar50_250m <- aggregate(ar50_50m, fact = 5, fun = "mean", na.rm = FALSE)

# Resample to ensure exact alignment with template_250m
ar50_250m_stack <- resample(ar50_250m, template_250m, method = "bilinear")
names(ar50_250m_stack) <- "artype_60"

### Read Norway mainland land mask from disk ####
land_mask_ar50 <- rast("output/ar50_250m_land_EPSG3035.tif")

# Apply land mask to AR50 stack: keep 0's on land, set non-land to NA
ar50_250m_stack <- mask(ar50_250m_stack, land_mask_ar50)

# Write multiband AR50 output
writeRaster(
  ar50_250m_stack,
  filename = "output/ar50_250m_cover_EPSG3035.tif",
  overwrite = TRUE,
  names = names(ar50_250m_stack)
)

cat("Processed AR50 artype 60 layer\n")

## Process terrain variables ####

# Project terrain stack from UTM33 to EPSG:3035
# All derivatives already calculated at 250m in original projection
terrain_stack_3035 <- project(
  x = terrain_stack_utm33,
  y = template_250m,
  method = "bilinear"
)

# Apply land mask to all terrain layers
terrain_stack_masked <- mask(terrain_stack_3035, land_mask_ar50)

# Final terrain stack with proper names
terrain_stack <- terrain_stack_masked
names(terrain_stack) <- names(terrain_stack_utm33)

# Clean up large objects
rm(
  dtm_mosaic_utm33,
  dtm_250m_utm33,
  terrain_layers_utm33,
  terrain_stack_utm33,
  terrain_stack_3035,
  terrain_stack_masked
)
gc()

## Process CHELSA variables ####

# Create Norway extent in WGS84 for cropping before reprojection
# Use raster extent instead of vector union for efficiency
norway_ext_wgs84 <- project(
  norway_ext_3035,
  from = "EPSG:3035",
  to = "EPSG:4326"
) # Transform to WGS84

# Crop to Norway extent first (in WGS84) to reduce memory requirements
chelsa_past_cropped_wgs84 <- crop(chelsa_past_stack, norway_ext_wgs84)

# Then transform to EPSG:3035 and resample to 250m
chelsa_past_3035 <- project(
  x = chelsa_past_cropped_wgs84,
  y = template_250m,
  method = "bilinear"
)

# Final mask to Norway boundaries using raster mask
chelsa_past_masked <- mask(chelsa_past_3035, land_mask_ar50)

# Crop to Norway extent first (in WGS84) to reduce memory requirements
chelsa_future_cropped_wgs84 <- crop(chelsa_future_stack, norway_ext_wgs84)

# Then transform to EPSG:3035 and resample to 250m
chelsa_future_3035 <- project(
  x = chelsa_future_cropped_wgs84,
  y = template_250m,
  method = "bilinear"
)

# Final mask to Norway boundaries using raster mask
chelsa_future_masked <- mask(chelsa_future_3035, land_mask_ar50)

## Process paleo-derived predictors for regional model ####

# Crop paleo stack to Norway extent (in WGS84) to reduce memory requirements
paleo_regional_cropped_wgs84 <- crop(paleo_stack, norway_ext_wgs84)

# Then transform to EPSG:3035 and resample to 250m
paleo_regional_3035 <- project(
  x = paleo_regional_cropped_wgs84,
  y = template_250m,
  method = "bilinear"
)

# Final mask to Norway boundaries using raster mask
paleo_regional_masked <- mask(paleo_regional_3035, land_mask_ar50)

cat(
  "Processed",
  nlyr(paleo_regional_masked),
  "paleo-derived predictors for regional model\n"
)

## Combine all predictors ####

# Combine all predictors for current conditions (19 CHELSA + 2 terrain + AR50 + 9 paleo)
predictors_current <- c(
  chelsa_past_masked,
  terrain_stack,
  ar50_250m_stack,
  paleo_regional_masked
)

# Combine predictors for future conditions (only CHELSA changes, paleo stays same)
predictors_future <- c(
  chelsa_future_masked,
  terrain_stack,
  ar50_250m_stack,
  paleo_regional_masked
)

stopifnot(nlyr(predictors_current) == nlyr(predictors_future))

## Fix threshold-based variables with truncated ranges ####

# Pre-correction coverage check - current scenario
cat("\nChecking coverage before correction (regional model - current):\n")
all_counts_current <- global(predictors_current, "notNA")
reference_count_current <- max(all_counts_current$notNA)
cat(
  "  Reference coverage (maximum across layers):",
  reference_count_current,
  "cells\n"
)

reduced_coverage_current <- all_counts_current$notNA < reference_count_current
if (any(reduced_coverage_current)) {
  cat("  Variables with reduced coverage compared to maximum:\n")
  for (i in which(reduced_coverage_current)) {
    missing <- reference_count_current - all_counts_current$notNA[i]
    cat(sprintf(
      "    - %s: %d cells (%d cells missing)\n",
      names(predictors_current)[i],
      all_counts_current$notNA[i],
      missing
    ))
  }
}

# Pre-correction coverage check - future scenario
cat("\nChecking coverage before correction (regional model - future):\n")
all_counts_future <- global(predictors_future, "notNA")
reference_count_future <- max(all_counts_future$notNA)
cat(
  "  Reference coverage (maximum across layers):",
  reference_count_future,
  "cells\n"
)

reduced_coverage_future <- all_counts_future$notNA < reference_count_future
if (any(reduced_coverage_future)) {
  cat("  Variables with reduced coverage comparied to maximum:\n")
  for (i in which(reduced_coverage_future)) {
    missing <- reference_count_future - all_counts_future$notNA[i]
    cat(sprintf(
      "    - %s: %d cells (%d cells missing)\n",
      names(predictors_future)[i],
      all_counts_future$notNA[i],
      missing
    ))
  }
}

# Variables like gdd5, gdd10, swe, gst are truncated to positive values,
# creating NAs in cold regions where they should logically be 0.
# Convert NA to 0 only in cells where all other layers have valid data.

threshold_vars <- c("gdd10", "gdd5", "gst", "swe")

# Process current scenario
vars_present_current <- names(predictors_current)[
  names(predictors_current) %in% threshold_vars
]

if (length(vars_present_current) > 0) {
  # Create reference mask: cells where ALL non-threshold layers have data
  non_threshold_idx <- which(!names(predictors_current) %in% threshold_vars)
  reference_layers <- predictors_current[[non_threshold_idx]]
  all_valid_mask <- all(!is.na(reference_layers))

  # Convert NA to 0 for threshold variables where reference layers have data
  cat(
    "\nConverting NA to 0 for threshold variables where reference layers have data (current)...\n"
  )
  for (var in vars_present_current) {
    var_layer <- predictors_current[[var]]
    na_mask <- is.na(var_layer)
    cells_to_fix <- na_mask & all_valid_mask
    n_changed <- global(cells_to_fix, "sum", na.rm = TRUE)$sum

    if (n_changed > 0) {
      predictors_current[[var]] <- ifel(cells_to_fix, 0, var_layer)
      cat(sprintf("  - %s: %d cells converted from NA to 0\n", var, n_changed))
    } else {
      cat(sprintf("  - %s: no cells needed conversion\n", var))
    }
  }
  cat("\n")
}

# Process future scenario
vars_present_future <- names(predictors_future)[
  names(predictors_future) %in% threshold_vars
]

if (length(vars_present_future) > 0) {
  # Create reference mask: cells where ALL non-threshold layers have data
  non_threshold_idx <- which(!names(predictors_future) %in% threshold_vars)
  reference_layers <- predictors_future[[non_threshold_idx]]
  all_valid_mask <- all(!is.na(reference_layers))

  # Convert NA to 0 for threshold variables where reference layers have data
  cat(
    "\nConverting NA to 0 for threshold variables where reference layers have data (future)...\n"
  )
  for (var in vars_present_future) {
    var_layer <- predictors_future[[var]]
    na_mask <- is.na(var_layer)
    cells_to_fix <- na_mask & all_valid_mask
    n_changed <- global(cells_to_fix, "sum", na.rm = TRUE)$sum

    if (n_changed > 0) {
      predictors_future[[var]] <- ifel(cells_to_fix, 0, var_layer)
      cat(sprintf("  - %s: %d cells converted from NA to 0\n", var, n_changed))
    } else {
      cat(sprintf("  - %s: no cells needed conversion\n", var))
    }
  }
  cat("\n")
}

## Write outputs ####
writeRaster(
  predictors_current,
  filename = "output/predictors_regional_250m_Norway_current_EPSG3035.tif",
  overwrite = TRUE,
  names = names(predictors_current)
)

writeRaster(
  predictors_future,
  filename = "output/predictors_regional_250m_Norway_future_EPSG3035.tif",
  overwrite = TRUE,
  names = names(predictors_future)
)

# Quick visualization check
if (interactive()) {
  # Plot bio1, elevation, and paleo predictor
  plot(
    predictors_current[[c(
      which(names(predictors_current) == "bio01"),
      which(names(predictors_current) == "elevation"),
      which(names(predictors_current) == "paleo_years_icefreeland")
    )]]
  )
}

## Spatial coverage validation ####

# Check spatial coverage across both current and future stacks
# Count non-NA cells for each layer
current_counts <- global(predictors_current, "notNA")
current_counts$layer <- names(predictors_current)

future_counts <- global(predictors_future, "notNA")
future_counts$layer <- names(predictors_future)

# Join the two tables
coverage_comparison <- left_join(
  current_counts,
  future_counts,
  by = "layer",
  suffix = c("_current", "_future")
) |>
  mutate(equal = notNA_current == notNA_future) |>
  arrange(notNA_current) |> # Sort by current coverage (ascending)
  select(layer, notNA_current, notNA_future, equal)

cat("\nSpatial coverage comparison (sorted by current coverage, ascending):\n")
print(coverage_comparison, row.names = FALSE)

# Check for inconsistencies within each stack
if (
  !all(
    coverage_comparison$notNA_current == coverage_comparison$notNA_current[1]
  )
) {
  cat("\nWARNING: Current predictors have inconsistent spatial coverage\n")
}
if (
  !all(coverage_comparison$notNA_future == coverage_comparison$notNA_future[1])
) {
  cat("\nWARNING: Future predictors have inconsistent spatial coverage\n")
}
if (!all(coverage_comparison$equal)) {
  cat("\nWARNING: Coverage differs between current and future scenarios\n")
}

# Clean up terra temporary files
terra::tmpFiles(remove = TRUE)

# sessionInfo ####

sessioninfo::session_info()
