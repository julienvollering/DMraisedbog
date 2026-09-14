# Build the predictor stacks ####

# PURPOSE: Assembles the 5 km EU and 250 m Norway predictor stacks (CHELSA, paleo, DTM50 terrain, AR50 land cover) on a common EPSG:3035 grid.

library(terra)
library(sf)
library(rnaturalearth)
library(tidyverse)

# Clean up old terra temporary files from previous runs
terra::tmpFiles(remove = TRUE)

# Threshold variables (gdd*, gst, swe, and gsp once its sentinel is masked) are truncated
# to positive values, so a cold cell where the quantity is really 0 arrives as NA. This
# sets NA to 0 only where every OTHER layer of the stack has data, so a genuine no-data
# cell stays NA. One implementation for the global stack and both regional scenarios;
# it reports per-layer coverage before the fix and the cells changed per variable.
fill_threshold_na <- function(s, vars, label) {
  cat("\nChecking coverage before correction (", label, "):\n", sep = "")
  counts <- global(s, "notNA")$notNA
  ref <- max(counts)
  cat("  Reference coverage (maximum across layers):", ref, "cells\n")
  for (i in which(counts < ref)) {
    cat(sprintf(
      "    - %s: %d cells (%d cells missing)\n", names(s)[i], counts[i], ref - counts[i]
    ))
  }
  present <- intersect(vars, names(s))
  if (length(present) == 0) {
    return(s)
  }
  all_valid <- all(!is.na(s[[setdiff(names(s), vars)]]))
  cat("  Converting NA to 0 where every non-threshold layer has data:\n")
  for (v in present) {
    fix <- is.na(s[[v]]) & all_valid
    n_changed <- global(fix, "sum", na.rm = TRUE)$sum
    if (n_changed > 0) {
      s[[v]] <- ifel(fix, 0, s[[v]])
    }
    cat(sprintf("    - %s: %d cells\n", v, n_changed))
  }
  s
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
# CHELSA CMIP6 scenario; must match files downloaded by pl0_downloadCHELSA.R.
# Filtering by scenario means leftover files from other scenarios (e.g. ssp585)
# in the same directory will not be stacked.
future_scenario <- "ssp370"
chelsa_future_files <- list.files(
  "data/CHELSA/2071-2100",
  pattern = paste0("CHELSA_gfdl-esm4_", future_scenario, "_.*\\.tif$"),
  full.names = TRUE
) %>%
  sort()

chelsa_future_stack <- rast(chelsa_future_files)
names(chelsa_future_stack) <- stringr::str_extract(
  names(chelsa_future_stack),
  paste0("(?<=CHELSA_gfdl-esm4_", future_scenario, "_).+(?=_2071)")
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

# A layer whose future median is more than five times its current one (or less than a
# fifth) has almost certainly been stored with a different scale factor between the two
# CHELSA releases. None has so far; stop rather than stack it if one ever does.
if (length(checkvars) > 0) {
  stop("Possible scaling mismatch between current and future CHELSA layers: ",
       paste(checkvars, collapse = ", "))
}
cat("Current and future CHELSA medians agree to within a factor of 5 for every layer\n")

### Unflagged no-data in CHELSA integer layers ####

# gsp is stored as an unsigned 32-bit integer with a 0.1 scale factor and NO NoData tag,
# so where the growing season has zero length the sentinel 4294967295 is read as a real
# value and scaled to 4.29e8 mm. The gdd/gst/swe layers carry a NoData tag and arrive as
# NA, which the NA-to-zero step below handles; gsp did not, and the sentinel reached the
# modelling frame in 8,987 training rows (notebook 2026-09-12). It is masked HERE, on the
# cropped native grid and before project(), because bilinear resampling blends a sentinel
# into its neighbours and those blends cannot be recognised afterwards. Once NA, gsp goes
# through the same NA-to-zero step as gdd/gst/swe: no growing season, no growing-season
# precipitation.
SENTINEL_VARS <- c("gsp")
MAX_PLAUSIBLE <- c(gsp = 1e5) # mm; the wettest cells in the frame are ~7,000

mask_sentinels <- function(s) {
  for (nm in intersect(SENTINEL_VARS, names(s))) {
    s[[nm]] <- classify(s[[nm]], cbind(MAX_PLAUSIBLE[[nm]], Inf, NA))
  }
  s
}

# Fail loudly if any CHELSA layer still carries a value no climate variable can take, so
# a sentinel in another layer or another CHELSA release cannot slip through silently again.
assert_plausible <- function(s, label, bound = 1e5) {
  mx <- global(s, "max", na.rm = TRUE)$max
  names(mx) <- names(s)
  bad <- names(mx)[is.finite(mx) & mx > bound]
  if (length(bad) > 0) {
    stop(
      label, ": implausible maxima in ", paste(bad, collapse = ", "),
      " (", paste(signif(mx[bad], 3), collapse = ", "), ")"
    )
  }
  cat(label, "- layer maxima all below", bound, "\n")
  invisible(mx)
}

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

# Crop to European extent first (in WGS84) to reduce memory requirements, then mask the
# unflagged no-data sentinel before any resampling (see above).
chelsa_cropped_wgs84 <- crop(chelsa_past_stack, europe_ext_wgs84) |>
  mask_sentinels()
assert_plausible(chelsa_cropped_wgs84, "CHELSA current, Europe crop")

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

# See fill_threshold_na() at the top. gsp is in the list because its no-data sentinel
# was turned into NA above (no growing season, no growing-season precipitation).
THRESHOLD_VARS_GLOBAL <- c("gdd10", "gst", "swe", "gsp")
chelsa_masked <- fill_threshold_na(chelsa_masked, THRESHOLD_VARS_GLOBAL, "global model")

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

# Crop to Norway extent first (in WGS84) to reduce memory requirements, then mask the
# unflagged no-data sentinel before any resampling (see the global model section).
chelsa_past_cropped_wgs84 <- crop(chelsa_past_stack, norway_ext_wgs84) |>
  mask_sentinels()
assert_plausible(chelsa_past_cropped_wgs84, "CHELSA current, Norway crop")

# Then transform to EPSG:3035 and resample to 250m
chelsa_past_3035 <- project(
  x = chelsa_past_cropped_wgs84,
  y = template_250m,
  method = "bilinear"
)

# Final mask to Norway boundaries using raster mask
chelsa_past_masked <- mask(chelsa_past_3035, land_mask_ar50)

# Crop to Norway extent first (in WGS84) to reduce memory requirements, then mask the
# unflagged no-data sentinel before any resampling.
chelsa_future_cropped_wgs84 <- crop(chelsa_future_stack, norway_ext_wgs84) |>
  mask_sentinels()
assert_plausible(chelsa_future_cropped_wgs84, "CHELSA future, Norway crop")

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

# See fill_threshold_na() at the top. gdd5 is in this list and NOT in the global one; the
# asymmetry is inherited from the two separate implementations this replaces and is kept
# as it was rather than resolved here.
THRESHOLD_VARS_REGIONAL <- c("gdd10", "gdd5", "gst", "swe", "gsp")
predictors_current <- fill_threshold_na(
  predictors_current, THRESHOLD_VARS_REGIONAL, "regional model - current"
)
predictors_future <- fill_threshold_na(
  predictors_future, THRESHOLD_VARS_REGIONAL, "regional model - future"
)

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
