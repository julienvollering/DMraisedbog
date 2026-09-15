# Build unified predictor stack (EU + Norway, 250m, EPSG:3035) ####

# PURPOSE: Assembles a single unified predictor stack at 250m resolution (EPSG:3035)
# covering EU + Norway, with CHELSA climate, paleo-derived, and terrain (elevation/slope)
# variables. This replaces the previous separate-stack approach (global 5km EU +
# regional 250m Norway), which created scale mismatches for slope and NA-to-zero
# inconsistencies.
#
# Elevation: elevatr for all EU + Norway (uniform quality, no separate DTM50 warp).
# Terrain: elevation + slope calculated at 250m in EPSG:3035.
# NA-to-zero: applied consistently to gdd5, gdd10, gst, swe, gsp for all domains.

library(terra)
library(sf)
library(rnaturalearth)
library(elevatr)
library(tidyverse)

source("R/config.R")

# Clean up old terra temporary files from previous runs
terra::tmpFiles(remove = TRUE)

# Threshold variables (gdd*, gst, swe, and gsp once its sentinel is masked) are truncated
# to positive values, so a cold cell where the quantity is really 0 arrives as NA. This
# sets NA to 0 only where every OTHER layer of the stack has data, so a genuine no-data
# cell stays NA. One implementation for all scenarios; it reports per-layer coverage before
# the fix and the cells changed per variable.
fill_threshold_na <- function(s, vars, label) {
  cat("\nChecking coverage before correction (", label, "):\n", sep = "")
  counts <- global(s, "notNA")$notNA
  ref <- max(counts)
  cat("  Reference coverage (maximum across layers):", ref, "cells\n")
  for (i in which(counts < ref)) {
    cat(sprintf(
      "    - %s: %d cells (%d cells missing)\n",
      names(s)[i],
      counts[i],
      ref - counts[i]
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

# Unflagged no-data in CHELSA integer layers ####

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
      label,
      ": implausible maxima in ",
      paste(bad, collapse = ", "),
      " (",
      paste(signif(mx[bad], 3), collapse = ", "),
      ")"
    )
  }
  cat(label, "- layer maxima all below", bound, "\n")
  invisible(mx)
}

# Load CHELSA and paleo data ####

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
future_scenario <- FUTURE_SCENARIO # R/config.R, shared with pl0_downloadCHELSA.R
chelsa_future_files <- list.files(
  "data/CHELSA/2071-2100",
  pattern = paste0(
    "CHELSA_",
    tolower(FUTURE_GCM),
    "_",
    future_scenario,
    "_.*\\.tif$"
  ),
  full.names = TRUE
) %>%
  sort()

chelsa_future_stack <- rast(chelsa_future_files)
names(chelsa_future_stack) <- stringr::str_extract(
  names(chelsa_future_stack),
  paste0(
    "(?<=CHELSA_",
    tolower(FUTURE_GCM),
    "_",
    future_scenario,
    "_).+(?=_2071)"
  )
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
  filter(median_ratio > 5 | median_ratio < 0.2) |>
  pull(variable)

if (length(checkvars) > 0) {
  stop(
    "Possible scaling mismatch between current and future CHELSA layers: ",
    paste(checkvars, collapse = ", ")
  )
}
cat(
  "Current and future CHELSA medians agree to within a factor of 5 for every layer\n"
)

# Create unified extent and template grid ####

## Define domain extent ####

# Get European countries
europe_countries <- ne_countries(
  continent = "europe",
  scale = 10,
  returnclass = "sf"
)

# EU member states
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

# Transform to EPSG:3035
target_countries <- st_transform(target_countries, crs = "EPSG:3035")

# Create combined polygon for masking
mask_polygon <- target_countries %>%
  st_union() %>%
  st_sf()

## Create unified 250m grid ####

# Define unified European extent in EPSG:3035 covering EU + Norway
europe_ext_3035 <- ext(c(
  xmin = 2000000,
  xmax = 7000000,
  ymin = 1000000,
  ymax = 5500000
))

# Create 250m resolution template grid in EPSG:3035
template_250m <- rast(europe_ext_3035, resolution = 250, crs = "EPSG:3035")

# Define extent for cropping before reprojection (WGS84)
europe_ext_wgs84 <- ext(c(xmin = -10, xmax = 35, ymin = 35, ymax = 72))

# Process CHELSA data ####

## Current climate ####

# Crop to European extent first (in WGS84) to reduce memory, then mask sentinel before resampling
chelsa_past_cropped <- crop(chelsa_past_stack, europe_ext_wgs84) |>
  mask_sentinels()
assert_plausible(chelsa_past_cropped, "CHELSA current, Europe crop")

# Project to EPSG:3035 at 250m resolution
chelsa_past_3035 <- project(
  x = chelsa_past_cropped,
  y = template_250m,
  method = "bilinear"
)

# Mask to EU + Norway boundaries
chelsa_past_masked <- mask(chelsa_past_3035, mask_polygon)

## Future climate ####

# Crop to European extent first (in WGS84) to reduce memory, then mask sentinel
chelsa_future_cropped <- crop(chelsa_future_stack, europe_ext_wgs84) |>
  mask_sentinels()
assert_plausible(chelsa_future_cropped, "CHELSA future, Europe crop")

# Project to EPSG:3035 at 250m resolution
chelsa_future_3035 <- project(
  x = chelsa_future_cropped,
  y = template_250m,
  method = "bilinear"
)

# Mask to EU + Norway boundaries
chelsa_future_masked <- mask(chelsa_future_3035, mask_polygon)

# Process paleo-derived predictors ####

# Crop paleo stack to European extent (in WGS84)
paleo_cropped <- crop(paleo_stack, europe_ext_wgs84)

# Project to EPSG:3035 at 250m resolution
paleo_3035 <- project(
  x = paleo_cropped,
  y = template_250m,
  method = "bilinear"
)

# Mask to EU + Norway boundaries
paleo_masked <- mask(paleo_3035, mask_polygon)

cat("Processed", nlyr(paleo_masked), "paleo-derived predictors\n")

# Build terrain variables (elevation + slope) using elevatr ####

# Note: This step may take significant time as elevatr downloads elevation data.
# Elevation is downloaded at the native resolution available (~1 arc-second ~30m),
# then aggregated to 250m.

cat("Fetching elevation data from elevatr for EU + Norway...\n")

# Create a simple raster template at ~500m for elevatr query (to get reasonable file sizes)
# elevatr will return its native resolution, then we aggregate
query_res <- 500
template_query <- rast(
  europe_ext_3035,
  resolution = query_res,
  crs = "EPSG:3035"
)

# Convert template to sf for elevatr
query_bbox_sf <- as.polygons(template_query) %>%
  st_as_sf() %>%
  st_transform(crs = "EPSG:4326")

# Fetch elevation from elevatr (z=9 is ~30m resolution)
# This uses the default GEBCO source; adjust z parameter if needed for different resolution
elev_raw <- get_elev_raster(query_bbox_sf, z = 9, clip = FALSE)

# Reproject to EPSG:3035
elev_3035 <- project(rast(elev_raw), template_250m, method = "bilinear")

# Mask to domain
elev_masked <- mask(elev_3035, mask_polygon)
names(elev_masked) <- "elevation"

cat("Elevation fetched and reprojected to 250m\n")

# Calculate slope from elevation at 250m
slope_250m <- terrain(elev_masked, v = "slope", unit = "degrees")
names(slope_250m) <- "slope"

cat("Slope calculated at 250m\n")

# Combine all predictors ####

predictors_current <- c(
  chelsa_past_masked,
  paleo_masked,
  elev_masked,
  slope_250m
)

predictors_future <- c(
  chelsa_future_masked,
  paleo_masked,
  elev_masked,
  slope_250m
)

stopifnot(nlyr(predictors_current) == nlyr(predictors_future))

cat("Combined predictor stack has", nlyr(predictors_current), "layers\n")

# Fix threshold-based variables with truncated ranges ####

# Consistent NA-to-zero handling for all domains: gdd5, gdd10, gst, swe, gsp
THRESHOLD_VARS <- c("gdd5", "gdd10", "gst", "swe", "gsp")

record_settings(
  "R/pl0_collatePredictors.R",
  future_scenario = future_scenario,
  future_gcm = FUTURE_GCM,
  threshold_vars = THRESHOLD_VARS,
  gsp_sentinel_bound_mm = MAX_PLAUSIBLE[["gsp"]],
  terrain_resolution_m = 250,
  elevation_source = "elevatr"
)

predictors_current <- fill_threshold_na(
  predictors_current,
  THRESHOLD_VARS,
  "unified model - current"
)
predictors_future <- fill_threshold_na(
  predictors_future,
  THRESHOLD_VARS,
  "unified model - future"
)

# Write outputs ####

writeRaster(
  predictors_current,
  filename = "output/predictors_unified_250m_EUNorway_current_EPSG3035.tif",
  overwrite = TRUE,
  names = names(predictors_current)
)

writeRaster(
  predictors_future,
  filename = "output/predictors_unified_250m_EUNorway_future_EPSG3035.tif",
  overwrite = TRUE,
  names = names(predictors_future)
)

cat("Unified predictor stacks written to:\n")
cat("  - output/predictors_unified_250m_EUNorway_current_EPSG3035.tif\n")
cat("  - output/predictors_unified_250m_EUNorway_future_EPSG3035.tif\n")

# Quick visualization check
if (interactive()) {
  plot(
    predictors_current[[c(
      which(names(predictors_current) == "bio01"),
      which(names(predictors_current) == "elevation"),
      which(names(predictors_current) == "paleo_years_icefreeland")
    )]]
  )
}

# Spatial coverage validation ####

# Check spatial coverage across both current and future stacks
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
  arrange(notNA_current) |>
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

# Clean up
rm(
  chelsa_past_stack,
  chelsa_past_cropped,
  chelsa_past_3035,
  chelsa_past_masked,
  chelsa_future_stack,
  chelsa_future_cropped,
  chelsa_future_3035,
  chelsa_future_masked,
  paleo_stack,
  paleo_cropped,
  paleo_3035,
  paleo_masked,
  elev_raw,
  elev_3035,
  predictors_current,
  predictors_future
)
terra::tmpFiles(remove = TRUE)

# sessionInfo ####

sessioninfo::session_info()
