# Build the predictor stacks ####

# PURPOSE: Assembles the 5 km EU working grid and the 250 m Norway predictor stack
# (CHELSA, paleo, DTM50 terrain, AR50 land cover) on a common EPSG:3035 origin.

# TWO GRIDS, AND WHY THEY ARE DIFFERENT RESOLUTIONS. The 250 m Norway stack is a
# modelling surface: the projection runs over it and the Norwegian training rows are read
# off it. The 5 km EU grid is NOT a modelling surface and nothing is ever trained on it --
# it is a coarse working grid, and its only consumers want it that way:
#
#   pl0_buildEUdomain.R           rasterize template + extent for delimiting the domain
#   pl0_collateEuropeanRaisedBog.R  bio10 at Natura 2000 polygon centroids
#   pl0_sampleEUabsences.R        bio10 as the absence stratification axis
#
# All three read bio10 and nothing else, at a grain where 5 km is the right answer: a
# domain polygon does not become more correct when rasterized 400x finer.
#
# WHAT THE EU BLOCK'S TRAINING ROWS USE INSTEAD. Not this grid. They are point-extracted
# from the native CHELSA and paleo files by extract_native_predictors() in R/functions.R,
# called from pl2_createModelingFrame.R. That is what closes the resolution asymmetry
# logged as notebook open issue 2 on 2026-09-14 -- EU rows used to be sampled off this
# 5 km grid, a genuinely 25x-smoothed field, while Norwegian rows came off 250 m.
#
# WHY NOT ONE UNIFIED 250 m EU+Norway STACK. That was tried on 2026-09-15 (commit
# 207dcab) and is what filled the disk mid-run: 20000 x 18000 cells x 43 layers is 58 GB
# per scenario, 115 GB across the two, before terra's scratch and before
# pl2_createModelingFrame.R copies both to scenario_*.tif -- against 57 GB free. Only 20%
# of that box is land and the EU part of it was only ever sampled at ~36k points. The
# same commit also replaced DTM50 with a single unbounded elevatr call over the whole box
# (12,561 tiles at z=9); the chunked, cached fetch in pl0_buildEUterrain.R is the working
# implementation of that idea and is where EU terrain still comes from.
#
# TERRAIN SOURCES DIFFER BY BLOCK, GRAIN DOES NOT. Norway uses DTM50 (the 50 m national
# model, mean-aggregated to 250 m); the EU uses elevatr at the same 250 m via
# pl0_buildEUterrain.R. Both derive slope AFTER aggregating, so neither is systematically
# steeper. The single cell that belongs to both blocks disagrees by 0.08 degrees of slope,
# which identify_dynamic_predictors() in R/functions.R already accounts for.

library(terra)
library(sf)
library(rnaturalearth)
library(tidyverse)

source("R/functions.R")
source("R/config.R")

# Clean up old terra temporary files from previous runs
terra::tmpFiles(remove = TRUE)

record_settings(
  "R/pl0_collatePredictors.R",
  future_scenario = FUTURE_SCENARIO,
  future_gcm = FUTURE_GCM,
  threshold_vars = THRESHOLD_VARS,
  gsp_sentinel_bound_mm = MAX_PLAUSIBLE[["gsp"]],
  eu_working_grid_m = 5000,
  norway_grid_m = 250,
  norway_elevation_source = "DTM50"
)

# Load native predictor sources ####

# Opened through the shared accessors so this script and extract_native_predictors()
# cannot disagree about layer names or about which paleo file is excluded.
chelsa_past_stack <- chelsa_native_stack("current")
chelsa_future_stack <- chelsa_native_stack("future")
paleo_stack <- paleo_native_stack()

cat("Loaded", nlyr(chelsa_past_stack), "CHELSA layers and", nlyr(paleo_stack), "paleo\n")

## Check scaling across variables ####

# Compare the scale factor and offset stored in each current and future CHELSA file.
# This replaces a comparison of medians over a fixed block of global cells: that block
# lay at ~84 N, where values near zero turned real warming (bio10, gdd0, gsl) and the
# unmasked gsp sentinel into false alarms. terra applies scoff on read, so matching
# metadata means the two stacks are in the same units.
if (!setequal(names(chelsa_past_stack), names(chelsa_future_stack))) {
  stop(
    "Current and future CHELSA stacks hold different layers: ",
    paste(
      setdiff(
        union(names(chelsa_past_stack), names(chelsa_future_stack)),
        intersect(names(chelsa_past_stack), names(chelsa_future_stack))
      ),
      collapse = ", "
    )
  )
}
scoff_past <- scoff(chelsa_past_stack)
rownames(scoff_past) <- names(chelsa_past_stack)
scoff_future <- scoff(chelsa_future_stack)
rownames(scoff_future) <- names(chelsa_future_stack)
scoff_future <- scoff_future[rownames(scoff_past), , drop = FALSE]

mismatch <- rownames(scoff_past)[
  rowSums(abs(scoff_past - scoff_future) > 1e-9) > 0
]
if (length(mismatch) > 0) {
  stop(
    "Scale/offset mismatch between current and future CHELSA layers: ",
    paste(mismatch, collapse = ", ")
  )
}
cat("Current and future CHELSA layers share scale factor and offset for every layer\n")

# EU working grid, 5 km ####

## Domain mask ####

europe_countries <- ne_countries(
  continent = "europe",
  scale = 10,
  returnclass = "sf"
)

eu_countries <- c(
  "Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark",
  "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Ireland",
  "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland",
  "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden"
)

target_countries <- europe_countries %>%
  filter(name %in% c(eu_countries, "Norway")) %>%
  st_transform(crs = "EPSG:3035")

mask_polygon <- target_countries %>%
  st_union() %>%
  st_sf()

## Template and projection ####

europe_ext_3035 <- ext(c(
  xmin = 2000000,
  xmax = 7000000,
  ymin = 1000000,
  ymax = 5500000
))

template_5km <- rast(europe_ext_3035, resolution = 5000, crs = "EPSG:3035")

europe_ext_wgs84 <- ext(c(xmin = -10, xmax = 35, ymin = 35, ymax = 72))

# Crop on the native grid first, then mask the sentinel, then project -- in that order,
# because bilinear resampling blends a sentinel into its neighbours and those blends
# cannot be recognised afterwards.
chelsa_eu_cropped <- crop(chelsa_past_stack, europe_ext_wgs84) |>
  mask_sentinels()
assert_plausible(chelsa_eu_cropped, "CHELSA current, Europe crop")

chelsa_eu_3035 <- project(
  x = chelsa_eu_cropped,
  y = template_5km,
  method = "bilinear"
)

paleo_eu_3035 <- crop(paleo_stack, europe_ext_wgs84) |>
  project(y = template_5km, method = "bilinear")

predictors_eu <- c(
  mask(chelsa_eu_3035, mask_polygon),
  mask(paleo_eu_3035, mask_polygon)
)

cat("EU working grid has", nlyr(predictors_eu), "layers\n")

predictors_eu <- fill_threshold_na(predictors_eu, THRESHOLD_VARS, "EU working grid")

writeRaster(
  predictors_eu,
  filename = "output/predictors_global_5km_EUNorway_EPSG3035.tif",
  overwrite = TRUE,
  names = names(predictors_eu),
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

# Norway stack, 250 m ####

## Terrain from DTM50 ####

# Aggregate to 250 m in the native UTM33 projection and derive slope there, BEFORE
# reprojecting. Deriving slope from the 50 m model and aggregating afterwards would give
# a systematically steeper answer; reprojecting first would distort the derivative.
# pl0_buildEUterrain.R mirrors this order for the EU side.
dtm_files <- list.files(
  "data/DTM50_UTM33_20250613",
  pattern = "\\.tif$",
  full.names = TRUE
)
stopifnot(length(dtm_files) > 0)

dtm_mosaic_utm33 <- do.call(mosaic, map(dtm_files, rast))
dtm_250m_utm33 <- aggregate(dtm_mosaic_utm33, fact = 5, fun = "mean")

terrain_stack_utm33 <- c(
  dtm_250m_utm33,
  terrain(dtm_250m_utm33, v = "slope", unit = "degrees")
)
names(terrain_stack_utm33) <- c("elevation", "slope")

cat("Terrain derived at 250 m in UTM33\n")

## Template aligned to the EU grid origin ####

# Same origin as template_5km so the two grids nest exactly: a 5 km cell boundary is
# always also a 250 m cell boundary, which is what lets pl0_sampleEUabsences.R resample
# 5 km bio10 onto the 250 m EU grid without introducing a half-cell shift.
norway_ext_3035 <- project(
  ext(dtm_mosaic_utm33),
  from = crs(dtm_mosaic_utm33),
  to = "EPSG:3035"
)

origin_x <- xmin(template_5km)
origin_y <- ymin(template_5km)

norway_ext_aligned <- ext(c(
  xmin = origin_x + floor((norway_ext_3035$xmin - origin_x) / 250) * 250,
  xmax = origin_x + ceiling((norway_ext_3035$xmax - origin_x) / 250) * 250,
  ymin = origin_y + floor((norway_ext_3035$ymin - origin_y) / 250) * 250,
  ymax = origin_y + ceiling((norway_ext_3035$ymax - origin_y) / 250) * 250
))

template_250m <- rast(norway_ext_aligned, resolution = 250, crs = "EPSG:3035")

## AR50 land cover ####

# artype_60 is the Norwegian LABEL source (pl0_labelNorwayBlock.R), not a predictor;
# pl2_createModelingFrame.R drops it before fitting so the response cannot leak.
land_mask_ar50 <- rast("output/ar50_250m_land_EPSG3035.tif")

ar50_250m <- rast("output/ar50_artype_layers/ar50_50m_EPSG3035_artype60.tif") |>
  aggregate(fact = 5, fun = "mean", na.rm = FALSE) |>
  resample(template_250m, method = "bilinear") |>
  mask(land_mask_ar50)
names(ar50_250m) <- "artype_60"

writeRaster(
  ar50_250m,
  filename = "output/ar50_250m_cover_EPSG3035.tif",
  overwrite = TRUE,
  names = names(ar50_250m),
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

## Climate and paleo over Norway ####

norway_ext_wgs84 <- project(norway_ext_3035, from = "EPSG:3035", to = "EPSG:4326")

terrain_250m <- project(terrain_stack_utm33, template_250m, method = "bilinear") |>
  mask(land_mask_ar50)

project_to_norway <- function(s, label) {
  cropped <- crop(s, norway_ext_wgs84) |> mask_sentinels()
  assert_plausible(cropped, label)
  project(cropped, y = template_250m, method = "bilinear") |>
    mask(land_mask_ar50)
}

chelsa_no_current <- project_to_norway(chelsa_past_stack, "CHELSA current, Norway crop")
chelsa_no_future <- project_to_norway(chelsa_future_stack, "CHELSA future, Norway crop")

paleo_no <- crop(paleo_stack, norway_ext_wgs84) |>
  project(y = template_250m, method = "bilinear") |>
  mask(land_mask_ar50)

## Combine and correct ####

predictors_current <- c(chelsa_no_current, terrain_250m, ar50_250m, paleo_no)
predictors_future <- c(chelsa_no_future, terrain_250m, ar50_250m, paleo_no)

stopifnot(nlyr(predictors_current) == nlyr(predictors_future))
cat("Norway stack has", nlyr(predictors_current), "layers\n")

predictors_current <- fill_threshold_na(
  predictors_current, THRESHOLD_VARS, "Norway - current"
)
predictors_future <- fill_threshold_na(
  predictors_future, THRESHOLD_VARS, "Norway - future"
)

writeRaster(
  predictors_current,
  filename = "output/predictors_regional_250m_Norway_current_EPSG3035.tif",
  overwrite = TRUE,
  names = names(predictors_current),
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

writeRaster(
  predictors_future,
  filename = "output/predictors_regional_250m_Norway_future_EPSG3035.tif",
  overwrite = TRUE,
  names = names(predictors_future),
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

cat("Predictor stacks written:\n")
cat("  - output/predictors_global_5km_EUNorway_EPSG3035.tif\n")
cat("  - output/predictors_regional_250m_Norway_current_EPSG3035.tif\n")
cat("  - output/predictors_regional_250m_Norway_future_EPSG3035.tif\n")

# Spatial coverage validation ####

coverage_comparison <- left_join(
  global(predictors_current, "notNA") |>
    mutate(layer = names(predictors_current)),
  global(predictors_future, "notNA") |>
    mutate(layer = names(predictors_future)),
  by = "layer",
  suffix = c("_current", "_future")
) |>
  mutate(equal = notNA_current == notNA_future) |>
  arrange(notNA_current) |>
  select(layer, notNA_current, notNA_future, equal)

cat("\nNorway coverage (sorted by current coverage, ascending):\n")
print(coverage_comparison, row.names = FALSE)

if (!all(coverage_comparison$notNA_current == coverage_comparison$notNA_current[1])) {
  cat("\nWARNING: Current predictors have inconsistent spatial coverage\n")
}
if (!all(coverage_comparison$notNA_future == coverage_comparison$notNA_future[1])) {
  cat("\nWARNING: Future predictors have inconsistent spatial coverage\n")
}
if (!all(coverage_comparison$equal)) {
  cat("\nWARNING: Coverage differs between current and future scenarios\n")
}

# Clean up
rm(
  chelsa_past_stack, chelsa_future_stack, paleo_stack,
  chelsa_eu_cropped, chelsa_eu_3035, paleo_eu_3035, predictors_eu,
  dtm_mosaic_utm33, dtm_250m_utm33, terrain_stack_utm33, terrain_250m,
  chelsa_no_current, chelsa_no_future, paleo_no,
  predictors_current, predictors_future
)
terra::tmpFiles(remove = TRUE)

# sessionInfo ####

sessioninfo::session_info()
