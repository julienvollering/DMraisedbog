# Extract the 43 local-model predictors at EU presence locations ####

# PURPOSE: Reconstructs the 43 local-model predictors at EU presence locations, the shared input for the pl3 reliability work.

# pl3 shared-input step. The European raised-bog presences are a separately-
# coloured domain in the weighted-PCA diagnostic (pl3_exploreWeightedPCA.R) and,
# under the plan's section 1.3 sweep, the high-novelty end of the training axis.
# They live natively in the global 5 km feature space, so here we reconstruct the
# 43 local-model predictors at their locations:
#   - 41 CHELSA climate + paleo features: read from the global 5 km EU+Norway
#     predictor stack (same variables as the 250 m regional stack, coarser grid)
#   - elevation + slope: not available outside Norway (derived from DTM50), so
#     fetched from a global DEM via {elevatr} and processed to mirror the
#     regional methodology in pl0_collatePredictors.R, i.e. mean-aggregate the
#     DEM to ~250 m, then terra::terrain() slope in degrees.
#
# rf_global is deliberately absent (plan section 1.2). The header of this script
# used to claim the anchors were leakage-free while reading their rf_global from
# pl0_modelGlobalScale.R's *in-sample* prediction surface -- that model trains on
# these very presences. The covariate is retired with the hierarchical cascade,
# so the leak is removed rather than patched. Until pl2_createModelingFrame.R is
# re-run without rf_global the guard below will fail loudly, which is correct.
#
# Output columns mirror the predictor columns of modeling_frame_regional.csv
# (plus x, y) so the EU rows can be row-bound onto the current-scenario training
# frame downstream (PCA diagnostic; pl3_assembleAnchors.R).
#
# Requires {elevatr} (install.packages("elevatr")); it fetches DEM tiles from AWS
# Terrain Tiles, so this step needs an internet connection.

library(readr)
library(dplyr)
library(tibble)
library(terra)
library(sf)
library(elevatr)

dir.create("output/pl3", showWarnings = FALSE, recursive = TRUE)

seed <- 42
set.seed(seed)

crs_3035 <- "EPSG:3035"

## Canonical predictor order (from the modeling frame) ####

mf_cols <- names(read_csv(
  "output/pl2/modeling_frame_regional.csv",
  n_max = 0,
  show_col_types = FALSE
))
# `dataset` is bookkeeping in the pooled frame, not a predictor, so it is dropped here
# alongside the response and coordinates. Without it the guard below fails on `dataset`,
# which is the correct behaviour but the wrong diagnosis.
feat <- setdiff(
  mf_cols,
  c("scenario", "response", "dataset", "x", "y")
) # 43 predictors

## EU presence coordinates (EPSG:3035, matching the global stacks) ####

# Use ALL EU presences, not the spatially-thinned set. Which side of the
# training/evaluation line these rows fall on is set per cut by the section 1.3
# sweep, not here, so this step reconstructs predictors for every one of them and
# leaves the cut to downstream scripts.
eu <- read_csv(
  "output/presence_coords_global_eu.csv",
  show_col_types = FALSE
)

## 41 climate + paleo features from the global 5 km stack ####

global_5km <- rast("output/predictors_global_5km_EUNorway_EPSG3035.tif")
feat_climate <- terra::extract(global_5km, eu[c("x", "y")], ID = FALSE)

## elevation + slope from a global DEM (elevatr) ####

# Per-point: fetch a small DEM patch around the presence, project to EPSG:3035,
# mean-aggregate to ~250 m (mirroring DTM50 50 m -> 250 m), derive slope in
# degrees, then read elevation & slope at the point. A ~1 km buffer guarantees a
# full 3x3 neighbourhood for a stable slope estimate at the central 250 m cell.

target_res <- 250 # m, matches the regional predictor grid
buffer_m <- 1000 # m, DEM neighbourhood fetched per point
dem_zoom <- 11 # elevatr zoom (~40-50 m cells at EU latitudes)

extract_elev_slope <- function(x, y) {
  pt <- st_sf(geometry = st_sfc(st_point(c(x, y)), crs = 3035))
  aoi <- st_buffer(pt, buffer_m)

  dem <- tryCatch(
    suppressWarnings(suppressMessages(
      rast(get_elev_raster(aoi, z = dem_zoom, clip = "bbox", verbose = FALSE))
    )),
    error = function(e) NULL
  )
  if (is.null(dem)) {
    return(c(elevation = NA_real_, slope = NA_real_))
  }

  if (!same.crs(dem, crs_3035)) {
    dem <- project(dem, crs_3035)
  }

  # Mean-aggregate to ~250 m, as in pl0_collatePredictors.R
  fact <- max(1, round(target_res / mean(res(dem))))
  dem250 <- aggregate(dem, fact = fact, fun = "mean")
  slp <- terrain(dem250, v = "slope", unit = "degrees")

  xy <- cbind(x, y)
  c(
    elevation = terra::extract(dem250, xy)[, 1],
    slope = terra::extract(slp, xy)[, 1]
  )
}

terrain_vals <- mapply(extract_elev_slope, eu$x, eu$y)
terrain_df <- as.data.frame(t(terrain_vals))

n_failed <- sum(is.na(terrain_df$elevation))
cat("EU presences:", nrow(eu), " elevation/slope failures:", n_failed, "\n")

## Assemble in modeling-frame predictor order ####

eu_predictors <- bind_cols(
  tibble(x = eu$x, y = eu$y),
  feat_climate,
  terrain_df
)

# Guard: every modeling-frame predictor must be reconstructed
missing <- setdiff(feat, names(eu_predictors))
if (length(missing) > 0) {
  stop(
    "Predictors missing from EU extraction: ",
    paste(missing, collapse = ", ")
  )
}

eu_predictors <- eu_predictors |>
  select(x, y, all_of(feat))

## Save ####

eu_predictors |>
  write_csv("output/pl3/eu_presence_predictors.csv", append = FALSE)

cat(
  "Wrote output/pl3/eu_presence_predictors.csv (",
  nrow(eu_predictors),
  "rows,",
  length(feat),
  "predictors )\n"
)

# sessionInfo ####

sessioninfo::session_info()
