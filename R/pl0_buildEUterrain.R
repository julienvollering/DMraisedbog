# Elevation and slope over the EU domain ####

# PURPOSE: Fetches and derives EU elevation and slope at 250 m via elevatr, mirroring the aggregate-then-slope recipe used for Norway DTM50.

# DTM50 is Norway-only, so the two terrain predictors have no EU source. This script
# builds them once, on the same 250 m grid the EU block is drawn from, mirroring the
# Norwegian methodology in pl0_collatePredictors.R: aggregate the DEM to ~250 m by mean,
# then terra::terrain() slope in degrees.
#
# WHY A RASTER RATHER THAN PER-POINT. The archived R/archive/pl3_extractEUpredictors.R
# fetched a small DEM patch around each EU presence -- fine for ~700 points, hopeless for the ~35k absences the EU
# block now carries, since every point is a separate network round-trip. Fetching the
# domain once and reading both blocks off it is the same methodology at a workable cost,
# and it makes the terrain reproducible instead of re-downloaded per script.
#
# Chunked and cached because the fetch is the slow, failure-prone part: {elevatr} pulls
# AWS Terrain Tiles over the network, and a single request covering a 1,710 x 2,268 km
# bounding box would be both enormous and all-or-nothing. Chunks that already exist on
# disk are skipped, so an interrupted run resumes instead of restarting.

library(sf)
library(terra)
library(dplyr)
library(elevatr)

dir.create("output/pl0/eu_terrain_parts", showWarnings = FALSE, recursive = TRUE)

crs_3035 <- "EPSG:3035"
target_res <- 250 # m, matches the regional predictor grid
chunk_m <- 100000 # 100 km tiles
halo_m <- 3000 # fetched beyond each chunk so slope is valid at the seams
dem_zoom <- 9 # ~87 m at these latitudes; aggregates cleanly toward 250 m
max_tries <- 3

force_rebuild <- FALSE

## Target grid and chunks ####

grid_250m <- rast("output/pl0/epm_map_cat_250m.tif")

domain <- st_read("data/DMraisedbog.gpkg", layer = "EU_domain", quiet = TRUE)
stopifnot(st_crs(domain)$epsg == 3035)

bb <- st_bbox(domain)
xs <- seq(floor(bb$xmin / chunk_m) * chunk_m, bb$xmax, by = chunk_m)
ys <- seq(floor(bb$ymin / chunk_m) * chunk_m, bb$ymax, by = chunk_m)

chunks <- expand.grid(x = xs, y = ys) |>
  as_tibble() |>
  mutate(id = sprintf("x%07.0f_y%07.0f", x, y))

chunk_poly <- function(x, y) {
  st_as_sfc(st_bbox(
    c(xmin = x, xmax = x + chunk_m, ymin = y, ymax = y + chunk_m),
    crs = st_crs(3035)
  ))
}

# Only chunks that actually meet the domain: the bounding box is nearly eight times the
# domain's area, and fetching the rest would be most of the runtime for no rows.
keep <- vapply(
  seq_len(nrow(chunks)),
  function(i) {
    length(st_intersects(chunk_poly(chunks$x[i], chunks$y[i]), domain)[[1]]) > 0
  },
  logical(1)
)
chunks <- chunks[keep, ]

cat(
  "Chunks intersecting the domain:", nrow(chunks),
  "of", length(xs) * length(ys), "in the bounding box\n"
)

## Fetch, aggregate, derive ####

build_chunk <- function(x, y, id) {
  part <- file.path("output/pl0/eu_terrain_parts", paste0(id, ".tif"))
  if (file.exists(part) && !force_rebuild) {
    return(part)
  }

  aoi <- st_as_sfc(st_bbox(
    c(
      xmin = x - halo_m, xmax = x + chunk_m + halo_m,
      ymin = y - halo_m, ymax = y + chunk_m + halo_m
    ),
    crs = st_crs(3035)
  )) |>
    st_sf(geometry = _)

  dem <- NULL
  for (try_i in seq_len(max_tries)) {
    dem <- tryCatch(
      suppressWarnings(suppressMessages(
        rast(get_elev_raster(aoi, z = dem_zoom, clip = "bbox", verbose = FALSE))
      )),
      error = function(e) NULL
    )
    if (!is.null(dem)) break
    cat("  retry", try_i, "for chunk", id, "\n")
  }
  if (is.null(dem)) {
    cat("  FAILED chunk", id, "\n")
    return(NA_character_)
  }

  if (!same.crs(dem, crs_3035)) {
    dem <- project(dem, crs_3035)
  }

  # Mean-aggregate toward 250 m, then slope -- the order used for DTM50 in
  # pl0_collatePredictors.R. Deriving slope from the native DEM and aggregating afterwards
  # would give a systematically steeper answer and would not match the Norway block.
  fact <- max(1, round(target_res / mean(res(dem))))
  dem_agg <- aggregate(dem, fact = fact, fun = "mean", na.rm = TRUE)
  slope_agg <- terrain(dem_agg, v = "slope", unit = "degrees")

  terr <- c(dem_agg, slope_agg)
  names(terr) <- c("elevation", "slope")

  # Snap onto the shared grid, then drop the halo so parts tile without overlap.
  target <- crop(grid_250m, ext(x, x + chunk_m, y, y + chunk_m), snap = "out")
  out <- resample(terr, target, method = "bilinear")

  writeRaster(
    out, part,
    overwrite = TRUE, gdal = c("COMPRESS=LZW", "TILED=YES")
  )
  part
}

parts <- character(nrow(chunks))
for (i in seq_len(nrow(chunks))) {
  parts[i] <- build_chunk(chunks$x[i], chunks$y[i], chunks$id[i])
  if (i %% 10 == 0) cat("  chunk", i, "/", nrow(chunks), "\n")
}

failed <- sum(is.na(parts))
parts <- parts[!is.na(parts)]
cat("Chunks built:", length(parts), " failed:", failed, "\n")
stopifnot(length(parts) > 0)

## Merge ####

terrain_250m <- merge(sprc(lapply(parts, rast)))
names(terrain_250m) <- c("elevation", "slope")

# AWS Terrain Tiles carry two things that are not "land elevation in this domain":
# bathymetry (the Baltic reads to about -260 m, which is correct and simply not land) and
# occasional decode artifacts (six cells came back above 5,000 m in the Oetztal Alps,
# where the true maximum is under 4,000). Neither touches a sampled row, but a future
# re-draw could land on one, so implausible values are made NA -- a row that hits one is
# then dropped by drop_na() in pl2_createModelingFrame.R rather than trained on. Sea cells
# are left alone: they are outside the domain and never sampled.
elev_plausible <- c(-500, 5000)
n_implausible <- as.numeric(global(
  terrain_250m[["elevation"]] < elev_plausible[1] |
    terrain_250m[["elevation"]] > elev_plausible[2],
  "sum", na.rm = TRUE
))
cat(
  "Elevation cells outside", paste(elev_plausible, collapse = " to "),
  "m, set NA:", n_implausible, "
"
)

terrain_250m[["elevation"]] <- clamp(
  terrain_250m[["elevation"]],
  lower = elev_plausible[1], upper = elev_plausible[2], values = FALSE
)
terrain_250m[["slope"]] <- mask(terrain_250m[["slope"]], terrain_250m[["elevation"]])

writeRaster(
  terrain_250m, "output/pl0/eu_terrain_250m.tif",
  overwrite = TRUE, gdal = c("COMPRESS=LZW", "TILED=YES")
)

## Report ####

cat("Elevation (m)  min/med/max:", paste(
  round(as.numeric(global(terrain_250m[["elevation"]], quantile,
    probs = c(0, 0.5, 1), na.rm = TRUE
  )), 1),
  collapse = " / "
), "\n")
cat("Slope (deg)    min/med/max:", paste(
  round(as.numeric(global(terrain_250m[["slope"]], quantile,
    probs = c(0, 0.5, 1), na.rm = TRUE
  )), 2),
  collapse = " / "
), "\n")

plot(terrain_250m[["elevation"]], main = "EU domain elevation, 250 m")

cat("Wrote output/pl0/eu_terrain_250m.tif\n")

# sessionInfo ####

sessioninfo::session_info()
