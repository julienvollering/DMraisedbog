# CHELSA-TraCE21k centennial time series over Norway ####

# PURPOSE: Builds per-variable multi-layer stacks of the centennial TraCE21k climate over Norway for 10-0 ka, keeping the time dimension that pl0_prepareCHELSATrace.R collapses.

library(tidyverse)
library(terra)

source("R/config.R")

# Read-only to the model. pl2_exploreClimateHistory.R asks whether the future climate of
# today's bogs lies inside the climate those bogs have lived through since they formed
# (plan_climateHistory.md). That needs the TIME dimension, which pl0_prepareCHELSATrace.R
# averages away into the static paleo_* predictors. This script keeps it.
#
# WHY A REMOTE WINDOW READ, NOT A DOWNLOAD. Every centennial slice for 10-0 ka is 101 files
# per variable at ~140 MB each (global, 30 arcsec), ~56 GB for four variables. The files
# are cloud-optimised GeoTIFFs (512 px tiles, LZW) on a server that honours byte ranges,
# so GDAL's /vsicurl/ fetches only the tiles inside the Norway window: ~4 s per slice.
#
# WHY NOT data/CHELSA/paleo/. pl0_prepareCHELSATrace.R globs every bio01/bio12/glz/orog
# file in that folder and stops unless all four share the same time slices. Centennial
# bio12 files there would either break it or, with matching glz/orog, silently change the
# paleo_* predictors the production model was fitted on. The series therefore never touches
# that folder, and pl0_downloadCHELSA.R is left unchanged.
#
# Output: data/CHELSA/paleo_derived/trace_series_<var>_Norway.tif, one layer per slice,
# layer names "t<time>" with TraCE time in hundred-year units (20 = the last slice;
# age b2k = (20 - time) * 100). Values in CHELSA v2.1 units: degC for temperatures (TraCE
# stores Kelvin), kg m-2 (= mm) for precipitation.

## Configuration ####

# bio10 and bio12 are the primary axes; bio18 (warmest-quarter precipitation) and bio06
# (coldest-month minimum) are the plan's sensitivity axes. None of gsp, gdd*, gsl, gst, npp
# exist in TraCE (checked 2026-10-06).
SERIES_VARS <- c("bio10", "bio12", "bio18", "bio06")
TEMPERATURE_VARS <- c("bio10", "bio06")

# Every centennial slice from 10 ka (time -80) to the last slice (time 20). Five-hundred-year
# spacing samples one century in five and can miss the warmest HTM century (plan, gap 3).
SERIES_TIMES <- seq(-80, 20, by = 1)

# Mainland Norway with a margin; the same box pl0_downloadCHELSA.R documents.
NORWAY_EXT <- ext(4, 32, 57, 72)

TRACE_URL <- "https://os.zhdk.cloud.switch.ch/chelsa01/chelsa_trace21k/global/bioclim"
OUT_DIR <- "data/CHELSA/paleo_derived"
CACHE_DIR <- "data/CHELSA/trace_series_cache"

record_settings(
  "R/pl0_prepareCHELSATraceSeries.R",
  series_vars = SERIES_VARS,
  series_times = range(SERIES_TIMES),
  series_step_hundreds = 1,
  norway_ext = as.vector(NORWAY_EXT)
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

## Remote read ####

# Same time formatting as pl0_downloadCHELSA.R: "-080" for negatives, "0020" otherwise.
format_trace_time <- function(t) {
  if (t < 0) sprintf("-%03d", abs(t)) else sprintf("%04d", t)
}

trace_url <- function(var, t) {
  paste0(
    "/vsicurl/", TRACE_URL, "/", var, "/CHELSA_TraCE21k_", var, "_",
    format_trace_time(t), "_V.1.0.tif"
  )
}

# One slice, cropped and cached. A cached file is reused, so an interrupted run resumes;
# a failed read is retried once before stopping (the server occasionally drops a range
# request, and a silent gap in the series would bias the per-pixel maximum).
read_slice <- function(var, t) {
  cache_file <- file.path(CACHE_DIR, paste0(var, "_", t, ".tif"))
  if (file.exists(cache_file)) {
    return(cache_file)
  }
  for (attempt in 1:2) {
    ok <- tryCatch(
      {
        r <- crop(rast(trace_url(var, t)), NORWAY_EXT)
        if (var %in% TEMPERATURE_VARS) {
          r <- r - 273.15
        }
        writeRaster(r, cache_file, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
        TRUE
      },
      error = function(e) {
        message("  read failed (", var, ", ", t, ", attempt ", attempt, "): ",
          conditionMessage(e))
        FALSE
      }
    )
    if (ok) {
      return(cache_file)
    }
  }
  stop("Could not read TraCE ", var, " at time ", t)
}

## Build the stacks ####

for (var in SERIES_VARS) {
  out_file <- file.path(OUT_DIR, paste0("trace_series_", var, "_Norway.tif"))
  if (file.exists(out_file)) {
    cat("Exists, skipping:", out_file, "\n")
    next
  }
  cat("Reading", var, "over", length(SERIES_TIMES), "slices...\n")
  files <- map_chr(SERIES_TIMES, \(t) read_slice(var, t))

  s <- rast(files)
  names(s) <- paste0("t", SERIES_TIMES)
  writeRaster(s, out_file, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  cat("  Written", out_file, "\n")
}

## Checks ####

# Units: a Kelvin temperature left unconverted, or a precipitation read at the wrong scale,
# would pass every later step and give a confident wrong answer. Ranges over the window,
# last slice.
for (var in SERIES_VARS) {
  s <- rast(file.path(OUT_DIR, paste0("trace_series_", var, "_Norway.tif")))
  stopifnot(nlyr(s) == length(SERIES_TIMES))
  rng <- as.vector(unlist(global(s[["t20"]], "range", na.rm = TRUE)))
  cat(var, "range at t20:", round(rng, 1), "\n")
  if (var %in% TEMPERATURE_VARS) {
    stopifnot(rng[1] > -60, rng[2] < 40)
  } else {
    stopifnot(rng[1] >= 0, rng[2] < 10000)
  }
}

# The per-slice cache is only an aid to resuming; the stacks are the product.
unlink(CACHE_DIR, recursive = TRUE)

# sessionInfo ####

sessioninfo::session_info()
