# EPOCH-2 pollen-based Holocene climate over Norway ####

# PURPOSE: Fetches the EPOCH-2 gridded pollen-based Holocene reconstruction (Mauri et al. 2015) and writes its summer-temperature and annual-precipitation anomalies over Norway, the second record of the bogs' past next to TraCE21k.

library(tidyverse)
library(terra)

source("R/config.R")

# Read-only to the model. pl2_exploreClimateHistory.R runs the bogs' climate history on two
# independent records in parallel (plan_climateHistory.md 3.9): TraCE21k (a transient climate
# simulation, centennial) and EPOCH-2 (pollen-based, millennial). They are never blended:
# TraCE's mid-Holocene error is not systematic in time (Karger et al. 2023, Fig. 5).
#
# EPOCH-2: 1 degree grid, 28.5-71.5 N, 11.5 W-50.5 E, 13 slices at 12000, 11000, ..., 1000 and
# 100 cal BP. Values are anomalies; the 100 BP slice is exactly 0 everywhere over Norway
# (checked 2026-10-07), so they are relative to ~100 cal BP (~1850 CE). Temperature in degC,
# precipitation in mm/day.
#
# Output: data/EPOCH2/epoch2_<var>_Norway.tif, one layer per slice, layer names "a<age>" with
# age in cal BP (b1950).

## Configuration ####

EPOCH2_URL <- "https://www.ncei.noaa.gov/pub/data/paleo/pollen/epoch2"
EPOCH2_DIR <- "data/EPOCH2"
EPOCH2_AGES <- c(12000, 11000, 10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, 100)

# Same window as the TraCE series, padded by one 1-degree cell so bilinear reads at the edge
# pixels have their neighbours.
NORWAY_EXT <- ext(3, 33, 56, 73)

record_settings(
  "R/pl0_prepareEPOCH2.R",
  epoch2_url = EPOCH2_URL,
  epoch2_ages_bp = range(EPOCH2_AGES),
  epoch2_vars = c("tanom_jja", "panom_ann", "tanom_jja_se")
)

dir.create(EPOCH2_DIR, showWarnings = FALSE, recursive = TRUE)

## Download ####

for (f in c("epoch-2.nc", "epoch-2_uncertainties.nc")) {
  dest <- file.path(EPOCH2_DIR, f)
  if (!file.exists(dest)) {
    download.file(file.path(EPOCH2_URL, f), dest, mode = "wb", quiet = TRUE)
  }
}

## Crop and write ####

write_epoch2 <- function(r, var) {
  stopifnot(nlyr(r) == length(EPOCH2_AGES))
  r <- crop(r, NORWAY_EXT)
  names(r) <- paste0("a", EPOCH2_AGES)
  out <- file.path(EPOCH2_DIR, paste0("epoch2_", var, "_Norway.tif"))
  writeRaster(r, out, overwrite = TRUE)
  cat("Written", out, "\n")
  r
}

tjja <- write_epoch2(rast(file.path(EPOCH2_DIR, "epoch-2.nc"), subds = "tanom_jja"), "tanom_jja")
pann <- write_epoch2(rast(file.path(EPOCH2_DIR, "epoch-2.nc"), subds = "panom_ann"), "panom_ann")
unc <- rast(file.path(EPOCH2_DIR, "epoch-2_uncertainties.nc"))
tse <- write_epoch2(unc[[grep("^tanom_jja_se_", names(unc))]], "tanom_jja_se")

## Checks ####

# The baseline: the 100 BP slice must be 0, or the anomalies are not relative to ~1850 and the
# preindustrial offset in pl2_exploreClimateHistory.R is wrong.
base <- values(tjja[["a100"]], na.rm = TRUE)
cat("tanom_jja at 100 BP over the window: range", range(base), "\n")
stopifnot(all(abs(base) < 1e-6))

# Orientation and magnitude: summer anomalies are a few degC, never tens.
stopifnot(all(abs(values(tjja, na.rm = TRUE)) < 15))

# sessionInfo ####

sessioninfo::session_info()
