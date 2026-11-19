# Observed summer warming to the 1981-2010 baseline, over Norway ####

# PURPOSE: Computes, from Berkeley Earth gridded observations, the summer (JJA) warming from each palaeoclimate record's own "present" to 1981-2010, so both records can be placed on the CHELSA 1981-2010 baseline.

library(tidyverse)
library(terra)

source("R/config.R")

# Read-only to the model. Both records of the bogs' past arrive as anomalies against their own
# present, while the future (CHELSA 2071-2100) is an absolute climate whose baseline is
# CHELSA 1981-2010 (plan_climateHistory.md 3.9, decision 1):
#   - EPOCH-2 is relative to ~100 cal BP (~1850 CE). Its offset uses 1850-1900, the IPCC AR6
#     preindustrial reference period.
#   - TraCE21k is relative to its last slice, the 1900-1990 mean (CHELSA-TraCE21k technical
#     specification, timeID 20).
# Each offset D = JJA(1981-2010) - JJA(record's present), from ONE observational dataset so
# the two records are treated identically. Adding (anomaly - D) to CHELSA 1981-2010 puts a
# past slice on the scale of the future.
#
# Berkeley Earth, land only, 1 degree, monthly anomalies from January 1750
# (Complete_TAVG_LatLong1.nc, ~200 MB). Layer k is month ((k - 1) %% 12) + 1 of year
# 1750 + (k - 1) %/% 12; checked 2026-10-07 against known Oslo extremes (July 2018 +4.7 C,
# January 1987 -6.7 C). The anomalies share one climatology, so differences of period means
# need no absolute temperatures.
#
# Output: data/BerkeleyEarth/baseline_offsets_jja_Norway.tif, layers D_epoch2 and D_trace (C).

## Configuration ####

BE_URL <- "https://berkeley-earth-temperature.s3.us-west-1.amazonaws.com/Global/Gridded/Complete_TAVG_LatLong1.nc"
BE_DIR <- "data/BerkeleyEarth"
BE_FILE <- file.path(BE_DIR, "Complete_TAVG_LatLong1.nc")

PERIOD_BASELINE <- c(1981, 2010) # CHELSA v2.1 current
PERIOD_EPOCH2 <- c(1850, 1900) # IPCC AR6 preindustrial, standing in for ~100 cal BP
PERIOD_TRACE <- c(1900, 1990) # TraCE21k last slice
JJA <- 6:8

NORWAY_EXT <- ext(3, 33, 56, 73)

record_settings(
  "R/pl0_prepareBaselineOffsets.R",
  source = BE_URL,
  period_baseline = PERIOD_BASELINE,
  period_epoch2 = PERIOD_EPOCH2,
  period_trace = PERIOD_TRACE,
  months = JJA
)

dir.create(BE_DIR, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(BE_FILE)) {
  download.file(BE_URL, BE_FILE, mode = "wb", quiet = TRUE)
}

## Period means ####

be <- rast(BE_FILE, subds = "temperature")
n <- nlyr(be)
months <- tibble(k = seq_len(n), year = 1750 + (k - 1) %/% 12, month = (k - 1) %% 12 + 1)

jja_mean <- function(period) {
  k <- months |>
    filter(month %in% JJA, year >= period[1], year <= period[2]) |>
    pull(k)
  s <- crop(be[[k]], NORWAY_EXT)
  # Berkeley Earth "Complete" is interpolated, so a missing month over Norway would be a
  # file problem rather than a data gap: stop rather than average around it.
  land <- !is.na(s[[length(k)]])
  stopifnot(global(is.na(s) & land, "sum")$sum == 0)
  mean(s)
}

baseline <- jja_mean(PERIOD_BASELINE)
offsets <- c(
  baseline - jja_mean(PERIOD_EPOCH2),
  baseline - jja_mean(PERIOD_TRACE)
)
names(offsets) <- c("D_epoch2", "D_trace")

out <- file.path(BE_DIR, "baseline_offsets_jja_Norway.tif")
writeRaster(offsets, out, overwrite = TRUE)
cat("Written", out, "\n")

## Checks ####

# Observed Norwegian summer warming is of order 1 C since the late 19th century and smaller
# since the mid 20th; anything far outside that is a period or indexing error.
print(global(offsets, c("min", "mean", "max"), na.rm = TRUE))
stopifnot(
  global(offsets[["D_epoch2"]], "mean", na.rm = TRUE)$mean > 0,
  global(offsets[["D_epoch2"]], "mean", na.rm = TRUE)$mean < 2.5,
  global(abs(offsets[["D_trace"]]), "max", na.rm = TRUE)$max < 2
)

# sessionInfo ####

sessioninfo::session_info()
