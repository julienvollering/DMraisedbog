# Pipeline configuration ####

# Two things live here, and every script that has either sources this file:
#
#  1. Constants that more than one script must agree on. The future scenario and GCM
#     were defined independently in pl0_downloadCHELSA.R and pl0_collatePredictors.R;
#     had they drifted apart, the filename filter in the second would have matched
#     nothing and the future stack would have been silently empty. The response levels
#     were typed out in six scripts.
#
#  2. A registry for the design constants each script sets (absences per bog, partition
#     count, tree count, bins, bootstrap reps, the calibrating arm, ...). They stay
#     defined in the script they belong to, next to the reasoning; record_settings()
#     copies them into output/pipeline_settings.csv so a reviewer can confirm the
#     configuration of a run from one file. RUNALL prints that table at the end of a run.

## Shared constants ####

# CHELSA CMIP6 scenario and GCM for the 2071-2100 stack. pl0_downloadCHELSA.R fetches
# these files; pl0_collatePredictors.R stacks only files whose names carry both.
FUTURE_SCENARIO <- "ssp370"
FUTURE_GCM <- "GFDL-ESM4"

# The 3-class response, in the order every probability matrix, confusion table and
# prediction layer uses.
RESPONSE_LEVELS <- c("nonpeat", "otherpeat", "bog")

## CHELSA no-data conventions ####

# These three constants describe how CHELSA encodes "the quantity does not occur here",
# and they live here because TWO code paths now have to agree on them: the Norwegian
# block reads its predictors off a projected 250 m raster, the EU block point-extracts
# from the native grid (see extract_native_predictors() in R/functions.R). They were
# previously written out in pl0_collatePredictors.R alone, with the raster path filling
# gdd5 in the Norway stack but not in the EU one -- the asymmetry logged as notebook
# open issue 1 on 2026-09-14. One list, sourced by both paths, is what stops that
# recurring.

# Threshold quantities are truncated to positive values, so a cell where the quantity is
# really 0 arrives as NA rather than as zero. gsp joins them only after its sentinel is
# masked (below): no growing season, no growing-season precipitation.
THRESHOLD_VARS <- c("gdd5", "gdd10", "gst", "swe", "gsp")

# gsp is stored as an unsigned 32-bit integer with a 0.1 scale factor and NO NoData tag,
# so where the growing season has zero length the sentinel 4294967295 is read as a real
# value and scaled to 4.29e8 mm. The gdd/gst/swe layers carry a NoData tag and arrive as
# NA; gsp did not, and the sentinel reached the modelling frame in 8,987 training rows
# (notebook 2026-09-12). It must be masked on the NATIVE grid, before any bilinear step,
# because interpolation blends a sentinel into its neighbours and those blends cannot be
# recognised afterwards -- which is equally true of project() and of a bilinear extract.
SENTINEL_VARS <- c("gsp")
MAX_PLAUSIBLE <- c(gsp = 1e5) # mm; the wettest cells in the frame are ~7,000

## Settings registry ####

SETTINGS_PATH <- "output/pipeline_settings.csv"

# record_settings("R/pl2_evaluate.R", NTREE = 1000, ...) replaces that script's rows in
# the registry and leaves every other script's rows alone, so the file always reflects
# the last run of each script. One call per script; vectors are stored ";"-joined.
record_settings <- function(script, ...) {
  vals <- list(...)
  stopifnot(length(vals) > 0, !is.null(names(vals)), all(nzchar(names(vals))))
  rows <- data.frame(
    script = script,
    setting = names(vals),
    value = vapply(
      vals,
      function(v) paste(format(v, trim = TRUE), collapse = ";"),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(SETTINGS_PATH), showWarnings = FALSE, recursive = TRUE)
  old <- if (file.exists(SETTINGS_PATH)) {
    read.csv(SETTINGS_PATH, stringsAsFactors = FALSE)
  } else {
    rows[0, ]
  }
  old <- old[old$script != script, , drop = FALSE]
  write.csv(rbind(old, rows), SETTINGS_PATH, row.names = FALSE)
  invisible(rows)
}
