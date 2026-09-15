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
