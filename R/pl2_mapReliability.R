# Map reliability over the projection domain ####

# PURPOSE: Computes DI for every projection cell on the frozen ruler and reads the error profiles onto it, giving per-pixel expected skill and an area-of-applicability mask.

# The last step, and the one the whole DI apparatus exists for. pl2_predict.R says WHAT the
# model predicts at each cell; this says HOW MUCH THAT PREDICTION CAN BE TRUSTED there, by
# measuring the cell's novelty on the frozen ruler and reading the cross-validated skill at
# that novelty off the curves fitted in pl2_fitErrorProfiles.R.
#
# The framing to keep in mind is forecast lead time. Skill-versus-lead-time is measured on
# past forecasts exactly as skill-versus-DI is measured on held-out folds, and it is why a
# single accuracy number for the map is the wrong thing to ask for. The sentence these
# layers are built to support:
#
#   "This pixel's climate is THIS far outside anything the model was trained on. In
#    cross-validation at that degree of novelty, the model correctly identified X% of known
#    bogs while falsely flagging Y% of known non-bogs."
#
# Per-pixel, in units ecologists use, no prevalence assumption, and no claim to be a
# posterior probability of bog.
#
# WHAT THESE LAYERS ARE NOT. `exp_recall_bog` is not P(bog) at the cell, and must never be
# read as one. It is the expected TRUE-POSITIVE RATE of the classifier among cells that
# really are bog, at this cell's novelty -- a property of the model's discrimination, not of
# the cell's ecology. Combining it with the class probability from rf_local_pred_brf.tif to
# get a posterior would require the prevalence of bog in the future domain, which is exactly
# the unknown being estimated.
#
# WHICH AXIS THE CURVES ARE READ ON (notebook 2026-09-12). Skill is interpolated on SIGNED
# bio10 OFFSET from the training bogs, not on DI. DI still gets a layer and still sets the
# AOA mask -- applicability is a multivariate question -- but it could not carry the skill
# read: its bog-count bins all ended below DI 0.187 while the whole projection sat above
# that, so every cell read one bin; and because DI is unsigned it scored cold-novel and
# warm-novel cells alike when the model treats them oppositely. The offset axis is
# one-dimensional, signed and bog-referenced, so a pixel position on it means "this much
# warmer than the typical bog the model learned from", and the projection spreads across the
# curve instead of collapsing into its end. The origin is `axis_ref_production` from
# error_profiles.rds, so map offsets and CV offsets share it.
#
# BEYOND THE LAST BIN. Interpolation is clamped (`rule = 2`), so cells more novel than the
# most novel CV bin inherit that bin's skill rather than an extrapolated trend. That is
# deliberate and optimistic-leaning, which is why the AOA mask matters: past the threshold
# the honest statement is "no estimate", not a clamped one. Report the two together.

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(sf)
library(terra)

source("R/functions.R")

## Configuration ####

SCENARIOS <- c(current = "output/pl2/scenario_current.tif",
               future = "output/pl2/scenario_future.tif")

# Layers written per scenario, and the profile column each is read from.
PROFILE_LAYERS <- c(
  exp_auc_bog = "auc_bog",
  exp_recall_bog = "recall_bog",
  exp_fpr_bog = "fpr_bog",
  exp_macro_auc = "macro_auc"
)

## Inputs ####

# The profiles must postdate the predictions they were fitted on, and those the frame they
# were scored on; a map read from stale curves is the failure this guards against.
assert_fresher(
  "output/pl2/error_profiles.rds",
  "output/pl2/predictions_cv_topfeature.csv"
)
assert_fresher(
  "output/pl2/predictions_cv_topfeature.csv",
  "output/pl2/modeling_frame_regional_partitioned_topfeature.csv"
)
assert_fresher(
  "output/pl2/di_ruler_production.rds",
  "output/pl2/weights_feature_data_partitioning.csv"
)

ruler <- readRDS("output/pl2/di_ruler_production.rds")
profiles_obj <- readRDS("output/pl2/error_profiles.rds")

profiles <- profiles_obj$profiles
profiles_ci <- profiles_obj$profiles_ci
aoa_threshold <- profiles_obj$aoa_threshold
mn <- ruler$metric_names

# The skill axis, frozen upstream. The axis feature must be one of the ruler variables,
# otherwise the per-block scenario subset read below would not carry it.
axis_feature <- profiles_obj$axis_feature
axis_ref <- profiles_obj$axis_ref_production
stopifnot(profiles_obj$axis == "offset", axis_feature %in% mn)

mf <- read_csv("output/pl2/modeling_frame_regional.csv", show_col_types = FALSE)
train <- mf |> filter(scenario == "current")

cat("Metric variables:", length(mn), "\n")
cat("Training rows:", nrow(train), "\n")
cat("AOA threshold (DI):", round(aoa_threshold, 4), "\n")
cat("Skill axis:", axis_feature, "offset, reference", round(axis_ref, 3), "\n")
cat("Profile bins:", nrow(profiles), "spanning offset",
    paste(round(range(profiles$axis_median), 2), collapse = " to "), "\n\n")

## Reference matrix ####

# The training set in the ruler's own space: standardised on the frozen centre/scale, then
# multiplied by the frozen weights. Every distance below is measured against this.
apply_ruler <- function(m) {
  sweep(
    scale(m, center = ruler$scaling$center[mn], scale = ruler$scaling$scale[mn]),
    2, pmax(ruler$weights[mn], 0), "*"
  )
}

train_ruled <- apply_ruler(as.matrix(train[, mn]))

## Interpolators ####

# One linear interpolator per reported quantity, over the bins' median offset. Clamped at
# both ends (see header).
make_interp <- function(col) {
  approxfun(profiles$axis_median, profiles[[col]], rule = 2)
}
interp <- lapply(PROFILE_LAYERS, make_interp)

# The band exists only where the bin had two or more independently fitted models behind it
# (pl2_fitErrorProfiles.R suppresses the rest rather than reporting zero width). With the
# partitions cut on this very axis that is a minority of bins, so the interpolator is built
# from what survives and returns NA outright when there is not enough to interpolate between.
# An NA band on the map is the correct statement here; a fabricated one is not.
ci_band <- profiles_ci |>
  filter(metric == "recall_bog", !is.na(lower), !is.na(upper))

cat("Bins carrying an estimable recall band:", nrow(ci_band), "of",
    nrow(filter(profiles_ci, metric == "recall_bog")), "\n\n")

# Clamped like the point estimates ONLY when every bin carries an interval. If the CI-bearing
# bins are a subset, clamping would pin the band to a different terminal bin than the estimate
# and produce a band entirely above or below its own point estimate -- under the LOPO arm, at
# the Lyngstad bogs' future offset of +3.15, that gave lower = 0.946 against an estimate of
# 0.793. So a partial band is not clamped (`rule = 1`): past its last bin, NA. Under the
# pairwise arm every bin is replicated and the two share one bin set, so both clamp together.
na_fun <- function(x) rep(NA_real_, length(x))
band_rule <- if (nrow(ci_band) == nrow(profiles)) 2 else 1
if (nrow(ci_band) >= 2) {
  interp_lo <- approxfun(ci_band$axis_median, ci_band$lower, rule = band_rule)
  interp_hi <- approxfun(ci_band$axis_median, ci_band$upper, rule = band_rule)
} else {
  interp_lo <- na_fun
  interp_hi <- na_fun
}

cat("Calibrating arm:", profiles_obj$cv_arm, "| band defined over offset",
    paste(round(range(ci_band$axis_median), 2), collapse = " to "),
    if (band_rule == 2) "(clamped with the estimates)\n\n" else "(NA outside)\n\n")

## Per-scenario reliability ####

# Processed in raster blocks rather than all at once: the projection grid holds ~5.3M
# non-NA cells and pulling 32 layers for all of them would need well over a gigabyte of
# doubles at once, for no gain.
map_scenario <- function(name, path) {
  cat("=== scenario:", name, "===\n")
  r <- rast(path)
  stopifnot(all(mn %in% names(r)))
  rsub <- r[[mn]]

  out_names <- c("DI", "offset", "aoa", names(PROFILE_LAYERS),
                 "exp_recall_bog_lo", "exp_recall_bog_hi")
  out <- rast(r[[1]], nlyrs = length(out_names))
  names(out) <- out_names

  out_path <- sprintf("output/pl2/reliability_%s.tif", name)
  b <- blocks(rsub)
  readStart(rsub)
  on.exit(readStop(rsub), add = TRUE)
  writeStart(out, out_path, overwrite = TRUE, datatype = "FLT4S")

  t0 <- Sys.time()
  for (i in seq_len(b$n)) {
    v <- readValues(rsub, b$row[i], b$nrows[i], 1, ncol(rsub), mat = TRUE)
    colnames(v) <- mn

    di <- rep(NA_real_, nrow(v))
    ok <- stats::complete.cases(v)
    if (any(ok)) {
      di[ok] <- FNN::get.knnx(
        train_ruled, apply_ruler(v[ok, , drop = FALSE]),
        k = 1, algorithm = "kd_tree"
      )$nn.dist[, 1] / ruler$train_avg_dist
    }

    # The skill read runs on the signed offset; DI only masks. NA propagates from `v`, so
    # cells outside the land mask stay NA on both.
    off <- v[, axis_feature] - axis_ref
    off[!ok] <- NA_real_

    block_out <- cbind(
      DI = di,
      offset = off,
      aoa = as.numeric(di <= aoa_threshold),
      vapply(interp, function(f) f(off), numeric(length(off))),
      exp_recall_bog_lo = interp_lo(off),
      exp_recall_bog_hi = interp_hi(off)
    )
    writeValues(out, block_out, b$row[i], b$nrows[i])
    cat("  block", i, "of", b$n, "\n")
  }
  writeStop(out)

  cat(sprintf(
    "  done in %.1f min -> %s\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")), out_path
  ))
  rast(out_path)
}

maps <- imap(SCENARIOS, ~ map_scenario(.y, .x))

## Summaries ####

summarise_scenario <- function(r, name) {
  di <- values(r[["DI"]])
  di <- di[!is.na(di)]
  off <- values(r[["offset"]])
  off <- off[!is.na(off)]
  rec <- values(r[["exp_recall_bog"]])
  rec <- rec[!is.na(rec)]
  tibble(
    scenario = name,
    n_cells = length(di),
    DI_median = median(di),
    DI_q90 = quantile(di, 0.90),
    DI_max = max(di),
    pct_outside_aoa = 100 * mean(di > aoa_threshold),
    offset_median = median(off),
    offset_q90 = quantile(off, 0.90),
    offset_max = max(off),
    exp_recall_bog_median = median(rec),
    exp_recall_bog_q10 = quantile(rec, 0.10)
  )
}

summary_tbl <- imap_dfr(maps, ~ summarise_scenario(.x, .y))
cat("\nReliability summary:\n")
print(as.data.frame(summary_tbl), row.names = FALSE, digits = 3)
write_csv(summary_tbl, "output/pl2/reliability_summary.csv", append = FALSE)

# The comparison that matters is not the level but the SHIFT: how much further the future
# domain sits from the training bogs than the present one does, on one axis.
shift <- function(col) {
  round(summary_tbl[[col]][summary_tbl$scenario == "future"] -
          summary_tbl[[col]][summary_tbl$scenario == "current"], 4)
}
cat("\nShift, current -> future: median DI", shift("DI_median"),
    "| median offset", shift("offset_median"), "C\n")

## Reliability at the Lyngstad bogs ####

# The deliverable is about known raised bogs, so report the layers where they are rather
# than only over the whole grid.
lyngstad <- st_read("data/DMraisedbog.gpkg", layer = "lyngstad-MTYPE_A", quiet = TRUE) |>
  st_transform(3035) |>
  vect()

# The band layers are NA wherever the offset falls outside the bins that carry an interval,
# and `fun = mean, na.rm = TRUE` would average whatever few cells remain and present it as if
# it described the bogs -- at the future scenario that was 7 cells of 1,689 reporting
# lower = 0.92 against a point estimate of 0.79. So coverage is measured cell-by-cell first,
# and the band is reported only if most cells actually have one.
MIN_BAND_COVERAGE <- 0.5

band_coverage <- function(r) {
  z <- terra::extract(r[[c("offset", "exp_recall_bog_lo")]], lyngstad)
  ok <- !is.na(z$offset)
  mean(!is.na(z$exp_recall_bog_lo[ok]))
}

lyng_tbl <- imap_dfr(maps, function(r, name) {
  z <- terra::extract(r, lyngstad, fun = mean, na.rm = TRUE)
  cover <- band_coverage(r)
  blank <- function(x) if (cover >= MIN_BAND_COVERAGE) x else NA_real_
  tibble(
    scenario = name,
    n_polygons = sum(!is.na(z$DI)),
    DI_median = median(z$DI, na.rm = TRUE),
    pct_outside_aoa = 100 * mean(z$DI > aoa_threshold, na.rm = TRUE),
    offset_median = median(z$offset, na.rm = TRUE),
    offset_q90 = quantile(z$offset, 0.90, na.rm = TRUE),
    exp_recall_bog = median(z$exp_recall_bog, na.rm = TRUE),
    band_coverage = cover,
    exp_recall_bog_lo = blank(median(z$exp_recall_bog_lo, na.rm = TRUE)),
    exp_recall_bog_hi = blank(median(z$exp_recall_bog_hi, na.rm = TRUE)),
    exp_fpr_bog = median(z$exp_fpr_bog, na.rm = TRUE)
  )
})

cat("\nAt the Lyngstad raised bogs:\n")
print(as.data.frame(lyng_tbl), row.names = FALSE, digits = 3)
write_csv(lyng_tbl, "output/pl2/reliability_at_lyngstad.csv", append = FALSE)

## Figures ####

for (nm in names(maps)) {
  plot(maps[[nm]][["DI"]], main = paste("Dissimilarity index —", nm))
  plot(maps[[nm]][["offset"]], main = paste("bio10 offset from training bogs (C) —", nm))
  plot(maps[[nm]][["exp_recall_bog"]], main = paste("Expected bog recall —", nm))
  plot(maps[[nm]][["aoa"]], main = paste("Inside area of applicability —", nm))
}

png("output/pl2/reliability_di_future.png", width = 900, height = 900, res = 120)
plot(maps[["future"]][["DI"]], main = "Dissimilarity index, future")
dev.off()

png("output/pl2/reliability_offset_future.png", width = 900, height = 900, res = 120)
plot(maps[["future"]][["offset"]], main = "bio10 offset from training bogs (C), future")
dev.off()

png("output/pl2/reliability_recall_future.png", width = 900, height = 900, res = 120)
plot(maps[["future"]][["exp_recall_bog"]], main = "Expected bog recall, future")
dev.off()

cat("\nWrote reliability rasters, summaries and figures to output/pl2/\n")

# sessionInfo ####

sessioninfo::session_info()
