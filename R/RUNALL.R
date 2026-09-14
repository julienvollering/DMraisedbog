# RUNALL: the analysis from start to finish ####

# Knits every pipeline script to HTML, in dependency order. This file is the table of
# contents for the project: if a script is not listed here it is not part of the
# analysis, and it belongs in R/archive/.
#
# STAGES
#   pl0  shared upstream — raw data, predictor stacks, and the two labelled training
#        blocks (Norway and Europe) that the pooled frame is built from.
#   pl2  everything else — pooled 3-class modelling frame, the frozen DI ruler, feature
#        weights, partitioning, cross-validation, production fit, projection,
#        interpretation, and reliability mapping.
#
# WHY THERE IS NO pl1. The pl1 scripts were the hierarchical two-scale cascade, in which
# a global EU model's prediction entered the local model as the `rf_global` covariate.
# That architecture was retired (plan_EUintegration.md section 1.1): the global model
# trained on the same predictor stack at coarser grain, so the covariate added no
# predictor axis, and it put the local model's DI on two training sets at once, which
# makes the reliability map uninterpretable. The scripts are kept in R/archive/ for
# reference. The pl2 name is left as it is so that the numbering in notebook.md and
# plan_EUintegration.md keeps pointing at the right things.
#
# WHY THERE IS NO pl3 EITHER. pl3 was to be a second, independent cross-validation whose
# purpose was to generate a (novelty, skill) calibration dataset — because at the time it
# was designed, CV partitions could not reach the novelty of the future projection, so
# ordinary CV could say nothing about how the projection would behave. Pooling the EU
# block into training and cutting partitions along bio10 removed that constraint. Measured
# on one axis, pl2's cross-validation now spans the future: its 90th percentile of DI
# (0.799) covers 96.2% of future cells and 100% of future bog cells, and 107k CV test rows
# (478 of them bog) sit above the future's median bog novelty of 0.398.
# `predictions_cv_topfeature.csv` IS the calibration dataset, so pl3_assembleAnchors.R,
# pl3_designCV.R and pl3_crossvalidate.R were dropped rather than written; the two
# reliability scripts survive as the tail of pl2 below. pl3_extractEUpredictors.R is in
# R/archive/ — it existed to score EU presences as held-out anchors, and they are training
# data now.
#
# R/functions.R holds the shared helpers (stratified absence draw, land-use screen,
# weighted DI, multiclass metrics) and is sourced by the scripts that need it.
# R/config.R holds the constants more than one script must agree on (future scenario and
# GCM, response levels) and the settings registry every script reports its constants to.

library(rmarkdown)

## Stage pl0 — shared upstream ####

# Ordered by dependency, not by name. Three constraints worth knowing:
#  - the EU domain runs before the EPM mask and the land-use screen, because it is what
#    bounds their warps;
#  - the two block scripts must follow pl0_rasterizePresenceAbsence.R, which establishes
#    the Norwegian survey footprint, and pl0_buildLandUseScreen.R, whose screen decides
#    absence eligibility in both;
#  - Norway runs before Europe, because it publishes the warm-flank threshold (the maximum
#    bio10 its *training* data reach) against which the EU draw is defined.
scripts_pl0 <- c(
  # Survey data
  "R/pl0_collateLyngstad.R",            # merge the raw Lyngstad raised-bog polygons
  "R/pl0_collateLyngstadExtent.R",      # reconstruct the area actually surveyed
  # Climate and paleo predictors
  "R/pl0_downloadCHELSA.R",             # fetch CHELSA current + future bioclim
  "R/pl0_prepareCHELSATrace.R",         # derive paleo_* predictors from TraCE21k
  "R/pl0_collatePredictors.R",          # assemble the 5 km EU and 250 m Norway stacks
  # EU domain and its layers
  "R/pl0_buildEUdomain.R",              # delimit the EU absence domain
  "R/pl0_buildEPMmask.R",               # EPM2025 map_cat -> EU absence label split
  "R/pl0_buildEUterrain.R",             # EU elevation and slope at 250 m (elevatr)
  "R/pl0_buildLandUseScreen.R",         # WorldCover human/water screen, both blocks
  "R/pl0_collateEuropeanRaisedBog.R",   # EU presences + exclusion mask M
  # Presence/absence cells and the two training blocks
  "R/pl0_rasterizePresenceAbsence.R",   # bog polygons -> presence/absence cells
  "R/pl0_labelNorwayBlock.R",           # Norway 3-class labels + stratified absence draw
  "R/pl0_sampleEUabsences.R"            # EU absence draw, same shared rule
)

## Stage pl2 — production pipeline ####

# pl2_freezeDIRuler.R sits third on purpose. Everything downstream that measures a
# distance — the CV, the projection's novelty, the reliability curves — must use ONE
# ruler, and deriving it inside pl2_predict.R (which runs last) left pl2_evaluate.R
# normalising each fold by its own constant. Those differ by a factor of 1.75, which is
# what previously made CV DI incomparable with the DI of the future map.
scripts_pl2 <- c(
  "R/pl2_createModelingFrame.R",           # pool the two blocks; write scenario stacks
  "R/pl2_weightFeaturesDataPartitioning.R",# feature weights for every weighted distance
  "R/pl2_freezeDIRuler.R",                 # freeze weights + scaling + normalising constant
  "R/pl2_exploreOccupancy.R",              # model-free: where the future goes, what is observed there
  "R/pl2_partitionDataByTopFeature.R",     # cut CV partitions along the top shifting feature
  "R/pl2_evaluate.R",                      # pairwise CV: skill, confusion, per-row DI
  "R/pl2_predict.R",                       # production fit and projection
  "R/pl2_interpret.R",                     # Lyngstad class transitions and change in P(bog)
  # Reliability mapping. Both read the CV output above rather than generating their own,
  # which is the whole reason pl3 collapsed into this stage.
  "R/pl2_fitErrorProfiles.R",              # skill vs DI, cluster-bootstrap CIs, AOA cut
  "R/pl2_mapReliability.R"                 # apply those curves to the projection domain
)

scripts <- c(scripts_pl0, scripts_pl2)

## Run ####

# Set to a stage (scripts_pl2) or a slice for a targeted run; `scripts` runs everything.
# Named rather than indexed so inserting a script upstream cannot silently repoint it.
to_run <- scripts

## Snapshot of the previous run's summaries ####

# The small tables that carry the results are copied aside before anything runs and
# diffed against the new ones at the end. Unchanged code should give an empty diff (the
# fold fits are seeded); a changed number then points at exactly the script to re-read.
SUMMARY_FILES <- c(
  "output/pl0/no_block_summary.csv",
  "output/pl0/eu_mask_summary.csv",
  "output/pl0/eu_domain_summary.csv",
  "output/pl2/modeling_frame_summary.csv",
  "output/pl2/modeling_frame_predictor_ranges.csv",
  "output/pl2/weights_feature_data_partitioning.csv",
  "output/pl2/di_ruler_summary.csv",
  "output/pl2/occupancy_grid.csv",
  "output/pl2/metrics_cv_topfeature.csv",
  "output/pl2/cv_novelty_coverage.csv",
  "output/pl2/oob_metrics_production.csv",
  "output/pl2/lyngstad_transitions.csv",
  "output/pl2/error_profiles_by_bin.csv",
  "output/pl2/error_profiles_ci.csv",
  "output/pl2/reliability_summary.csv",
  "output/pl2/reliability_at_lyngstad.csv"
)
PREVIOUS_DIR <- "output/_previous"
dir.create(PREVIOUS_DIR, showWarnings = FALSE, recursive = TRUE)
for (f in SUMMARY_FILES) {
  if (file.exists(f)) {
    file.copy(f, file.path(PREVIOUS_DIR, basename(f)), overwrite = TRUE)
  }
}

# One line per summary file: how many cells changed, the largest absolute numeric change,
# and which columns moved. Non-numeric columns count as changed on any inequality.
diff_summary <- function(f, tol = 1e-9) {
  prev <- file.path(PREVIOUS_DIR, basename(f))
  if (!file.exists(f)) {
    return(sprintf("%-52s not written by this run", f))
  }
  if (!file.exists(prev)) {
    return(sprintf("%-52s no previous copy", f))
  }
  a <- read.csv(prev, stringsAsFactors = FALSE, check.names = FALSE)
  b <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) {
    return(sprintf(
      "%-52s shape changed: %d x %d -> %d x %d", f, nrow(a), ncol(a), nrow(b), ncol(b)
    ))
  }
  changed <- vapply(names(b), function(k) {
    x <- a[[k]]
    y <- b[[k]]
    if (is.numeric(x) && is.numeric(y)) {
      sum(xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & abs(x - y) > tol))
    } else {
      sum(xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y))
    }
  }, integer(1))
  max_abs <- max(c(0, unlist(lapply(names(b), function(k) {
    x <- a[[k]]
    y <- b[[k]]
    if (is.numeric(x) && is.numeric(y)) abs(x - y)[!is.na(x) & !is.na(y)] else numeric(0)
  }))))
  if (sum(changed) == 0) {
    return(sprintf("%-52s unchanged", f))
  }
  sprintf(
    "%-52s %d cells changed, max |diff| %.4g, in: %s",
    f, sum(changed), max_abs, paste(names(changed)[changed > 0], collapse = ", ")
  )
}

## Run manifest ####

# One row per script per run: when it started and finished, how long it took, whether it
# succeeded, and every file under output/ whose modification time falls inside that
# window -- which is what the script produced. Appended after every script, so a run that
# crashes still leaves the manifest of everything before the crash. It is the map from
# script to output file that otherwise needs grepping, and the end-of-run report below
# uses it to list files under output/ that nothing in the run wrote.
MANIFEST_PATH <- "output/run_manifest.csv"
run_id <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

output_files <- function() {
  f <- list.files("output", recursive = TRUE, full.names = TRUE)
  f <- f[!grepl("^output/(_previous/|run_manifest)", f)]
  data.frame(path = sub("^output/", "", f), mtime = file.mtime(f), stringsAsFactors = FALSE)
}

record_run <- function(script, t0, t1, status) {
  files <- output_files()
  touched <- files$path[files$mtime >= t0 & files$mtime <= t1]
  row <- data.frame(
    run = run_id,
    script = script,
    start = format(t0, "%Y-%m-%d %H:%M:%S"),
    end = format(t1, "%Y-%m-%d %H:%M:%S"),
    minutes = round(as.numeric(difftime(t1, t0, units = "mins")), 2),
    status = status,
    n_files = length(touched),
    files = paste(touched, collapse = ";"),
    stringsAsFactors = FALSE
  )
  new_file <- !file.exists(MANIFEST_PATH)
  write.table(
    row, MANIFEST_PATH,
    sep = ",", row.names = FALSE, col.names = new_file, append = !new_file
  )
  cat(sprintf("  %.1f min, %d output files touched\n", row$minutes, row$n_files))
}

# Every render() runs in THIS R process, so all the scripts share one tempdir and one
# terra scratch pool. terra only clears its spill files when a session ends, so across a
# full run they accumulate instead of turning over: on 2026-08-21 pl0_collatePredictors.R
# alone reached 15 GB of `spat_*.tif` and filled the disk mid-pipeline. Clearing after
# each script holds the high-water mark at whatever a single script needs rather than the
# sum over all of them. `tmpFiles()` only ever removes terra's own scratch files, and only
# those orphaned by this process -- outputs already written are untouched.
for (script in to_run) {
  t0 <- Sys.time()
  cat("Executing:", script, format(t0, "%H:%M:%S"), "\n")
  err <- tryCatch(
    {
      render(script, output_format = "html_document", knit_root_dir = "../")
      NULL
    },
    error = function(e) e,
    finally = {
      # gc() first so SpatRasters left behind by the render environment are finalised
      # and release their files; only then are the scratch files safe to delete.
      gc(verbose = FALSE)
      terra::tmpFiles(current = TRUE, orphan = TRUE, remove = TRUE)
    }
  )
  t1 <- Sys.time()
  record_run(script, t0, t1, if (is.null(err)) "ok" else "error")
  if (!is.null(err)) {
    cat("✗ ERROR in", script, ":\n", conditionMessage(err), "\n\n")
    stop("Script execution failed at: ", script)
  }
  cat("✓ Completed:", script, format(t1, "%H:%M:%S"), "\n\n")
}

## Manifest report: what this run did not write ####

# Only meaningful after a full run. The cached warp tiles (*_parts/) are skipped on
# purpose: they are rebuilt only when force_rebuild is set in their scripts, so an
# untouched part is expected. Anything else listed here is either stale (written by a
# script that no longer exists or no longer writes it) or written outside the pipeline.
if (identical(to_run, scripts)) {
  manifest <- read.csv(MANIFEST_PATH, stringsAsFactors = FALSE)
  manifest <- manifest[manifest$run == run_id, ]
  claimed <- unique(unlist(strsplit(manifest$files[manifest$files != ""], ";")))
  untouched <- setdiff(output_files()$path, claimed)
  untouched <- untouched[!grepl("_parts/", untouched)]
  cat("Files under output/ that no script of this run wrote:", length(untouched), "\n")
  if (length(untouched) > 0) {
    cat(paste0("- ", untouched), sep = "\n")
  }
}

## Settings of this run ####

# Every script with design constants registers them through record_settings() in
# R/config.R; this is the one table that says how the run was configured.
if (file.exists("output/pipeline_settings.csv")) {
  settings <- read.csv("output/pipeline_settings.csv", stringsAsFactors = FALSE)
  cat("\nPipeline settings (output/pipeline_settings.csv):\n")
  print(settings[order(settings$script), ], row.names = FALSE, right = FALSE)
}

## Diff against the previous run ####

cat("\nSummary tables against the previous run (output/_previous/):\n")
cat(vapply(SUMMARY_FILES, diff_summary, character(1)), sep = "\n")
cat("\n")

## Guard: nothing in R/ should be outside the pipeline ####

# R/archive/ is excluded on purpose -- that is where retired and exploratory scripts go.
# Anything else showing up here is either a new script that has not been wired in, or one
# that should have been archived.
all_r_scripts <- list.files("R", pattern = "[.]R$", full.names = TRUE)
unused_scripts <- setdiff(
  all_r_scripts, c(scripts, "R/RUNALL.R", "R/functions.R", "R/config.R")
)

if (length(unused_scripts) > 0) {
  cat("Scripts in R/ that are not part of any stage:\n")
  for (script in unused_scripts) {
    cat("- ", script, "\n")
  }
  cat("Wire them into a stage above, or move them to R/archive/.\n")
}
