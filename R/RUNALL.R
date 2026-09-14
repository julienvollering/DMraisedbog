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
#        interpretation, and (to come) reliability mapping.
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

# Every render() runs in THIS R process, so all 25 scripts share one tempdir and one
# terra scratch pool. terra only clears its spill files when a session ends, so across a
# full run they accumulate instead of turning over: on 2026-08-21 pl0_collatePredictors.R
# alone reached 15 GB of `spat_*.tif` and filled the disk mid-pipeline. Clearing after
# each script holds the high-water mark at whatever a single script needs rather than the
# sum over all of them. `tmpFiles()` only ever removes terra's own scratch files, and only
# those orphaned by this process -- outputs already written are untouched.
for (script in to_run) {
  cat("Executing:", script, format(Sys.time(), "%H:%M:%S"), "\n")
  tryCatch(
    {
      render(script, output_format = "html_document", knit_root_dir = "../")
      cat("✓ Completed:", script, format(Sys.time(), "%H:%M:%S"), "\n\n")
    },
    error = function(e) {
      cat("✗ ERROR in", script, ":\n", e$message, "\n\n")
      stop("Script execution failed at: ", script)
    },
    finally = {
      # gc() first so SpatRasters left behind by the render environment are finalised
      # and release their files; only then are the scratch files safe to delete.
      gc(verbose = FALSE)
      terra::tmpFiles(current = TRUE, orphan = TRUE, remove = TRUE)
    }
  )
}

## Guard: nothing in R/ should be outside the pipeline ####

# R/archive/ is excluded on purpose -- that is where retired and exploratory scripts go.
# Anything else showing up here is either a new script that has not been wired in, or one
# that should have been archived.
all_r_scripts <- list.files("R", pattern = "[.]R$", full.names = TRUE)
unused_scripts <- setdiff(all_r_scripts, c(scripts, "R/RUNALL.R", "R/functions.R"))

if (length(unused_scripts) > 0) {
  cat("Scripts in R/ that are not part of any stage:\n")
  for (script in unused_scripts) {
    cat("- ", script, "\n")
  }
  cat("Wire them into a stage above, or move them to R/archive/.\n")
}
