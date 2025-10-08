library(rmarkdown)

# Define script execution order (maintains dependencies)
scripts_shared <- c(
  "R/pl0_collateLyngstad.R",
  "R/pl0_collateLyngstadExtent.R", 
  "R/pl0_collateEuropeanRaisedBog.R",
  "R/pl0_collatePredictors.R",
  "R/pl0_rasterizePresenceAbsence.R",
  "R/pl0_spatiallyThinCells.R",
  "R/pl0_modelGlobalScale.R")
scripts_pipeline1 <- c(
  "R/pl1_weightFeaturesDataPartitioning.R",
  "R/pl1_partitionData.R", # Avoid open GIS locking TIFF-files
  "R/pl1_tuneHyperparametersLocal.R",
  "R/pl1_modelLocalScale.R",
  "R/pl1_evaluateLocalScale.R",
  "R/pl1_modelLocalScaleFuture.R", # Avoid open GIS locking TIFF-files
  "R/pl1_mapExtrapolation.R",
  "R/pl1_interpretLocalScale.R"
)
scripts_pipeline2 <- c(
  "R/pl2_createModelingFrame.R",
  "R/pl2_weightFeaturesDataPartitioning.R",
  "R/pl2_exploreFeatureSpaceDistances.R",
  # Three parallel partitioning schemes for performance bounds:
  "R/pl2_partitionDataByImportance.R",  # Optimistic bound
  "R/pl2_partitionDataByPCA.R",         # Best estimate
  "R/pl2_partitionDataByMaxShift.R",    # Pessimistic bound
  "R/pl2_tuneHyperparameters.R",
  "R/pl2_evaluate.R"
)

# Combine all scripts
scripts <- c(scripts_shared, scripts_pipeline1, scripts_pipeline2)

# Execute scripts with error handling
for (script in scripts) {
  cat("Executing:", script, "\n")
  tryCatch({
    render(script, output_format = "html_document", knit_root_dir = "../")
    cat("✓ Completed:", script, "\n\n")
  }, error = function(e) {
    cat("✗ ERROR in", script, ":\n", e$message, "\n\n")
    stop("Script execution failed at: ", script)
  })
}

# Detect scripts that are not used in any pipeline
all_r_scripts <- list.files("R", pattern = "\\.R$", full.names = TRUE)
used_scripts <- scripts
unused_scripts <- setdiff(all_r_scripts, used_scripts)

if (length(unused_scripts) > 0) {
  cat("Scripts found in R/ directory that are not included in any pipeline:\n")
  for (script in unused_scripts) {
    cat("- ", script, "\n")
  }
}

# for (script in unused_scripts) {
#   cat("Executing:", script, "\n")
#   tryCatch({
#     render(script, output_format = "html_document", knit_root_dir = "../")
#     cat("✓ Completed:", script, "\n\n")
#   }, error = function(e) {
#     cat("✗ ERROR in", script, ":\n", e$message, "\n\n")
#     stop("Script execution failed at: ", script)
#   })
# }
