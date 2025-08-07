library(rmarkdown)

# Create output directory if it doesn't exist
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
  cat("Created output/ directory\n")
}

# Define script execution order (maintains dependencies)
scripts <- c(
  "R/collateLyngstad.R",
  "R/collateLyngstadExtent.R", 
  "R/collateEuropeanRaisedBog.R",
  "R/collatePredictors.R",
  # "R/spatiallyThinGeometries.R",
  "R/rasterizePresenceAbsence.R",
  "R/spatiallyThinCells.R",
  "R/pilotStructuredCV.R"
)

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
