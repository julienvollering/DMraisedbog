# Download CHELSA climate data ####

# PURPOSE: Fetches CHELSA v2.1 bioclim+ layers for the current (1981-2010) and future (2071-2100) periods used as climate predictors.

library(tidyverse)
library(terra)

source("R/config.R")

# Download CHELSA Climate Data
#
# Objective: Download CHELSA climate data for raised bog distribution modeling
# - Current period bioclim (1981-2010, v2.1)
# - Near future bioclim (2071-2100, SSP370, CMIP6)
# - Optional: Paleo data from CHELSA-TRACE
# - Data source: Direct download from CHELSA bucket (>https://os.unil.cloud.switch.ch/chelsa02/>; <https://os.zhdk.cloud.switch.ch/chelsa01/chelsa_trace21k/>)
# - Output: GeoTIFF files organized by time period in data/CHELSA/

## Setup ####

### Output directories ####

# Create directory structure for downloaded data
dir.create("data/CHELSA", showWarnings = FALSE)
dir.create("data/CHELSA/1981-2010", showWarnings = FALSE)
dir.create("data/CHELSA/2071-2100", showWarnings = FALSE)
dir.create("data/CHELSA/paleo", showWarnings = FALSE)

## Configuration ####

### Time periods ####

current_period <- "1981-2010"
future_period <- "2071-2100"

### Climate scenarios ####

# Climate scenario for future projections (CHELSA CMIP6: ssp126, ssp370, ssp585)
scenario <- FUTURE_SCENARIO # R/config.R, shared with pl0_collatePredictors.R
ssp_scenarios <- c(scenario)

### Spatial extent ####

# Approximate Norway bounding box in WGS84
# TODO: Fine-tune extent if needed for specific study area
norway_extent <- ext(c(
  xmin = 4,
  xmax = 32,
  ymin = 57,
  ymax = 72
))

## Variable discovery ####

# Note: Rchelsa package does not provide access to the data we need.
# Proceeding with direct download from CHELSA bucket.

### Base URLs ####

# CHELSA data bucket
chelsa_base_url <- "https://os.unil.cloud.switch.ch/chelsa02"

# Subdirectories for different data types
# Current period structure: chelsa/global/bioclim/{variable}/{timeperiod}/
current_bioclim_url <- file.path(chelsa_base_url, "chelsa/global/bioclim")

# Future period structure: chelsa/global/bioclim/{variable}/{timeperiod}/{model}/{scenario}/
future_bioclim_url <- file.path(chelsa_base_url, "chelsa/global/bioclim")

### Climate model selection ####

# For future scenarios, we need model-specific paths
# Common GCMs available: GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL
gcm_model <- FUTURE_GCM # R/config.R, shared with pl0_collatePredictors.R

### Automatic variable discovery ####

cat("Discovering available variables from CHELSA bucket...\n")

# Known CHELSA bioclim variables based on documentation
# See: <https://www.chelsa-climate.org/datasets/chelsa_bioclim>
candidate_vars <- c(
  sprintf("bio%02d", 1:19), # Standard 19 bioclimatic variables (bio01-bio19)
  "fcf", # Frost change frequency
  "gdd0", # Growing degree days (base 0°C)
  "gdd5", # Growing degree days (base 5°C)
  "gdd10", # Growing degree days (base 10°C)
  "gsl", # Growing season length
  "gsp", # Growing season precipitation
  "gst", # Growing season temperature
  "ngd0", # Number of growing days (base 0°C)
  "ngd5", # Number of growing days (base 5°C)
  "ngd10", # Number of growing days (base 10°C)
  "npp", # Net primary productivity
  "scd", # Snow cover days
  "swe" # Snow water equivalent
)

# Helper function to check if a file exists at URL
url_exists <- function(url) {
  tryCatch(
    {
      con <- url(url, "rb")
      close(con)
      TRUE
    },
    error = function(e) FALSE,
    warning = function(w) FALSE
  )
}

# Test which variables exist in both current and future periods
cat("  Testing variable availability (this may take a minute)...\n")

available_vars <- character(0)
for (var in candidate_vars) {
  # Construct test URLs for current and future
  current_file <- paste0("CHELSA_", var, "_1981-2010_V.2.1.tif")
  current_url <- file.path(current_bioclim_url, var, "1981-2010", current_file)

  future_file <- paste0(
    "CHELSA_",
    tolower(gcm_model),
    "_", scenario, "_",
    var,
    "_2071-2100_V.2.1.tif"
  )
  future_url <- file.path(
    future_bioclim_url,
    var,
    "2071-2100",
    gcm_model,
    scenario,
    future_file
  )

  # Check if both files exist
  if (url_exists(current_url) && url_exists(future_url)) {
    available_vars <- c(available_vars, var)
    cat("    Found:", var, "\n")
  }
}

cat("\n  Variables available in both periods:", length(available_vars), "\n")

common_vars <- available_vars

### Variable selection ####

# Use all discovered common variables
variables_to_download <- common_vars

## Data download ####

### Current period ####

cat(
  "Downloading current period (",
  current_period,
  ") variables...\n",
  sep = ""
)

for (var in variables_to_download) {
  cat("  Downloading", var, "...\n")

  # Construct filename following CHELSA naming convention
  # Format: CHELSA_bio01_1981-2010_V.2.1.tif or CHELSA_fcf_1981-2010_V.2.1.tif
  chelsa_filename <- paste0("CHELSA_", var, "_1981-2010_V.2.1.tif")

  # Full URL to the file
  # Structure: chelsa/global/bioclim/{variable}/{timeperiod}/filename
  file_url <- file.path(current_bioclim_url, var, "1981-2010", chelsa_filename)

  # Local output path
  output_file <- file.path("data/CHELSA", current_period, chelsa_filename)

  # Skip if file already exists
  if (file.exists(output_file)) {
    cat("    Already exists, skipping\n")
    next
  }

  # Download the file with extended timeout (5 minutes for large files)
  download_status <- tryCatch(
    {
      download.file(
        file_url,
        output_file,
        mode = "wb",
        quiet = TRUE,
        timeout = 300
      )
    },
    error = function(e) {
      return(1)
    }
  )

  # Check if download was successful by verifying file exists and has content
  if (file.exists(output_file) && file.size(output_file) > 0) {
    cat("    Saved to", output_file, "\n")
  } else {
    cat("    ERROR: Download failed or file is empty\n")
    if (file.exists(output_file)) {
      file.remove(output_file) # Remove empty file
    }
  }
}

cat("\n")

### Future period ####

cat(
  "Downloading future period (",
  future_period,
  ", ", toupper(scenario), ") variables...\n",
  sep = ""
)

for (var in variables_to_download) {
  cat("  Downloading", var, "...\n")

  # Construct filename following CHELSA CMIP6 naming convention
  # Format: CHELSA_gfdl-esm4_<scenario>_bio01_2071-2100_V.2.1.tif
  chelsa_filename <- paste0(
    "CHELSA_",
    tolower(gcm_model),
    "_", scenario, "_",
    var,
    "_2071-2100_V.2.1.tif"
  )

  # Full URL to the file
  # Structure: chelsa/global/bioclim/{variable}/{timeperiod}/{MODEL}/{scenario}/filename
  file_url <- file.path(
    future_bioclim_url,
    var,
    "2071-2100",
    gcm_model,
    scenario,
    chelsa_filename
  )

  # Local output path
  output_file <- file.path("data/CHELSA", future_period, chelsa_filename)

  # Skip if file already exists
  if (file.exists(output_file)) {
    cat("    Already exists, skipping\n")
    next
  }

  # Download the file with extended timeout (5 minutes for large files)
  download_status <- tryCatch(
    {
      download.file(
        file_url,
        output_file,
        mode = "wb",
        quiet = TRUE,
        timeout = 300
      )
    },
    error = function(e) {
      return(1)
    }
  )

  # Check if download was successful by verifying file exists and has content
  if (file.exists(output_file) && file.size(output_file) > 0) {
    cat("    Saved to", output_file, "\n")
  } else {
    cat("    ERROR: Download failed or file is empty\n")
    if (file.exists(output_file)) {
      file.remove(output_file) # Remove empty file
    }
  }
}

cat("\n")

### Paleo period ####

#### Time slices ####

# Specify paleo time periods of interest
# CHELSA-TRACE timeline file naming:
# -200 = 20 ka (Last Glacial Maximum)
# -120 = 12 ka (Younger Dryas)
# -60 = 6 ka (Mid-Holocene)
# 0000 = 2 ka
# 0020 = 0 ka (present)

# Example: Download every 500 years from 16 ka to present
# See <https://doi.org/10.5194/cp-19-439-2023> fig. S1
interval_years <- 500
paleo_time_slices_hundreds <- c(
  seq(from = -140, to = 0, by = interval_years / 100), # 16 ka to 3 ka
  seq(from = 0, to = 20, by = interval_years / 100) # 2 ka to 0 ka
)

#### Variables ####

# CHELSA-TRACE Bioclim variables
# See: <https://www.chelsa-climate.org/datasets/chelsa-trace21k-centennial-bioclim>
paleo_bioclim_vars <- c(
  "bio01", # Annual mean temperature
  "bio12", # Annual precipitation
  "glz", # Ice Sheet Surface Altitude
  "orog" # Surface Altitude. Altitude is the (geometric) height above the geoid which is the reference geopotential surface. The geoid is similar to mean sea level.
)

#### Download loop ####

# Base URL for CHELSA-TRACE (note: different server than CHELSA v2!)
record_settings(
  "R/pl0_downloadCHELSA.R",
  scenario = scenario,
  gcm_model = gcm_model,
  current_period = current_period,
  future_period = future_period,
  paleo_interval_years = interval_years,
  paleo_slices_hundreds = range(paleo_time_slices_hundreds),
  paleo_bioclim_vars = paleo_bioclim_vars
)

trace_base_url <- "https://os.zhdk.cloud.switch.ch/chelsa01/chelsa_trace21k/global/bioclim"

# Helper function to format time values for CHELSA-TRACE
format_trace_time <- function(time_hundreds) {
  if (time_hundreds < 0) {
    # Negative values: sign + 3 digits with leading zeros (e.g., "-090")
    return(sprintf("-%03d", abs(time_hundreds)))
  } else {
    # Zero and positive: 4 digits with leading zeros (e.g., "0000", "0020")
    return(sprintf("%04d", time_hundreds))
  }
}

if (length(paleo_time_slices_hundreds) > 0 && length(paleo_bioclim_vars) > 0) {
  cat("Downloading CHELSA-TRACE paleoclimate data...\n")

  for (time_hundreds in paleo_time_slices_hundreds) {
    for (var in paleo_bioclim_vars) {
      # Convert to actual calendar year for display
      time_year <- time_hundreds * 100
      if (time_year < 0) {
        time_label <- paste0(abs(time_year), " BCE")
      } else {
        time_label <- paste0(time_year, " CE")
      }
      cat("  Downloading", var, "for year", time_label, "...\n")

      # Format time string for filename
      time_str <- format_trace_time(time_hundreds)

      # Construct filename following CHELSA-TRACE naming convention
      # Format: CHELSA_TraCE21k_bio01_-190_V.1.0.tif
      chelsa_filename <- paste0(
        "CHELSA_TraCE21k_",
        var,
        "_",
        time_str,
        "_V.1.0.tif"
      )

      # Full URL to the file
      # Structure: chelsa_trace21k/global/bioclim/{variable}/filename
      file_url <- file.path(
        trace_base_url,
        var,
        chelsa_filename
      )

      # Local output path
      output_file <- file.path("data/CHELSA/paleo", chelsa_filename)

      # Skip if file already exists
      if (file.exists(output_file)) {
        cat("    Already exists, skipping\n")
        next
      }

      # Download the file with extended timeout (5 minutes for large files)
      download_status <- tryCatch(
        {
          download.file(
            file_url,
            output_file,
            mode = "wb",
            quiet = TRUE,
            timeout = 300
          )
        },
        error = function(e) {
          return(1)
        }
      )

      # Check if download was successful by verifying file exists and has content
      if (file.exists(output_file) && file.size(output_file) > 0) {
        cat("    Saved to", output_file, "\n")
      } else {
        cat("    ERROR: Download failed or file is empty\n")
        if (file.exists(output_file)) {
          file.remove(output_file) # Remove empty file
        }
      }
    }
  }
} else {
  cat("Paleo download skipped (no variables specified)\n")
}

cat("\n")

## File integrity check ####

# READ THE WINDOW THE PIPELINE ACTUALLY USES. This check previously read a 21 x 21 block
# at the centre of each global raster -- which is 0 N, 0 E, in the Atlantic off Africa,
# roughly 5,000 km from the study area and 441 cells out of ~900 million. Every file
# passed it. Three future files (bio19, gdd0, gdd5) were nonetheless truncated mid-tile
# over Europe by an interrupted download, and the failure surfaced 130 minutes into
# pl0_collatePredictors.R instead of here (notebook 2026-09-16).
#
# Forcing a full read of the European window is the only test that finds that class of
# corruption, because a partial download leaves a valid header and valid tiles everywhere
# it did finish. It costs a couple of minutes against the hours a late failure costs.
#
# This STOPS rather than repairing: deleting a data file is the user's call, and the
# download loops above skip any file that already exists, so removing the named files and
# re-running this script is the repair.
check_ext <- ext(-10, 35, 35, 72) # the crop pl0_collatePredictors.R applies

chelsa_files <- list.files(
  "data/CHELSA/",
  recursive = TRUE,
  pattern = "CHELSA_.*\\.tif$",
  full.names = TRUE
)

cat("\nIntegrity check: reading the European window of", length(chelsa_files), "files\n")

issues <- list()
for (file_path in chelsa_files) {
  result <- tryCatch(
    {
      r <- rast(file_path)
      # Paleo TraCE files are regional and may not span the window; intersect first.
      win <- intersect(ext(r), check_ext)
      if (!is.null(win)) {
        invisible(global(crop(r, win), "max", na.rm = TRUE))
      }
      NULL
    },
    error = function(e) list(type = "ERROR", msg = conditionMessage(e)),
    warning = function(w) list(type = "WARNING", msg = conditionMessage(w))
  )
  if (!is.null(result)) {
    issues[[basename(file_path)]] <- result
  }
}

if (length(issues) > 0) {
  cat("\nFile integrity issues found:\n")
  for (fname in names(issues)) {
    cat("  ", issues[[fname]]$type, ": ", fname, "\n", sep = "")
    cat("    Message: ", issues[[fname]]$msg, "\n", sep = "")
  }
  stop(
    "Corrupt CHELSA file(s): ", paste(names(issues), collapse = ", "),
    ". Delete them from data/CHELSA/ and re-run this script to re-download."
  )
}
cat("All files passed integrity check over the European window\n")

# Cleanup
gc()
terra::tmpFiles(remove = TRUE)

# sessionInfo ####

sessioninfo::session_info()
