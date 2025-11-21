library(tidyverse)
library(terra)

# Prepare CHELSA-TraCE21k Paleoclimate Data
#
# Objective: Process CHELSA-TraCE21k paleoclimate data to derive predictors
# - Consecutive years as ice-free land (calculated from present -- ice-free and above sea level)
# - Mean bio01 and bio12 during ice-free land period
# - Sum bio01 and bio12 during ice-free land period
# - Mean bio01 and bio12 during entire 16ka period
# - Sum bio01 and bio12 during entire 16ka period
#
# Ice-free definition: glz <= 10m AND orog > 0
# Output: Native CHELSA resolution (~1km), cropped to European extent

## Setup ####

dir.create("output/pl0", showWarnings = FALSE, recursive = TRUE)
dir.create("output/pl0/chelsa_trace", showWarnings = FALSE, recursive = TRUE)

# WGS84 extent for cropping CHELSA data
# Matches extent used in pl0_collatePredictors.R
eu_extent_wgs84 <- ext(c(xmin = -10, xmax = 35, ymin = 35, ymax = 72))

# Ice-free: ice sheet surface altitude <= 10m (tolerance for data uncertainty)
ice_threshold <- 10
# Above sea level: surface altitude > 0m
sea_threshold <- 0

## Data loading ####

cat("Loading CHELSA-TraCE21k data...\n")

### List paleoclimate files ####

paleo_dir <- "data/CHELSA/paleo"

# Get all downloaded files
bio01_files <- list.files(
  paleo_dir,
  pattern = "CHELSA_TraCE21k_bio01_.*\\.tif$",
  full.names = TRUE
)
bio12_files <- list.files(
  paleo_dir,
  pattern = "CHELSA_TraCE21k_bio12_.*\\.tif$",
  full.names = TRUE
)
glz_files <- list.files(
  paleo_dir,
  pattern = "CHELSA_TraCE21k_glz_.*\\.tif$",
  full.names = TRUE
)
orog_files <- list.files(
  paleo_dir,
  pattern = "CHELSA_TraCE21k_orog_.*\\.tif$",
  full.names = TRUE
)

cat("  Found", length(bio01_files), "bio01 files\n")
cat("  Found", length(bio12_files), "bio12 files\n")
cat("  Found", length(glz_files), "glz files\n")
cat("  Found", length(orog_files), "orog files\n")

### Parse time values from filenames ####

# Extract time values from filenames
# Format: CHELSA_TraCE21k_{var}_{time}_V.1.0.tif
# Time: -140 (16ka), 0000 (2ka), 0020 (0ka)
parse_trace_time <- function(filepath) {
  basename(filepath) |>
    str_extract("(?<=_)[-0-9]+(?=_V\\.1\\.0\\.tif$)") |>
    as.integer()
}

# Create data frame with file paths and time values
bio01_df <- tibble(
  time = parse_trace_time(bio01_files),
  file = bio01_files
) |>
  arrange(time)

bio12_df <- tibble(
  time = parse_trace_time(bio12_files),
  file = bio12_files
) |>
  arrange(time)

glz_df <- tibble(
  time = parse_trace_time(glz_files),
  file = glz_files
) |>
  arrange(time)

orog_df <- tibble(
  time = parse_trace_time(orog_files),
  file = orog_files
) |>
  arrange(time)

# Verify all variables have same time slices
if (
  !identical(bio01_df$time, bio12_df$time) ||
    !identical(bio01_df$time, glz_df$time) ||
    !identical(bio01_df$time, orog_df$time)
) {
  stop("Time slices do not match across variables")
}

time_slices <- bio01_df$time
n_times <- length(time_slices)

cat(
  "  Time range:",
  min(time_slices),
  "to",
  max(time_slices),
  "(hundred-year units)\n"
)
cat("  Number of time slices:", n_times, "\n\n")

### Load raster data ####

cat("Loading and cropping rasters to European extent...\n")

# Create temporary directory for intermediate cropped files
temp_dir <- file.path("data/CHELSA/", "temp_cropped")
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

# Helper function to crop and write individual rasters to avoid memory issues
crop_and_stack <- function(files, var_name, time_values) {
  # Build expected temporary file paths
  temp_files <- file.path(
    temp_dir,
    paste0(var_name, "_", time_values, ".tif")
  )

  # Check if all cropped files already exist
  all_exist <- all(file.exists(temp_files))

  if (all_exist) {
    cat("    Reusing existing cropped files from", temp_dir, "\n")
  } else {
    cat("    Cropping rasters and writing to", temp_dir, "\n")

    for (i in seq_along(files)) {
      # Crop individual raster and write to temp file
      r <- rast(files[i])
      r_crop <- crop(r, eu_extent_wgs84, filename = "", overwrite = TRUE)

      # Write to temporary file
      writeRaster(
        r_crop,
        temp_files[i],
        overwrite = TRUE,
        gdal = c("COMPRESS=LZW")
      )

      # Clean up
      rm(r, r_crop)
      gc(verbose = FALSE)
    }
  }

  # Load all cropped files as a stack (now much smaller)
  stack <- rast(temp_files)
  names(stack) <- paste0(var_name, "_", time_values)

  return(stack)
}

cat("  Processing bio01 (temperature)...\n")
bio01_stack <- crop_and_stack(bio01_df$file, "bio01", time_slices)

cat("  Processing bio12 (precipitation)...\n")
bio12_stack <- crop_and_stack(bio12_df$file, "bio12", time_slices)

cat("  Processing glz (ice sheet altitude)...\n")
glz_stack <- crop_and_stack(glz_df$file, "glz", time_slices)

cat("  Processing orog (surface altitude)...\n")
orog_stack <- crop_and_stack(orog_df$file, "orog", time_slices)

cat("  Loaded", n_times, "time slices\n")
cat(
  "  Resolution:",
  paste(res(bio01_stack)[1:2], collapse = " x "),
  "degrees\n"
)
cat(
  "  Dimensions:",
  paste(dim(bio01_stack)[1:2], collapse = " x "),
  "cells\n\n"
)

## Calculate ice-free land masks ####

output_dir <- "data/CHELSA/paleo_derived"

# Define output file paths
years_icefreeland_file <- file.path(
  output_dir,
  "paleo_years_icefreeland_EUextent_EPSG4326.tif"
)
n_consecutive_icefree_file <- file.path(
  output_dir,
  "paleo_icefree_n_consecutive_EUextent_EPSG4326.tif"
)

# Check if both output files already exist
if (
  file.exists(years_icefreeland_file) && file.exists(n_consecutive_icefree_file)
) {
  cat("Loading existing ice-free land calculations from", output_dir, "\n")

  years_icefreeland <- rast(years_icefreeland_file)
  n_consecutive_icefree <- rast(n_consecutive_icefree_file)

  cat("  Loaded years_icefreeland and n_consecutive_icefree\n")
  minmax(years_icefreeland)
  minmax(n_consecutive_icefree)
  plot(years_icefreeland)
  plot(n_consecutive_icefree)
} else {
  cat("Calculating ice-free land masks...\n")

  # For each time slice: ice-free land = glz <= 10 AND orog > 0
  ice_free_land_stack <- ((glz_stack <= ice_threshold) | (is.na(glz_stack))) &
    (orog_stack > sea_threshold)
  names(ice_free_land_stack) <- paste0("icefreeland_", time_slices)

  plot(ice_free_land_stack)

  ## Calculate years since deglaciation ####

  cat("Calculating years since deglaciation...\n")

  # Convert time slices to actual years (hundred-year units to years)
  # Present is time = 20
  years_bp <- (max(time_slices) - time_slices) * 100

  # Create function to calculate years since deglaciation for each cell
  # Working backward from present to past
  # Count consecutive years where cell is ice-free
  calc_years_since_deglac <- function(ice_free_values) {
    # ice_free_values: vector of TRUE/FALSE from oldest to most recent
    # Return: c(years since deglaciation, number of consecutive ice-free slices)

    n <- length(ice_free_values)

    # If currently has missing data, return NA for both
    if (is.na(ice_free_values[n])) {
      return(c(NA_real_, NA_real_))
    }

    # If currently ice-covered or below sea level, return 0 for both
    if (!ice_free_values[n]) {
      return(c(0, 0))
    }

    # Count consecutive ice-free periods from present backward
    years_free <- 0
    n_consecutive <- 0
    for (i in n:1) {
      # Handle NA values as break points (treat as ice-covered)
      if (is.na(ice_free_values[i])) {
        break
      }

      if (ice_free_values[i]) {
        n_consecutive <- n_consecutive + 1
        # Add the interval this observation represents
        if (i > 1) {
          # Interval from previous (older) time to current time
          time_step <- years_bp[i - 1] - years_bp[i]
          years_free <- years_free + time_step
        }
        # Note: for i=1 (oldest observation), we don't add extra time
        # since we've already counted all intervals between observations
      } else {
        break # Stop at first ice-covered period
      }
    }

    return(c(years_free, n_consecutive))
  }

  # Apply function across all cells (returns 2 layers: years and count)
  deglac_results <- app(ice_free_land_stack, calc_years_since_deglac)

  # Split into separate rasters
  years_icefreeland <- deglac_results[[1]]
  names(years_icefreeland) <- "paleo_years_icefreeland"

  n_consecutive_icefree <- deglac_results[[2]]
  names(n_consecutive_icefree) <- "paleo_icefree_n_consecutive"

  cat("  Years since deglaciation calculated\n")
  cat("  Number of consecutive ice-free slices calculated\n")
  minmax(years_icefreeland)
  minmax(n_consecutive_icefree)
  plot(years_icefreeland)
  plot(n_consecutive_icefree)

  # Years since deglaciation
  writeRaster(
    years_icefreeland,
    years_icefreeland_file,
    overwrite = TRUE,
    names = "paleo_years_icefreeland"
  )

  # Number of consecutive ice-free time slices
  writeRaster(
    n_consecutive_icefree,
    n_consecutive_icefree_file,
    overwrite = TRUE,
    names = "paleo_icefree_n_consecutive"
  )
}

## Prepare climate data ####

# Ice-free period (bio01, bio12):
#   - Mean: average conditions during ice-free land period
#   - Sum: cumulative energy/water during ice-free land period
# Entire 16ka period (bio01, bio12):
#   - Mean: average conditions across all time slices
#   - Sum: cumulative energy/water across all time slices
# Total: 8 derived paleoclimate variables (2 vars × 2 periods x 2 stats)

# Function to calculate mean during consecutive ice-free period (matrix version for app)
calc_mean_consecutive_rows <- function(x) {
  # x: matrix where rows=cells, columns=layers
  # Last column is n_consecutive, other columns are climate time slices
  # Return: vector of mean values (one per cell/row)

  n_consecutive <- x[, ncol(x)] # Extract count column
  climate_data <- x[, 1:(ncol(x) - 1)] # All climate columns
  n_times <- ncol(climate_data)

  # Process each row (cell)
  result <- sapply(1:nrow(x), function(i) {
    n <- n_consecutive[i]

    # If n_consecutive is NA, return NA
    if (is.na(n)) {
      return(NA_real_)
    }
    # If n_consecutive is 0, return 0
    if (n == 0) {
      return(0)
    }

    # Extract the most recent n_consecutive slices for this cell
    start_idx <- n_times - n + 1
    consecutive_climate <- climate_data[i, start_idx:n_times]

    return(mean(consecutive_climate, na.rm = TRUE))
  })

  return(result)
}

# Function to calculate sum during consecutive ice-free period (matrix version for app)
calc_sum_consecutive_rows <- function(x) {
  # x: matrix where rows=cells, columns=layers
  # Last column is n_consecutive, other columns are climate time slices
  # Return: vector of sum values (one per cell/row)

  n_consecutive <- x[, ncol(x)] # Extract count column
  climate_data <- x[, 1:(ncol(x) - 1)] # All climate columns
  n_times <- ncol(climate_data)

  # Process each row (cell)
  result <- sapply(1:nrow(x), function(i) {
    n <- n_consecutive[i]

    # If n_consecutive is NA, return NA
    if (is.na(n)) {
      return(NA_real_)
    }
    # If n_consecutive is 0, return 0
    if (n == 0) {
      return(0)
    }

    # Extract the most recent n_consecutive slices for this cell
    start_idx <- n_times - n + 1
    consecutive_climate <- climate_data[i, start_idx:n_times]

    return(sum(consecutive_climate, na.rm = TRUE))
  })

  return(result)
}

# Function to calculate mean across all time slices
calc_mean_all <- function(climate_values) {
  # climate_values: vector of climate values
  # Return: mean across all time slices (no ice-free filter)

  if (all(is.na(climate_values))) {
    return(NA_real_)
  }

  return(mean(climate_values, na.rm = TRUE))
}

# Function to calculate sum across all time slices
calc_sum_all <- function(climate_values) {
  # climate_values: vector of climate values
  # Return: sum across all time slices (no ice-free filter)

  if (all(is.na(climate_values))) {
    return(NA_real_)
  }

  return(sum(climate_values, na.rm = TRUE))
}

## During consecutive ice-free land period ####

cat("Calculating climate during consecutive ice-free periods...\n")

# Mean bio01

# Combine bio01 stack and consecutive count
bio01_consecutive_stack <- c(bio01_stack, n_consecutive_icefree)
mean_bio01_icefree <- app(bio01_consecutive_stack, calc_mean_consecutive_rows)
names(mean_bio01_icefree) <- "paleo_bio01_mean_icefree"

cat("  Mean bio01 during consecutive ice-free period calculated\n")

# Mean bio12

# Same for bio12
bio12_consecutive_stack <- c(bio12_stack, n_consecutive_icefree)
mean_bio12_icefree <- app(bio12_consecutive_stack, calc_mean_consecutive_rows)
names(mean_bio12_icefree) <- "paleo_bio12_mean_icefree"

cat("  Mean bio12 during consecutive ice-free period calculated\n")

# Sum bio01

# Sum of bio01 values during consecutive ice-free period
sum_bio01_icefree <- app(bio01_consecutive_stack, calc_sum_consecutive_rows)
names(sum_bio01_icefree) <- "paleo_bio01_sum_icefree"

cat("  Sum bio01 during consecutive ice-free period calculated\n")

# Sum bio12

# Sum of bio12 values during consecutive ice-free period
sum_bio12_icefree <- app(bio12_consecutive_stack, calc_sum_consecutive_rows)
names(sum_bio12_icefree) <- "paleo_bio12_sum_icefree"

cat("  Sum bio12 during consecutive ice-free period calculated\n\n")

## During entire 16ka period ####

cat("Calculating climate over entire 16ka period...\n")

# Mean of bio01 across all time slices
mean_bio01_16ka <- app(bio01_stack, calc_mean_all)
names(mean_bio01_16ka) <- "paleo_bio01_mean_16ka"

cat("  Mean bio01 over 16ka period calculated\n")

# Mean of bio12 across all time slices
mean_bio12_16ka <- app(bio12_stack, calc_mean_all)
names(mean_bio12_16ka) <- "paleo_bio12_mean_16ka"

cat("  Mean bio12 over 16ka period calculated\n")

# Sum of bio01 across all time slices
sum_bio01_16ka <- app(bio01_stack, calc_sum_all)
names(sum_bio01_16ka) <- "paleo_bio01_sum_16ka"

cat("  Sum bio01 over 16ka period calculated\n")

# Sum of bio12 across all time slices
sum_bio12_16ka <- app(bio12_stack, calc_sum_all)
names(sum_bio12_16ka) <- "paleo_bio12_sum_16ka"

cat("  Sum bio12 over 16ka period calculated\n\n")

## Write outputs ####

cat("Writing output files...\n")

output_dir <- "data/CHELSA/paleo_derived"

# Mean bio01 during ice-free period
writeRaster(
  mean_bio01_icefree,
  file.path(output_dir, "paleo_bio01_mean_icefree_EUextent_EPSG4326.tif"),
  overwrite = TRUE,
  names = "paleo_bio01_mean_icefree"
)

# Mean bio12 during ice-free period
writeRaster(
  mean_bio12_icefree,
  file.path(output_dir, "paleo_bio12_mean_icefree_EUextent_EPSG4326.tif"),
  overwrite = TRUE,
  names = "paleo_bio12_mean_icefree"
)

# Sum bio01 during ice-free period
writeRaster(
  sum_bio01_icefree,
  file.path(output_dir, "paleo_bio01_sum_icefree_EUextent_EPSG4326.tif"),
  overwrite = TRUE,
  names = "paleo_bio01_sum_icefree"
)

# Sum bio12 during ice-free period
writeRaster(
  sum_bio12_icefree,
  file.path(output_dir, "paleo_bio12_sum_icefree_EUextent_EPSG4326.tif"),
  overwrite = TRUE,
  names = "paleo_bio12_sum_icefree"
)

# Mean bio01 over entire 16ka period
writeRaster(
  mean_bio01_16ka,
  file.path(output_dir, "paleo_bio01_mean_16ka_EUextent_EPSG4326.tif"),
  overwrite = TRUE,
  names = "paleo_bio01_mean_16ka"
)

# Mean bio12 over entire 16ka period
writeRaster(
  mean_bio12_16ka,
  file.path(output_dir, "paleo_bio12_mean_16ka_EUextent_EPSG4326.tif"),
  overwrite = TRUE,
  names = "paleo_bio12_mean_16ka"
)

# Sum bio01 over entire 16ka period
writeRaster(
  sum_bio01_16ka,
  file.path(output_dir, "paleo_bio01_sum_16ka_EUextent_EPSG4326.tif"),
  overwrite = TRUE,
  names = "paleo_bio01_sum_16ka"
)

# Sum bio12 over entire 16ka period
writeRaster(
  sum_bio12_16ka,
  file.path(output_dir, "paleo_bio12_sum_16ka_EUextent_EPSG4326.tif"),
  overwrite = TRUE,
  names = "paleo_bio12_sum_16ka"
)

# Clean up temporary directory
cat("\nCleaning up temporary files...\n")
unlink(temp_dir, recursive = TRUE)

cat("\nProcessing complete!\n")
cat("Output directory:", output_dir, "\n\n")

# sessionInfo ####

sessioninfo::session_info()
