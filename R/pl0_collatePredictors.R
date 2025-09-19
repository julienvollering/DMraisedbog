library(terra)
library(sf)
library(rnaturalearth)
library(tidyverse)

# Function to get memory usage of all objects
memory_usage <- function() {
  obj_names <- ls(envir = .GlobalEnv)
  obj_sizes <- sapply(obj_names, function(x) {
    object.size(get(x, envir = .GlobalEnv))
  })
  # Convert to MB and sort
  obj_sizes_mb <- round(obj_sizes / 1024^2, 2)
  sort(obj_sizes_mb, decreasing = TRUE)
}

# Global model ####

# Plan for Predictor Stack Preparation for Global (EU) Model
# 
# Objective: Prepare CHELSA bioclimatic variables for nested SDM global model
# - Input: 19 CHELSA variables at 30 arcsec (~1km) resolution in WGS84
# - Output: Single GeoTIFF with 19 bands at 5000m resolution in EPSG:3035, masked to EU+Norway

## Step 1: Load CHELSA data ####
chelsa_files <- list.files("data/CHELSA/past",
                          pattern = "CHELSA_bio.*\\.tif$",
                          full.names = TRUE) %>%
  sort() # Ensure bio1, bio10, bio11, ..., bio19, bio2, bio3, etc.

# Reorder to get bio1-bio19 in correct sequence
bio_numbers <- str_extract(chelsa_files, "bio\\d+") %>%
  str_extract("\\d+") %>%
  as.numeric()
chelsa_files <- chelsa_files[order(bio_numbers)]

# Load as SpatRaster stack
chelsa_stack <- rast(chelsa_files)
names(chelsa_stack) <- paste0("bio", 1:19)
chelsa_stack

## Step 2: Create EU + Norway mask ####

# Get European countries
europe_countries <- ne_countries(continent = "europe", scale = 10, returnclass = "sf")

# Get EU member states (approximate list - may need updating)
eu_countries <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia",
                  "Denmark", "Estonia", "Finland", "France", "Germany", "Greece",
                  "Hungary", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg",
                  "Malta", "Netherlands", "Poland", "Portugal", "Romania", "Slovakia",
                  "Slovenia", "Spain", "Sweden")

# Filter to EU countries + Norway
target_countries <- europe_countries %>%
  filter(name %in% c(eu_countries, "Norway"))

# Transform to EPSG:3035 (European Equal Area)
target_countries <- st_transform(target_countries, crs = "EPSG:3035")

# Create combined polygon for masking
mask_polygon <- target_countries %>%
  st_union() %>%
  st_sf()

## Step 3: Create template grid and project to EPSG:3035 ####

# Define European extent in EPSG:3035 coordinates
# Approximate bounding box covering EU + Norway in EPSG:3035
europe_ext_3035 <- ext(c(xmin = 2000000, xmax = 7000000,
                         ymin = 1000000, ymax = 5500000))

# Create 5000m resolution template grid in EPSG:3035
template_5km <- rast(europe_ext_3035, resolution = 5000, crs = "EPSG:3035")

# Define European extent in WGS84 for cropping before reprojection
# Approximate bounding box covering EU + Norway in WGS84
europe_ext_wgs84 <- ext(c(xmin = -10, xmax = 35,
                         ymin = 35, ymax = 72))

# Crop to European extent first (in WGS84) to reduce memory requirements
chelsa_cropped_wgs84 <- crop(chelsa_stack, europe_ext_wgs84)

# Project CHELSA data from WGS84 to EPSG:3035 at 5000m resolution
chelsa_3035 <- project(x = chelsa_cropped_wgs84,
                       y = template_5km,
                       method = "bilinear")

## Step 4: Mask to EU + Norway ####

# Mask to country boundaries (removes sea and non-target countries)
chelsa_masked <- mask(chelsa_3035, mask_polygon)

## Step 5: Write output ####
writeRaster(chelsa_masked,
            filename = "output/predictors_global_5km_EUNorway_EPSG3035.tif",
            overwrite = TRUE,
            names = names(chelsa_masked))

# Quick visualization check
if(interactive()) {
  plot(chelsa_masked[[1]], main = "Bio1: Annual Mean Temperature")
  plot(st_geometry(target_countries), add = TRUE, border = "black", lwd = 0.5)
}

# Regional model ####

# Plan for Predictor Stack Preparation for Regional (Norway) Model
# 
# Objective: Prepare 28 predictors for nested SDM regional model at 250m resolution
# - 19 CHELSA bioclimatic variables (current + future scenarios)
# - 2 terrain variables (elevation, slope) from DTM50
# - 7 AR50 land cover classes (artype: 10,20,30,50,60,70,81)
# - Output: EPSG:3035, 250m resolution, mainland Norway extent

## Step 1: Create template grid in EPSG:3035 ####

# List all DTM tiles to define Norway extent
dtm_files <- list.files("data/DTM50_UTM33_20250613", 
                       pattern = "\\.tif$", 
                       full.names = TRUE)

# Load and mosaic DTM tiles in original UTM33 CRS
dtm_tiles <- map(dtm_files, rast)
dtm_mosaic_utm33 <- do.call(mosaic, dtm_tiles)

# Get Norway extent in UTM33, then transform to EPSG:3035
norway_ext_utm33 <- ext(dtm_mosaic_utm33)
norway_ext_3035 <- project(norway_ext_utm33, 
                           from = crs(dtm_mosaic_utm33), 
                           to = "EPSG:3035")

# Create 250m template using same origin as global model (template_5km)
# Get origin from global template (5km grid)
global_origin_x <- xmin(template_5km)
global_origin_y <- ymin(template_5km)

# Calculate grid-aligned extent using the same origin
# Align Norway extent to 250m grid with same origin as 5km grid
xmin_aligned <- global_origin_x + floor((norway_ext_3035$xmin - global_origin_x) / 250) * 250
xmax_aligned <- global_origin_x + ceiling((norway_ext_3035$xmax - global_origin_x) / 250) * 250
ymin_aligned <- global_origin_y + floor((norway_ext_3035$ymin - global_origin_y) / 250) * 250
ymax_aligned <- global_origin_y + ceiling((norway_ext_3035$ymax - global_origin_y) / 250) * 250

norway_ext_aligned <- ext(c(xmin = xmin_aligned, xmax = xmax_aligned,
                           ymin = ymin_aligned, ymax = ymax_aligned))

# Create 250m resolution template grid with aligned extent and shared origin
template_250m <- rast(norway_ext_aligned, resolution = 250, crs = "EPSG:3035")

## Step 2: Process AR50 land cover ####

# Load individual AR50 artype rasters created in QGIS (50m resolution, 0/1 values)
ar50_files <- list.files("output/ar50_artype_layers", 
                         pattern = "ar50_50m_EPSG3035_artype.*\\.tif$", 
                         full.names = TRUE) |>
  sort() # Ensure consistent ordering

# Extract artype numbers for naming
artype_numbers <- str_extract(basename(ar50_files), "artype\\d+") |>
  str_extract("\\d+") |>
  as.numeric()

cat("Processing", length(ar50_files), "AR50 artype layers:", paste(artype_numbers, collapse = ", "), "\n")

# Process each artype layer: aggregate from 50m to 250m using mean for continuous [0,1] values
ar50_layers <- list()
for(i in seq_along(ar50_files)) {
  # Load 50m raster
  ar50_50m <- rast(ar50_files[i])

  # Aggregate to 250m using mean (5x5 cells -> continuous [0,1] values)
  ar50_250m <- aggregate(ar50_50m, fact = 5, fun = "mean", na.rm = FALSE)
  
  # Resample to ensure exact alignment with template_250m
  ar50_250m_aligned <- resample(ar50_250m, template_250m, method = "bilinear")
  
  # Store with artype name
  ar50_layers[[paste0("artype_", artype_numbers[i])]] <- ar50_250m_aligned
  
  cat("Processed artype", artype_numbers[i], "\n")
}

# Combine all artype layers into final stack
ar50_250m_stack <- rast(ar50_layers)
names(ar50_250m_stack) <- names(ar50_layers)

### Create Norway mainland land mask from AR50 ####
# Sum all layers to identify areas with any AR50 coverage
ar50_sum <- sum(ar50_250m_stack, na.rm = TRUE)
land_mask_ar50 <- ifel(ar50_sum > 0, 0, NA)  # 0 = land, NA = non-land 
names(land_mask_ar50) <- "land_mask"

# Write land mask to disk
writeRaster(land_mask_ar50, 
            filename = "output/ar50_250m_land_EPSG3035.tif",
            overwrite = TRUE)

# Apply land mask to AR50 stack: keep 0's on land, set non-land to NA
ar50_250m_masked <- mask(ar50_250m_stack, land_mask_ar50)

# Update the stack to the masked version
ar50_250m_stack <- ar50_250m_masked

# Write multiband AR50 output
writeRaster(ar50_250m_stack, 
            filename = "output/ar50_250m_cover_EPSG3035.tif", 
            overwrite = TRUE,
            names = names(ar50_250m_stack))

cat("Created AR50 stack with", nlyr(ar50_250m_stack), "layers:", paste(names(ar50_250m_stack), collapse = ", "), "\n")

# Clean up
rm(ar50_layers)
gc()

## Step 3: Process terrain variables ####

# Project DTM from UTM33 to EPSG:3035 and resample to 250m
dtm_3035 <- project(x = dtm_mosaic_utm33,
                    y = template_250m,
                    method = "bilinear")

# Calculate elevation (already at 250m after projection)
elevation_250m <- dtm_3035
names(elevation_250m) <- "elevation"

# Calculate slope at 250m resolution in EPSG:3035
slope_250m <- terrain(dtm_3035, v = "slope", unit = "degrees")
names(slope_250m) <- "slope"

# Apply mask to both elevation and slope
elevation <- mask(elevation_250m, land_mask_ar50)
slope <- mask(slope_250m, land_mask_ar50)

# Create terrain stack
terrain_stack <- c(elevation, slope)
names(terrain_stack) <- c("elevation", "slope")

# Clean up large objects
rm(dtm_mosaic_utm33, dtm_3035, elevation_250m, slope_250m)
gc()

## Step 4: Process CHELSA variables ####

### Current climate ####
chelsa_past_files <- list.files("data/CHELSA/past", 
                                pattern = "CHELSA_bio.*\\.tif$", 
                                full.names = TRUE) %>%
  sort()

# Reorder to get bio1-bio19 in correct sequence
bio_numbers_past <- str_extract(chelsa_past_files, "bio\\d+") %>%
  str_extract("\\d+") %>%
  as.numeric()
chelsa_past_files <- chelsa_past_files[order(bio_numbers_past)]

chelsa_past_stack <- rast(chelsa_past_files)
names(chelsa_past_stack) <- paste0("bio", 1:19)

# Create Norway extent in WGS84 for cropping before reprojection
# Use raster extent instead of vector union for efficiency
norway_ext_wgs84 <- project(norway_ext_3035, 
                            from = "EPSG:3035",
                            to = "EPSG:4326") # Transform to WGS84

# Crop to Norway extent first (in WGS84) to reduce memory requirements
chelsa_past_cropped_wgs84 <- crop(chelsa_past_stack, norway_ext_wgs84)

# Then transform to EPSG:3035 and resample to 250m
chelsa_past_3035 <- project(x = chelsa_past_cropped_wgs84, 
                            y = template_250m,
                            method = "bilinear")

# Final mask to Norway boundaries using raster mask
chelsa_past_masked <- mask(chelsa_past_3035, land_mask_ar50)

### Future climate ####
chelsa_future_files <- list.files("data/CHELSA/future", 
                                  pattern = "CHELSA_bio.*\\.tif$", 
                                  full.names = TRUE) %>%
  sort()

# Reorder future files
bio_numbers_future <- str_extract(chelsa_future_files, "bio\\d+") %>%
  str_extract("\\d+") %>%
  as.numeric()
chelsa_future_files <- chelsa_future_files[order(bio_numbers_future)]

chelsa_future_stack <- rast(chelsa_future_files)
names(chelsa_future_stack) <- paste0("bio", 1:19)

# Crop to Norway extent first (in WGS84) to reduce memory requirements
chelsa_future_cropped_wgs84 <- crop(chelsa_future_stack, norway_ext_wgs84)

# Then transform to EPSG:3035 and resample to 250m
chelsa_future_3035 <- project(x = chelsa_future_cropped_wgs84, 
                              y = template_250m,
                              method = "bilinear")

# Final mask to Norway boundaries using raster mask
chelsa_future_masked <- mask(chelsa_future_3035, land_mask_ar50)

## Step 5: Combine all predictors ####

# Combine all predictors for current conditions (19 CHELSA + 2 terrain + multiple AR50 layers)
predictors_current <- c(chelsa_past_masked, terrain_stack, ar50_250m_stack)

# Combine predictors for future conditions (terrain and landcover remain same)
predictors_future <- c(chelsa_future_masked, terrain_stack, ar50_250m_stack)

# Verify we have expected number of predictors (19 + 2 + number of AR50 classes)
n_ar50_classes <- nlyr(ar50_250m_stack)
expected_layers <- 19 + 2 + n_ar50_classes
cat("Expected layers:", expected_layers, "(19 CHELSA + 2 terrain +", n_ar50_classes, "AR50 classes)\n")
stopifnot(nlyr(predictors_current) == expected_layers)
stopifnot(nlyr(predictors_future) == expected_layers)

# Write outputs
writeRaster(predictors_current, 
            filename = "output/predictors_regional_250m_Norway_current_EPSG3035.tif",
            overwrite = TRUE,
            names = names(predictors_current))

writeRaster(predictors_future, 
            filename = "output/predictors_regional_250m_Norway_future_EPSG3035.tif",
            overwrite = TRUE,
            names = names(predictors_future))

# Quick visualization check
if(interactive()) {
  # Plot bio1, elevation, and first AR50 layer
  last_layer <- nlyr(predictors_current)
  plot(predictors_current[[c(1, 20, last_layer)]], 
       main = c("Bio1: Annual Mean Temperature", "Elevation", paste("AR50:", names(predictors_current)[last_layer])))
}

## Spatial coverage validation ####

# Check that all predictors have identical spatial coverage
check_coverage <- function(raster_stack) {
  # Count non-NA cells for each layer using terra's efficient global() function
  cell_counts <- global(raster_stack, "notNA")
  # Check if all layers have the same number of non-NA cells
  all(cell_counts$notNA == cell_counts$notNA[1])
}

# Validate both stacks
stopifnot("Current predictors have inconsistent spatial coverage" = check_coverage(predictors_current))
stopifnot("Future predictors have inconsistent spatial coverage" = check_coverage(predictors_future))

# sessionInfo ####

sessioninfo::session_info()
