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
# 
# Steps:
# 1. Load all 19 CHELSA bioclimatic variables from data/CHELSA/past/
# 2. Create EU + Norway mask from rnaturalearth data in EPSG:3035
# 3. Create 5000m resolution template grid in EPSG:3035 for European extent
# 4. Project CHELSA data from WGS84 to EPSG:3035 at 5000m resolution
# 5. Mask to EU + Norway extent (remove sea and non-target countries)
# 6. Write final stack as multi-band GeoTIFF

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
# - Output: EPSG:25833, 250m resolution, mainland Norway extent
# 
# Steps:
# 1. Process AR50 land cover (rasterize vector to 50m, create binary layers, aggregate to 250m)
# 2. Create Norway mainland mask (derived from AR50)
# 3. Process CHELSA variables (transform to UTM33, resample to 250m)
# 4. Process terrain variables (mosaic DTM tiles, calculate slope)
# 5. Combine into final multi-band outputs (current + future)

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

# Load AR50 land cover data
ar50 <- st_read("data/0000_25833_ar50_gdb.gdb")
ar50_land <- filter(ar50, artype != 82 & artype != 99) # Exclude sea (82) and not mapped (99)

# Transform AR50 polygons from UTM33 to EPSG:3035
ar50_land_3035 <- st_transform(ar50_land, crs = "EPSG:3035")
rm(ar50, ar50_land) # Free memory
gc()

# Rasterize AR50 directly to 250m EPSG:3035 grid
# Create raster::RasterLayer for fasterize compatibility
raster_template <- raster::raster(template_250m)

# Rasterize classes: 10 - Bebygd og samferdsel; 20 - Jordbruksareal; 30 - Skog; 50 - Snaumark; 60 - Myr; 70 - SnøIsbre; 81 - Ferskvann
ar50_250m <- fasterize::fasterize(
  sf = ar50_land_3035,
  raster = raster_template,
  field = "artype"
)
rm(ar50_land_3035, raster_template)
gc()

# Convert to terra SpatRaster and write to file
ar50_250m_terra <- rast(ar50_250m)
crs(ar50_250m_terra) <- "EPSG:3035"
writeRaster(ar50_250m_terra, "output/ar50_250m_EPSG3035.tif", overwrite = TRUE)
rm(ar50_250m)
gc()

### Create multi-layer from single layer ####
# ar50_250m <- rast("output/ar50_250m_EPSG3035.tif")
# 
# # Create SpatRaster with 0 for all non-NA cells (land mask)
# ar50_land_250m <- ifel(is.na(ar50_250m), NA, 0)
# writeRaster(ar50_land_250m, 
#             filename = "output/ar50_250m_land_EPSG3035.tif",
#             overwrite = TRUE)
# 
# # Create individual SpatRasters for each class in one multi-layer SpatRaster
# class_values <- c(10, 20, 30, 50, 60, 70, 81)
# ar50_250m_stack <- c()
# for (val in class_values) {
#   ar50_250m_stack <- c(ar50_250m_stack, (ar50_250m == val) * 1)
# }
# ar50_250m_stack <- do.call(c, ar50_250m_stack)
# names(ar50_250m_stack) <- paste0("class_", class_values)
# 
# rm(ar50_250m)
# gc()
# 
# writeRaster(ar50_250m_stack, 
#             filename = "output/ar50_250m_layers_EPSG3035.tif",
#             overwrite = TRUE)

## Step 3: Create Norway mainland mask ####

# Load AR50 raster and create land mask: any non-NA value = land
ar50_250m <- rast("output/ar50_250m_EPSG3035.tif")
ar50_land_250m <- ifel(is.na(ar50_250m), NA, 0)

# Use ar50_land_250m as land mask: 0 = land, NA = all else
land_mask <- ar50_land_250m
names(land_mask) <- "land_mask"

## Step 4: Process terrain variables ####

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
elevation <- mask(elevation_250m, land_mask)
slope <- mask(slope_250m, land_mask)

# Create terrain stack
terrain_stack <- c(elevation, slope)
names(terrain_stack) <- c("elevation", "slope")

# Clean up large objects
rm(dtm_mosaic_utm33, dtm_3035, elevation_250m, slope_250m)
gc()

## Step 5: Process CHELSA variables ####

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
chelsa_past_masked <- mask(chelsa_past_3035, land_mask)

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
chelsa_future_masked <- mask(chelsa_future_3035, land_mask)

## Step 6: Combine all predictors ####

# Add proper name to AR50 layer
names(ar50_250m) <- "artype"

# Combine all predictors for current conditions
predictors_current <- c(chelsa_past_masked, terrain_stack, ar50_250m)

# Combine predictors for future conditions (terrain and landcover remain same)
predictors_future <- c(chelsa_future_masked, terrain_stack, ar50_250m)

# Verify we have 22 predictors
stopifnot(nlyr(predictors_current) == 22)
stopifnot(nlyr(predictors_future) == 22)

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
  plot(predictors_current[[c(1, 20, 22)]], 
       main = c("Bio1: Annual Mean Temperature", "Elevation", "AR50 Land Cover"))
}

## Spatial coverage validation ####

# Check that all 22 predictors have identical spatial coverage
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
