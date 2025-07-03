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
# - Output: Single GeoTIFF with 19 bands at ~5km resolution, masked to EU+Norway
# 
# Steps:
# 1. Load all 19 CHELSA bioclimatic variables from data/CHELSA/past/
# 2. Create EU + Norway mask from rnaturalearth data
# 3. Aggregate from 30 arcsec to approximately 5km (factor ~5-6) without CRS transformation
# 4. Mask to EU + Norway extent (remove sea and non-target countries)
# 5. Write final stack as multi-band GeoTIFF

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

# Transform to same CRS as CHELSA data (WGS84)
target_countries <- st_transform(target_countries, crs = crs(chelsa_stack))

# Create combined polygon for masking
mask_polygon <- target_countries %>%
  st_union() %>%
  st_sf()

## Step 3: Aggregate to ~5km resolution ####

# Calculate aggregation factor
# 30 arcsec ≈ 0.008333 degrees, 5km ≈ 0.045 degrees at 60°N latitude
# Aggregation factor of 5-6 should work
agg_factor <- 5

# Aggregate using mean (appropriate for climate variables)
chelsa_agg <- aggregate(chelsa_stack, fact = agg_factor, fun = "mean", na.rm = TRUE)

## Step 4: Mask to EU + Norway ####

# Crop to approximate European extent first (for efficiency)
europe_ext <- ext(c(xmin = -10, xmax = 32, 
                    ymin = 36, ymax = 72)) # Approximate bounding box for Europe
chelsa_cropped <- crop(chelsa_agg, europe_ext)

# Mask to country boundaries (removes sea and non-target countries)
chelsa_masked <- mask(chelsa_cropped, mask_polygon)

## Step 5: Write output ####
writeRaster(chelsa_masked, 
            filename = "output/predictors_global_5km_EUNorway.tif",
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

## Step 1: Process AR50 land cover ####

### Rasterize to 50m ####
# Load AR50 land cover data
ar50 <- st_read("data/0000_25833_ar50_gdb.gdb")
ar50_land <- filter(ar50, artype != 82 & artype != 99) # Exclude sea (82) and not mapped (99)
rm(ar50) # Free memory

# Create DTM-based template grid for rasterization
# List all DTM tiles
dtm_files <- list.files("data/DTM50_UTM33_20250613", 
                       pattern = "\\.tif$", 
                       full.names = TRUE)

# Load and mosaic DTM tiles
dtm_tiles <- map(dtm_files, rast)
dtm_mosaic <- do.call(mosaic, dtm_tiles)

# Create 50m template grid for AR50 fasterization
grid_template_50m <- raster::raster(dtm_mosaic)

# Rasterize classes: 10 - Bebygd og samferdsel; 20 - Jordbruksareal; 30 - Skog; 50 - Snaumark; 60 - Myr; 70 - SnøIsbre; 81 - Ferskvann
ar50_50m <- fasterize::fasterize(
  sf = ar50_land,
  raster = grid_template_50m,
  field = "artype"
) # 120 sec
rm(ar50_land)
gc()

# Write the RasterLayer to GeoTIFF first
raster::writeRaster(ar50_50m, "output/ar50_50m.tif", overwrite = TRUE)
rm(ar50_50m)
gc()

### Multi-layer from single layer ####
ar50_50m <- rast("output/ar50_50m.tif")
crs(ar50_50m) <- "EPSG:25833" # Set CRS to UTM33, to match grid_template_50m/dtm_mosaic

# Create SpatRaster with 0 for all non-NA cells
ar50_land_50m <- ifel(is.na(ar50_50m), NA, 0)
writeRaster(ar50_land_50m, 
            filename = "output/ar50_50m_land.tif",
            overwrite = TRUE)

# Create individual SpatRasters for each class in one multi-layer SpatRaster
class_values <- c(10, 20, 30, 50, 60, 70, 81)
ar50_50m_stack <- c()
for (val in class_values) {
  ar50_50m_stack <- c(ar50_50m_stack, (ar50_50m == val) * 1)
}
ar50_50m_stack <- do.call(c, ar50_50m_stack)
names(ar50_50m_stack) <- paste0("class_", class_values)
plot(ar50_50m_stack)

rm(ar50_50m)
gc()

writeRaster(ar50_50m_stack, 
            filename = "output/ar50_50m_layers.tif",
            overwrite = TRUE)

### Aggregate to 250m ####
# Use mean aggregation to get proportional coverage within each 250m cell
ar50_250m_stack <- aggregate(ar50_50m_stack, fact = 5, fun = "mean", na.rm = TRUE)
plot(ar50_250m_stack)

rm(ar50_50m_stack)
gc()

## Step 2: Create Norway mainland mask ####

# Create land mask by aggregating ar50_land_50m: 0 = land, NA = all else
land_mask <- aggregate(ar50_land_50m, fact = 5, fun = "max", na.rm = TRUE)
names(land_mask) <- "land_mask"

## Step 3: Process CHELSA variables ####

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
norway_bbox_utm33 <- ext(land_mask)
norway_ext_wgs84 <- norway_bbox_utm33 |> 
  project(from = crs(land_mask),
          to="EPSG:4326") # Transform to WGS84

# Crop to Norway extent first (in WGS84) to reduce memory requirements
chelsa_past_cropped_wgs84 <- crop(chelsa_past_stack, norway_ext_wgs84)

# Then transform to UTM33 and resample to 250m
chelsa_past_utm33 <- project(x = chelsa_past_cropped_wgs84, 
                             y = land_mask,
                             method = "bilinear")

# Final mask to Norway boundaries using raster mask
chelsa_past_masked <- mask(chelsa_past_utm33, land_mask)

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

# Then transform to UTM33 and resample to 250m
chelsa_future_utm33 <- project(x = chelsa_future_cropped_wgs84, 
                               y = land_mask,
                               method = "bilinear")

# Final mask to Norway boundaries using raster mask
chelsa_future_masked <- mask(chelsa_future_utm33, land_mask)

## Step 4: Process terrain variables ####

# Calculate elevation as mean of 50m resolution
elevation_250m <- aggregate(dtm_mosaic, fact = 5, fun = "mean", na.rm = TRUE)
names(elevation_250m) <- "elevation"

# Calculate slope at 50m resolution first, then aggregate
slope_50m <- terrain(dtm_mosaic, v = "slope", unit = "degrees")

# Aggregate slope from 50m to 250m using mean
slope_250m <- aggregate(slope_50m, fact = 5, fun = "mean", na.rm = TRUE)

# Apply mask to both elevation and slope
elevation <- mask(elevation_250m, land_mask)
slope <- mask(slope_250m, land_mask)

# Create terrain stack
terrain_stack <- c(elevation, slope)
names(terrain_stack) <- c("elevation", "slope")

## Step 5: Combine all predictors ####

# Apply land mask to landcover stack
landcover_stack <- mask(ar50_250m_stack, land_mask)

# Combine all predictors for current conditions
predictors_current <- c(chelsa_past_masked, terrain_stack, landcover_stack)

# Combine predictors for future conditions (terrain and landcover remain same)
predictors_future <- c(chelsa_future_masked, terrain_stack, landcover_stack)

# Verify we have 28 predictors
stopifnot(nlyr(predictors_current) == 28)
stopifnot(nlyr(predictors_future) == 28)

# Write outputs
writeRaster(predictors_current, 
            filename = "output/predictors_regional_250m_Norway_current.tif",
            overwrite = TRUE,
            names = names(predictors_current))

writeRaster(predictors_future, 
            filename = "output/predictors_regional_250m_Norway_future.tif",
            overwrite = TRUE,
            names = names(predictors_future))

# Quick visualization check
if(interactive()) {
  plot(predictors_current[[c(1, 20, 22)]])
}

## Spatial coverage validation ####

# Check that all 28 predictors have identical spatial coverage
check_coverage <- function(raster_stack) {
  # Count non-NA cells for each layer using terra's efficient global() function
  cell_counts <- global(raster_stack, "notNA")
  # Check if all layers have the same number of non-NA cells
  all(cell_counts$notNA == cell_counts$notNA[1])
}

# Validate both stacks
stopifnot("Current predictors have inconsistent spatial coverage" = check_coverage(predictors_current))
stopifnot("Future predictors have inconsistent spatial coverage" = check_coverage(predictors_future))
