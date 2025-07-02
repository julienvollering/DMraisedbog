library(terra)
library(sf)
library(rnaturalearth)
library(tidyverse)

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
            filename = "output/CHELSA_global_predictors_5km_EU_Norway.tif",
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
# Objective: Prepare 29 predictors for nested SDM regional model at 250m resolution
# - 19 CHELSA bioclimatic variables (current + future scenarios)
# - 2 terrain variables (elevation, slope) from DTM50
# - 7 AR50 land cover classes (artype: 10,20,30,50,60,70,81)
# - Output: EPSG:25833, 250m resolution, mainland Norway extent
# 
# Steps:
# 1. Create Norway mainland mask
# 2. Process CHELSA variables (transform to UTM33, resample to 250m)
# 3. Process terrain variables (mosaic DTM tiles, calculate slope)
# 4. Process AR50 land cover (rasterize vector to 250m binary layers)
# 5. Combine into final multi-band outputs (current + future)

## Step 1: Create Norway mainland mask from DTM50 ####

# List all DTM tiles
dtm_files <- list.files("data/DTM50_UTM33_20250613", 
                       pattern = "\\.tif$", 
                       full.names = TRUE)

# Load and mosaic DTM tiles
dtm_tiles <- map(dtm_files, rast)
dtm_mosaic <- do.call(mosaic, dtm_tiles)

# Define target grid parameters
target_crs <- st_crs(25833)
target_res <- 250 # meters

# Aggregate DTM from 50m to 250m to create template grid
dtm_250m <- aggregate(dtm_mosaic, fact = 5, fun = "mean", na.rm = TRUE)
elevation_250m <- dtm_250m
names(elevation_250m) <- "elevation"

# Create 250m template grid for QGIS rasterization
grid_template <- dtm_250m
values(grid_template) <- NA
names(grid_template) <- "template_250m"

# Write template grid to file for QGIS
writeRaster(grid_template, 
            filename = "output/rasterized_ar50_250m_UTM33.tif",
            overwrite = FALSE)

# Skip R rasterization - will be done in QGIS
# After QGIS rasterization, read the result back:
ar50_rasterized_250m <- rast("output/rasterized_ar50_250m_UTM33.tif")  # TODO: Create this file in QGIS
names(ar50_rasterized_250m) <- "ar50_artype"

# Create land mask by reclassifying ar50_rasterized_250m
# Combine all land cells (exclude sea/ocean - artype 82, and NA values)
# Create binary mask: 1 = land, 0 = sea/water, NA = no data
land_mask <- classify(ar50_rasterized_250m, 
                      matrix(c(82, NA,     # sea/ocean -> 0 (water)
                               81, 1,     # ferskvann -> 1 ("land")
                               10, 1,      # bebygd_samferdsel -> 1 (land)
                               20, 1,      # jordbruk -> 1 (land)
                               30, 1,      # skog -> 1 (land)
                               50, 1,      # snaumark -> 1 (land)
                               60, 1,      # myr -> 1 (land)
                               70, 1),     # sno_isbre -> 1 (land)
                             ncol = 2, byrow = TRUE),
                      others = NA)
names(land_mask) <- "land_mask"

## Step 2: Process CHELSA variables ####

# Load CHELSA past data (current climate)
chelsa_past_files <- list.files("data/CHELSA/past", 
                                pattern = "CHELSA_bio.*\\.tif$", 
                                full.names = TRUE) %>%
  sort()

# Reorder to get bio1-bio19 in correct sequence
bio_numbers_past <- str_extract(chelsa_past_files, "bio\\d+") %>%
  str_extract("\\d+") %>%
  as.numeric()
chelsa_past_files <- chelsa_past_files[order(bio_numbers_past)]

# Load CHELSA future data
chelsa_future_files <- list.files("data/CHELSA/future", 
                                  pattern = "CHELSA_bio.*\\.tif$", 
                                  full.names = TRUE) %>%
  sort()

# Reorder future files
bio_numbers_future <- str_extract(chelsa_future_files, "bio\\d+") %>%
  str_extract("\\d+") %>%
  as.numeric()
chelsa_future_files <- chelsa_future_files[order(bio_numbers_future)]

# Process past climate data
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

# Process future climate data
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

## Step 3: Process terrain variables ####

# DTM processing already done in Step 1 for elevation
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

## Step 4: Process AR50 land cover ####

# ar50_rasterized_250m already contains all artype data at correct CRS and resolution
# Simply create binary layers for each artype class using raster operations

# Define artype classes (excluding 99 - Ikke kartlagt)
artype_classes <- c(10, 20, 30, 50, 60, 70, 81)
artype_names <- c("bebygd_samferdsel", "jordbruk", "skog", "snaumark", 
                  "myr", "sno_isbre", "ferskvann")

# Create binary rasters for each artype class
artype_rasters <- map2(artype_classes, artype_names, function(class_code, class_name) {
  
  # Create binary layer: 1 where artype matches class_code, 0 elsewhere, NA where original is NA
  class_raster <- ifel(ar50_rasterized_250m == class_code, 1, 0)
  names(class_raster) <- class_name
  return(class_raster)
})

# Combine into single stack
landcover_stack <- rast(artype_rasters)
landcover_stack <- mask(landcover_stack, land_mask)

## Step 5: Combine all predictors ####

# Combine all predictors for current conditions
predictors_current <- c(chelsa_past_masked, terrain_stack, landcover_stack)

# Combine predictors for future conditions (terrain and landcover remain same)
predictors_future <- c(chelsa_future_masked, terrain_stack, landcover_stack)

# Verify we have 29 predictors
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
