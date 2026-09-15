# Rasterize presences and absences onto the modelling grid ####

# PURPOSE: Converts the Lyngstad bog polygons and the survey footprint into presence and absence cell coordinate tables on the 250 m Norway grid.

library(tidyverse)
library(sf)
library(terra)

# Load raster data ####
regional_raster <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")
crs_common <- st_crs(regional_raster)

# Load vector data ####
lyngstad <- st_read("data/DMraisedbog.gpkg", layer = "lyngstad-MTYPE_A") |>
  st_transform(crs_common)
nib_footprint <- st_read("data/DMraisedbog.gpkg", layer = "nib-lyngstad-footprint") |>
  st_transform(crs_common)

# Table 1: Regional presence coordinates (lyngstad-MTYPE_A × regional raster) ####
presence_raster <- rasterize(
  lyngstad,
  regional_raster,
  background = 0,
  touches = FALSE) # only cells with covered center, sensitivity down specificity up
presence_cells <- which(values(presence_raster) == 1, arr.ind = FALSE)
table1 <- xyFromCell(presence_raster, presence_cells) |>
  as_tibble()

# Table 2: Regional absence coordinates (footprint - presence) ####
footprint_raster <- rasterize(
  nib_footprint,
  regional_raster,
  background = 0,
  touches = FALSE) # only cells with covered center
absence_raster <- footprint_raster - presence_raster # including some cells with partial lyngstad coverage
absence_cells <- which(values(absence_raster) == 1, arr.ind = FALSE)
table2 <- xyFromCell(absence_raster, absence_cells) |>
  as_tibble()

# Write output tables ####
write_csv(table1, "output/presence_coords_regional.csv", append = FALSE)
write_csv(table2, "output/absence_coords_regional.csv", append = FALSE)

# Summary ####
cat("Regional presences:", nrow(table1), "coordinates\n")
cat("Regional absences:", nrow(table2), "coordinates\n")

# sessionInfo ####

sessioninfo::session_info()
