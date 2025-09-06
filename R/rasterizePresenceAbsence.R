library(tidyverse)
library(sf)
library(terra)

# Load raster data ####
regional_raster <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")
global_raster <- rast("output/predictors_global_5km_EUNorway_EPSG3035.tif")

crs_common <- st_crs(global_raster)

# Load vector data ####
lyngstad <- st_read("data/DMraisedbog.gpkg", layer = "lyngstad-MTYPE_A") |>
  st_transform(crs_common)
nib_footprint <- st_read("data/DMraisedbog.gpkg", layer = "nib-lyngstad-footprint") |> 
  st_transform(crs_common)
presences_eu <- st_read("data/DMraisedbog.gpkg", layer = "Presences_EU") |> 
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

# Table 3: Global presence coordinates, Norway (lyngstad-MTYPE_A × global raster) ####
presence_raster <- rasterize(
  lyngstad, 
  global_raster, 
  background = 0, 
  touches = TRUE) # any touched cell, sensitivity up specificity down
presence_cells <- which(values(presence_raster) == 1, arr.ind = FALSE)
table3 <- xyFromCell(presence_raster, presence_cells) |> 
  as_tibble()

# Table 4: Global presence coordinates, EU (Presences_EU × global raster) ####
presence_raster_1 <- rasterize(
  svc(presences_eu)[[1]], 
  global_raster, 
  background = 0, 
  touches = TRUE) # any touched cell, sensitivity up specificity down
presence_raster_2 <- rasterize(
  svc(presences_eu)[[2]], 
  global_raster, 
  background = 0, 
  touches = TRUE) # any touched cell, sensitivity up specificity down
presence_raster <- max(presence_raster_1, presence_raster_2)
presence_cells <- which(values(presence_raster) == 1, arr.ind = FALSE)
table4 <- xyFromCell(presence_raster, presence_cells) |> 
  as_tibble()

# Write output tables ####
write_csv(table1, "output/presence_coords_regional.csv", append = FALSE)
write_csv(table2, "output/absence_coords_regional.csv", append = FALSE) 
write_csv(table3, "output/presence_coords_global_no.csv", append = FALSE)
write_csv(table4, "output/presence_coords_global_eu.csv", append = FALSE)

# Summary ####
cat("Regional presences:", nrow(table1), "coordinates\n")
cat("Regional absences:", nrow(table2), "coordinates\n") 
cat("Global presences, Norway:", nrow(table3), "coordinates\n")
cat("Global presences, EU:", nrow(table4), "coordinates\n")

# sessionInfo ####

sessioninfo::session_info()
