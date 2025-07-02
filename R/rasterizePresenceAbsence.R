library(tidyverse)
library(sf)
library(terra)

# Load raster data ####
regional_raster <- rast("output/predictors_regional_250m_Norway_current.tif")
global_raster <- rast("output/CHELSA_global_predictors_5km_EU_Norway.tif")

# Load vector data ####
lyngstad_mtype_a <- st_read("data/DMraisedbog.gpkg", layer = "lyngstad-MTYPE_A")
nib_footprint <- st_read("data/DMraisedbog.gpkg", layer = "nib-lyngstad-footprint") 
presences_eu <- st_read("data/DMraisedbog.gpkg", layer = "Presences_EU")

# Table 1: Presence coordinates (lyngstad-MTYPE_A × regional raster) ####
st_crs(lyngstad_mtype_a) == st_crs(regional_raster) # Check CRS compatibility
presence_raster <- rasterize(lyngstad_mtype_a, regional_raster, background = 0)
presence_cells <- which(values(presence_raster) == 1, arr.ind = FALSE)
table1 <- xyFromCell(presence_raster, presence_cells) |> 
  as_tibble()

# Table 2: Absence coordinates (footprint - presence) ####
st_crs(nib_footprint) == st_crs(regional_raster) # Check CRS compatibility
footprint_raster <- rasterize(nib_footprint, regional_raster, background = 0)
absence_raster <- footprint_raster - presence_raster
absence_cells <- which(values(absence_raster) == 1, arr.ind = FALSE)
table2 <- xyFromCell(absence_raster, absence_cells) |> 
  as_tibble()

# Table 3: EU presence coordinates (Presences_EU × global raster) ####
st_crs(presences_eu) == st_crs(global_raster) # Check CRS compatibility
presences_eu <- st_transform(presences_eu, crs(global_raster)) # Ensure CRS match
eu_presence_raster_polygon <- rasterize(svc(presences_eu)[[1]], global_raster, 
                                        background = 0, touches = TRUE)
eu_presence_raster <- rasterize(svc(presences_eu)[[2]], eu_presence_raster_polygon, 
                                update = TRUE, touches = TRUE)
eu_presence_cells <- which(values(eu_presence_raster) == 1, arr.ind = FALSE)
table3 <- xyFromCell(eu_presence_raster, eu_presence_cells) |> 
  as_tibble()

# Write output tables ####
write_csv(table1, "output/presence_coords_regional.csv")
write_csv(table2, "output/absence_coords_regional.csv") 
write_csv(table3, "output/presence_coords_eu_global.csv")

# Summary ####
cat("Table 1 (Regional presences):", nrow(table1), "coordinates\n")
cat("Table 2 (Regional absences):", nrow(table2), "coordinates\n") 
cat("Table 3 (EU presences):", nrow(table3), "coordinates\n")
