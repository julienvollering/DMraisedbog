library(tidyverse)
library(sf)
library(GeoThinneR)

# Load presence data ####
presences_eu <- st_read("data/DMraisedbog.gpkg", layer = "Presences_EU")
lyngstad_mtype_a <- st_read("data/DMraisedbog.gpkg", layer = "lyngstad-MTYPE_A")

cat("Initial data counts:\n")
cat("- Presences_EU:", nrow(presences_eu), "features\n")
cat("- lyngstad-MTYPE_A:", nrow(lyngstad_mtype_a), "features\n\n")

## Transform all data to common CRS (EPSG:3035 - LAEA Europe) ####
lyngstad_mtype_a_3035 <- lyngstad_mtype_a |> 
  st_transform(crs = st_crs(presences_eu))

## Convert all geometries to centroids for distance calculations ####
presences_eu <- presences_eu |> 
  mutate(source = "EU", 
         presence_id = row_number())

lyngstad <- lyngstad_mtype_a_3035 |> 
  mutate(source = "Lyngstad",
         presence_id = row_number())

## Combine all presences as points ####
all_presences <- bind_rows(
  presences_eu |> st_centroid() |> select(source, presence_id, geom = geom),
  lyngstad_points |> st_centroid() |> select(source, presence_id, geom = geom)
)

# Calculate nearest neighbor distances from EU presences ####

## Remove spatial duplicates ####
duplicates = st_equals(presences_eu_points, retain_unique = TRUE) |> 
  unlist() |> 
  unique()

if(length(duplicates) > 0) {
  presences_eu_points_unique <- presences_eu_points[-duplicates, ]
  cat("Removed", length(duplicates), "spatial duplicates, keeping", nrow(presences_eu_points_unique), "unique locations\n")
} else {
  presences_eu_points_unique <- presences_eu_points
  cat("No spatial duplicates found\n")
}

## Calculate nearest neighbor distances ####
nearest <- st_nearest_feature(presences_eu_points_unique)
nn_dist <- st_distance(presences_eu_points_unique, 
                       presences_eu_points_unique[nearest,], 
                       by_element=TRUE)

# Calculate summary statistics
mean_nn_dist <- mean(nn_dist)
median_nn_dist <- median(nn_dist)

cat("EU nearest neighbor distance statistics (meters):\n")
cat("- Mean:", round(mean_nn_dist, 0), "m\n")
cat("- Median:", round(median_nn_dist, 0), "m\n") 

## Set thinning distance ####
thinning_distance <- mean_nn_dist |> 
  units::set_units("km")

cat("Using thinning distance:", round(thinning_distance, 1), "km\n\n")

# Apply spatial thinning using GeothinneR ####
presences_WGS84_df <- all_presences |> 
  st_drop_geometry() |> 
  bind_cols(st_coordinates(st_transform(all_presences, "epsg:4326")))

thinned_points <- thin_points(
  data = presences_WGS84_df,
  lon_col = "X",
  lat_col = "Y",
  method = "distance",
  thin_dist = as.numeric(thinning_distance), # Distance in km
  trials = 10,
  all_trials = FALSE,
  seed = 123
)

summary(thinned_points)

# Identify which original features to keep based on thinned points ####
# Match thinned centroids back to original geometries
all_presences_thinned <- all_presences[thinned_points$retained[[1]],] |> 
  st_drop_geometry()
presences_eu_thinned <- inner_join(presences_eu, all_presences_thinned,
                                   by = join_by(source, presence_id))
lyngstad_thinned <- inner_join(lyngstad, all_presences_thinned,
                               by = join_by(source, presence_id))

# Write thinned datasets to GeoPackage ####

all_presences_thinned <-
  bind_rows(presences_eu_thinned, lyngstad_thinned)
all_presences_thinned |> 
  st_write("data/DMraisedbog.gpkg", 
           layer = "Presence_Geometries_Thinned", 
           append = FALSE)

# Quality control: Check final nearest neighbor distances ####
final_centroids <- all_presences_thinned |> 
  st_centroid()
nearest <- st_nearest_feature(final_centroids)
nn_dist <- st_distance(final_centroids, 
                       final_centroids[nearest,], 
                       by_element=TRUE)
min(nn_dist)

# sessionInfo ####

sessioninfo::session_info()
