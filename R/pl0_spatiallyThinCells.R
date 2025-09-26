library(tidyverse)
library(sf)
library(GeoThinneR)

# Load coordinate data ####
coords_no <- read_csv("output/presence_coords_global_no.csv") |> 
  mutate(source = "NO",
         presence_id = row_number())

coords_eu <- read_csv("output/presence_coords_global_eu.csv") |> 
  mutate(source = "EU", 
         presence_id = row_number())

# Combine datasets
all_coords <- bind_rows(coords_no, coords_eu)

cat("Initial data counts:\n")
cat("- NO coordinates:", nrow(coords_no), "points\n")
cat("- EU coordinates:", nrow(coords_eu), "points\n")
cat("- Total:", nrow(all_coords), "points\n\n")

# Nearest neighbor distances in EU data ####
coords_eu_sf <- coords_eu |> 
  st_as_sf(coords = c("x", "y"), crs = 3035) 

# Calculate nearest neighbor distances within EU dataset
nearest <- st_nearest_feature(coords_eu_sf)
nn_dist <- st_distance(coords_eu_sf, 
                       coords_eu_sf[nearest,], 
                       by_element = TRUE)
quantile(nn_dist, probs = seq(0.5, 1, 0.05))

# Apply thinning to combined dataset ####

# Thinning distance threshold
thinning_distance <- units::set_units(min(nn_dist)*5, "km") 

# Convert to WGS84 for GeoThinneR compatibility
all_coords_sf <- all_coords |> 
  st_as_sf(coords = c("x", "y"), crs = 3035) |> 
  st_transform(crs = 4326)

all_coords_wgs84 <- all_coords_sf |>
  st_drop_geometry() |> 
  bind_cols(st_coordinates(all_coords_sf))

thinned_points <- thin_points(
  data = all_coords_wgs84,
  lon_col = "X",
  lat_col = "Y", 
  method = "distance",
  thin_dist = as.numeric(thinning_distance),
  trials = 10,
  all_trials = FALSE,
  seed = 123
)

summary(thinned_points)

all_coords_thinned <- all_coords[thinned_points$retained[[1]], ]

# Split back into separate datasets using inner_join
coords_no_thinned <- coords_no |> 
  select(source, presence_id) |> 
  inner_join(all_coords_thinned, by = join_by(source, presence_id))

coords_eu_thinned <- coords_eu |> 
  select(source, presence_id) |> 
  inner_join(all_coords_thinned, by = join_by(source, presence_id))

# Write thinned coordinates to CSV files
coords_no_thinned |> 
  select(x, y) |> 
  write_csv("output/presence_coords_global_no_thinned.csv", append = FALSE)
coords_eu_thinned |> 
  select(x, y) |> 
  write_csv("output/presence_coords_global_eu_thinned.csv", append = FALSE)

# Quality control ####
# Verify minimum separation in final dataset
final_coords_sf <- 
  bind_rows(coords_no_thinned, coords_eu_thinned) |> 
  st_as_sf(coords = c("x", "y"), crs = 3035)

nearest_final <- st_nearest_feature(final_coords_sf)
nn_dist_final <- st_distance(final_coords_sf,
                             final_coords_sf[nearest_final,],
                             by_element = TRUE)
min_separation <- min(nn_dist_final)

cat("Quality control:\nminimum separation in thinned dataset:", round(min_separation, 0), "m\n")

# sessionInfo ####

sessioninfo::session_info()
