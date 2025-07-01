library(tidyverse)
library(sf)
library(rnaturalearth)

# Load background map ####
physical_regions <- ne_download(type = "geography_regions_polys", 
                                category = "physical", scale = 10)
europe_physical <- filter(physical_regions, REGION == "Europe")
europe_countries <- ne_countries(continent = "europe", scale = 10)
europe_countries |> 
  st_geometry() |>
  plot(axes = TRUE, border = 'grey', 
       extent = st_bbox(europe_physical))

# Read in mire regions ####
mire_regions <- st_read("data/Tanneberger/mire_region_general/mire_region_general.shp")
mire_regions |> 
  arrange(desc(region_num)) |> # NB: spatial overlap between Region 5 and 7!  
  select(region) |> 
  st_transform(crs = st_crs(europe_countries)) |>
  plot(add = TRUE)

# Read in Natura 2000 data ####

natura <- st_read("data/Natura2000_end2021_rev1_gpkg/Natura2000_end2021_rev1.gpkg", 
                  layer = "NaturaSite_polygon")
habitats <- st_read("data/Natura2000_end2021_rev1_gpkg/Natura2000_end2021_rev1.gpkg", 
                  layer = "HABITATS")
raisedbogs <- habitats |> 
  filter(HABITATCODE %in% c("7110", "7120"))
# "Sites judged to be still capable of natural regeneration will include those areas where the hydrology can be repaired and where, with appropriate rehabilitation management, there is a reasonable expectation of re-establishing vegetation with peat-forming capability within 30 years." 

# Explore Natura 2000 data ####

raisedbogs |>
  group_by(DATAQUALITY) |> 
  count()

raisedbogs |> 
  group_by(COVER_HA > 1) |> # how much of the site is covered by the habitat?
  count()

naturaRB <- semi_join(natura, raisedbogs, by = "SITECODE") |> 
  left_join(raisedbogs, by = "SITECODE")

naturaRB |> 
  st_write("data/DMraisedbog.gpkg", 
            layer = "Natura2000_raisedbogs", append = FALSE)

naturaRB |> 
  mutate(area = st_area(geom)) |> 
  st_drop_geometry() |> 
  group_by(area < units::set_units(1, "km^2")) |>
  count() # Polygons finer than the modeling resolution may be included without intersecting EUNIS

# Read in EUNIS data ####
eunis <- st_read("data/EUNIS/Wetland_Plot.gpkg", layer = "Wetland")
eunisQ11 <- eunis |>
  filter(HabitatCod == "Q11")

# Filtered presences by intersection ####

naturaRBG <- naturaRB |> 
  filter(DATAQUALITY == "G")
regionIV <- mire_regions |> 
  filter(region == "IV") |> 
  st_transform(crs = st_crs(naturaRBG))
naturaRBGIV <- st_intersection(naturaRBG, regionIV) |> 
  select(SITECODE, HABITATCODE, NON_PRESENCE_IN_SITE, COVER_HA, geom)

## Small polygons ####

presences1 <- naturaRBGIV |> 
  mutate(area = st_area(geom)) |>
  filter(area < units::set_units(1, "km^2")) # Polygons finer than the modeling resolution may be included without intersecting EUNIS

europe_countries |> 
  st_geometry() |>
  plot(axes = TRUE, border = 'grey', 
       extent = st_bbox(europe_physical))
presences1 |> 
  st_geometry() |>
  st_transform(crs = st_crs(europe_countries)) |>
  plot(add = TRUE, col = 'blue', pch = 20, cex = 1)
plot(st_geometry(presences1), add=TRUE)

## Points inside large polygons ####

largepolygons <- naturaRBGIV |> 
  mutate(area = st_area(geom)) |>
  filter(area >= units::set_units(1, "km^2")) # Polygons finer than the modeling resolution may be included without intersecting EUNIS

presences2 <- eunisQ11[largepolygons, ] |> 
  select(PLOTOBSID, HabitatCod, Shape) |> 
  st_set_geometry("geom")

presences2 |> 
  st_geometry() |>
  st_transform(crs = st_crs(europe_countries)) |>
  plot(add = TRUE, col = 'red', pch = 20, cex = 1)
plot(st_geometry(presences1), add=TRUE)

# Combine presences ####

bind_rows(presences1, presences2) |> 
  st_write("data/DMraisedbog.gpkg", 
            layer = "Presences_EU", append = FALSE)
  
