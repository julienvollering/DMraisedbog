library(tidyverse)
library(foreign)
library(sf)
library(rnaturalearth)

# Read in raw data #### 
shpfiles <- list.files("data/Lyngstad/rawdata20240212", pattern = '.shp', full.names = TRUE)
polygons <- shpfiles |> 
  map(st_read)

polygons |> 
  map(st_crs)
shpfilecrs <- tribble(~file, ~crs,
                      "1_Hedmark 2012_terratec_20120613f.shp", 25832,
                      "2_Hedmark 2012_blom_13062012f.shp", 25832,
                      "3_Akershus Østfold 2014 Vestfold 2016_høgmyrvestfoldmmsosi20170110f.shp", 25832,
                      "4_Hedmark 2015_hoegmyr_nordlige_oestlandet_20170110f.shp", 25832,         
                      "5_Hedmark 2015_detaljert2_20150416f.shp", 25832,                                
                      "6_Sørlandet 2016_høgmyrsørlandet20170110f.shp", 25832,                          
                      "7_Telemark Buskerud 2016_hoegmyrtelemarkmm20170217f.shp", 25832,                
                      "8_Nordland_200312f.shp", 25833,                                                 
                      "9_Namdalen og fosen_230201f.shp", 25832)

polygons <- map2(polygons, shpfilecrs$crs, \(x,y) st_set_crs(x, y))
polygons <- polygons |> 
  map(\(x) st_transform(x, crs=25833)) |> 
  bind_rows(.id = 'fileindex')

norway <- ne_countries(scale = 10, country = "Norway", returnclass = 'sf')
norway <- norway |> 
  st_crop(xmin=4, xmax = 32, ymin = 57, ymax = 72) |> 
  st_transform(crs = 25833)

norway |> 
  st_geometry() |> 
  plot(border = 'grey', axes = TRUE)
plot(st_geometry(st_centroid(polygons)), pch = 3, col = 'red', add = TRUE)

# Check if different spatial formats contain the same attributes ####
dbffiles <- list.files("data/Lyngstad/rawdata20240212", pattern = '.dbf', full.names = TRUE)
dbfs <- dbffiles |> 
  map(read.dbf) |> 
  bind_rows(.id = "fileindex")
colnames(validpolygons)
colnames(dbfs)
all(colnames(dbfs) %in% colnames(validpolygons))

# dbf <- read.dbf('data/Lyngstad/rawdata20240212/1_Hedmark 2012_terratec_20120613f.dbf')
# Hedmark20120613 <- st_as_sf(dbf, coords =c('X', 'Y'), crs = 25832)
# plot(Hedmark20120613)

# Check for spatial overlap ####
all(st_geometry_type(polygons) == "POLYGON")
all(st_is_valid(polygons))
validpolygons <- st_make_valid(polygons)
overlap <- st_intersection(validpolygons)
plot(overlap$n.overlaps)

# Check for duplicates ####
nrow(validpolygons)
distinct(validpolygons, fileindex, POLY_ID, MYRID) |> 
  nrow()
distinct(validpolygons, fileindex, POLY_ID) |> 
  nrow()

# Filter raised bogs ####
rb <- validpolygons |> 
  filter(grepl("A", M_TYPE))

norway |> 
  st_geometry() |> 
  plot(border = 'grey', axes = TRUE)
rb |> 
  group_by(MYRID) |> 
  summarize() |> 
  st_centroid() |> 
  plot(pch = 3, col = 'red', add = TRUE)

nrow(rb)
rb |> 
  distinct(POLY_ID, MYRID) |> 
  nrow()

rb <- rb |> 
  st_geometry() |> 
  st_union() |> 
  st_cast("POLYGON") |> 
  st_as_sf()
nrow(rb)

rb |>
  st_area() |> 
  units::set_units("ha") |> 
  hist(xlab = "Area", main = "Area of raised bogs (contiguous polygons)")
rb |> 
  st_area() |>
  units::set_units("ha") |> 
  median() |> 
  units::set_units("m^2") |> 
  sqrt() # median side length of raised bogs (assuming square shape)

# Write to file ####

validpolygons |> 
  select(-AREA, -PERIMETER, -X, -Y, -KOORDH) |> 
  st_write("data/DMraisedbog.gpkg", 
           layer = "lyngstad-polygons",
           append = FALSE)
rb |> 
  st_write("data/DMraisedbog.gpkg", 
           layer = "lyngstad-MTYPE_A",
           append = TRUE)

# sessionInfo ####

sessioninfo::session_info()
