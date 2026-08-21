# Build the EU training domain ####

# PURPOSE: Delimits the EU absence domain (Tanneberger mire region IV, intersected with EU-27, EPM coverage and the global stack, minus Norway).

# Plan section 3 ("Domain") and section 5. The EU block is drawn from
#
#   Tanneberger mire region IV  n  EU-27  n  EPM2025 country coverage
#                               n  global 5 km stack extent  -  Norway
#
# The EU-27 clip is not optional: Natura 2000 exists only in member states, so an
# unclipped Region IV would treat Russia, Belarus and Norway as absence domain where
# raised bogs exist but are unmapped -- systematic false absences concentrated in the
# continental east. Norway is dropped separately because it is the Norway block; under
# the EU-27 clip that subtraction should be a no-op, and the script asserts it is.
#
# Runs BEFORE pl0_buildEPMmask.R even though the plan lists the mask first: the domain
# is what bounds the EPM warp, and warping only the domain countries avoids reading
# ~10 Gcell of 1-arcsec source over countries that can never contribute a row. EPM
# *coverage* is a country-level fact taken from the raster filenames here, so there is
# no circular dependency on the warped rasters themselves.

library(sf)
library(terra)
library(dplyr)
library(readr)
library(ggplot2)

dir.create("output/pl0", showWarnings = FALSE, recursive = TRUE)

crs_3035 <- "EPSG:3035"

## Reference grid ####

grid_5km <- rast("output/predictors_global_5km_EUNorway_EPSG3035.tif")
stopifnot(crs(grid_5km, describe = TRUE)$code == "3035")

stack_extent <- st_as_sfc(st_bbox(
  c(
    xmin = xmin(grid_5km),
    xmax = xmax(grid_5km),
    ymin = ymin(grid_5km),
    ymax = ymax(grid_5km)
  ),
  crs = st_crs(3035)
))

## Mire region IV ####

# The shapefile has a .prj but no EPSG code, so st_crs()$epsg is NA. Assert the
# projection we expect (World Mollweide) before transforming rather than trusting it.
mire_regions <- st_read(
  "data/Tanneberger/mire_region_general/mire_region_general.shp",
  quiet = TRUE
)
stopifnot(
  is.na(st_crs(mire_regions)$epsg),
  grepl("proj=moll", st_crs(mire_regions)$proj4string, fixed = TRUE)
)

region_iv <- mire_regions |>
  filter(region == "IV") |>
  st_transform(crs_3035) |>
  st_make_valid() |>
  st_union()

## Countries ####

countries <- st_read(
  "data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp",
  quiet = TRUE
) |>
  st_transform(crs_3035) |>
  st_make_valid()

# ADM0_A3, never ISO_A2/ISO_A3: those are -99 for France and Norway, and France holds
# more raised-bog polygons (3,207) than any other country.
eu27_iso3 <- c(
  "AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA",
  "DEU", "GRC", "HUN", "IRL", "ITA", "LVA", "LTU", "LUX", "MLT", "NLD",
  "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE"
)
stopifnot(all(eu27_iso3 %in% countries$ADM0_A3))

eu27 <- countries |>
  filter(ADM0_A3 %in% eu27_iso3) |>
  st_union()

norway <- countries |>
  filter(ADM0_A3 == "NOR") |>
  st_union()

## EPM2025 coverage ####

# Coverage is a country-level fact: EPM ships one raster per country, and inside a
# covered country an absent (NoData) cell means "no peatland", while outside coverage
# it means "unknown". Conflating the two would manufacture absences. The per-country
# rasters are bounding boxes, not country outlines, so coverage comes from the country
# polygons, not from raster footprints.
epm_raster_dir <- "data/EPM2025/Tegetmeyer_etal_EPM2025_geodata_/EPM2025_raster"
epm_files <- list.files(epm_raster_dir, pattern = "tif.zip$")

# Filenames are not perfectly regular: AUT/BEL carry a trailing underscore, Slovakia is
# EPM24 not EPM25, and the Netherlands is filed as NDL rather than the ISO3 NLD.
epm_iso3 <- substr(epm_files, 1, 3) |>
  recode(NDL = "NLD") |>
  sort()

stopifnot(
  length(epm_iso3) == length(unique(epm_iso3)),
  all(epm_iso3 %in% countries$ADM0_A3)
)

epm_coverage <- countries |>
  filter(ADM0_A3 %in% epm_iso3) |>
  st_union()

cat("EPM2025 country rasters:", length(epm_iso3), "\n")
cat(
  "EU-27 states without EPM coverage:",
  paste(setdiff(eu27_iso3, epm_iso3), collapse = ", "),
  "\n"
)

## Domain ####

domain <- region_iv |>
  st_intersection(eu27) |>
  st_intersection(epm_coverage) |>
  st_intersection(stack_extent) |>
  st_make_valid()

# Norway subtraction: expected to be a no-op under the EU-27 clip. Report, do not assume.
area_before <- as.numeric(sum(st_area(domain))) / 1e6
domain <- st_difference(domain, norway) |>
  st_make_valid()
area_after <- as.numeric(sum(st_area(domain))) / 1e6
cat(
  "Norway subtraction removed",
  round(area_before - area_after, 1),
  "km2 of",
  round(area_before, 1),
  "km2\n"
)

domain_sf <- st_sf(domain_id = 1L, geometry = st_sfc(domain, crs = 3035))

## Verification against plan section 4.1 ####

# Two load-bearing checks. Region IV is partly *defined* by mire occurrence, and the
# EU-27 clip removes most of its area -- either could have cut off the warm flank that
# is the entire reason for recruiting EU data. Recompute the table rather than trust it.
#
# Cell rule: touches = TRUE. A 5 km cell that is only partly inside the domain still
# supplies EU information, and the centre rule would silently drop the coastal and
# border cells that carry much of the warm flank. This is the rule behind the plan's
# section 4.1 table.

domain_mask_5km <- rasterize(
  vect(domain_sf), grid_5km[[1]],
  field = 1L, touches = TRUE
)
names(domain_mask_5km) <- "eu_domain"
writeRaster(
  domain_mask_5km,
  "output/pl0/eu_domain_5km.tif",
  overwrite = TRUE,
  datatype = "INT1U"
)

bio10 <- grid_5km[["bio10"]]
bio10_domain <- as.numeric(values(mask(bio10, domain_mask_5km), na.rm = TRUE))

# Norway's envelope here comes from the 250 m regional stacks. This is a DIAGNOSTIC
# figure only: the stacks cover terrain Lyngstad never surveyed, so their maximum
# overstates what Norway's training data actually cover.
#
# The OPERATIVE warm-flank threshold is published by pl0_labelNorwayBlock.R, computed
# over the Norway training population (the survey footprint), and is what
# pl0_sampleEUabsences.R consumes. That script runs later in the pipeline because it
# needs the presence/absence rasterization; this one cannot see it, and must not pretend
# otherwise.
bio10_norway_cur <- as.numeric(values(
  rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")[["bio10"]],
  na.rm = TRUE
))
bio10_norway_fut <- as.numeric(values(
  rast("output/predictors_regional_250m_Norway_future_EPSG3035.tif")[["bio10"]],
  na.rm = TRUE
))
warm_threshold <- max(bio10_norway_cur)

envelope <- function(v) {
  q <- quantile(v, c(0, 0.5, 0.95, 0.99, 1), names = FALSE)
  sprintf("%.1f / %.1f / %.1f / %.1f / %.1f", q[1], q[2], q[3], q[4], q[5])
}

# The domain earns its place only if it spans where Norway is going.
domain_q99 <- quantile(bio10_domain, 0.99, names = FALSE)

domain_summary <- tibble(
  quantity = c(
    "domain area (km2)",
    "domain cells at 5km (touches rule)",
    "domain bio10 min/med/q95/q99/max",
    "Norway current bio10 min/med/q95/q99/max (250 m grid)",
    "Norway future bio10 min/med/q95/q99/max (250 m grid)",
    "warm-flank threshold (Norway current max bio10)",
    "warm-flank cells in domain",
    "warm-flank cells above the legacy 17.0 threshold",
    "Norway future cells inside domain bio10 range (%)",
    "Norway future cells above domain bio10 q99 (%)"
  ),
  value = c(
    sprintf("%.0f", area_after),
    format(length(bio10_domain)),
    envelope(bio10_domain),
    envelope(bio10_norway_cur),
    envelope(bio10_norway_fut),
    sprintf("%.1f", warm_threshold),
    format(sum(bio10_domain > warm_threshold)),
    format(sum(bio10_domain > 17.0)),
    sprintf("%.1f", 100 * mean(
      bio10_norway_fut >= min(bio10_domain) & bio10_norway_fut <= max(bio10_domain)
    )),
    sprintf("%.1f", 100 * mean(bio10_norway_fut > domain_q99))
  )
)
print(as.data.frame(domain_summary))
write_csv(domain_summary, "output/pl0/eu_domain_summary.csv", append = FALSE)

# Area by country, for the record (plan section 4.1 reports this for full Region IV).
area_by_country <- countries |>
  filter(ADM0_A3 %in% epm_iso3) |>
  select(ADM0_A3) |>
  st_intersection(domain_sf) |>
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) |>
  st_drop_geometry() |>
  group_by(ADM0_A3) |>
  summarise(area_km2 = sum(area_km2), .groups = "drop") |>
  filter(area_km2 > 50) |>
  arrange(desc(area_km2))
print(as.data.frame(area_by_country))
write_csv(
  area_by_country,
  "output/pl0/eu_domain_area_by_country.csv",
  append = FALSE
)

## Plot ####

ggplot() +
  geom_sf(
    data = st_crop(st_geometry(countries), st_bbox(grid_5km)),
    fill = "grey95",
    colour = "grey70",
    linewidth = 0.2
  ) +
  geom_sf(data = domain_sf, fill = "#2c7fb8", colour = NA, alpha = 0.7) +
  coord_sf(expand = FALSE) +
  labs(
    title = "EU training domain",
    subtitle = "Region IV, EU-27, EPM coverage and stack extent; Norway removed"
  ) +
  theme_minimal()

## Save ####

st_write(
  domain_sf,
  "data/DMraisedbog.gpkg",
  layer = "EU_domain",
  append = FALSE,
  quiet = TRUE
)
st_write(
  st_sf(geometry = st_sfc(epm_coverage, crs = 3035)),
  "data/DMraisedbog.gpkg",
  layer = "EPM_coverage",
  append = FALSE,
  quiet = TRUE
)

cat("Wrote data/DMraisedbog.gpkg layers EU_domain, EPM_coverage\n")
cat("Wrote output/pl0/eu_domain_5km.tif\n")

# sessionInfo ####

sessioninfo::session_info()
