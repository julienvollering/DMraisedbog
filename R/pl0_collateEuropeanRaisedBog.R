# Collate European raised-bog presences and the exclusion mask ####

# PURPOSE: Selects EU raised-bog presences from Natura 2000 and EUNIS under stated evidence bars, and builds the exclusion mask M.

# Plan section 3, rules 1-8 and 10. Two products, deliberately asymmetric:
#
#   Presences P  -- narrow. Natura 2000 7110/7120 only, G-quality, representativity D
#                   excluded, inside the unbuffered domain.
#   Exclusion M  -- maximally inclusive, and deliberately BROADER than P. A broader mask
#                   removes false absences at no cost to the presence definition, so it
#                   takes 7130 blanket bog, every DATAQUALITY and REPRESENTATIVITY grade,
#                   EUNIS Q11 and Q12 plots, and EPM bog-like polygons.
#
# Rule 10 is the point of the asymmetry: a large Natura polygon that fails the presence
# rules stays in M, so it is excluded from BOTH absence classes and becomes unlabelled
# rather than being quietly relabelled other-peat or non-peat. That is structural here,
# not a filter that can silently fail.
#
# EUNIS Q-codes are used for LOCATION only, never as labels (plan section 2). Poor fen
# and quaking mire are constituents of raised-bog complexes -- lagg, hollows -- so a Q22
# relevé inside a bog complex would be a spurious absence, systematically rather than
# randomly. And the response is active raised-bog *function*, a landform and hydrology
# concept; vegetation-type labels would silently redefine it.
#
# Site attributes are not used. COVER_HA, PERCENTAGE_COVER, CONSERVATION and RELSURFACE
# are not trusted and play no part in any rule (plan section 3, section 4.3).

library(sf)
library(terra)
library(dplyr)
library(readr)
library(ggplot2)

source("R/config.R")

dir.create("output/pl0", showWarnings = FALSE, recursive = TRUE)

crs_3035 <- "EPSG:3035"

# Buffers. The plan left "at least 1 cell" open (section 6); it is fixed here at one
# 250 m cell, applied twice for different reasons:
#   - plot_buffer_m: EUNIS relevés are ~1 m2 point records with positional error, so a
#     plot marks its cell, not its coordinate.
#   - the final one-cell dilation of M, so a bog that clips the corner of a cell does not
#     leave the neighbouring cell available as an absence.
plot_buffer_m <- 250
domain_buffer_m <- 5000 # matches pl0_buildEPMmask.R, so M can reach the domain edge

## Grids and domain ####

# The 250 m EU grid is the one pl0_buildEPMmask.R already established, so M lands on
# exactly the grid the absence sampler will draw from.
grid_250m <- rast("output/pl0/epm_map_cat_250m.tif")

domain <- st_read("data/DMraisedbog.gpkg", layer = "EU_domain", quiet = TRUE)
stopifnot(st_crs(domain)$epsg == 3035)
domain_buf <- st_buffer(domain, domain_buffer_m) |> st_make_valid()

## Natura 2000 ####

natura_gpkg <- "data/Natura2000_end2021_rev1_gpkg/Natura2000_end2021_rev1.gpkg"
sites <- st_read(natura_gpkg, layer = "NaturaSite_polygon", quiet = TRUE)
habitats <- st_read(natura_gpkg, layer = "HABITATS", quiet = TRUE)

# Mask side: every raised-bog and blanket-bog record, no quality filter at all.
hab_mask <- habitats |>
  filter(HABITATCODE %in% c("7110", "7120", "7130"))

# Presence side: rule 5 (7110/7120, G quality) and rule 7 (drop representativity D,
# keep NA). Excluding on missing metadata would induce bias correlated with national
# reporting completeness -- the same error pattern avoided with PRECISION -- and the NA
# rows are immaterial either way.
hab_pres <- habitats |>
  filter(
    HABITATCODE %in% c("7110", "7120"),
    DATAQUALITY == "G",
    is.na(REPRESENTATIVITY) | REPRESENTATIVITY != "D"
  )

cat("Natura habitat records (Europe-wide):\n")
print(count(habitats |> filter(HABITATCODE %in% c("7110", "7120", "7130")), HABITATCODE))
cat(
  "  presence-eligible after quality and representativity rules:",
  nrow(hab_pres), "\n"
)
cat(
  "  dropped by REPRESENTATIVITY == D:",
  sum(habitats$HABITATCODE %in% c("7110", "7120") &
    habitats$DATAQUALITY == "G" &
    habitats$REPRESENTATIVITY %in% "D", na.rm = TRUE),
  "; retained with REPRESENTATIVITY NA:",
  sum(habitats$HABITATCODE %in% c("7110", "7120") &
    habitats$DATAQUALITY == "G" &
    is.na(habitats$REPRESENTATIVITY)),
  "\n"
)

natura_polygons <- function(hab, clip_to) {
  sites |>
    semi_join(st_drop_geometry(hab), by = "SITECODE") |>
    st_cast("POLYGON", warn = FALSE) |>
    st_make_valid() |>
    st_filter(clip_to, .predicate = st_intersects)
}

sites_mask <- natura_polygons(hab_mask, domain_buf)
sites_pres <- natura_polygons(hab_pres, domain) |>
  left_join(
    hab_pres |>
      st_drop_geometry() |>
      group_by(SITECODE) |>
      summarise(
        HABITATCODE = paste(sort(unique(HABITATCODE)), collapse = "/"),
        REPRESENTATIVITY = paste(sort(unique(REPRESENTATIVITY)), collapse = "/"),
        .groups = "drop"
      ),
    by = "SITECODE"
  ) |>
  mutate(area_km2 = as.numeric(st_area(geom)) / 1e6)

cat("Natura polygons in domain -- mask:", nrow(sites_mask),
    " presence-eligible:", nrow(sites_pres), "\n")

## EUNIS plots ####

# Location only. Q11 raised bog drives presence recovery inside large polygons; Q12
# blanket bog is mask-only, exactly as 7130 is on the Natura side.
eunis <- st_read("data/EUNIS/Wetland_Plot.gpkg", layer = "Wetland", quiet = TRUE) |>
  st_transform(crs_3035)
eunis_q11 <- filter(eunis, HabitatCod == "Q11")
eunis_mask <- filter(eunis, HabitatCod %in% c("Q11", "Q12"))

cat("EUNIS plots -- Q11:", nrow(eunis_q11), " Q11+Q12:", nrow(eunis_mask), "\n")

## Presences (rules 5-7) ####

# Rule 6. Polygons finer than the modelling resolution are taken directly; larger ones
# are only usable where a Q11 relevé pins the bog inside the designated complex.
small <- sites_pres |> filter(area_km2 < 1)
large <- sites_pres |> filter(area_km2 >= 1)

presences_small <- small |>
  transmute(
    source = "natura_small",
    SITECODE, HABITATCODE, REPRESENTATIVITY,
    geom
  )

presences_q11 <- eunis_q11[large, ] |>
  transmute(
    source = "eunis_q11_in_large_natura",
    SITECODE = NA_character_,
    HABITATCODE = "7110/7120",
    REPRESENTATIVITY = NA_character_,
    geom = Shape
  ) |>
  st_set_geometry("geom")

presences <- bind_rows(presences_small, presences_q11)

# Plan section 4.2: the large-polygon loss is severe but climate-neutral. Recompute it
# rather than trust the table -- it is the evidence that the loss costs power, not
# coverage, and it is the first thing a reviewer will probe.
large_with_q11 <- lengths(st_intersects(large, eunis_q11)) > 0
bio10 <- rast("output/predictors_global_5km_EUNorway_EPSG3035.tif")[["bio10"]]

bio10_at <- function(g) {
  if (nrow(g) == 0) {
    return(numeric(0))
  }
  v <- terra::extract(bio10, vect(st_point_on_surface(st_geometry(g))), ID = FALSE)[[1]]
  v[!is.na(v)]
}

envelope <- function(v) {
  if (length(v) == 0) {
    return(NA_character_)
  }
  q <- quantile(v, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
  sprintf("%.1f / %.1f / %.1f / %.1f / %.1f", q[1], q[2], q[3], q[4], q[5])
}

presence_audit <- tibble(
  set = c("small, retained", "large with Q11, retained", "large, discarded"),
  polygons = c(nrow(small), sum(large_with_q11), sum(!large_with_q11)),
  area_km2 = c(
    sum(small$area_km2),
    sum(large$area_km2[large_with_q11]),
    sum(large$area_km2[!large_with_q11])
  ),
  bio10_min_q25_med_q75_max = c(
    envelope(bio10_at(small)),
    envelope(bio10_at(large[large_with_q11, ])),
    envelope(bio10_at(large[!large_with_q11, ]))
  )
)
print(as.data.frame(presence_audit))
write_csv(presence_audit, "output/pl0/eu_presence_audit.csv", append = FALSE)

cat("Presences -- small polygons:", nrow(presences_small),
    " Q11 points in large polygons:", nrow(presences_q11),
    " total:", nrow(presences), "\n")

## Exclusion mask M (rules 1-4) ####

# Built as a raster on the 250 m EU grid, because that is the form the absence sampler
# consumes. The vector components are written alongside it for inspection.
mask_vector <- bind_rows(
  sites_mask |> transmute(component = "natura_7110_7120_7130", geom),
  eunis_mask |>
    st_buffer(plot_buffer_m) |>
    transmute(component = "eunis_q11_q12_buffered", geom = Shape) |>
    st_set_geometry("geom")
)

mask_rast <- rasterize(vect(mask_vector), grid_250m, field = 1L, touches = TRUE)

# Rule 3: EPM bog-like polygons, where available. This is the component that stops a
# mapped European bog from being labelled other-peat by the EPM raster, so it is not
# optional decoration -- without it the EU block would manufacture false negatives for
# the bog class exactly where bogs are known.
#
# Rasterized straight out of the FileGDB by GDAL: the bog-like classes run to millions
# of polygons and materializing their geometry in R would take hours (see
# R/archive/scratch_exploreEPM2025.R). Note st_layers() on this gdb is itself very slow, so the
# layer names are constructed rather than listed; Finland ships three, and Austria ships
# none.
Sys.setenv(OGR_ORGANIZE_POLYGONS = "SKIP")
epm_gdb <- "data/EPM2025/Tegetmeyer_etal_EPM2025_geodata_/EPM2025_vector/EPM_2025.gdb"

domain_iso3 <- read_csv(
  "output/pl0/eu_domain_area_by_country.csv",
  show_col_types = FALSE
)$ADM0_A3
epm_layers <- c(
  paste0(setdiff(domain_iso3, "FIN"), "_peat"),
  paste0("FIN_peat_", c("d", "e", "u"))
)

# One SQL predicate, stated rather than hidden: peatl_type carries 60+ unharmonised
# values, and every bog-like one contains the substring (inventory in
# R/archive/scratch_exploreEPM2025.R; the full scan takes too long to repeat here).
bog_where <- "peatl_type LIKE '%bog%'"

record_settings(
  "R/pl0_collateEuropeanRaisedBog.R",
  plot_buffer_m = plot_buffer_m,
  domain_buffer_m = domain_buffer_m,
  natura_presence_codes = c("7110", "7120"),
  natura_mask_codes = c("7110", "7120", "7130"),
  large_polygon_km2 = 1,
  epm_bog_predicate = bog_where
)

# One output file per layer, combined once at the end. Two constraints force this shape:
#   - gdal_rasterize requires -te/-tr here and then CREATES the target, so successive
#     calls against a single path overwrite each other instead of accumulating; burning
#     all twelve layers into one file would keep only the last.
#   - terra's max() returns a lazy raster that still references its inputs, so deleting
#     a part inside the loop silently turns the accumulator into NA. Parts are therefore
#     kept on disk until the combine, and cached across runs -- Finland's three layers
#     alone are ~29M polygons and take the better part of an hour to scan.
# Austria ships no vector layer, so a missing layer is expected and reported, not fatal.
te <- c(xmin(grid_250m), ymin(grid_250m), xmax(grid_250m), ymax(grid_250m))
force_rebuild_bog <- FALSE

burn_layer <- function(lyr) {
  part <- file.path("output/pl0/epm_parts", paste0("bog_", lyr, ".tif"))
  if (file.exists(part) && !force_rebuild_bog) {
    cat("  EPM layer", lyr, "- cached\n")
    return(part)
  }
  ok <- tryCatch(
    {
      sf::gdal_utils(
        "rasterize", epm_gdb, part,
        options = c(
          "-l", lyr, "-where", bog_where, "-burn", "1", "-init", "0",
          "-ot", "Byte", "-a_srs", crs_3035,
          "-te", te[1], te[2], te[3], te[4], "-tr", "250", "250",
          "-co", "COMPRESS=LZW", "-co", "TILED=YES"
        ),
        quiet = TRUE
      )
      TRUE
    },
    error = function(e) {
      cat("  EPM layer", lyr, "unavailable:", conditionMessage(e), "\n")
      FALSE
    }
  )
  if (!ok) {
    unlink(part)
    return(NA_character_)
  }
  # An all-NoData part means the country carries no bog-like peatl_type at all, which
  # is a real answer, not a failure: Sweden is exactly that case. Report it as 0 rather
  # than letting terra's na.rm sum over an empty set surface as NaN.
  burned_cells <- as.numeric(global(rast(part), "sum", na.rm = TRUE))
  if (is.na(burned_cells)) burned_cells <- 0
  cat("  EPM layer", lyr, "burned; bog-like cells:", burned_cells, "\n")
  part
}

dir.create("output/pl0/epm_parts", showWarnings = FALSE, recursive = TRUE)
bog_parts <- vapply(epm_layers, burn_layer, character(1), USE.NAMES = FALSE)
bog_parts <- bog_parts[!is.na(bog_parts)]
stopifnot(length(bog_parts) > 0)

# gdal_rasterize leaves un-burned cells as NoData rather than 0, so "any layer burned
# here" is a na.rm sum across the parts.
# A cell where every part is NoData must come out 0, not NA: terra's na.rm sum returns
# NA when all inputs are NA, and an NA here would propagate into the mask union below.
epm_bog <- sum(rast(bog_parts), na.rm = TRUE)
epm_bog <- ifel(is.na(epm_bog) | epm_bog == 0, 0L, 1L)

epm_bog_path <- "output/pl0/epm_boglike_250m.tif"
writeRaster(
  epm_bog, epm_bog_path,
  overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=LZW", "TILED=YES")
)
epm_bog <- rast(epm_bog_path)
cat("EPM bog-like cells in mask:", as.numeric(global(epm_bog, "sum", na.rm = TRUE)), "\n")

# Rule 4: union, then dilate by one cell. focal max over a 3x3 window is exactly a
# one-cell dilation on this grid.
mask_M <- (!is.na(mask_rast)) | (epm_bog == 1)
mask_M <- focal(mask_M, w = 3, fun = "max", na.policy = "omit", na.rm = TRUE)
mask_M <- ifel(mask_M >= 1, 1L, 0L)
mask_M <- mask(mask_M, grid_250m)
names(mask_M) <- "exclusion_mask"

writeRaster(
  mask_M, "output/pl0/eu_exclusion_mask_250m.tif",
  overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=LZW", "TILED=YES")
)

## Report ####

cell_km2 <- (res(grid_250m)[1] / 1000)^2
mask_cells <- as.numeric(global(mask_M, "sum", na.rm = TRUE))
domain_cells <- as.numeric(global(!is.na(grid_250m), "sum", na.rm = TRUE))

mask_summary <- tibble(
  quantity = c(
    "domain cells at 250 m",
    "cells in exclusion mask M",
    "mask share of domain",
    "mask area (km2)",
    "presences (small Natura polygons)",
    "presences (EUNIS Q11 in large Natura polygons)",
    "presences total"
  ),
  value = c(
    format(domain_cells),
    format(mask_cells),
    sprintf("%.3f", mask_cells / domain_cells),
    sprintf("%.0f", mask_cells * cell_km2),
    format(nrow(presences_small)),
    format(nrow(presences_q11)),
    format(nrow(presences))
  )
)
print(as.data.frame(mask_summary))
write_csv(mask_summary, "output/pl0/eu_mask_summary.csv", append = FALSE)

plot(mask_M, main = "EU exclusion mask M, 250 m")

## Diagnostic: the EU presence set, drawn ####

# Provenance is the reviewable question here -- these presences are the only EU signal the
# pooled model gets, and they are assembled from two sources under several filters. The map
# shows what survived and where, so the geography of the evidence is visible rather than
# implied by counts. The warm south is what the EU block was recruited for, so a presence
# set that turned out to be all Baltic would undercut the whole architecture.

countries_bg <- st_read(
  "data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp", quiet = TRUE
) |>
  st_transform(3035) |>
  st_geometry()

pres_xy <- presences |>
  st_point_on_surface() |>
  st_coordinates() |>
  as_tibble() |>
  rename(x = X, y = Y) |>
  mutate(source = presences$source)

bb <- st_bbox(domain)

p_map <- ggplot() +
  geom_sf(data = countries_bg, fill = "grey96", colour = "grey80", linewidth = 0.2) +
  geom_sf(data = st_geometry(domain), fill = "#dbe7d8", colour = "#2b6a3f",
          linewidth = 0.3, alpha = 0.5) +
  geom_point(data = pres_xy, aes(x = x, y = y, colour = source), size = 0.9, alpha = 0.85) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"])) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    title = "EU raised-bog presences entering the pooled training frame",
    subtitle = paste0(
      nrow(presences), " presence geometries inside the domain (green), by source"
    ),
    x = NULL, y = NULL, colour = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_map)
ggsave("output/pl0/eu_presences_map.png", p_map, width = 8, height = 8, dpi = 150)

# Presences by latitude band: the EU block exists to supply the WARM flank, so where the
# presences sit north-to-south is the thing to check, not just how many there are.
p_lat <- pres_xy |>
  mutate(northing_km = y / 1000) |>
  ggplot(aes(x = northing_km, fill = source)) +
  geom_histogram(bins = 40, alpha = 0.85) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Where the EU presences sit, north to south",
    x = "EPSG:3035 northing (km)", y = "presence geometries", fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_lat)
ggsave("output/pl0/eu_presences_latitude.png", p_lat, width = 8, height = 4, dpi = 150)

## Save ####

st_write(presences, "data/DMraisedbog.gpkg", layer = "Presences_EU",
         append = FALSE, quiet = TRUE)
st_write(mask_vector, "data/DMraisedbog.gpkg", layer = "ExclusionMask_EU_vector",
         append = FALSE, quiet = TRUE)

cat("Wrote Presences_EU and ExclusionMask_EU_vector to data/DMraisedbog.gpkg\n")
cat("Wrote output/pl0/eu_exclusion_mask_250m.tif and epm_boglike_250m.tif\n")

# sessionInfo ####

sessioninfo::session_info()
