# Resample EPM2025 map_cat onto the modelling grids ####

# PURPOSE: Warps EPM2025 map_cat to the 250 m grid per country, giving the EU absence label split (0 = non-peat, 1 = other-peat, 2 = excluded) plus a coverage mask.

# Plan section 3 (rules 8-9) and section 5. EPM2025 is the absence-label source for the
# EU block: map_cat 0 = no peatland (non-peat), 1 = peat dominated (other-peat),
# 2 = peat in soil mosaic (excluded -- that class *is* the mixed-cell problem, already
# labelled by the data producers). EPM is a compilation of national peatland maps, not a
# climate-driven model, so using it as an absence source is not circular with respect to
# the CHELSA predictors.
#
# Runs AFTER pl0_buildEUdomain.R: the domain bounds the warp. Only the ten countries
# that contribute domain area are read, which is the difference between minutes and
# hours over ~10 Gcell of 1-arcsec source.
#
# THREE THINGS THE SOURCE DOES NOT TELL YOU UP FRONT:
#
# 1. The NoData value is not harmonised across countries. Luxembourg stores NBITS=1 with
#    NoData=0; Czechia stores NBITS=2 with NoData=3. "No peatland" is therefore encoded
#    as whatever that country's NoData happens to be. Warping with the nodata mask left
#    on would resample only the peat pixels and lose the 0 class entirely, so the mask is
#    switched off (-srcnodata None), the majority is taken over raw codes, and the
#    country's own NoData code is remapped to 0 afterwards. Taking the majority over raw
#    codes is safe because the remap is a relabelling, not a recount.
#
# 2. The per-country rasters are bounding boxes, not country outlines, so a neighbour's
#    "no peatland" background would clobber real data across the border once the nodata
#    mask is off. Each country is therefore warped under a cutline of its own territory
#    intersected with the buffered domain, and the parts are merged afterwards.
#
# 3. Reprojecting a lon/lat bounding box into LAEA leaves uncovered corners. Those must
#    stay NA rather than becoming a manufactured "no peatland", hence -dstnodata 255.
#
# Rule 9 (NA handling) is the reason for the separate coverage mask: inside a covered
# country NA means 0, outside coverage it means unknown, and uncovered cells are
# excluded rather than defaulted to non-peat.

library(sf)
library(terra)
library(dplyr)
library(readr)

dir.create("output/pl0/epm_parts", showWarnings = FALSE, recursive = TRUE)

crs_3035 <- "EPSG:3035"
target_res <- 250 # m, matches the regional predictor grid
buffer_m <- 5000 # >= 1 cell of the coarse grid, so mask/presence work near the edge

# Warping is the expensive step; set TRUE to discard cached country parts and redo it.
force_rebuild <- FALSE

epm_raster_dir <- "data/EPM2025/Tegetmeyer_etal_EPM2025_geodata_/EPM2025_raster"

## Domain and countries ####

domain <- st_read("data/DMraisedbog.gpkg", layer = "EU_domain", quiet = TRUE)
stopifnot(st_crs(domain)$epsg == 3035)
domain_buf <- st_buffer(domain, buffer_m) |> st_make_valid()

countries <- st_read(
  "data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp",
  quiet = TRUE
) |>
  st_transform(crs_3035) |>
  st_make_valid()

epm_files <- list.files(epm_raster_dir, pattern = "tif.zip$")
epm_lookup <- tibble(
  iso3 = recode(substr(epm_files, 1, 3), NDL = "NLD"),
  file = epm_files
)

# Cutlines: country territory clipped to the buffered domain. This both masks the
# neighbouring-bbox problem and, via -crop_to_cutline, keeps gdalwarp from reading
# source blocks that can never contribute a row (most of France, Sweden, Finland).
# Norway is dropped even though the buffer reaches across the border: it is the other
# block, labelled from artype + Lyngstad, and must never enter the EU absence pool. The
# resulting uncovered ring along the border is the honest answer, not a gap to fill.
cutlines <- countries |>
  select(ADM0_A3) |>
  filter(ADM0_A3 %in% epm_lookup$iso3, ADM0_A3 != "NOR") |>
  st_intersection(st_geometry(domain_buf)) |>
  st_make_valid() |>
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) |>
  filter(area_km2 > 1)

cutline_path <- "output/pl0/epm_parts/cutlines.gpkg"
st_write(cutlines, cutline_path, layer = "cut", append = FALSE, quiet = TRUE)

work <- epm_lookup |>
  filter(iso3 %in% cutlines$ADM0_A3) |>
  left_join(st_drop_geometry(cutlines), by = c("iso3" = "ADM0_A3")) |>
  arrange(desc(area_km2))
print(as.data.frame(work))

## Warp each country to the 250 m grid ####

# Read a band's NoData value from GDAL rather than terra: terra reports NAflag as NaN
# for these files while still masking on the real value, so it cannot be trusted here.
epm_nodata <- function(src) {
  info <- sf::gdal_utils("info", src, quiet = TRUE)
  m <- regmatches(info, regexpr("NoData Value=[-0-9.]+", info))
  if (length(m) == 0) {
    return(NA_real_)
  }
  as.numeric(sub("NoData Value=", "", m))
}

warp_country <- function(iso3, file) {
  part <- file.path("output/pl0/epm_parts", paste0(iso3, "_map_cat_250m.tif"))
  if (file.exists(part) && !force_rebuild) {
    cat(iso3, "- cached\n")
    return(part)
  }

  src <- file.path(
    "/vsizip", epm_raster_dir, file, sub(".zip", "", file, fixed = TRUE)
  )
  nodata <- epm_nodata(src)

  # If a country ever encoded NoData as 1 or 2 the remap below would destroy a real
  # class, so refuse rather than silently mislabel.
  stopifnot(!is.na(nodata), !(nodata %in% c(1, 2)))

  raw <- tempfile(fileext = ".tif")
  on.exit(unlink(raw), add = TRUE)

  sf::gdal_utils(
    "warp", src, raw,
    options = c(
      "-t_srs", crs_3035,
      "-tr", target_res, target_res, "-tap",
      "-r", "mode",
      "-srcnodata", "None", "-dstnodata", "255", "-ot", "Byte",
      "-cutline", cutline_path, "-cl", "cut",
      "-cwhere", sprintf("ADM0_A3 = '%s'", iso3),
      "-crop_to_cutline",
      "-co", "COMPRESS=LZW", "-co", "TILED=YES",
      "-multi", "-wo", "NUM_THREADS=ALL_CPUS", "-overwrite"
    ),
    quiet = TRUE
  )

  # Relabel this country's NoData code to 0 = no peatland. Cells left at 255 are the
  # uncovered reprojection corners and stay NA.
  r <- rast(raw)
  NAflag(r) <- 255
  r <- subst(r, nodata, 0L)
  writeRaster(
    r, part,
    overwrite = TRUE, datatype = "INT1U",
    gdal = c("COMPRESS=LZW", "TILED=YES")
  )

  cat(iso3, "- warped (source NoData =", nodata, ")\n")
  part
}

parts <- mapply(warp_country, work$iso3, work$file, USE.NAMES = FALSE)

# The source band types are as unharmonised as the NoData values -- Byte, Int8 and
# Float32 all occur -- so -ot Byte silently clamps out-of-range codes. Finland's Int8
# -128 lands on 0, which happens to be the label we want. Do not rely on that: assert
# the class set on every part, cached ones included, so a country with an encoding this
# script has not seen fails loudly instead of contributing mislabelled absences.
for (p in parts) {
  observed <- freq(rast(p))$value
  if (!all(observed %in% c(0L, 1L, 2L))) {
    stop(
      basename(p), ": unexpected map_cat codes ",
      paste(setdiff(observed, c(0L, 1L, 2L)), collapse = ", ")
    )
  }
}

## Merge onto the common 250 m grid ####

map_cat_250m <- merge(sprc(lapply(parts, rast)))
names(map_cat_250m) <- "map_cat"

## Coverage, and rule 9 ####

# Coverage is the union of countries EPM actually ships, clipped to the buffered domain.
# Inside it, an NA left by the merge is a real "no peatland"; outside it, NA is unknown
# and the cell must be excluded from the absence pool rather than defaulted to non-peat.
coverage_vect <- cutlines |>
  st_union() |>
  st_sf(geometry = _) |>
  vect()

coverage_250m <- rasterize(coverage_vect, map_cat_250m, field = 1L, touches = TRUE)
names(coverage_250m) <- "coverage"

map_cat_250m <- ifel(is.na(map_cat_250m) & coverage_250m == 1, 0L, map_cat_250m)
map_cat_250m <- mask(map_cat_250m, coverage_250m)

writeRaster(
  map_cat_250m, "output/pl0/epm_map_cat_250m.tif",
  overwrite = TRUE, datatype = "INT1U",
  gdal = c("COMPRESS=LZW", "TILED=YES")
)
writeRaster(
  coverage_250m, "output/pl0/epm_coverage_250m.tif",
  overwrite = TRUE, datatype = "INT1U",
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

## 5 km version ####

# Derived from the 250 m product rather than from a second pass over the source: a
# second warp would double the read cost for a layer used only for domain-level
# diagnostics and climate stratification. Majority-of-majorities is an approximation and
# is not used for any label.
grid_5km <- rast("output/predictors_global_5km_EUNorway_EPSG3035.tif")[[1]]
target_5km <- crop(grid_5km, ext(map_cat_250m), snap = "out")

map_cat_5km <- resample(map_cat_250m, target_5km, method = "mode")
names(map_cat_5km) <- "map_cat"
coverage_5km <- resample(coverage_250m, target_5km, method = "near")
names(coverage_5km) <- "coverage"
map_cat_5km <- mask(map_cat_5km, coverage_5km)

writeRaster(
  map_cat_5km, "output/pl0/epm_map_cat_5km.tif",
  overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=LZW")
)
writeRaster(
  coverage_5km, "output/pl0/epm_coverage_5km.tif",
  overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=LZW")
)

## Report ####

label <- c("0" = "non-peat", "1" = "other-peat", "2" = "excluded (soil mosaic)")

summarise_grid <- function(r, res_label) {
  f <- freq(r) |>
    as_tibble() |>
    transmute(
      grid = res_label,
      map_cat = value,
      label = label[as.character(value)],
      cells = count,
      area_km2 = count * (res(r)[1] / 1000)^2
    )
  f |> mutate(share = cells / sum(cells))
}

epm_summary <- bind_rows(
  summarise_grid(map_cat_250m, "250 m"),
  summarise_grid(map_cat_5km, "5 km")
)
print(as.data.frame(epm_summary))
write_csv(epm_summary, "output/pl0/epm_map_cat_summary.csv", append = FALSE)

plot(map_cat_250m, main = "EPM2025 map_cat, 250 m, EU domain")

cat("Wrote output/pl0/epm_map_cat_{250m,5km}.tif and epm_coverage_{250m,5km}.tif\n")

# sessionInfo ####

sessioninfo::session_info()
