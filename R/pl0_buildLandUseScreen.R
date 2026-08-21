# Human land cover and water fraction on the modelling grid ####

# PURPOSE: Builds the ESA WorldCover fractional cropland/built-up/water/no-data screen at 250 m -- one absence-eligibility rule for both blocks.

# One screen, ONE rule, applied to BOTH blocks (plan_EUintegration.md section 3's
# "one rule, stated once" principle, extended from the absence draw to the absence
# *eligibility*). It exists because the first full projection (notebook 2026-08-19) put
# all 564 Lyngstad raised bogs into `nonpeat`, and the cells anchoring that answer turned
# out to be coastal, urban-fringe land: 48.5% water and 23.8% built-up or cultivated on
# average, against 5.1% and 7.3% for ordinary Norwegian non-peat.
#
# The problem those numbers describe is that `nonpeat` currently conflates three things:
# land whose climate cannot support peat, land that is half sea, and land that was drained
# and built on. Only the first is evidence about climate, and the error runs in the unsafe
# direction -- it charges anthropogenic conversion and open ocean to warming.
#
# WHY THE SAME SCREEN ON BOTH BLOCKS, NOT JUST NORWAY. If Norway's non-peat is screened
# and Europe's is not, then AT MATCHED CLIMATE the EU block retains converted land that
# the Norway block has lost. The EU block is what supplies the warm flank, so the
# asymmetry would land exactly where the projection reads its answer, inflating apparent
# climatic hostility at the warm end -- the same error as before, moved to the other block.
#
# WHY ESA WORLDCOVER AND NOT CORINE. This is not a preference, it is arithmetic. CORINE's
# minimum mapping unit is 25 ha; a 250 m cell is 6.25 ha. Every cell is therefore SMALLER
# than the smallest thing CLC maps, so CLC cannot resolve within-cell composition at all --
# it can only stamp each cell with one class inherited from a polygon at least four times
# the cell's area. A fractional screen is impossible from it. WorldCover is 10 m, giving
# 625 subpixels per cell, and being global it applies identically to both blocks, which is
# the symmetry argument above satisfied by construction rather than by assertion.
#
# WHY FULL RESOLUTION (-ovr NONE). The source COGs carry overviews down to 562x562, and
# gdalwarp will happily use one: a 25x downscale picks an overview around 148 m, leaving
# 2-3 source pixels per target cell and quantising every "fraction" to 0, 1/3, 1/2 or 1.
# Worse, overviews of a categorical palette raster are built by subsampling, not by
# majority, so they are not even a fair sample of composition. Forcing full resolution
# costs ~2 minutes per tile instead of ~2 seconds, and is the whole point of the exercise.
#
# Emits FOUR bands, deliberately keeping water and no-data apart (the same care
# plan section 3 rule 9 takes over EPM's NA): class 80 is mapped permanent water, while
# no-data is territory outside the WorldCover footprint, i.e. open ocean. Both mean "not
# land" for our purposes but they are not the same measurement, and a reader should be
# able to see which is which.

library(sf)
library(terra)
library(dplyr)
library(readr)

dir.create("output/pl0/landuse_parts", showWarnings = FALSE, recursive = TRUE)
dir.create("data/ESA_WorldCover", showWarnings = FALSE, recursive = TRUE)

crs_3035 <- "EPSG:3035"
target_res <- 250 # m, matches both blocks' grids (verified aligned: both on the 250 m lattice)

# Warping at full resolution is the expensive step; TRUE discards cached parts and redoes it.
force_rebuild <- FALSE

wc_version <- "v200"
wc_year <- "2021"
wc_base <- sprintf(
  "https://esa-worldcover.s3.eu-central-1.amazonaws.com/%s/%s/map",
  wc_version, wc_year
)
wc_dir <- "data/ESA_WorldCover"

# WorldCover classes kept, and the band each becomes. Grassland (30) is deliberately NOT
# included: improved pasture lives there but so does semi-natural grassland, and Norway's
# AR50 counts innmarksbeite as agriculture while WorldCover would often call it 30. Rather
# than guess, the mismatch is left out of the screen and reported as a caveat -- and since
# the same rule runs on both blocks, whatever it misses, it misses symmetrically.
wc_bands <- list(
  cropland = 40,
  builtup = 50,
  water = 80
)

## Domain: BOTH blocks ####

domain_eu <- st_read("data/DMraisedbog.gpkg", layer = "EU_domain", quiet = TRUE)
stopifnot(st_crs(domain_eu)$epsg == 3035)

# Norway enters via its own grid's footprint, not via the EU domain, which explicitly
# excludes it (it is the other block).
grid_no <- rast("output/ar50_250m_land_EPSG3035.tif")
domain_no <- st_as_sfc(st_bbox(grid_no))

domain <- st_union(st_union(st_geometry(domain_eu)), domain_no) |> st_make_valid()

## Tiles ####

# WorldCover ships 3 x 3 degree tiles named by their south-west corner. Enumerate those
# that intersect the domain, then keep the ones that exist -- all-ocean tiles are simply
# absent from the bucket, and a HEAD request is cheaper than discovering that mid-warp.
sf_use_s2(FALSE)
domain_ll <- st_transform(domain, 4326)
bb <- st_bbox(domain_ll)

tiles <- expand.grid(
  lon = seq(floor(bb[["xmin"]] / 3) * 3, floor(bb[["xmax"]] / 3) * 3, by = 3),
  lat = seq(floor(bb[["ymin"]] / 3) * 3, floor(bb[["ymax"]] / 3) * 3, by = 3)
) |>
  mutate(
    tile = sprintf(
      "%s%02d%s%03d",
      ifelse(lat >= 0, "N", "S"), abs(lat),
      ifelse(lon >= 0, "E", "W"), abs(lon)
    )
  )

intersects_domain <- vapply(
  seq_len(nrow(tiles)),
  function(i) {
    cell <- st_as_sfc(st_bbox(
      c(
        xmin = tiles$lon[i], xmax = tiles$lon[i] + 3,
        ymin = tiles$lat[i], ymax = tiles$lat[i] + 3
      ),
      crs = st_crs(4326)
    ))
    length(st_intersects(cell, domain_ll)[[1]]) > 0
  },
  logical(1)
)
tiles <- tiles[intersects_domain, ]

tile_url <- function(tile) {
  sprintf(
    "%s/ESA_WorldCover_10m_%s_%s_%s_Map.tif",
    wc_base, wc_year, wc_version, tile
  )
}
tile_path <- function(tile) {
  file.path(
    wc_dir,
    sprintf("ESA_WorldCover_10m_%s_%s_%s_Map.tif", wc_year, wc_version, tile)
  )
}

tiles$available <- file.exists(tile_path(tiles$tile))
cat(
  "Tiles intersecting the domain:", nrow(tiles),
  "| present locally:", sum(tiles$available), "\n"
)
if (!all(tiles$available)) {
  cat(
    "Missing tiles (download from", wc_base, "):\n",
    paste(tiles$tile[!tiles$available], collapse = ", "), "\n"
  )
}
tiles <- tiles[tiles$available, ]

## Per-tile warp to 250 m fractional cover ####

# A VRT with one band per class, each carrying a LUT that turns the categorical source
# into a 0/1 indicator, so a single `-r average` warp yields all fractions in one read of
# the source rather than one read per class.
#
# GDAL's LUT interpolates linearly between the points given, so each class is bracketed by
# its neighbours (39:0, 40:1, 41:0) to make a step function rather than a ramp. Only
# integers occur in the source, so the interpolated interior is never sampled.
#
# The fourth band inverts the source NoData (0): WorldCover's footprint stops offshore, so
# no-data IS open ocean here, and a cell that is half ocean must be visible as such.
build_vrt <- function(tile, vrt_path) {
  src <- tile_path(tile)
  r <- rast(src)
  n_x <- ncol(r)
  n_y <- nrow(r)
  e <- as.vector(ext(r))
  res_x <- (e[2] - e[1]) / n_x
  res_y <- (e[4] - e[3]) / n_y

  luts <- c(
    vapply(
      wc_bands,
      function(code) sprintf("0:0,%d:0,%d:1,%d:0,100:0", code - 1, code, code + 1),
      character(1)
    ),
    nodata = "0:1,1:0,100:0"
  )

  band_xml <- vapply(
    seq_along(luts),
    function(i) {
      sprintf(
        paste0(
          '  <VRTRasterBand dataType="Float32" band="%d">\n',
          "    <Description>%s</Description>\n",
          "    <ComplexSource>\n",
          '      <SourceFilename relativeToVRT="0">%s</SourceFilename>\n',
          "      <SourceBand>1</SourceBand>\n",
          '      <SrcRect xOff="0" yOff="0" xSize="%d" ySize="%d"/>\n',
          '      <DstRect xOff="0" yOff="0" xSize="%d" ySize="%d"/>\n',
          "      <LUT>%s</LUT>\n",
          "    </ComplexSource>\n",
          "  </VRTRasterBand>"
        ),
        i, names(luts)[i], src, n_x, n_y, n_x, n_y, luts[i]
      )
    },
    character(1)
  )

  writeLines(
    sprintf(
      paste0(
        '<VRTDataset rasterXSize="%d" rasterYSize="%d">\n',
        "  <SRS>EPSG:4326</SRS>\n",
        "  <GeoTransform>%.12f, %.14f, 0, %.12f, 0, %.14f</GeoTransform>\n",
        "%s\n</VRTDataset>"
      ),
      n_x, n_y, e[1], res_x, e[4], -res_y,
      paste(band_xml, collapse = "\n")
    ),
    vrt_path
  )

  names(luts)
}

warp_tile <- function(tile) {
  part <- file.path("output/pl0/landuse_parts", paste0(tile, "_250m.tif"))
  if (file.exists(part) && !force_rebuild) {
    cat(tile, "- cached\n")
    return(part)
  }

  vrt_path <- tempfile(fileext = ".vrt")
  on.exit(unlink(vrt_path), add = TRUE)
  band_names <- build_vrt(tile, vrt_path)

  raw <- tempfile(fileext = ".tif")
  on.exit(unlink(raw), add = TRUE)

  t0 <- Sys.time()
  sf::gdal_utils(
    "warp", vrt_path, raw,
    options = c(
      "-t_srs", crs_3035,
      "-tr", target_res, target_res, "-tap",
      "-r", "average",
      # See the header: overviews would silently substitute ~148 m subsampled data.
      "-ovr", "NONE",
      # 0 is the source NoData but it is meaningful here (band 4 counts it), so the mask
      # is switched off and the class is read as an ordinary value.
      "-srcnodata", "None",
      "-ot", "Float32",
      "-co", "COMPRESS=LZW", "-co", "TILED=YES",
      "-multi", "-wo", "NUM_THREADS=ALL_CPUS", "-overwrite"
    ),
    quiet = TRUE
  )

  # Stored as integer percent: the thresholds this feeds are at 10% and 50%, so 1% is
  # ample precision and it keeps a 84 Mcell x 4 band grid at a manageable size.
  r <- rast(raw)
  names(r) <- band_names
  r <- round(r * 100)
  writeRaster(
    r, part,
    overwrite = TRUE, datatype = "INT1U",
    gdal = c("COMPRESS=LZW", "TILED=YES")
  )

  cat(
    tile, "- warped in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n"
  )
  part
}

parts <- vapply(tiles$tile, warp_tile, character(1), USE.NAMES = FALSE)

## Merge onto the common grid ####

screen <- merge(sprc(lapply(parts, rast)))
names(screen) <- c(names(wc_bands), "nodata")

out_path <- "output/pl0/landuse_screen_250m.tif"
writeRaster(
  screen, out_path,
  overwrite = TRUE, datatype = "INT1U",
  gdal = c("COMPRESS=LZW", "TILED=YES")
)
cat("\nWrote", out_path, "\n")
print(screen)

## Summary ####

# `human` and `nonland` are the two quantities the labelling scripts will threshold. They
# are reported here, not applied here: this script builds the measurement, and
# pl0_labelNorwayBlock.R / pl0_sampleEUabsences.R decide what to do with it.
summary_tbl <- global(
  c(
    screen[["cropland"]] + screen[["builtup"]],
    screen[["water"]] + screen[["nodata"]]
  ),
  fun = c("mean"), na.rm = TRUE
)
rownames(summary_tbl) <- c("human_pct", "nonland_pct")
print(summary_tbl)

write_csv(
  tibble(tile = tiles$tile, part = parts),
  "output/pl0/landuse_parts_index.csv",
  append = FALSE
)

# sessionInfo ####

sessioninfo::session_info()
