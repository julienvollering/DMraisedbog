# Scratch: explore the information content of EPM2025 `peatl_type` ####

# Exploratory, NOT part of any RUNALL pipeline. Two questions about the European
# Peatland Map 2025 (Tegetmeyer et al.), toward judging whether classified fens
# can serve as high-DI hard-negative anchors for pl3 reliability mapping:
#
#   Q1. Breakdown of `peatl_type` per country and in total.
#   Q2. Breakdown of `peatl_type` for the EPM polygons that intersect our EU
#       presences (data/DMraisedbog.gpkg | Presences_EU).
#
# EFFICIENCY NOTES (the gdb is huge: Sweden 9.9M, Finland ~55M polygons):
#  - The sf/OGR FileGDB reader ALWAYS materializes geometry; assembling tens of
#    millions of complex multipolygons takes ~tripple-digit minutes. Avoid it.
#  - Q1 aggregates INSIDE gdal via the SQLITE dialect (`GROUP BY`), so no geometry
#    is built in R (Finland's 29M layer counts in ~8 s).
#  - Q2 needs geometry, but only for the few polygons near each presence: a terra
#    proxy per layer + an extent query hits the FileGDB spatial index and reads
#    just those features.
#
# CAVEAT surfaced while building this: `peatl_type` is NOT harmonized across
# countries (values like "fen", "fen peat", "transition or other", "swamp",
# "undefined", plus blank/NA). Treat the raw values as-is here; a crosswalk is a
# separate decision before any of these become hard negatives.

library(sf)
library(terra)
library(dplyr)
library(readr)
library(tidyr)

Sys.setenv(OGR_ORGANIZE_POLYGONS = "SKIP") # we never need assembled geometry for Q1

gdb <- "data/EPM2025/Tegetmeyer_etal_EPM2025_geodata_/EPM2025_vector/EPM_2025.gdb"
dir.create("output/pl3", showWarnings = FALSE, recursive = TRUE)

lyrs <- st_layers(gdb)$name
# ISO3 country from layer name (Finland is split into FIN_peat_e/u/d -> all "FIN")
iso3_of <- function(layer) sub("_peat.*$", "", layer)

## Q1: peatl_type breakdown per country (layer) and total ####

# GROUP BY runs in gdal's SQLITE dialect; only a tiny count table returns to R.
group_count <- function(layer) {
  dst <- tempfile(fileext = ".csv")
  on.exit(unlink(dst), add = TRUE)
  sql <- sprintf("SELECT peatl_type, COUNT(*) AS n FROM %s GROUP BY peatl_type", layer)
  gdal_utils(
    "vectortranslate", gdb, dst,
    options = c("-f", "CSV", "-dialect", "SQLITE", "-sql", sql)
  )
  read_csv(dst, show_col_types = FALSE) |>
    mutate(peatl_type = if_else(is.na(peatl_type) | peatl_type == "", "<NA>", peatl_type))
}

q1_list <- vector("list", length(lyrs))
for (i in seq_along(lyrs)) {
  cat("Q1", i, "/", length(lyrs), ":", lyrs[i], "\n")
  q1_list[[i]] <- tryCatch(
    group_count(lyrs[i]) |> mutate(layer = lyrs[i], iso3 = iso3_of(lyrs[i]), .before = 1),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
}
q1 <- bind_rows(q1_list)

# Long table by layer, plus a wide country x type matrix and an overall total.
by_country <- q1 |>
  group_by(iso3, peatl_type) |>
  summarise(n = sum(n), .groups = "drop")

wide_country <- by_country |>
  pivot_wider(names_from = peatl_type, values_from = n, values_fill = 0) |>
  arrange(iso3)

total_by_type <- by_country |>
  group_by(peatl_type) |>
  summarise(n = sum(n), .groups = "drop") |>
  arrange(desc(n))

cat("\n=== Q1: peatl_type totals across all countries ===\n")
print(as.data.frame(total_by_type))
cat("\n=== Q1: peatl_type by country (wide) ===\n")
print(as.data.frame(wide_country))

write_csv(q1, "output/pl3/epm_peatl_type_by_layer.csv")
write_csv(wide_country, "output/pl3/epm_peatl_type_by_country_wide.csv")
write_csv(total_by_type, "output/pl3/epm_peatl_type_total.csv")

## Q2: peatl_type of EPM polygons intersecting our EU presences ####

presences <- st_read("data/DMraisedbog.gpkg", layer = "Presences_EU", quiet = TRUE) |>
  st_zm() |> # drop stray Z
  st_make_valid()
stopifnot(st_crs(presences)$epsg == 3035) # must match the gdb (LAEA Europe)

presences$pid <- seq_len(nrow(presences))
presences$geom_type <- as.character(st_geometry_type(presences))

# One representative point per presence (handles mixed POINT/POLYGON geometry),
# used only to shortlist candidate country layers by extent.
rep_xy <- suppressWarnings(st_coordinates(st_point_on_surface(st_geometry(presences))))

# Cheap per-layer extents via proxy (no geometry read).
exts <- lapply(lyrs, function(l) as.vector(ext(vect(gdb, layer = l, proxy = TRUE))))
names(exts) <- lyrs
in_ext <- function(xy, e) xy[1] >= e[1] && xy[1] <= e[2] && xy[2] >= e[3] && xy[2] <= e[4]

# Outer loop over layers so each FileGDB layer is opened once (proxy reused).
hits <- list()
for (l in lyrs) {
  e <- exts[[l]]
  members <- which(apply(rep_xy, 1, in_ext, e = e))
  if (!length(members)) next
  cat("Q2 layer", l, "-- candidate presences:", length(members), "\n")
  v <- vect(gdb, layer = l, proxy = TRUE)
  for (j in members) {
    fi <- vect(st_geometry(presences)[j]) # single-type feature
    sub <- tryCatch(query(v, extent = ext(fi), vars = "peatl_type"),
                    error = function(e) NULL)
    if (is.null(sub) || nrow(sub) == 0) next
    keep <- is.related(sub, fi, "intersects")
    if (!any(keep)) next
    pt <- as.character(sub$peatl_type[keep])
    pt[is.na(pt)] <- "<NA>"
    hits[[length(hits) + 1]] <- tibble(
      pid = presences$pid[j], iso3 = iso3_of(l), layer = l, peatl_type = pt
    )
  }
}
hits <- bind_rows(hits)

# Overall: every EPM polygon that intersects any presence.
overall_q2 <- hits |>
  count(peatl_type, name = "n_polygons", sort = TRUE)

# Per presence: how many polygons, and which types (incl. presences with none).
per_presence <- presences |>
  st_drop_geometry() |>
  select(pid, geom_type) |>
  left_join(
    hits |>
      group_by(pid) |>
      summarise(
        n_polygons = n(),
        types = paste(sort(unique(peatl_type)), collapse = "; "),
        .groups = "drop"
      ),
    by = "pid"
  ) |>
  mutate(n_polygons = coalesce(n_polygons, 0L), types = coalesce(types, "<none>"))

n_zero <- sum(per_presence$n_polygons == 0)

cat("\n=== Q2: peatl_type across all intersecting EPM polygons ===\n")
print(as.data.frame(overall_q2))
cat("\nPresences with >=1 intersecting EPM polygon:",
    sum(per_presence$n_polygons > 0), "of", nrow(presences),
    "(", n_zero, "with none)\n")
cat("\n=== Q2: presence count by set of intersecting types ===\n")
print(as.data.frame(count(per_presence, types, sort = TRUE)))

write_csv(hits, "output/pl3/epm_peatl_type_at_presences_polygons.csv")
write_csv(per_presence, "output/pl3/epm_peatl_type_at_presences_summary.csv")
write_csv(overall_q2, "output/pl3/epm_peatl_type_at_presences_overall.csv")

# sessionInfo ####

sessioninfo::session_info()
