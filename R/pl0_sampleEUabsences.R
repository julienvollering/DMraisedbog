# Draw the EU absence block ####

# PURPOSE: Draws the EU absence block from domain minus M by the same shared stratified rule and land-use screen as the Norway block.

# Plan section 3, rules 9-12. Absences A = domain \ M, split categorically by the EPM
# map_cat label:
#
#   map_cat 0  no peatland          -> nonpeat
#   map_cat 1  peat dominated       -> otherpeat
#   map_cat 2  peat in soil mosaic  -> EXCLUDED. This class *is* the mixed-cell problem,
#                                     already labelled as such by the data producers.
#
# Rule 9 (NA handling) is settled upstream: pl0_buildEPMmask.R masks map_cat to EPM
# coverage, so an NA here is out-of-coverage and unknown, never a defaulted non-peat.
#
# Rules 11 and 12 are carried by draw_stratified_absences() in R/functions.R, which
# pl0_labelNorwayBlock.R calls with the same bins and the same budget constant. Keeping
# the rule in one place is what makes the two blocks comparable: an "other-peat cell at
# 18 degrees" means the same thing on both sides, which is what the dataset-blocked
# comparison of section 1.5 relies on.

library(sf)
library(terra)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)

source("R/functions.R")
source("R/config.R")

dir.create("output/pl0", showWarnings = FALSE, recursive = TRUE)

seed <- 42
set.seed(seed)

record_settings(
  "R/pl0_sampleEUabsences.R",
  seed = seed,
  absences_per_bog_per_class = ABSENCES_PER_BOG_PER_CLASS,
  landuse_max_human_pct = LANDUSE_MAX_HUMAN_PCT,
  landuse_max_nonland_pct = LANDUSE_MAX_NONLAND_PCT
)

## Inputs ####

map_cat <- rast("output/pl0/epm_map_cat_250m.tif")
mask_M <- rast("output/pl0/eu_exclusion_mask_250m.tif")

# The UNBUFFERED domain is the sampling frame. pl0_buildEPMmask.R and the mask work on a
# buffered domain so that exclusions can reach across the edge; drawing absences from
# that buffer would put rows outside the domain the plan defines.
domain <- st_read("data/DMraisedbog.gpkg", layer = "EU_domain", quiet = TRUE)

# Centre rule here, unlike the domain diagnostics: a drawn cell must actually lie in the
# domain, not merely touch it.
domain_mask <- rasterize(vect(domain), map_cat, field = 1L)

# bio10 is the stratification axis because it is also the axis the section 1.3 sweep cuts
# along; stratifying on one axis and cutting on another would leave the warm end of the
# sweep to chance.
bio10_5km <- rast("output/predictors_global_5km_EUNorway_EPSG3035.tif")[["bio10"]]
bio10 <- resample(crop(bio10_5km, map_cat, snap = "out"), map_cat, method = "near")
names(bio10) <- "bio10"

## Candidate pool ####

# map_cat 2 dropped here, not later: it is an exclusion, not a class.
candidate <- ifel(
  !is.na(domain_mask) & mask_M == 0 & (map_cat == 0 | map_cat == 1),
  map_cat, NA
)
names(candidate) <- "map_cat"

pool <- as.data.frame(c(candidate, bio10), xy = TRUE, na.rm = TRUE) |>
  as_tibble() |>
  filter(!is.na(bio10)) |>
  mutate(
    response = if_else(map_cat == 0, "nonpeat", "otherpeat"),
    stratum = findInterval(bio10, BIO10_BREAKS, rightmost.closed = TRUE)
  )

rm(candidate)
gc()

## Land-use and water screen ####

# Same rule, same thresholds, same function as the Norway block -- see
# landuse_screen_keep() in R/functions.R. Applied to the candidate pool before the draw so
# the stratified targets are filled from eligible cells rather than trimmed afterwards.
keep_pool <- landuse_screen_keep(pool, label = "EU absence pool")
screen_pool <- tibble(
  class = as.character(pool$response),
  keep = keep_pool
) |>
  count(class, keep) |>
  tidyr::pivot_wider(names_from = keep, values_from = n, values_fill = 0) |>
  rename(kept = `TRUE`, dropped = `FALSE`)
pool <- pool[keep_pool, ]

cat("Candidate absence cells:", nrow(pool), "\n")
print(count(pool, response))

## Warm flank ####

# Read from the Norway block, not from eu_domain_summary.csv: the threshold is a
# property of the Norway TRAINING POPULATION, and the domain script can only see the full
# 250 m raster, which includes terrain Lyngstad never surveyed and so overstates what
# Norway actually covers.
warm_threshold <- read_csv(
  "output/pl0/no_block_summary.csv",
  show_col_types = FALSE
) |>
  filter(quantity == "warm-flank threshold (Norway training max bio10)") |>
  pull(value) |>
  as.numeric()

warm_strata <- which(BIO10_BREAKS[-1] > warm_threshold)
cat(
  "Warm-flank threshold:", warm_threshold,
  "-- strata above it:", paste(warm_strata, collapse = ", "), "\n"
)
cat(
  "EU candidate cells beyond Norway's envelope:",
  sum(pool$bio10 > warm_threshold), "of", nrow(pool), "\n"
)

## Budget ####

# Scaled to the block's own presence support, so neither block's absences swamp the
# other and the pooled prior is not reset by whichever domain happens to be larger.
# Counted as unique 250 m CELLS, not as geometries. Several EUNIS Q11 releves can fall in
# one cell of one Natura polygon, and pl2_createModelingFrame.R collapses them to a single
# training row -- 701 geometries become 524 cells. Budgeting off the geometry count would
# make "per bog" mean something different on the two blocks and quietly break the symmetry
# the shared rule exists to provide.
eu_presence_xy <- st_read(
  "data/DMraisedbog.gpkg",
  layer = "Presences_EU", quiet = TRUE
) |>
  st_geometry() |>
  st_point_on_surface() |>
  st_coordinates()

# Screened on the same rule as everything else, so the budget counts the bog cells that
# survive rather than the ones that were mapped. pl2_createModelingFrame.R applies the
# identical screen when it builds the presence ROWS, so the count here and the rows there
# cannot drift apart.
eu_presence_cells <- tibble(
  x = eu_presence_xy[, "X"], y = eu_presence_xy[, "Y"]
) |>
  mutate(cell = cellFromXY(map_cat, cbind(x, y))) |>
  filter(!is.na(cell)) |>
  distinct(cell, .keep_all = TRUE)

keep_pres <- landuse_screen_keep(eu_presence_cells, label = "EU presences")
n_presence_eu <- sum(keep_pres)

SCREEN_EFFECT <- bind_rows(
  screen_pool,
  tibble(class = "bog", kept = sum(keep_pres), dropped = sum(!keep_pres))
)

target <- n_presence_eu * ABSENCES_PER_BOG_PER_CLASS
targets <- c(nonpeat = target, otherpeat = target)

cat(
  "EU presence cells: ", n_presence_eu, " (from ", nrow(eu_presence_xy),
  " geometries) -> ", ABSENCES_PER_BOG_PER_CLASS,
  " per bog per class = ", target, " per absence class\n",
  sep = ""
)

## Draw ####

drawn <- draw_stratified_absences(as.data.frame(pool), targets, warm_strata)

absences <- drawn$rows |>
  as_tibble() |>
  mutate(dataset = "EU") |>
  select(x, y, response, dataset, map_cat, bio10, stratum)

## Report ####

print(drawn$allocation[drawn$allocation$available > 0, ], row.names = FALSE)
write_csv(
  as_tibble(drawn$allocation),
  "output/pl0/eu_absence_allocation.csv",
  append = FALSE
)

cat("\nDrawn absences by class:\n")
print(as.data.frame(count(absences, response)))
cat(
  "On the warm flank (bio10 >", warm_threshold, "):",
  sum(absences$bio10 > warm_threshold), "of", nrow(absences), "\n"
)

shortfall <- drawn$allocation |>
  as_tibble() |>
  group_by(response) |>
  summarise(drawn = sum(drawn), .groups = "drop") |>
  mutate(target = target, shortfall = target - drawn)
print(as.data.frame(shortfall))

## Save ####

absences |>
  write_csv("output/pl0/eu_absence_coords.csv", append = FALSE)

cat("Wrote output/pl0/eu_absence_coords.csv (", nrow(absences), "rows )\n")

BLOCK_LABEL <- "EU block"
BLOCK_ALLOC_PNG <- "output/pl0/eu_absence_allocation.png"
BLOCK_SCREEN_PNG <- "output/pl0/eu_landuse_screen.png"

## Diagnostic: did the stratified draw reach the warm flank? ####

# The single most consequential thing this script can get silently wrong. The absence draw
# exists to put labelled data on the warm climate flank, so a draw that fills its budget
# from the cool core while quietly under-filling the warm strata would leave the projection
# unsupported exactly where it is read. The allocation table says whether that happened;
# this shows it.

alloc_plot <- as_tibble(drawn$allocation) |>
  mutate(
    bio10_lower = BIO10_BREAKS[stratum],
    warm = stratum %in% warm_strata
  ) |>
  filter(available > 0)

p_alloc <- alloc_plot |>
  tidyr::pivot_longer(c(available, drawn), names_to = "what", values_to = "cells") |>
  ggplot(aes(x = bio10_lower, y = cells, fill = what)) +
  geom_col(position = "identity", alpha = 0.75, width = 0.9) +
  facet_wrap(~response, ncol = 1, scales = "free_y") +
  scale_y_log10() +
  scale_fill_manual(values = c(available = "#c8c8c8", drawn = "#2b6a3f")) +
  labs(
    title = paste(BLOCK_LABEL, "-- stratified absence draw by bio10 stratum"),
    subtitle = paste(
      "Grey = cells available in the candidate pool, green = cells actually drawn.",
      "\nLog scale: the point is whether the warm strata are filled at all, not their size."
    ),
    x = "bio10 stratum lower edge (degrees C)", y = "cells (log10)", fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_alloc)
ggsave(BLOCK_ALLOC_PNG, p_alloc, width = 8, height = 6, dpi = 150)

cat("\nFill rate by stratum (drawn / available), warmest strata last:\n")
alloc_plot |>
  arrange(desc(stratum)) |>
  head(12) |>
  transmute(stratum, bio10_lower, response, available, drawn,
            fill_pct = round(100 * drawn / available, 1), warm) |>
  as.data.frame() |>
  print(row.names = FALSE)

## Diagnostic: what the land-use screen removed, by class ####

# Applied to every class, so its effect on the PRESENCES has to be visible next to its
# effect on the absences -- an absences-only reading would hide a reshaped presence set.
p_screen <- SCREEN_EFFECT |>
  tidyr::pivot_longer(c(kept, dropped), names_to = "what", values_to = "cells") |>
  ggplot(aes(x = class, y = cells, fill = what)) +
  geom_col() +
  geom_text(
    data = SCREEN_EFFECT,
    aes(x = class, y = kept + dropped,
        label = sprintf("%.1f%% dropped", 100 * dropped / (kept + dropped))),
    vjust = -0.3, size = 3, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c(kept = "#2b6a3f", dropped = "#c1462c")) +
  labs(
    title = paste(BLOCK_LABEL, "-- effect of the land-use / water screen"),
    subtitle = "Cropland+built-up > 50% or water+no-data > 50%, applied to all classes alike",
    x = NULL, y = "candidate cells", fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_screen)
ggsave(BLOCK_SCREEN_PNG, p_screen, width = 7, height = 4.5, dpi = 150)

# sessionInfo ####

sessioninfo::session_info()
