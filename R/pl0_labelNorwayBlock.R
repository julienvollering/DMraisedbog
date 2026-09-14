# Label and draw the Norway block, 3-class ####

# PURPOSE: Labels the Norwegian block into bog / other-peat / non-peat using artype_60 as a LABEL, applies the land-use screen, and draws the stratified absences.

# The plan's remaining blocker (section 6): Norway's response has to mean the same thing
# as the EU block's before the two can be row-bound. This script builds it.
#
# THE CHANGE IN ONE SENTENCE: artype_60 stops being a MASK and becomes a LABEL.
#
# pl2_createModelingFrame.R currently masks the population to `artype_60 >= 0.5 |
# presence` and fits a binary response. That mask is what deleted the warm flank -- every
# cell that is warm enough to have no peat at all was removed from the population, so the
# model never saw where bogs stop on the warm side, only where they stop on the cold
# side. Under the 3-class reframing (notebook 2026-06-05) the same layer instead splits
# the absences into the two ecologically opposite flanks that raised bog sits between:
#
#   Lyngstad presence         -> bog
#   artype_60 >= 0.5          -> otherpeat   (cold flank: other mire)
#   0.1 <= artype_60 < 0.5    -> EXCLUDED    (mixed cell)
#   artype_60 <  0.1          -> nonpeat     (warm flank: forest, heath, agriculture)
#   artype_60 NA              -> EXCLUDED    (no label available)
#
# The excluded middle band mirrors EPM map_cat 2 on the EU side, which the data producers
# already label as the mixed-cell problem and which rule 9 drops. Without it Norway's
# mixed cells would land in nonpeat while Europe's equivalents are discarded, and
# "other-peat" would not mean the same thing on the two blocks -- which is precisely the
# coherence that pooling depends on.
#
# The absence draw uses the shared sampler in R/functions.R, with the same bins and the
# same budget constant as pl0_sampleEUabsences.R. Retiring the mask grows the Norway
# candidate population from ~63k rows to the full ~1.76M-cell survey footprint, nearly
# all of it nonpeat, so the same climate stratification that the EU block needs applies
# here for the same reason.
#
# The footprint itself is NOT relaxed: absences are only meaningful inside the area
# Lyngstad actually surveyed, so the population stays
# `output/absence_coords_regional.csv` and this script only relabels and subsamples it.

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

# Band edges for artype_60, the fractional peatland cover layer.
otherpeat_min <- 0.5
nonpeat_max <- 0.1

record_settings(
  "R/pl0_labelNorwayBlock.R",
  seed = seed,
  otherpeat_min_artype60 = otherpeat_min,
  nonpeat_max_artype60 = nonpeat_max,
  absences_per_bog_per_class = ABSENCES_PER_BOG_PER_CLASS,
  bio10_stratum_width = unique(diff(BIO10_BREAKS)),
  landuse_max_human_pct = LANDUSE_MAX_HUMAN_PCT,
  landuse_max_nonland_pct = LANDUSE_MAX_NONLAND_PCT
)

## Inputs ####

preds <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")

presence <- read_csv("output/presence_coords_regional.csv", show_col_types = FALSE)
absence <- read_csv("output/absence_coords_regional.csv", show_col_types = FALSE)

## Label ####

# Extract only the two layers needed over 1.76e6 points; pulling the full stack here
# would be wasteful and the predictors are attached downstream anyway.
absence <- absence |>
  mutate(
    artype_60 = terra::extract(
      preds[["artype_60"]], absence[c("x", "y")], ID = FALSE
    )[[1]],
    bio10 = terra::extract(
      preds[["bio10"]], absence[c("x", "y")], ID = FALSE
    )[[1]]
  )

banding <- absence |>
  mutate(
    band = case_when(
      is.na(artype_60) ~ "excluded (artype_60 NA)",
      artype_60 >= otherpeat_min ~ "otherpeat",
      artype_60 < nonpeat_max ~ "nonpeat",
      TRUE ~ "excluded (mixed cell)"
    )
  )

cat("Norway absence population banded by artype_60:\n")
print(as.data.frame(count(banding, band)))

pool <- banding |>
  filter(band %in% c("nonpeat", "otherpeat"), !is.na(bio10)) |>
  mutate(
    response = band,
    stratum = findInterval(bio10, BIO10_BREAKS, rightmost.closed = TRUE)
  )

## Land-use and water screen ####

# Applied to the CANDIDATE POOL, before the draw -- not to the drawn rows afterwards.
# Filtering after the draw would leave the stratified targets under-filled in exactly the
# strata the screen bites hardest, which is the coastal warm end. See
# landuse_screen_keep() in R/functions.R for why this exists and why it runs on both
# blocks and all three classes.
keep_pool <- landuse_screen_keep(pool, label = "NO absence pool")
screen_pool <- tibble(
  class = as.character(pool$response),
  keep = keep_pool
) |>
  count(class, keep) |>
  tidyr::pivot_wider(names_from = keep, values_from = n, values_fill = 0) |>
  rename(kept = `TRUE`, dropped = `FALSE`)
pool <- pool[keep_pool, ]

# Presences get the same rule. This changes n_bog and therefore the absence budget, which
# is intended: the budget is defined per bog cell, so it must count the bog cells that
# actually survive.
keep_pres <- landuse_screen_keep(presence, label = "NO presences")
SCREEN_EFFECT <- bind_rows(
  screen_pool,
  tibble(class = "bog", kept = sum(keep_pres), dropped = sum(!keep_pres))
)
presence <- presence[keep_pres, ]

## Warm flank ####

# Norway's own current maximum: the point beyond which the Norwegian training data say
# nothing, and the reason the EU block exists at all.
warm_threshold <- max(pool$bio10)
warm_strata <- which(BIO10_BREAKS[-1] > warm_threshold)
cat(
  "Norway current max bio10:", round(warm_threshold, 1),
  "-- warm strata above it:",
  if (length(warm_strata)) paste(warm_strata, collapse = ", ") else "none (by definition)",
  "\n"
)

## Budget ####

n_bog <- nrow(presence)
target <- n_bog * ABSENCES_PER_BOG_PER_CLASS
targets <- c(nonpeat = target, otherpeat = target)

cat(
  "Norway bog cells: ", n_bog, " -> ", ABSENCES_PER_BOG_PER_CLASS,
  " per bog per class = ", target, " per absence class\n",
  sep = ""
)

## Draw ####

# warm_strata is empty for Norway by construction (the threshold IS Norway's maximum),
# so the floor pass is inert here and the draw is plain equal-per-stratum water-filling.
drawn <- draw_stratified_absences(as.data.frame(pool), targets, warm_strata)

absences <- drawn$rows |>
  as_tibble() |>
  select(x, y, response, bio10, stratum)

norway_block <- bind_rows(
  presence |>
    select(x, y) |>
    mutate(
      response = "bog",
      bio10 = terra::extract(
        preds[["bio10"]], presence[c("x", "y")], ID = FALSE
      )[[1]],
      stratum = findInterval(bio10, BIO10_BREAKS, rightmost.closed = TRUE)
    ),
  absences
) |>
  mutate(dataset = "NO") |>
  select(x, y, response, dataset, bio10, stratum)

## Report ####

print(drawn$allocation[drawn$allocation$available > 0, ], row.names = FALSE)
write_csv(
  as_tibble(drawn$allocation),
  "output/pl0/no_absence_allocation.csv",
  append = FALSE
)

cat("\nNorway block by class:\n")
print(as.data.frame(count(norway_block, response)))

shortfall <- drawn$allocation |>
  as_tibble() |>
  group_by(response) |>
  summarise(drawn = sum(drawn), .groups = "drop") |>
  mutate(target = target, shortfall = target - drawn)
print(as.data.frame(shortfall))

## Save ####

# The warm-flank threshold is published HERE, not in pl0_buildEUdomain.R, because it is a
# property of the Norway TRAINING POPULATION: the point beyond which Norwegian data say
# nothing. Reading it off the full 250 m raster instead pulls in unsurveyed terrain and
# overstates Norway's coverage, which moves every cut in the section 1.3 sweep.
no_block_summary <- tibble(
  quantity = c(
    "bog cells",
    "absences drawn per class",
    "warm-flank threshold (Norway training max bio10)",
    "Norway training bio10 min/med/q95/q99/max"
  ),
  value = c(
    format(n_bog),
    format(target),
    sprintf("%.1f", warm_threshold),
    {
      q <- quantile(norway_block$bio10, c(0, 0.5, 0.95, 0.99, 1), names = FALSE)
      sprintf("%.1f / %.1f / %.1f / %.1f / %.1f", q[1], q[2], q[3], q[4], q[5])
    }
  )
)
print(as.data.frame(no_block_summary))
write_csv(no_block_summary, "output/pl0/no_block_summary.csv", append = FALSE)

norway_block |>
  write_csv("output/pl0/no_block_coords.csv", append = FALSE)

cat("Wrote output/pl0/no_block_coords.csv (", nrow(norway_block), "rows )\n")

BLOCK_LABEL <- "Norway block"
BLOCK_ALLOC_PNG <- "output/pl0/no_absence_allocation.png"
BLOCK_SCREEN_PNG <- "output/pl0/no_landuse_screen.png"

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
