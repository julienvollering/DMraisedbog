# Milder against harsher future: SSP1-2.6 and SSP3-7.0 compared ####

# PURPOSE: Compares the production projection under the milder scenario with the primary one: bog area by class, current-to-future transitions, and which present-day bog stays bog (the cells and Lyngstad polygons that "still stand a chance").

# The model is fitted once on current climate and never sees either scenario, so the two
# projections differ only in the future climate they are handed. This script reads both and
# asks how sensitive the predicted loss is to the degree of warming. Four pieces:
#
#  1. Climate. How much of the bio10 shift at present-day bog cells does the milder
#     scenario avoid? Everything else is read against this number.
#  2. Area. Argmax class area and expected bog area (sum of P(bog) x cell area) under
#     current / SSP1-2.6 / SSP3-7.0. The expected area uses the probability and so does not
#     hide small shifts behind a class boundary.
#  3. Transitions. The 3 x 3 current -> future argmax table per scenario, in km2.
#  4. Persistence. Cells predicted bog today, split by whether they are still bog under
#     both scenarios, under the milder one only (the cells that still stand a chance), under
#     the harsher one only, or under neither; and the same split over the Lyngstad polygons.
#     Each persistence class is also read against the area of applicability from
#     pl2_mapReliability.R, because a cell that "survives" outside the area where the model
#     is reliable is weaker evidence than one that survives inside it.
#
# Argmax is the operating point of the balanced forest, as in pl2_interpret.R.

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(sf)
library(terra)

source("R/functions.R")
source("R/config.R")

# Layer prefixes the prediction rasters use for each scenario (see pl2_predict.R).
SCEN <- tibble(
  scenario = c("current", FUTURE_SCENARIO, MILD_SCENARIO),
  prefix = c(
    "current_", "future_", paste0("future", scenario_suffix(MILD_SCENARIO), "_")
  ),
  stack = c(
    "output/pl2/scenario_current.tif",
    "output/pl2/scenario_future.tif",
    paste0("output/pl2/scenario_future", scenario_suffix(MILD_SCENARIO), ".tif")
  ),
  reliability = c(
    "output/pl2/reliability_current.tif",
    "output/pl2/reliability_future.tif",
    paste0("output/pl2/reliability_future", scenario_suffix(MILD_SCENARIO), ".tif")
  )
)
HARSH <- FUTURE_SCENARIO
MILD <- MILD_SCENARIO
AXIS <- PARTITION_AXIS

record_settings(
  "R/pl2_compareScenarios.R",
  scenarios = SCEN$scenario,
  persistence_reference = "argmax class at current climate"
)

## Load predictions ####

pred_path <- c(
  "output/pl2/rf_local_pred_brf.tif",
  paste0("output/pl2/rf_local_pred_brf", scenario_suffix(MILD), ".tif")
)
assert_fresher(pred_path[2], SCEN$stack[SCEN$scenario == MILD])
pred <- c(rast(pred_path[1]), rast(pred_path[2]))
stopifnot(all(unlist(map(SCEN$prefix, ~ paste0(.x, RESPONSE_LEVELS))) %in% names(pred)))

P <- values(pred)
ok <- complete.cases(P)
P <- P[ok, , drop = FALSE]
cell_km2 <- prod(res(pred)) / 1e6
cat("Cells with a prediction:", nrow(P), "| cell area:", cell_km2, "km2\n\n")

argmax_of <- function(prefix) {
  m <- P[, paste0(prefix, RESPONSE_LEVELS), drop = FALSE]
  factor(RESPONSE_LEVELS[max.col(m, ties.method = "first")], levels = RESPONSE_LEVELS)
}
cls <- set_names(map(SCEN$prefix, argmax_of), SCEN$scenario)
pbog <- set_names(map(SCEN$prefix, ~ P[, paste0(.x, "bog")]), SCEN$scenario)

## 1. Climate: how much warming does the milder scenario avoid? ####

# bio10 (the pinned partition axis) at the cells that are predicted bog today, read from
# the same scenario stacks the model projected on.
bog_now <- cls[["current"]] == "bog"
bio10_at <- function(path) {
  v <- values(rast(path)[[AXIS]])[ok]
  v
}
bio10 <- set_names(map(SCEN$stack, bio10_at), SCEN$scenario)

climate <- map_dfr(SCEN$scenario, function(s) {
  tibble(
    scenario = s,
    bio10_median_domain = median(bio10[[s]]),
    bio10_median_bog_cells = median(bio10[[s]][bog_now]),
    shift_domain = median(bio10[[s]] - bio10[["current"]]),
    shift_bog_cells = median(bio10[[s]][bog_now] - bio10[["current"]][bog_now])
  )
})
cat("bio10 at present-day predicted-bog cells (n =", sum(bog_now), "):\n")
print(as.data.frame(climate), row.names = FALSE, digits = 3)
write_csv(climate, "output/pl2/scenario_comparison_climate.csv", append = FALSE)

## 2. Area by class ####

area <- map_dfr(SCEN$scenario, function(s) {
  tibble(
    scenario = s,
    class = RESPONSE_LEVELS,
    argmax_km2 = as.numeric(table(cls[[s]])[RESPONSE_LEVELS]) * cell_km2,
    expected_km2 = map_dbl(
      RESPONSE_LEVELS,
      ~ sum(P[, paste0(SCEN$prefix[SCEN$scenario == s], .x)]) * cell_km2
    )
  )
}) |>
  group_by(class) |>
  mutate(
    argmax_change_pct = 100 * (argmax_km2 / argmax_km2[scenario == "current"] - 1),
    expected_change_pct = 100 * (expected_km2 / expected_km2[scenario == "current"] - 1)
  ) |>
  ungroup()
cat("\nArea by class (km2) and change against current:\n")
print(as.data.frame(area), row.names = FALSE, digits = 4)
write_csv(area, "output/pl2/scenario_comparison_area.csv", append = FALSE)

## 3. Transitions ####

transitions <- map_dfr(c(MILD, HARSH), function(s) {
  tibble(from = cls[["current"]], to = cls[[s]]) |>
    count(from, to, name = "cells") |>
    mutate(scenario = s, km2 = cells * cell_km2, .before = 1)
}) |>
  arrange(scenario, desc(cells))
cat("\nArgmax transitions, current -> scenario:\n")
print(as.data.frame(transitions), row.names = FALSE, digits = 4)
write_csv(transitions, "output/pl2/scenario_comparison_transitions.csv", append = FALSE)

## 4. Persistence of present-day bog ####

persist_levels <- c(
  "bog under both", paste("bog under", MILD, "only"),
  paste("bog under", HARSH, "only"), "bog under neither"
)
persistence <- function(mild, harsh) {
  factor(
    case_when(
      mild & harsh ~ persist_levels[1],
      mild & !harsh ~ persist_levels[2],
      !mild & harsh ~ persist_levels[3],
      TRUE ~ persist_levels[4]
    ),
    levels = persist_levels
  )
}

# Area of applicability, per scenario, aligned to the cells above. A bog cell is only
# counted as reliably persistent if the model is inside its applicability area there.
aoa_of <- function(path) values(rast(path)[["aoa"]])[ok] == 1
aoa <- set_names(map(SCEN$reliability, aoa_of), SCEN$scenario)

cell_persist <- persistence(cls[[MILD]] == "bog", cls[[HARSH]] == "bog")
cells_tbl <- tibble(persist = cell_persist[bog_now], aoa_mild = aoa[[MILD]][bog_now],
                    aoa_harsh = aoa[[HARSH]][bog_now]) |>
  group_by(persist) |>
  summarise(
    cells = n(),
    km2 = cells * cell_km2,
    pct_inside_aoa_mild = 100 * mean(aoa_mild, na.rm = TRUE),
    pct_inside_aoa_harsh = 100 * mean(aoa_harsh, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(pct_of_current_bog = 100 * cells / sum(cells))
cat("\nPresent-day predicted-bog cells, by what each scenario does to them:\n")
print(as.data.frame(cells_tbl), row.names = FALSE, digits = 4)
write_csv(cells_tbl, "output/pl2/scenario_comparison_persistence_cells.csv", append = FALSE)

# The same split over the surveyed raised bogs themselves. Polygon means of the class
# probabilities, then argmax, exactly as pl2_interpret.R does for the primary scenario.
lyngstad <- st_read("data/DMraisedbog.gpkg", layer = "lyngstad-MTYPE_A", quiet = TRUE) |>
  st_transform(3035) |>
  vect()
zonal <- terra::extract(
  pred[[unlist(map(SCEN$prefix, ~ paste0(.x, RESPONSE_LEVELS)))]],
  lyngstad, fun = mean, na.rm = TRUE
) |>
  as_tibble()
zone_class <- function(prefix) {
  m <- as.matrix(zonal[, paste0(prefix, RESPONSE_LEVELS)])
  out <- rep(NA_character_, nrow(m))
  good <- complete.cases(m)
  out[good] <- RESPONSE_LEVELS[max.col(m[good, , drop = FALSE], ties.method = "first")]
  out
}
zc <- set_names(map(SCEN$prefix, zone_class), SCEN$scenario)

# Polygon-mean bio10 under each scenario, to say what the polygons that keep the class
# have in common (the stacks are the ones the model projected on).
zone_bio10 <- map(SCEN$stack, function(p) {
  terra::extract(rast(p)[[AXIS]], lyngstad, fun = mean, na.rm = TRUE)[[AXIS]]
})
names(zone_bio10) <- SCEN$scenario

poly_all <- tibble(
  ID = zonal$ID,
  class_current = zc[["current"]],
  class_mild = zc[[MILD]],
  class_harsh = zc[[HARSH]],
  pbog_current = zonal[[paste0(SCEN$prefix[1], "bog")]],
  pbog_mild = zonal[[paste0(SCEN$prefix[SCEN$scenario == MILD], "bog")]],
  pbog_harsh = zonal[[paste0(SCEN$prefix[SCEN$scenario == HARSH], "bog")]],
  bio10_current = zone_bio10[["current"]],
  bio10_mild = zone_bio10[[MILD]],
  bio10_harsh = zone_bio10[[HARSH]]
)
poly <- poly_all |>
  filter(!is.na(class_current), !is.na(class_mild), !is.na(class_harsh))

# One row per classified polygon, ID = row of the lyngstad-MTYPE_A layer (the order
# pl2_interpret.R binds its predictions in). ms/make_figures.R draws the scenario figure
# and the transition table from this file.
write_csv(poly, "output/pl2/scenario_comparison_polygons.csv", append = FALSE)

polygons_tbl <- poly |>
  mutate(
    persist = persistence(class_mild == "bog", class_harsh == "bog"),
    now = if_else(class_current == "bog", "bog today", "not bog today")
  ) |>
  group_by(now, persist) |>
  summarise(polygons = n(), mean_pbog_current = mean(pbog_current),
            mean_pbog_mild = mean(pbog_mild), mean_pbog_harsh = mean(pbog_harsh),
            .groups = "drop")
cat("\nLyngstad polygons (", nrow(poly), " classified):\n", sep = "")
print(as.data.frame(polygons_tbl), row.names = FALSE, digits = 3)
write_csv(polygons_tbl, "output/pl2/scenario_comparison_persistence_polygons.csv", append = FALSE)

poly |>
  summarise(
    n = n(),
    mean_change_mild = mean(pbog_mild - pbog_current),
    mean_change_harsh = mean(pbog_harsh - pbog_current),
    prop_declining_mild = mean(pbog_mild < pbog_current),
    prop_declining_harsh = mean(pbog_harsh < pbog_current)
  ) |>
  write_csv("output/pl2/scenario_comparison_polygon_change.csv", append = FALSE)

## Map ####

# Present-day predicted-bog cells coloured by persistence. At 250 m over Norway the cells
# are small, so this is a diagnostic; the manuscript figure is built by ms/make_figures.R.
codes <- rep(NA_integer_, length(ok))
codes_ok <- rep(NA_integer_, nrow(P))
codes_ok[bog_now] <- as.integer(cell_persist[bog_now])
codes[ok] <- codes_ok
persist_map <- rast(pred[[1]])
values(persist_map) <- codes
levels(persist_map) <- data.frame(id = seq_along(persist_levels), persistence = persist_levels)

# Colours are tied to the class id: only the classes that occur are drawn, and a bare
# `col` vector would then be matched by position rather than by class.
coltab(persist_map) <- data.frame(
  value = seq_along(persist_levels),
  col = c("#2b6a3f", "#e0a030", "#8a8a8a", "#c1462c")
)

png("output/pl2/scenario_comparison_persistence.png", width = 1100, height = 1100, res = 120)
plot(
  persist_map,
  main = "Present-day predicted bog: persistence under the two scenarios",
  mar = c(3, 3, 3, 14), plg = list(cex = 0.8)
)
dev.off()

cat("\nWrote output/pl2/scenario_comparison_*.csv and the persistence map\n")

# sessionInfo ####

sessioninfo::session_info()
