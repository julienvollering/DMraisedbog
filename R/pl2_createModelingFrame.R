# Create the pooled 3-class modelling frame ####

# PURPOSE: Row-binds the Norway and EU blocks into the pooled 3-class training frame and writes the current/future scenario raster stacks.

# Plan section 5 step 6. Row-binds the two blocks built in pl0 into one training frame:
#
#   Norway  pl0_labelNorwayBlock.R  -> bog / otherpeat / nonpeat, dataset = "NO"
#   Europe  pl0_sampleEUabsences.R  -> otherpeat / nonpeat,       dataset = "EU"
#           Presences_EU            -> bog,                        dataset = "EU"
#
# THREE CHANGES FROM THE PREVIOUS VERSION, all of them deliberate:
#
# 1. `rf_global` is gone (plan section 1.2). The global model trains on the same 41-layer
#    stack at coarser grain, so the covariate added no predictor axis -- only EU label
#    information, through a channel redundant with bio10/gdd -- and it put the local
#    model's DI on two training sets at once, which makes the reliability map
#    uninterpretable. 43 predictors, not 44.
#
# 2. The `artype_60 >= 0.5 | presence` mask is gone. That mask deleted every cell warm
#    enough to hold no peat at all, which is precisely the warm flank the model needs in
#    order to learn where bogs stop on the warm side. artype_60 is now a LABEL source in
#    pl0_labelNorwayBlock.R, and it is dropped from the predictors here so the response
#    cannot leak into the features.
#
# 3. The response is 3-class and carries a `dataset` column. Downstream scripts that
#    tested `response == 1` / `response == 0` broke against it, which was the intent --
#    breaking loudly beats carrying a compatibility column that invites silent misuse.
#    All of pipeline 2 was converted on 2026-08-19 (plan section 5 step 8). `dataset` is
#    bookkeeping for the section 1.5 dataset-blocked comparison and must never enter a
#    model as a predictor; every downstream script excludes it explicitly.
#
# ONE ASYMMETRY TO STATE RATHER THAN HIDE: Norway's climate comes from the 250 m regional
# stack, Europe's from the global 5 km stack, because no 250 m EU climate stack exists.
# CHELSA is 30 arcsec native, so both are resamples of the same source and neither is
# "the" native resolution -- but the EU rows are the coarser of the two, and any feature
# that lives on fine topographic gradients is smoother on the EU side. Terrain is not
# affected: pl0_buildEUterrain.R builds EU elevation and slope at 250 m by the same
# aggregate-then-slope recipe used for DTM50.

library(tidyverse)
library(terra)
library(sf)

source("R/functions.R")

dir.create("output/pl2", showWarnings = FALSE, recursive = TRUE)

seed <- 42
set.seed(seed)

# Future rows exist for the DI and projected-shift diagnostics, not for fitting. They must
# cover the block's own coordinates exactly, because pl2_exploreOccupancy.R joins future
# conditions onto the Norwegian bog cells by (x, y); a purely random sample would miss
# them and silently produce NAs.
n_future_sample <- 200000

## Predictors ####

preds_no_cur <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")
preds_no_fut <- rast("output/predictors_regional_250m_Norway_future_EPSG3035.tif")

# artype_60 is the Norway label source; keeping it as a predictor would leak the response.
preds_no_cur <- preds_no_cur[[!grepl("artype", names(preds_no_cur))]]
preds_no_fut <- preds_no_fut[[!grepl("artype", names(preds_no_fut))]]

feat <- names(preds_no_cur)
cat("Predictors:", length(feat), "\n")

writeRaster(preds_no_cur, "output/pl2/scenario_current.tif", overwrite = TRUE)
writeRaster(preds_no_fut, "output/pl2/scenario_future.tif", overwrite = TRUE)

## Norway block ####

no_block <- read_csv("output/pl0/no_block_coords.csv", show_col_types = FALSE)

df_no <- terra::extract(preds_no_cur, no_block[c("x", "y")], ID = FALSE) |>
  bind_cols(no_block |> select(x, y, response, dataset)) |>
  drop_na() |>
  as_tibble()

cat("Norway rows:", nrow(df_no), "of", nrow(no_block), "after dropping incomplete\n")

## EU block ####

# Presences: one row per 250 m cell. Several EUNIS Q11 relevés can fall inside the same
# cell of the same Natura polygon, and they are the same observation as far as the model
# is concerned.
eu_presence <- st_read(
  "data/DMraisedbog.gpkg",
  layer = "Presences_EU", quiet = TRUE
) |>
  st_geometry() |>
  st_point_on_surface() |>
  st_coordinates() |>
  as_tibble() |>
  rename(x = X, y = Y) |>
  mutate(response = "bog", dataset = "EU")

grid_250m <- rast("output/pl0/epm_map_cat_250m.tif")
eu_presence <- eu_presence |>
  mutate(cell = cellFromXY(grid_250m, cbind(x, y))) |>
  filter(!is.na(cell)) |>
  distinct(cell, .keep_all = TRUE) |>
  select(-cell)

# The same land-use / water screen the absence pools get, on the same thresholds. It runs
# here for the EU presence ROWS and in pl0_sampleEUabsences.R for the presence COUNT that
# sets the absence budget; both call landuse_screen_keep() over the same cells, so they
# agree by construction. Norway's presences are screened upstream in
# pl0_labelNorwayBlock.R and arrive already filtered.
eu_presence <- eu_presence[
  landuse_screen_keep(eu_presence, label = "EU presence rows"),
]

eu_absence <- read_csv("output/pl0/eu_absence_coords.csv", show_col_types = FALSE) |>
  select(x, y, response, dataset)

eu_block <- bind_rows(eu_presence, eu_absence)

cat(
  "EU rows: ", nrow(eu_presence), " presences (from 701 geometries) + ",
  nrow(eu_absence), " absences\n",
  sep = ""
)

# 41 climate and paleo features from the global stack, terrain from the EU 250 m build.
global_5km <- rast("output/predictors_global_5km_EUNorway_EPSG3035.tif")
eu_terrain <- rast("output/pl0/eu_terrain_250m.tif")

df_eu <- bind_cols(
  eu_block |> select(x, y, response, dataset),
  terra::extract(global_5km, eu_block[c("x", "y")], ID = FALSE),
  terra::extract(eu_terrain, eu_block[c("x", "y")], ID = FALSE)
) |>
  as_tibble()

# Guard: the two blocks must carry exactly the same predictor set, or the row-bind below
# silently fills one block's missing columns with NA.
missing_eu <- setdiff(feat, names(df_eu))
extra_eu <- setdiff(setdiff(names(df_eu), c("x", "y", "response", "dataset")), feat)
if (length(missing_eu) > 0 || length(extra_eu) > 0) {
  stop(
    "EU predictor mismatch. Missing: ", paste(missing_eu, collapse = ", "),
    " | unexpected: ", paste(extra_eu, collapse = ", ")
  )
}

df_eu <- df_eu |>
  select(x, y, response, dataset, all_of(feat)) |>
  drop_na()

cat("EU rows:", nrow(df_eu), "of", nrow(eu_block), "after dropping incomplete\n")

## Current scenario ####

response_levels <- c("nonpeat", "otherpeat", "bog")

current <- bind_rows(df_no, df_eu) |>
  mutate(
    response = factor(response, levels = response_levels),
    dataset = factor(dataset, levels = c("NO", "EU"))
  ) |>
  select(response, dataset, all_of(feat), x, y)

stopifnot(!any(is.na(current$response)))

cat("\nPooled current frame:\n")
print(as.data.frame(count(current, dataset, response)))

## Future scenario ####

# Block coordinates first, so the (x, y) join downstream is exact; then a random sample of
# the wider projection domain, so the DI histogram still describes where Norway is going
# rather than only where the block happens to have drawn.
fut_block <- terra::extract(
  preds_no_fut, no_block[c("x", "y")],
  xy = TRUE, ID = FALSE
)
# Sample cells off ONE layer, then extract the stack at those coordinates. Sampling the
# 43-layer stack directly makes terra carry every layer through the rejection sampling
# that na.rm implies, which is what turns this step from seconds into many minutes.
fut_xy <- spatSample(
  preds_no_fut[[1]],
  size = n_future_sample, method = "random", na.rm = TRUE, xy = TRUE
) |>
  select(x, y)
fut_sample <- terra::extract(preds_no_fut, fut_xy, xy = TRUE, ID = FALSE)

future <- bind_rows(as_tibble(fut_block), as_tibble(fut_sample)) |>
  distinct(x, y, .keep_all = TRUE) |>
  drop_na() |>
  mutate(response = NA_character_, dataset = "NO") |>
  mutate(
    response = factor(response, levels = response_levels),
    dataset = factor(dataset, levels = c("NO", "EU"))
  ) |>
  select(response, dataset, all_of(feat), x, y)

cat("Future rows:", nrow(future),
    "( block coords +", n_future_sample, "sampled, deduplicated )\n")

rm(fut_block, fut_sample)
gc()

## Combine and save ####

mf <- bind_rows(current = current, future = future, .id = "scenario") |>
  mutate(scenario = factor(scenario, levels = c("current", "future")))

rm(current, future)
gc()

print(count(mf, scenario, dataset, response) |> as.data.frame())

frame_summary <- mf |>
  count(scenario, dataset, response, name = "rows")
write_csv(frame_summary, "output/pl2/modeling_frame_summary.csv", append = FALSE)

## Predictor ranges by block and scenario ####

# Written every run so a corrupt or mis-scaled layer is visible in the first table a
# reviewer opens rather than three scripts later. The gsp no-data sentinel (4.29e8 mm in
# 8,987 training rows, notebook 2026-09-12) would have shown here as a q99 in the hundreds
# of millions on the first run. The check warns when a predictor's Norwegian and European
# CURRENT ranges do not overlap at all: the two blocks share a climate source, so disjoint
# ranges mean a units or scaling fault in one of them, never a real contrast.
predictor_ranges <- mf |>
  group_by(scenario, dataset) |>
  summarise(
    n = n(),
    across(
      all_of(feat),
      list(
        min = ~ min(.x), q01 = ~ unname(quantile(.x, 0.01)), median = ~ median(.x),
        q99 = ~ unname(quantile(.x, 0.99)), max = ~ max(.x)
      ),
      .names = "{.col}__{.fn}"
    ),
    .groups = "drop"
  ) |>
  pivot_longer(
    -c(scenario, dataset, n),
    names_to = c("predictor", "stat"), names_sep = "__"
  ) |>
  pivot_wider(names_from = stat, values_from = value) |>
  arrange(predictor, scenario, dataset)
write_csv(
  predictor_ranges, "output/pl2/modeling_frame_predictor_ranges.csv", append = FALSE
)

cat("\nPredictor ranges, current scenario, by block:\n")
predictor_ranges |>
  filter(scenario == "current") |>
  transmute(
    predictor, dataset,
    min = signif(min, 4), median = signif(median, 4), max = signif(max, 4)
  ) |>
  as.data.frame() |>
  print(row.names = FALSE)

disjoint <- predictor_ranges |>
  filter(scenario == "current") |>
  select(predictor, dataset, min, max) |>
  pivot_wider(names_from = dataset, values_from = c(min, max)) |>
  filter(max_NO < min_EU | max_EU < min_NO)
if (nrow(disjoint) > 0) {
  warning(
    "Norwegian and European current ranges do not overlap for: ",
    paste(disjoint$predictor, collapse = ", ")
  )
} else {
  cat("\nNorwegian and European current ranges overlap for every predictor\n")
}

write_csv(mf, "output/pl2/modeling_frame_regional.csv", append = FALSE)

cat(
  "\nWrote output/pl2/modeling_frame_regional.csv (", nrow(mf), "rows,",
  length(feat), "predictors )\n"
)

# sessionInfo ####

sessioninfo::session_info()
