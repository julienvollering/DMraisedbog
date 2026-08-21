# Where the future goes, and what is observed there ####

# PURPOSE: Compares future bog cells against the three current classes along the most-shifting predictors, then reports observed bog prevalence in the climate space the future moves into.

# A model-free companion to the reliability machinery. `pl2_fitErrorProfiles.R` asks how
# well the CLASSIFIER performs at a given novelty; this asks what the DATA say about the
# climates the projection moves into. The second question turns out to carry more weight
# than the first (notebook 2026-08-21), because skill-versus-DI is flat mostly for a reason
# that has nothing to do with the warm edge: the CV's high-DI bog cells are cooler, not
# warmer, so that curve never probes the axis the projection travels.
#
# TWO PARTS.
#
#  1. Distribution comparison. Quantiles and box-whiskers of the three current classes plus
#     the future bog cells, over the most strongly shifting predictors. Answers "does the
#     future put existing bogs where other bogs are today, or where non-bogs are?"
#  2. Occupancy. Observed bog prevalence across a climate grid, with the future bog cells'
#     position on that same grid. Answers "and what do we actually observe there?"
#
# WHY PART 2 IS NOT BUILT ON THE TOP SHIFTING FEATURES. All 17 materially shifting
# predictors are temperature or degree-day measures -- npp, bio06, bio11, gdd*, ngd*, bio01,
# bio05, bio09, bio10, gsl, gst, scd, fcf -- and they are strongly collinear. A grid built
# from the top three would be nearly one-dimensional and would say only that the future is
# warmer, which is not in dispute. Every precipitation variable shifts less than 0.4 SD.
#
# That is exactly what makes moisture the right SECOND axis: it is approximately held
# constant between scenarios, so crossing it with temperature isolates the temperature
# effect instead of confounding it. Raised bogs are ombrotrophic, so precipitation surplus
# is the controlling variable and a warm-dry analogue would say little about a warm-wet
# future. The grid is therefore bio10 (highest-VI shifting predictor, and the partitioning
# axis) crossed with gsp (growing-season precipitation, highest-VI moisture predictor).
#
# The `dataset` split in the occupancy table is not decoration: the warm-and-wet evidence
# is almost entirely European, and a reader should see that rather than infer it.

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(sf)

## Configuration ####

RESPONSE_LEVELS <- c("nonpeat", "otherpeat", "bog")

# Predictors shown in part 1, taken as the most strongly shifting.
N_TOP_SHIFT <- 15

# Part 2 axes. Fixed, interpretable breaks rather than quantiles: this table is meant to be
# read directly, and the bands are chosen to bracket where the future cells land.
TEMP_VAR <- "bio10"
TEMP_BREAKS <- c(-Inf, 14, 16, 18, Inf)
MOIST_VAR <- "gsp"
MOIST_BREAKS <- c(-Inf, 500, 650, Inf)

## Inputs ####

mf <- read_csv("output/pl2/modeling_frame_regional.csv", show_col_types = FALSE)
weights <- read_csv(
  "output/pl2/weights_feature_data_partitioning.csv",
  show_col_types = FALSE
) |>
  filter(method == "Balanced Random Forest")

train <- mf |>
  filter(scenario == "current") |>
  mutate(response = factor(response, levels = RESPONSE_LEVELS))
future <- mf |> filter(scenario == "future")

# The future bog cells: the cells that are bog TODAY, carried forward under the scenario.
# This is a persistence comparison, so it must be the same cells, not a resample.
bog_xy <- train |> filter(response == "bog", dataset == "NO") |> select(x, y)
future_bog <- bog_xy |> inner_join(future, by = c("x", "y"))

cat("Current training rows:", nrow(train), "\n")
print(count(train, dataset, response))
cat("\nNorwegian bog cells carried to the future scenario:", nrow(future_bog), "\n\n")

## Part 1: where does the future put existing bogs? ####

top_shift <- weights |>
  filter(dynamic) |>
  mutate(abs_shift = abs(shift_sd)) |>
  arrange(desc(abs_shift)) |>
  head(N_TOP_SHIFT)

cat("Most strongly shifting predictors (all temperature-like -- see header):\n")
top_shift |>
  transmute(feature, VI = round(median, 4), shift_sd = round(shift_sd, 3)) |>
  as.data.frame() |>
  print(row.names = FALSE)

feats <- top_shift$feature

groups <- bind_rows(
  train |> select(all_of(feats), response) |> mutate(group = paste("current", response)),
  future_bog |> select(all_of(feats)) |> mutate(group = "future bog")
) |>
  select(-any_of("response")) |>
  mutate(group = factor(
    group,
    levels = c("current nonpeat", "current otherpeat", "current bog", "future bog")
  ))

long <- groups |>
  pivot_longer(all_of(feats), names_to = "feature", values_to = "value") |>
  mutate(feature = factor(feature, levels = feats))

quantiles_tbl <- long |>
  group_by(feature, group) |>
  summarise(
    n = n(),
    q05 = quantile(value, 0.05), q25 = quantile(value, 0.25),
    median = median(value),
    q75 = quantile(value, 0.75), q95 = quantile(value, 0.95),
    .groups = "drop"
  )

cat("\nQuantiles by group (first four predictors shown; full table written to CSV):\n")
quantiles_tbl |>
  filter(feature %in% feats[1:4]) |>
  mutate(across(q05:q95, ~ round(.x, 1))) |>
  as.data.frame() |>
  print(row.names = FALSE)

write_csv(quantiles_tbl, "output/pl2/occupancy_quantiles_by_group.csv", append = FALSE)

# Raw units, free y per facet: ecologically readable, and the only honest way to put
# degree-days and degrees Celsius on one figure.
p_box <- ggplot(long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  facet_wrap(~feature, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c(
    "current nonpeat" = "#b8b8b8", "current otherpeat" = "#7fa8c9",
    "current bog" = "#2b6a3f", "future bog" = "#c1462c"
  )) +
  labs(
    title = "Future bog cells against the three current classes",
    subtitle = paste("Top", N_TOP_SHIFT, "predictors by projected shift (all temperature-like)"),
    x = NULL, y = NULL, fill = NULL
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), legend.position = "bottom")
print(p_box)
ggsave("output/pl2/occupancy_boxplots.png", p_box, width = 12, height = 7.5, dpi = 150)

# Standardised on the current training distribution, so shift magnitudes are comparable
# ACROSS predictors on one panel -- the complement to the raw-unit facets above.
z <- long |>
  group_by(feature) |>
  mutate(z = (value - mean(value[group != "future bog"])) /
    sd(value[group != "future bog"])) |>
  ungroup()

p_z <- ggplot(z, aes(x = feature, y = z, fill = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
  coord_flip(ylim = quantile(z$z, c(0.005, 0.995))) +
  scale_fill_manual(values = c(
    "current nonpeat" = "#b8b8b8", "current otherpeat" = "#7fa8c9",
    "current bog" = "#2b6a3f", "future bog" = "#c1462c"
  )) +
  labs(
    title = "Same comparison, standardised on the current training distribution",
    x = NULL, y = "SD from current training mean", fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_z)
ggsave("output/pl2/occupancy_boxplots_z.png", p_z, width = 9, height = 8, dpi = 150)

## Part 2: observed bog prevalence in the climate space the future moves into ####

band <- function(d) {
  d |>
    mutate(
      temp_band = cut(.data[[TEMP_VAR]], TEMP_BREAKS,
        labels = c("<14", "14-16", "16-18", ">18")
      ),
      moist_band = cut(.data[[MOIST_VAR]], MOIST_BREAKS,
        labels = c("dry <500", "mid 500-650", "wet >650")
      )
    )
}

train_b <- band(train)
future_b <- band(future_bog)

occupancy <- train_b |>
  group_by(temp_band, moist_band) |>
  summarise(
    n_train = n(),
    n_bog = sum(response == "bog"),
    prevalence = n_bog / n_train,
    n_bog_EU = sum(response == "bog" & dataset == "EU"),
    n_train_EU = sum(dataset == "EU"),
    .groups = "drop"
  ) |>
  left_join(
    future_b |>
      count(temp_band, moist_band, name = "n_future_bog") |>
      mutate(pct_future_bog = 100 * n_future_bog / sum(n_future_bog)),
    by = c("temp_band", "moist_band")
  ) |>
  replace_na(list(n_future_bog = 0L, pct_future_bog = 0))

cat("\n=== Occupancy:", TEMP_VAR, "x", MOIST_VAR, "===\n")
occupancy |>
  mutate(prevalence = round(prevalence, 4), pct_future_bog = round(pct_future_bog, 1)) |>
  as.data.frame() |>
  print(row.names = FALSE)

write_csv(occupancy, "output/pl2/occupancy_grid.csv", append = FALSE)

# The headline reading: hold moisture constant in the wet band and walk up temperature.
cat("\nWithin the wet band, prevalence against temperature:\n")
occupancy |>
  filter(moist_band == "wet >650") |>
  transmute(
    temp_band, n_train, n_bog,
    prevalence = round(prevalence, 4),
    pct_future_bog_here = round(pct_future_bog, 1),
    pct_of_cells_from_EU = round(100 * n_train_EU / n_train, 1)
  ) |>
  as.data.frame() |>
  print(row.names = FALSE)

p_occ <- occupancy |>
  ggplot(aes(x = temp_band, y = moist_band)) +
  geom_tile(aes(fill = prevalence), colour = "white") +
  geom_text(
    aes(label = sprintf("%.2f%%\nn=%s\nbog=%d\nfuture %.0f%%",
      100 * prevalence, format(n_train, big.mark = ","), n_bog, pct_future_bog)),
    size = 2.9, lineheight = 0.95
  ) +
  scale_fill_gradient(low = "#f7f7f7", high = "#2b6a3f", labels = scales::percent) +
  labs(
    title = "Observed bog prevalence, and where the future bog cells land",
    subtitle = paste(
      "Rows: growing-season precipitation. Columns: warmest-quarter temperature.",
      "\n'future %' = share of the 1,116 current Norwegian bog cells falling here under the scenario"
    ),
    x = paste(TEMP_VAR, "(degrees C)"), y = paste(MOIST_VAR, "(mm)"), fill = "bog prevalence"
  ) +
  theme_minimal()
print(p_occ)
ggsave("output/pl2/occupancy_grid.png", p_occ, width = 9, height = 5.5, dpi = 150)

## How independent is the warm evidence, really? ####

# The occupancy table's authority comes from its counts, so the counts have to be honest
# about autocorrelation. At 250 m resolution, thousands of cells drawn from a few
# landscapes are not thousands of observations, and the warm band is the cell the whole
# argument leans on. This reports where those cells actually are.

warm <- train |> filter(.data[[TEMP_VAR]] > 18)

warm_ll <- warm |>
  st_as_sf(coords = c("x", "y"), crs = 3035) |>
  st_transform(4326) |>
  st_coordinates()

warm <- warm |>
  mutate(
    lon = warm_ll[, 1], lat = warm_ll[, 2],
    grid_lon = floor(lon), grid_lat = floor(lat)
  )

cells <- warm |>
  count(grid_lon, grid_lat, sort = TRUE) |>
  mutate(share = n / sum(n), cumshare = cumsum(share))

cat("\n=== Spatial concentration of the warm (", TEMP_VAR, "> 18 ) evidence ===\n")
cat("cells:", nrow(warm), "| bog:", sum(warm$response == "bog"),
    "| distinct 1-degree squares:", nrow(cells), "\n")
cat("longitude range:", round(range(warm$lon), 1),
    "| latitude range:", round(range(warm$lat), 1), "\n")
cat("share in the top 1 / 3 / 10 squares:",
    paste(round(100 * cells$cumshare[c(1, min(3, nrow(cells)), min(10, nrow(cells)))], 1),
          collapse = " / "), "%\n\n")
print(as.data.frame(head(cells, 8)), row.names = FALSE, digits = 3)

cat("\nElevation of the warm cells (m):\n")
print(round(quantile(warm$elevation, c(0, .25, .5, .75, 1))))

cat(paste0(
  "\nCAVEAT TO CARRY WITH THE OCCUPANCY TABLE: these cells come from a narrow southern\n",
  "arc (Massif Central, Jura, Alpine foothills, Bohemian belt) at 45-49 N, whereas the\n",
  "future Norwegian bogs sit at 58-66 N. Matching on ", TEMP_VAR, " and ", MOIST_VAR,
  " does not match\n",
  "photoperiod or continentality, and for Sphagnum systems growing-season light and\n",
  "evaporative demand are not incidental. Effective sample size is far below the raw n.\n"
))

p_warm <- warm |>
  ggplot(aes(x = lon, y = lat)) +
  geom_bin2d(bins = 40) +
  geom_point(
    data = warm |> filter(response == "bog"),
    aes(x = lon, y = lat), colour = "#c1462c", size = 2.5, shape = 4, stroke = 1.2
  ) +
  scale_fill_viridis_c(trans = "log10") +
  labs(
    title = paste("Where the warm evidence comes from:", TEMP_VAR, "> 18"),
    subtitle = paste0(
      nrow(warm), " training cells in ", nrow(cells),
      " one-degree squares; crosses = the ", sum(warm$response == "bog"), " bog cells.",
      "\nThe occupancy table's counts are not independent observations."
    ),
    x = "longitude", y = "latitude", fill = "cells"
  ) +
  theme_minimal()
print(p_warm)
ggsave("output/pl2/occupancy_warm_evidence_map.png", p_warm, width = 8, height = 6, dpi = 150)

write_csv(cells, "output/pl2/occupancy_warm_evidence_cells.csv", append = FALSE)

## Robustness: the same grid on annual precipitation ####

# bio12 has marginally higher VI than gsp and shifts even less, so it is the natural check
# that the pattern is about moisture rather than about this one moisture variable.
occ_bio12 <- train |>
  mutate(
    temp_band = cut(.data[[TEMP_VAR]], TEMP_BREAKS,
      labels = c("<14", "14-16", "16-18", ">18")
    ),
    moist_band = cut(bio12, c(-Inf, 700, 1100, Inf),
      labels = c("dry <700", "mid 700-1100", "wet >1100")
    )
  ) |>
  group_by(temp_band, moist_band) |>
  summarise(
    n_train = n(), n_bog = sum(response == "bog"),
    prevalence = round(n_bog / n_train, 4), .groups = "drop"
  )

cat("\n=== Robustness check: same grid using bio12 (annual precipitation) ===\n")
print(as.data.frame(occ_bio12), row.names = FALSE)
write_csv(occ_bio12, "output/pl2/occupancy_grid_bio12.csv", append = FALSE)

# sessionInfo ####

sessioninfo::session_info()
