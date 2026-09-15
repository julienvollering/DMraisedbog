# Build the manuscript figures and tables ####

# PURPOSE: Draws the static figures and writes the tables used in ms/ms.qmd, from pipeline
# outputs only. Fits nothing and recomputes nothing that the pipeline already reports.
#
# Run from the repository root:
#   Rscript --vanilla ms/make_figures.R
#
# Inputs: output/pl0, output/pl2 (see R/RUNALL.R) and data/DMraisedbog.gpkg.
# Outputs: ms/figures/*.png, ms/tables/*.csv.
#
# It lives under ms/ rather than R/ on purpose: it is part of the manuscript, not of the
# analysis, and R/RUNALL.R's guard would otherwise report it as an unwired script.

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(sf)
library(terra)

FIG_DIR <- "ms/figures"
TBL_DIR <- "ms/tables"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

GPKG <- "data/DMraisedbog.gpkg"
CRS <- 3035

# Aggregation factor for map rasters: 250 m -> 2 km is plenty for a page-width panel and
# keeps the data frames behind geom_raster to ~100k rows.
AGG <- 8

CLASS_LAB <- c(nonpeat = "non-peat", otherpeat = "other peatland", bog = "raised bog")
CLASS_COL <- c(
  "non-peat" = "#b8b8b8",
  "other peatland" = "#7fa8c9",
  "raised bog" = "#2b6a3f",
  "raised bog, future climate" = "#c1462c"
)

theme_ms <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 10),
    legend.key.size = unit(0.8, "lines")
  )
theme_map <- theme_ms +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

save_fig <- function(p, name, w, h) {
  ggsave(
    file.path(FIG_DIR, paste0(name, ".png")), p,
    width = w, height = h, dpi = 300, bg = "white"
  )
  cat("wrote", file.path(FIG_DIR, paste0(name, ".png")), "\n")
}
save_tbl <- function(d, name) {
  write_csv(d, file.path(TBL_DIR, paste0(name, ".csv")), append = FALSE)
  cat("wrote", file.path(TBL_DIR, paste0(name, ".csv")), "\n")
}

numbers <- list()

## Spatial layers ####

countries <- st_read(
  "data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp", quiet = TRUE
) |>
  st_transform(CRS)
norway <- countries |> filter(ADM0_A3 == "NOR")

domain <- st_read(GPKG, layer = "EU_domain", quiet = TRUE)
footprint <- st_read(GPKG, layer = "nib-lyngstad-footprint", quiet = TRUE) |>
  st_transform(CRS)
projects <- st_read(GPKG, layer = "nib-lyngstad-projects", quiet = TRUE)
bogs <- st_read(GPKG, layer = "lyngstad-MTYPE_A", quiet = TRUE) |>
  st_transform(CRS)
bog_pts <- st_centroid(st_geometry(bogs))
eu_pres <- st_read(GPKG, layer = "Presences_EU", quiet = TRUE) |>
  st_point_on_surface()

no_grid <- rast("output/ar50_250m_land_EPSG3035.tif")
no_bbox <- st_as_sfc(st_bbox(no_grid))

numbers$n_lyngstad_polygons <- nrow(bogs)
numbers$lyngstad_area_km2 <- as.numeric(sum(st_area(bogs))) / 1e6
numbers$footprint_km2 <- as.numeric(sum(st_area(footprint))) / 1e6
numbers$n_survey_projects <- nrow(projects)
numbers$n_eu_presence_geometries <- nrow(eu_pres)
numbers$eu_domain_km2 <- as.numeric(sum(st_area(domain))) / 1e6

## Fig 1: study area ####

bb_a <- st_bbox(st_union(st_union(st_geometry(domain)), no_bbox))
pad <- 150000
p1a <- ggplot() +
  geom_sf(data = countries, fill = "grey96", colour = "grey70", linewidth = 0.15) +
  geom_sf(data = domain, fill = "#dbe7d8", colour = "#2b6a3f", linewidth = 0.25) +
  geom_sf(data = footprint, fill = "#f6dccb", colour = "#c1462c", linewidth = 0.25) +
  geom_sf(data = eu_pres, colour = "#2b6a3f", size = 0.35) +
  geom_sf(data = bog_pts, colour = "#c1462c", size = 0.35) +
  coord_sf(
    xlim = c(bb_a["xmin"] - pad, bb_a["xmax"] + pad),
    ylim = c(bb_a["ymin"] - pad, bb_a["ymax"] + pad),
    expand = FALSE
  ) +
  labs(title = "A") +
  theme_map

bio10_cur <- rast("output/pl2/scenario_current.tif")[["bio10"]] |>
  aggregate(fact = AGG, fun = "mean", na.rm = TRUE)
bio10_df <- as.data.frame(bio10_cur, xy = TRUE)

p1b <- ggplot() +
  geom_raster(data = bio10_df, aes(x, y, fill = bio10)) +
  scale_fill_viridis_c(option = "magma", name = "bio10 (°C)\n1981-2010") +
  geom_sf(data = footprint, fill = NA, colour = "black", linewidth = 0.3) +
  geom_sf(data = bog_pts, colour = "#c1462c", size = 0.5) +
  coord_sf(expand = FALSE) +
  labs(title = "B") +
  theme_map +
  theme(legend.position = "right")

save_fig(p1a + p1b + plot_layout(widths = c(1.1, 1)), "fig_studyarea", 7.5, 4.6)

## Fig 2: climate space and occupancy ####

mf <- read_csv(
  "output/pl2/modeling_frame_regional.csv",
  col_select = c("scenario", "response", "dataset", "x", "y", "bio10", "gsp"),
  show_col_types = FALSE
)
train <- mf |> filter(scenario == "current")
future_xy <- mf |> filter(scenario == "future") |> select(x, y, bio10, gsp)

# DATA-QUALITY FLAG (found 2026-09-12 while drawing this figure). CHELSA's gsp layer
# stores no-data as 4294967295 (x 0.1 on read = 4.29e8) where the growing season has
# zero length, and pl0_collatePredictors.R only converts NA -> 0 for gdd10, gdd5, gst and
# swe -- gsp is not on that list, so the code is carried into the frame as a real value
# (and bilinear resampling smears it into intermediate values). Affected rows are cold,
# high-elevation non-peat cells (bio10 <= 16.8, median elevation ~1,800 m); no Norwegian
# bog row is affected. They are dropped from the scatter only; the pipeline outputs
# themselves are reported as they stand, with a TODO in the manuscript.
GSP_MAX_PLAUSIBLE <- 5000
numbers$gsp_nodata_rows_training <- sum(train$gsp > GSP_MAX_PLAUSIBLE)
numbers$gsp_nodata_rows_training_bog <- sum(train$gsp > GSP_MAX_PLAUSIBLE & train$response == "bog")
numbers$gsp_nodata_rows_future <- sum(future_xy$gsp > GSP_MAX_PLAUSIBLE)
numbers$gsp_nodata_max_bio10 <- max(train$bio10[train$gsp > GSP_MAX_PLAUSIBLE])
future_bog <- train |>
  filter(response == "bog", dataset == "NO") |>
  select(x, y) |>
  inner_join(future_xy, by = c("x", "y"))

set.seed(42)
bg <- train |>
  filter(response != "bog") |>
  group_by(response) |>
  slice_sample(n = 15000) |>
  ungroup() |>
  mutate(group = unname(CLASS_LAB[response]))

pts <- bind_rows(
  bg,
  train |> filter(response == "bog") |> mutate(group = "raised bog"),
  future_bog |> mutate(group = "raised bog, future climate")
) |>
  mutate(group = factor(group, levels = names(CLASS_COL))) |>
  filter(gsp <= GSP_MAX_PLAUSIBLE)

p2a <- ggplot(pts, aes(bio10, gsp, colour = group)) +
  geom_point(
    data = ~ filter(.x, group %in% c("non-peat", "other peatland")),
    size = 0.25, alpha = 0.25
  ) +
  geom_point(
    data = ~ filter(.x, !group %in% c("non-peat", "other peatland")),
    size = 0.45, alpha = 0.7
  ) +
  geom_vline(xintercept = c(14, 16, 18), linetype = 3, colour = "grey40") +
  geom_hline(yintercept = c(500, 650), linetype = 3, colour = "grey40") +
  scale_colour_manual(values = CLASS_COL, name = NULL) +
  coord_cartesian(ylim = c(0, quantile(pts$gsp, 0.999))) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 2)) +
  labs(
    title = "A",
    x = "Mean temperature of the warmest quarter, bio10 (°C)",
    y = "Growing-season precipitation, gsp (mm)"
  ) +
  theme_ms +
  theme(legend.position = "bottom")

occ <- read_csv("output/pl2/occupancy_grid.csv", show_col_types = FALSE) |>
  mutate(
    temp_band = factor(temp_band, levels = c("<14", "14-16", "16-18", ">18")),
    moist_band = factor(
      moist_band,
      levels = c("dry <500", "mid 500-650", "wet >650"),
      labels = c("< 500", "500-650", "> 650")
    )
  )

p2b <- ggplot(occ, aes(temp_band, moist_band)) +
  geom_tile(aes(fill = prevalence), colour = "white", linewidth = 0.8) +
  geom_text(
    aes(label = sprintf(
      "%.2f%% bog\n%d of %s\nfuture bogs: %.0f%%",
      100 * prevalence, n_bog, format(n_train, big.mark = ","), pct_future_bog
    )),
    size = 2.4, lineheight = 0.9
  ) +
  scale_fill_gradient(
    low = "#f7f7f7", high = "#2b6a3f", labels = scales::percent,
    name = "observed bog prevalence"
  ) +
  guides(fill = guide_colourbar(title.position = "top", barwidth = unit(4, "cm"))) +
  labs(
    title = "B",
    x = "Mean temperature of the warmest quarter, bio10 (°C)",
    y = "Growing-season precipitation, gsp (mm)"
  ) +
  theme_ms +
  theme(legend.position = "bottom")

save_fig(p2a + p2b + plot_layout(widths = c(1.15, 1)), "fig_occupancy", 7.5, 4.4)

save_tbl(
  occ |>
    transmute(
      bio10_band = temp_band, gsp_band = moist_band,
      n_cells = n_train, n_bog = n_bog,
      prevalence_pct = round(100 * prevalence, 2),
      share_eu_pct = round(100 * n_train_EU / n_train, 1),
      n_future_bog = n_future_bog,
      future_bog_pct = round(pct_future_bog, 1)
    ),
  "tbl_occupancy"
)

## Fig 3 and Table 4: skill against warming beyond the training bogs ####

ep <- readRDS("output/pl2/error_profiles.rds")
ci <- ep$profiles_ci
axis_ref <- ep$axis_ref_production
bog_edge <- max(train$bio10[train$response == "bog"]) - axis_ref

rel_fut <- rast("output/pl2/reliability_future.tif")
rel_cur <- rast("output/pl2/reliability_current.tif")
# Guard against reading a raster that pl2_mapReliability.R is still writing.
stopifnot(
  as.numeric(global(rel_fut[["offset"]], "notNA")) > 0,
  as.numeric(global(rel_cur[["offset"]], "notNA")) > 0
)
lyng_off_fut <- terra::extract(rel_fut[["offset"]], vect(bogs), fun = mean, na.rm = TRUE)$offset
lyng_off_cur <- terra::extract(rel_cur[["offset"]], vect(bogs), fun = mean, na.rm = TRUE)$offset

numbers$axis_ref_bio10 <- axis_ref
numbers$bog_edge_offset <- bog_edge
numbers$lyngstad_offset_future_q10 <- quantile(lyng_off_fut, 0.10, na.rm = TRUE)
numbers$lyngstad_offset_future_median <- median(lyng_off_fut, na.rm = TRUE)
numbers$lyngstad_offset_future_q90 <- quantile(lyng_off_fut, 0.90, na.rm = TRUE)
numbers$lyngstad_offset_current_median <- median(lyng_off_cur, na.rm = TRUE)
numbers$lyngstad_share_beyond_bog_edge <- mean(lyng_off_fut > bog_edge, na.rm = TRUE)
numbers$aoa_threshold <- ep$aoa_threshold

metric_lab <- c(
  recall_bog = "recall, raised bog",
  fpr_bog = "false-positive rate, raised bog",
  auc_bog = "AUC, raised bog vs rest",
  macro_auc = "macro AUC"
)

rug_df <- bind_rows(
  tibble(offset = lyng_off_cur, when = "current climate"),
  tibble(offset = lyng_off_fut, when = "2071-2100")
)

# The reference arm (leave-one-partition-out, production-sized training sets), point
# estimates only: drawn as thin dashed lines so the gap between the small-model curve the
# map reads and a production-sized model stays visible.
stopifnot(ep$cv_arm == "pairwise", ep$compare_arm == "lopo")
compare_long <- ep$profiles_compare |>
  select(axis_median, all_of(names(metric_lab))) |>
  pivot_longer(-axis_median, names_to = "metric", values_to = "estimate")

offset_panel <- function(metrics, title, ylab, ylim = c(0, 1)) {
  d <- ci |>
    filter(metric %in% metrics) |>
    mutate(metric = factor(metric_lab[metric], levels = metric_lab[metrics]))
  d_cmp <- compare_long |>
    filter(metric %in% metrics) |>
    mutate(metric = factor(metric_lab[metric], levels = metric_lab[metrics]))
  ggplot(d, aes(axis_median, estimate, colour = metric, fill = metric)) +
    annotate("rect",
      xmin = quantile(lyng_off_fut, 0.1, na.rm = TRUE),
      xmax = quantile(lyng_off_fut, 0.9, na.rm = TRUE),
      ymin = -Inf, ymax = Inf, fill = "#c1462c", alpha = 0.12
    ) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.6) +
    geom_line(data = d_cmp, linetype = 2, linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey30") +
    geom_vline(xintercept = bog_edge, linetype = 2, colour = "grey30") +
    geom_rug(
      data = rug_df, aes(x = offset, colour = NULL, fill = NULL),
      sides = "b", alpha = 0.3, length = unit(0.04, "npc"), colour = "#c1462c",
      inherit.aes = FALSE
    ) +
    scale_colour_brewer(palette = "Dark2", name = NULL) +
    scale_fill_brewer(palette = "Dark2", name = NULL) +
    coord_cartesian(ylim = ylim, xlim = c(-6, 6)) +
    labs(
      title = title,
      x = "bio10 offset from the typical training raised bog (°C)",
      y = ylab
    ) +
    theme_ms +
    theme(legend.position = "bottom")
}

p3a <- offset_panel(c("recall_bog", "fpr_bog"), "A", "rate at the argmax operating point")
p3b <- offset_panel(c("auc_bog", "macro_auc"), "B", "AUC", ylim = c(0.5, 1))
save_fig(p3a + p3b, "fig_skill_offset", 7.5, 3.8)

fmt_ci <- function(est, lo, hi) {
  ifelse(
    is.na(lo), sprintf("%.2f", est),
    sprintf("%.2f (%.2f-%.2f)", est, lo, hi)
  )
}
tbl_bins <- ci |>
  filter(metric %in% c("recall_bog", "fpr_bog", "auc_bog", "macro_auc")) |>
  mutate(cell = fmt_ci(estimate, lower, upper)) |>
  select(bin, n, n_bog, n_folds, axis_lo, axis_hi, axis_median, metric, cell) |>
  pivot_wider(names_from = metric, values_from = cell) |>
  mutate(
    offset_range = sprintf("%+.1f to %+.1f", axis_lo, axis_hi),
    offset_median = sprintf("%+.2f", axis_median)
  ) |>
  select(bin, offset_range, offset_median, n, n_bog, n_folds,
         recall_bog, fpr_bog, auc_bog, macro_auc)
save_tbl(tbl_bins, "tbl_offset_bins")

## S1 Fig: the retained DI axis ####

ci_di <- ep$profiles_ci_di |>
  filter(metric %in% c("recall_bog", "fpr_bog", "auc_bog", "macro_auc")) |>
  mutate(metric = factor(metric_lab[metric], levels = metric_lab))

di_fut <- values(rel_fut[["DI"]], na.rm = TRUE)[, 1]
di_fut_bog <- terra::extract(rel_fut[["DI"]], vect(bogs), fun = mean, na.rm = TRUE)$DI

ps1 <- ggplot(ci_di, aes(axis_median, estimate, colour = metric, fill = metric)) +
  annotate("rect",
    xmin = quantile(di_fut, 0.1), xmax = quantile(di_fut, 0.9),
    ymin = -Inf, ymax = Inf, fill = "grey60", alpha = 0.15
  ) +
  annotate("rect",
    xmin = quantile(di_fut_bog, 0.1, na.rm = TRUE),
    xmax = quantile(di_fut_bog, 0.9, na.rm = TRUE),
    ymin = -Inf, ymax = Inf, fill = "#c1462c", alpha = 0.15
  ) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  geom_vline(xintercept = ep$aoa_threshold, linetype = 2, colour = "grey30") +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  scale_fill_brewer(palette = "Dark2", name = NULL) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "dissimilarity index (log scale)",
    y = "rate or AUC"
  ) +
  theme_ms +
  theme(legend.position = "bottom")
save_fig(ps1, "figS_skill_di", 5.5, 3.8)

## Fig 4: projected class probabilities over Norway ####

pred <- rast("output/pl2/rf_local_pred_brf.tif") |>
  aggregate(fact = AGG, fun = "mean", na.rm = TRUE)
pred_df <- as.data.frame(pred, xy = TRUE) |>
  mutate(
    current_class = c("nonpeat", "otherpeat", "bog")[
      max.col(cbind(current_nonpeat, current_otherpeat, current_bog), ties.method = "first")
    ],
    future_class = c("nonpeat", "otherpeat", "bog")[
      max.col(cbind(future_nonpeat, future_otherpeat, future_bog), ties.method = "first")
    ],
    current_class = factor(CLASS_LAB[current_class], levels = CLASS_LAB),
    future_class = factor(CLASS_LAB[future_class], levels = CLASS_LAB)
  )

prob_map <- function(col, title) {
  ggplot() +
    geom_raster(data = pred_df, aes(x, y, fill = .data[[col]])) +
    scale_fill_viridis_c(limits = c(0, 1), name = "P(raised bog)") +
    geom_sf(data = footprint, fill = NA, colour = "white", linewidth = 0.2) +
    coord_sf(expand = FALSE) +
    labs(title = title) +
    theme_map
}
class_map <- function(col, title) {
  ggplot() +
    geom_raster(data = pred_df, aes(x, y, fill = .data[[col]])) +
    scale_fill_manual(values = CLASS_COL, name = "majority class", drop = FALSE) +
    geom_sf(data = footprint, fill = NA, colour = "black", linewidth = 0.2) +
    coord_sf(expand = FALSE) +
    labs(title = title) +
    theme_map
}

p4 <- (prob_map("current_bog", "A  1981-2010") | prob_map("future_bog", "B  2071-2100")) /
  (class_map("current_class", "C  1981-2010") | class_map("future_class", "D  2071-2100")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
save_fig(p4, "fig_projection", 7.5, 8.2)

# Share of the projection domain in each class, current and future, for the text.
class_share <- pred_df |>
  summarise(
    across(c(current_nonpeat, current_otherpeat, current_bog,
             future_nonpeat, future_otherpeat, future_bog), ~ mean(.x, na.rm = TRUE))
  ) |>
  pivot_longer(everything(), names_to = c("scenario", "class"), names_sep = "_",
               values_to = "mean_probability") |>
  left_join(
    bind_rows(
      pred_df |> count(scenario = "current", class = current_class, name = "n_cells"),
      pred_df |> count(scenario = "future", class = future_class, name = "n_cells")
    ) |>
      mutate(class = names(CLASS_LAB)[match(class, CLASS_LAB)]) |>
      group_by(scenario) |>
      mutate(majority_share = n_cells / sum(n_cells)) |>
      ungroup() |>
      select(-n_cells),
    by = c("scenario", "class")
  ) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))
save_tbl(class_share, "tbl_domain_class_share")

## Fig 5 and Table 5: change at the Lyngstad raised bogs ####

lyng_pred <- st_read("output/pl2/lyngstad_predictions.gpkg", layer = "brf", quiet = TRUE)
lyng_pts <- lyng_pred |>
  mutate(geom = st_centroid(st_geometry(lyng_pred))) |>
  st_set_geometry("geom")

p5a <- ggplot(lyng_pred, aes(prob_bog_current, prob_bog_future)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey40") +
  geom_point(size = 0.9, alpha = 0.6, colour = "#2b6a3f") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "A",
    x = "P(raised bog), 1981-2010",
    y = "P(raised bog), 2071-2100"
  ) +
  theme_ms

norway_crop <- st_crop(norway, st_bbox(no_grid))
p5b <- ggplot() +
  geom_sf(data = norway_crop, fill = "grey95", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = footprint, fill = NA, colour = "grey30", linewidth = 0.2) +
  geom_sf(data = lyng_pts, aes(colour = prob_bog_future), size = 0.8) +
  scale_colour_viridis_c(limits = c(0, 1), name = "P(raised bog)\n2071-2100") +
  coord_sf(
    xlim = c(st_bbox(footprint)["xmin"] - 50000, st_bbox(footprint)["xmax"] + 50000),
    ylim = c(st_bbox(footprint)["ymin"] - 50000, st_bbox(footprint)["ymax"] + 50000),
    expand = FALSE
  ) +
  labs(title = "B") +
  theme_map
save_fig(p5a + p5b + plot_layout(widths = c(1, 1.1)), "fig_lyngstad_change", 7.5, 4.2)

tbl_trans <- lyng_pred |>
  st_drop_geometry() |>
  filter(!is.na(transition)) |>
  count(class_current, class_future, name = "n_polygons") |>
  mutate(
    class_current = CLASS_LAB[class_current],
    class_future = CLASS_LAB[class_future],
    share_pct = round(100 * n_polygons / sum(n_polygons), 1)
  ) |>
  arrange(desc(n_polygons))
save_tbl(tbl_trans, "tbl_transitions")

pb <- lyng_pred |> st_drop_geometry()
numbers$lyngstad_pbog_current_mean <- mean(pb$prob_bog_current, na.rm = TRUE)
numbers$lyngstad_pbog_future_mean <- mean(pb$prob_bog_future, na.rm = TRUE)
numbers$lyngstad_pbog_change_median <- median(pb$change_bog, na.rm = TRUE)
numbers$lyngstad_pbog_change_q05 <- quantile(pb$change_bog, 0.05, na.rm = TRUE)
numbers$lyngstad_pbog_change_q95 <- quantile(pb$change_bog, 0.95, na.rm = TRUE)
numbers$lyngstad_share_declining <- mean(pb$change_bog < 0, na.rm = TRUE)
numbers$lyngstad_future_pbog_above_0.5 <- mean(pb$prob_bog_future > 0.5, na.rm = TRUE)

## Fig 6: reliability layers over the future domain ####

rel_agg <- aggregate(rel_fut[[c("offset", "DI", "exp_recall_bog", "exp_fpr_bog")]],
                     fact = AGG, fun = "mean", na.rm = TRUE)
rel_df <- as.data.frame(rel_agg, xy = TRUE)

p6a <- ggplot() +
  geom_raster(data = rel_df, aes(x, y, fill = offset)) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "grey95", high = "#b2182b", midpoint = 0,
    name = "bio10 offset (°C)"
  ) +
  geom_sf(data = footprint, fill = NA, colour = "black", linewidth = 0.2) +
  coord_sf(expand = FALSE) +
  labs(title = "A") +
  theme_map
p6b <- ggplot() +
  geom_raster(data = rel_df, aes(x, y, fill = DI)) +
  scale_fill_viridis_c(option = "cividis", name = "dissimilarity\nindex") +
  geom_sf(data = footprint, fill = NA, colour = "white", linewidth = 0.2) +
  coord_sf(expand = FALSE) +
  labs(title = "B") +
  theme_map
p6c <- ggplot() +
  geom_raster(data = rel_df, aes(x, y, fill = exp_recall_bog)) +
  scale_fill_viridis_c(limits = c(0, 1), name = "expected recall\nof raised bog") +
  geom_sf(data = footprint, fill = NA, colour = "white", linewidth = 0.2) +
  coord_sf(expand = FALSE) +
  labs(title = "C") +
  theme_map
save_fig(p6a + p6b + p6c, "fig_reliability", 7.5, 5.2)

## Table 1: training frame ####

frame <- read_csv("output/pl2/modeling_frame_summary.csv", show_col_types = FALSE) |>
  filter(scenario == "current") |>
  mutate(
    response = CLASS_LAB[response],
    dataset = c(NO = "Norway", EU = "Europe")[dataset]
  ) |>
  select(dataset, response, rows) |>
  pivot_wider(names_from = dataset, values_from = rows) |>
  mutate(total = Norway + Europe)
save_tbl(frame, "tbl_frame")

## Table: partitions ####

part <- read_csv(
  "output/pl2/modeling_frame_regional_partitioned_topfeature.csv",
  col_select = c("scenario", "partition", "envelope_side", "response", "dataset", "bio10"),
  show_col_types = FALSE
) |>
  filter(scenario == "current")

tbl_part <- part |>
  group_by(partition) |>
  summarise(
    bio10_min = round(min(bio10), 2),
    bio10_max = round(max(bio10), 2),
    bog_bio10_min = round(min(bio10[response == "bog"]), 2),
    bog_bio10_max = round(max(bio10[response == "bog"]), 2),
    n_bog = sum(response == "bog"),
    n_otherpeat = sum(response == "otherpeat"),
    n_nonpeat = sum(response == "nonpeat"),
    n_norway = sum(dataset == "NO"),
    n_europe = sum(dataset == "EU"),
    .groups = "drop"
  )
save_tbl(tbl_part, "tbl_partitions")

numbers$rows_above_envelope <- sum(part$envelope_side == "above")
numbers$rows_below_envelope <- sum(part$envelope_side == "below")

## Table 3: cross-validated skill ####

met <- read_csv("output/pl2/metrics_cv_topfeature.csv", show_col_types = FALSE)

tbl_cv <- met |>
  filter(slice %in% c("all", "NO", "EU", "inside", "below", "above")) |>
  group_by(arm, slice) |>
  summarise(
    n_folds = n(),
    n_test_rows = sum(n),
    median_n_train = median(n_train),
    macro_auc = mean(macro_auc, na.rm = TRUE),
    macro_gmean = mean(Gmean_macro, na.rm = TRUE),
    recall_bog = mean(recall_bog, na.rm = TRUE),
    recall_otherpeat = mean(recall_otherpeat, na.rm = TRUE),
    recall_nonpeat = mean(recall_nonpeat, na.rm = TRUE),
    fpr_bog = mean(fpr_bog, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    slice = factor(slice, levels = c("all", "NO", "EU", "inside", "below", "above")),
    arm = factor(arm, levels = c("lopo", "pairwise"))
  ) |>
  arrange(arm, slice) |>
  mutate(across(c(macro_auc:fpr_bog), ~ round(.x, 3)))
save_tbl(tbl_cv, "tbl_cv")

conf <- read_csv("output/pl2/confusion_cv_topfeature.csv", show_col_types = FALSE) |>
  filter(arm == "lopo") |>
  group_by(Actual, Predicted) |>
  summarise(n = sum(Freq), .groups = "drop") |>
  mutate(
    Actual = factor(CLASS_LAB[Actual], levels = CLASS_LAB),
    Predicted = factor(CLASS_LAB[Predicted], levels = CLASS_LAB)
  ) |>
  arrange(Actual, Predicted) |>
  pivot_wider(names_from = Predicted, values_from = n)
save_tbl(conf, "tbl_confusion_lopo")

## S1 Table: predictors ####

pred_desc <- c(
  bio01 = "Annual mean temperature", bio02 = "Mean diurnal temperature range",
  bio03 = "Isothermality", bio04 = "Temperature seasonality",
  bio05 = "Maximum temperature of the warmest month",
  bio06 = "Minimum temperature of the coldest month", bio07 = "Annual temperature range",
  bio08 = "Mean temperature of the wettest quarter",
  bio09 = "Mean temperature of the driest quarter",
  bio10 = "Mean temperature of the warmest quarter",
  bio11 = "Mean temperature of the coldest quarter", bio12 = "Annual precipitation",
  bio13 = "Precipitation of the wettest month", bio14 = "Precipitation of the driest month",
  bio15 = "Precipitation seasonality", bio16 = "Precipitation of the wettest quarter",
  bio17 = "Precipitation of the driest quarter",
  bio18 = "Precipitation of the warmest quarter",
  bio19 = "Precipitation of the coldest quarter",
  fcf = "Frost change frequency", gdd0 = "Growing degree days above 0 °C",
  gdd5 = "Growing degree days above 5 °C", gdd10 = "Growing degree days above 10 °C",
  gsl = "Growing season length", gsp = "Growing season precipitation",
  gst = "Growing season mean temperature", ngd0 = "Days with mean temperature above 0 °C",
  ngd5 = "Days with mean temperature above 5 °C",
  ngd10 = "Days with mean temperature above 10 °C",
  npp = "Potential net primary productivity", scd = "Snow cover days",
  swe = "Snow water equivalent",
  elevation = "Elevation", slope = "Slope",
  paleo_years_icefreeland = "Years as ice-free land (last 16 ka)",
  paleo_bio01_mean_icefree = "Mean annual temperature over the ice-free period",
  paleo_bio01_sum_icefree = "Summed annual temperature over the ice-free period",
  paleo_bio12_mean_icefree = "Mean annual precipitation over the ice-free period",
  paleo_bio12_sum_icefree = "Summed annual precipitation over the ice-free period",
  paleo_bio01_mean_16ka = "Mean annual temperature, 16 ka to present",
  paleo_bio01_sum_16ka = "Summed annual temperature, 16 ka to present",
  paleo_bio12_mean_16ka = "Mean annual precipitation, 16 ka to present",
  paleo_bio12_sum_16ka = "Summed annual precipitation, 16 ka to present"
)
pred_source <- function(f) {
  case_when(
    grepl("^paleo_", f) ~ "CHELSA-TraCE21k, derived",
    f %in% c("elevation", "slope") ~ "DTM50 (Norway); AWS terrain tiles via elevatr (Europe)",
    grepl("^bio", f) ~ "CHELSA V2.1 bioclim",
    TRUE ~ "CHELSA V2.1 bioclim+"
  )
}

w <- read_csv("output/pl2/weights_feature_data_partitioning.csv", show_col_types = FALSE) |>
  filter(method == "Balanced Random Forest")
vi <- read_csv("output/pl2/variable_importance_brf.csv", show_col_types = FALSE)

tbl_pred <- w |>
  left_join(vi, by = c("feature" = "variable")) |>
  transmute(
    predictor = feature,
    description = unname(pred_desc[feature]),
    source = pred_source(feature),
    changes_under_scenario = if_else(dynamic, "yes", "no"),
    in_distance_metric = if_else(dynamic, "yes", "no"),
    projected_shift_sd = round(shift_sd, 2),
    importance_all = round(all, 4),
    importance_bog = round(bog, 4)
  ) |>
  arrange(desc(importance_all))
save_tbl(tbl_pred, "tblS_predictors")

## Study numbers ####

numbers_tbl <- tibble(
  quantity = names(numbers),
  value = vapply(numbers, function(v) as.numeric(v)[1], numeric(1))
)
print(as.data.frame(numbers_tbl), row.names = FALSE)
save_tbl(numbers_tbl, "study_numbers")

cat("\nDone.\n")

# sessionInfo ####

sessioninfo::session_info()
