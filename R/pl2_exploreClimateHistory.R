# The climate history of today's bogs against their future climate ####

# PURPOSE: Asks whether the projected 2071-2100 climate of today's Norwegian raised bogs lies inside the climate those same bogs have lived through since an assumed transition age t0, using the TraCE21k centennial series.

# Model-free and read-only to the model, like pl2_exploreOccupancy.R. The design, its
# caveats and the order of work are in plan_climateHistory.md; section numbers below
# (3.3, 3.5, ...) refer to it.
#
# WHY. pl2_exploreOccupancy.R compares the future with the climates bogs occupy TODAY, across
# space. Raised bogs are millennial structures, so the relevant reference is the climate each
# bog has persisted through. A future that is novel against today's niche may lie inside
# the Holocene thermal maximum (HTM, ca. 8-6 ka) envelope.
#
# THE ONE UNKNOWN IS t0, the age at which a bog became a bog. The headline is therefore a
# coverage CURVE over t0, not a number: the share of future bog pixels inside the climate
# space experienced since t0, with t* = the smallest t0 at which coverage reaches
# COVERAGE_THRESHOLD.
#
# UNIT. The ~1 km TraCE/CHELSA pixel holding each Norwegian bog cell (several 250 m cells
# share a pixel). Pixels are NOT independent: after the delta correction each pixel's past
# is its own CHELSA baseline plus a CCSM3 anomaly from a ~3.75 degree grid, so all of Norway
# sits in a handful of GCM cells. Intervals are bootstrapped by GCM cell (3.2).
#
# THREE CONTAINMENT DEFINITIONS (3.5), reported together:
#   1  per-pixel headroom: future bio10 minus the warmest slice in that pixel's own window
#      ("this bog has been this warm before"); 1b adds the pixel's own bio12 range.
#   2  pooled convex hull of every bog pixel's experienced (bio10, bio12).
#   3  pooled density: inside the 95% highest-density region. Not monotone in t0 and trims
#      the rare warm slices by design, so t* is defined on 1 and 2 only.

library(tidyverse)
library(terra)

source("R/config.R")
source("R/functions.R")

## Configuration ####

SERIES_DIR <- "data/CHELSA/paleo_derived"

# t0 sweep, ka b2k. TraCE time t (hundred-year units) is age (20 - t) * 100 years b2k; the
# literature quotes BP (= b1950), 50 years younger. Stated once here, ignored below.
T0_GRID <- seq(0, 10, by = 0.5)

# Placeholder: JV's call, read against the bog-age literature (plan 3.7, 6).
COVERAGE_THRESHOLD <- 0.95

# Pooled density region (definition 3): the highest-density region holding this share of
# the experienced climate mass.
HDR_LEVEL <- 0.95

# Bootstrap over GCM cells. CCSM3 T31: 3.75 degrees in longitude, Gaussian latitudes ~3.71
# degrees apart. The grid below approximates it; that is enough to group pixels that share
# one anomaly field.
N_BOOT <- 200
GCM_DLON <- 3.75
GCM_DLAT <- 3.71

# Pixels whose mean bog-cell elevation is below this may have emerged from the sea later
# than the TraCE palaeotopography says (3.4): the deglaciation cap is a MAXIMUM age there.
# Flagged, never dropped. A rough stand-in for the Norwegian marine limit.
LOW_ELEVATION_FLAG_M <- 150

# Ages (ka) at which per-pixel headroom is written out and slices are drawn.
REPORT_AGES <- c(0, 2, 4, 6, 8, 10)

record_settings(
  "R/pl2_exploreClimateHistory.R",
  t0_grid = range(T0_GRID),
  t0_step = diff(T0_GRID)[1],
  coverage_threshold = COVERAGE_THRESHOLD,
  hdr_level = HDR_LEVEL,
  n_boot = N_BOOT,
  gcm_cell_deg = c(GCM_DLON, GCM_DLAT),
  low_elevation_flag_m = LOW_ELEVATION_FLAG_M,
  scenarios = FUTURE_SCENARIOS,
  axes_primary = c("bio10", "bio12")
)

set.seed(20261006)

## Bog pixels ####

# The same cells as future_bog in pl2_exploreOccupancy.R: Norwegian bog presences in the
# current training frame. Never a resample.
mf <- read_csv("output/pl2/modeling_frame_regional.csv", show_col_types = FALSE)
bog_cells <- mf |>
  filter(scenario == "current", response == "bog", dataset == "NO") |>
  select(x, y, elevation)

series <- c("bio10", "bio12", "bio18", "bio06") |>
  set_names() |>
  map(\(v) rast(file.path(SERIES_DIR, paste0("trace_series_", v, "_Norway.tif"))))
times <- as.integer(sub("^t", "", names(series$bio10)))
ages_ka <- (20 - times) / 10 # ka b2k

pts <- vect(bog_cells, geom = c("x", "y"), crs = "EPSG:3035") |> project("EPSG:4326")
bog_cells$pixel <- cellFromXY(series$bio10, crds(pts))

pixels <- bog_cells |>
  group_by(pixel) |>
  summarise(n_cells = n(), elevation = mean(elevation), .groups = "drop")
pixels <- bind_cols(pixels, as_tibble(xyFromCell(series$bio10, pixels$pixel))) |>
  rename(lon = x, lat = y)
pix_xy <- as.matrix(pixels[c("lon", "lat")])

cat("Norwegian bog cells:", nrow(bog_cells), "-> distinct TraCE pixels:", nrow(pixels), "\n")

## Climate at the pixels ####

# TraCE series: pixels x slices, oldest first (layer order is time ascending).
trace_mat <- map(series, \(s) as.matrix(terra::extract(s, pixels$pixel)))

# CHELSA v2.1 at the same pixel centres. Same 30 arcsec grid as TraCE, so a centre read is
# that pixel's value.
axis_vars <- names(series)
chelsa_at <- function(scenario, ssp = FUTURE_SCENARIO) {
  chelsa_native_stack(scenario, ssp)[[axis_vars]] |>
    terra::extract(pix_xy, method = "simple") |>
    as_tibble()
}
current <- chelsa_at("current")
future <- FUTURE_SCENARIOS |>
  set_names() |>
  map(\(ssp) chelsa_at("future", ssp))

# Deglaciation cap (3.4): a bog cannot be older than its ground.
icefree_ka <- rast(file.path(SERIES_DIR, "paleo_years_icefreeland_EUextent_EPSG4326.tif")) |>
  terra::extract(pix_xy, method = "simple") |>
  pull(1) / 1000

ok <- complete.cases(do.call(cbind, trace_mat), current, do.call(cbind, future), icefree_ka)
cat("Pixels dropped for missing climate (coast/edge):", sum(!ok), "of", length(ok), "\n")

pixels <- pixels[ok, ] |>
  mutate(
    icefree_ka = icefree_ka[ok],
    low_elevation = elevation < LOW_ELEVATION_FLAG_M,
    gcm_cell = paste(
      floor((lon + GCM_DLON / 2) / GCM_DLON), floor(lat / GCM_DLAT),
      sep = "_"
    )
  )
trace_mat <- map(trace_mat, \(m) m[ok, , drop = FALSE])
current <- current[ok, ]
future <- map(future, \(f) f[ok, ])
n_pix <- nrow(pixels)

cat("Pixels analysed:", n_pix, "in", n_distinct(pixels$gcm_cell), "approximate GCM cells\n")
print(count(pixels, gcm_cell, name = "n_pixels"))
cat("Pixels flagged low-elevation (<", LOW_ELEVATION_FLAG_M, "m):", sum(pixels$low_elevation), "\n")
cat("Deglaciation age at the pixels (ka), quantiles:\n")
print(round(quantile(pixels$icefree_ka, c(0, 0.05, 0.5, 0.95, 1)), 1))

## TraCE baseline and the CHELSA bridge (3.3) ####

# The last TraCE slice (t = 20) is a centennial mean ending with the TraCE run (1990 CE),
# i.e. mostly 20th-century climate; CHELSA v2.1 current is 1981-2010. The delta method
# anchors the TraCE anomaly to CHELSA 1981-2010, so any warming between the two periods
# makes the reconstructed past too warm by that amount (bias toward "inside"). CHELSA-TraCE
# is itself delta-downscaled onto a CHELSA climatology, so this is a plumbing check of the
# baseline, not a test of HTM fidelity.
last <- length(times)
offset <- map_dfr(axis_vars, \(v) {
  d <- trace_mat[[v]][, last] - current[[v]]
  r <- trace_mat[[v]][, last] / current[[v]]
  tibble(
    variable = v,
    median_trace_minus_chelsa = median(d),
    q05 = quantile(d, 0.05), q95 = quantile(d, 0.95),
    median_ratio_trace_over_chelsa = median(r)
  )
})
cat("\nTraCE last slice (t = 20) minus CHELSA 1981-2010 at the bog pixels:\n")
print(offset |> mutate(across(where(is.numeric), \(x) round(x, 3))))
write_csv(offset, "output/pl2/climatehistory_offset.csv", append = FALSE)

# The simulated HTM anomaly, for comparison with proxy summer-temperature reconstructions
# (3.8): per pixel, warmest 8-6 ka slice and mean of 8-6 ka, minus the last slice.
htm <- ages_ka >= 6 & ages_ka <= 8
htm_anomaly <- tibble(
  variable = c("bio10", "bio06"),
  median_htm_max_minus_t20 = map_dbl(variable, \(v)
    median(apply(trace_mat[[v]][, htm], 1, max) - trace_mat[[v]][, last])),
  median_htm_mean_minus_t20 = map_dbl(variable, \(v)
    median(rowMeans(trace_mat[[v]][, htm]) - trace_mat[[v]][, last]))
)
cat("\nSimulated HTM (8-6 ka) anomaly against the last slice, median over bog pixels:\n")
print(htm_anomaly |> mutate(across(where(is.numeric), \(x) round(x, 2))))
write_csv(htm_anomaly, "output/pl2/climatehistory_htm_anomaly.csv", append = FALSE)

## Experienced climate: delta and raw versions ####

# Temperatures additive; precipitation as a RATIO (an additive delta can go negative at dry
# sites). Raw TraCE is the sensitivity.
is_temperature <- function(v) v %in% c("bio10", "bio06")
history <- list(
  delta = map(set_names(axis_vars), \(v) {
    m <- trace_mat[[v]]
    if (is_temperature(v)) {
      current[[v]] + (m - m[, last])
    } else {
      current[[v]] * (m / m[, last])
    }
  }),
  raw = trace_mat
)

## Containment ####

# Window mask for t0: slice ages <= min(t0, deglaciation age). The present slice is always in.
window_mask <- function(t0) {
  outer(pmin(t0, pixels$icefree_ka), ages_ka, ">=") | col(matrix(0, n_pix, length(ages_ka))) == last
}

row_max_in <- function(m, w) apply(ifelse(w, m, -Inf), 1, max)
row_min_in <- function(m, w) apply(ifelse(w, m, Inf), 1, min)

# Point-in-convex-polygon by edge cross products; orientation-free.
in_hull <- function(px, py, ex, ey) {
  keep <- is.finite(ex) & is.finite(ey)
  ex <- ex[keep]
  ey <- ey[keep]
  h <- chull(ex, ey)
  hx <- ex[h]
  hy <- ey[h]
  if (length(h) < 3) {
    return(rep(FALSE, length(px)))
  }
  j <- c(2:length(h), 1)
  cr <- outer(py, hy, "-") * rep(hx[j] - hx, each = length(px)) -
    outer(px, hx, "-") * rep(hy[j] - hy, each = length(px))
  tol <- 1e-9 * max(abs(cr))
  apply(cr, 1, \(r) all(r >= -tol) || all(r <= tol))
}

# Weighted 2-D Gaussian KDE on axes standardised by the current bog pixels. Each pixel
# carries weight w_i spread evenly over its window, so a pixel with a long history does not
# outweigh one with a short history. Evaluated on a grid; points read off the nearest node.
std <- c(sd(current$bio10), sd(current$bio12))
in_hdr <- function(px, py, ex, ey, ew, level = HDR_LEVEL, n_grid = 120) {
  ex <- ex / std[1]
  ey <- ey / std[2]
  px <- px / std[1]
  py <- py / std[2]
  ew <- ew / sum(ew)
  bw <- c(MASS::bandwidth.nrd(ex), MASS::bandwidth.nrd(ey)) / 4 # kde2d convention
  gx <- seq(min(ex, px) - 3 * bw[1], max(ex, px) + 3 * bw[1], length.out = n_grid)
  gy <- seq(min(ey, py) - 3 * bw[2], max(ey, py) + 3 * bw[2], length.out = n_grid)
  kx <- dnorm(outer(gx, ex, "-") / bw[1])
  ky <- dnorm(outer(gy, ey, "-") / bw[2])
  dens <- kx %*% (ew * t(ky))
  node <- function(v, g) pmin(pmax(round((v - g[1]) / (g[2] - g[1])) + 1, 1), length(g))
  d_e <- dens[cbind(node(ex, gx), node(ey, gy))]
  o <- order(d_e)
  cut <- d_e[o][which(cumsum(ew[o]) >= 1 - level)[1]]
  dens[cbind(node(px, gx), node(py, gy))] >= cut
}

# Coverage of the future pixels for one t0, one history version, one scenario, one axis
# pair. `w` is a per-pixel multiplicity (bootstrap); 1 everywhere for the point estimate.
coverage <- function(t0, hist, fut, tvar, mvar, w = rep(1, n_pix), dens = TRUE) {
  win <- window_mask(t0)
  tx <- hist[[tvar]]
  my <- hist[[mvar]]
  fx <- fut[[tvar]]
  fy <- fut[[mvar]]

  headroom <- fx - row_max_in(tx, win)
  inside_1 <- headroom <= 0
  inside_1b <- inside_1 & fy >= row_min_in(my, win) & fy <= row_max_in(my, win)

  sel <- w > 0
  e_x <- tx[sel, ][win[sel, ]]
  e_y <- my[sel, ][win[sel, ]]
  inside_2 <- in_hull(fx, fy, e_x, e_y)

  out <- tibble(
    definition = c("1 headroom", "1b own range", "2 pooled hull"),
    coverage = c(
      weighted.mean(inside_1, w), weighted.mean(inside_1b, w), weighted.mean(inside_2, w)
    )
  )
  if (dens) {
    # Each pixel's weight spread over its window slices.
    e_w <- (w / rowSums(win))[sel] * win[sel, ]
    e_w <- e_w[win[sel, ]]
    inside_3 <- in_hdr(fx, fy, e_x, e_y, e_w)
    out <- add_row(out, definition = "3 pooled density", coverage = weighted.mean(inside_3, w))
  }
  out
}

# Headroom-only coverage on a single temperature axis (bio06 sensitivity).
coverage_headroom <- function(t0, hist, fut, tvar) {
  mean(fut[[tvar]] - row_max_in(hist[[tvar]], window_mask(t0)) <= 0)
}

## Coverage curves ####

AXES <- list(
  primary = c("bio10", "bio12"),
  bio18 = c("bio10", "bio18")
)

grid_runs <- expand_grid(
  scenario = FUTURE_SCENARIOS,
  version = names(history),
  axes = names(AXES),
  t0 = T0_GRID
)

curves <- grid_runs |>
  mutate(res = pmap(list(scenario, version, axes, t0), \(ssp, ver, ax, t0) {
    coverage(t0, history[[ver]], future[[ssp]], AXES[[ax]][1], AXES[[ax]][2])
  })) |>
  unnest(res)

winter <- expand_grid(scenario = FUTURE_SCENARIOS, version = names(history), t0 = T0_GRID) |>
  mutate(
    axes = "bio06",
    definition = "1 headroom",
    coverage = pmap_dbl(list(scenario, version, t0), \(ssp, ver, t0) {
      coverage_headroom(t0, history[[ver]], future[[ssp]], "bio06")
    })
  )
curves <- bind_rows(curves, winter)

## Bootstrap by GCM cell ####

# Primary axes, delta version only: the headline. Resample GCM cells with replacement; a
# pixel's multiplicity is the number of times its cell was drawn.
cells <- unique(pixels$gcm_cell)
boot_w <- map(seq_len(N_BOOT), \(b) {
  drawn <- table(factor(sample(cells, replace = TRUE), levels = cells))
  as.numeric(drawn[pixels$gcm_cell])
})

boot <- expand_grid(scenario = FUTURE_SCENARIOS, t0 = T0_GRID) |>
  mutate(res = map2(scenario, t0, \(ssp, t0) {
    map_dfr(boot_w, \(w) coverage(t0, history$delta, future[[ssp]], "bio10", "bio12", w = w))
  })) |>
  unnest(res) |>
  group_by(scenario, t0, definition) |>
  summarise(lo = quantile(coverage, 0.05), hi = quantile(coverage, 0.95), .groups = "drop") |>
  mutate(version = "delta", axes = "primary")

curves <- curves |> left_join(boot, by = c("scenario", "version", "axes", "t0", "definition"))
write_csv(curves, "output/pl2/climatehistory_coverage.csv", append = FALSE)

## Reference: today's pooled NO+EU bog climate ####

# fig_occupancy's reference is all training bogs, Norwegian and European, today. Drawn as a
# line on the coverage figure so the two can be read together. Read off the 250 m frame,
# not the pixels, because the EU bogs have no TraCE pixel here.
ref_bogs <- mf |> filter(scenario == "current", response == "bog")
reference <- map_dfr(FUTURE_SCENARIOS, \(ssp) {
  fut <- future[[ssp]]
  tibble(
    scenario = ssp,
    definition = c("2 pooled hull", "3 pooled density"),
    coverage = c(
      mean(in_hull(fut$bio10, fut$bio12, ref_bogs$bio10, ref_bogs$bio12)),
      mean(in_hdr(fut$bio10, fut$bio12, ref_bogs$bio10, ref_bogs$bio12, rep(1, nrow(ref_bogs))))
    )
  )
})
cat("\nReference: future NO bog pixels inside today's pooled NO+EU bog climate:\n")
print(reference)
write_csv(reference, "output/pl2/climatehistory_reference.csv", append = FALSE)

## t* ####

t_star <- curves |>
  filter(definition %in% c("1 headroom", "1b own range", "2 pooled hull")) |>
  group_by(scenario, version, axes, definition) |>
  summarise(
    coverage_t0_0 = coverage[t0 == 0],
    coverage_t0_max = coverage[t0 == max(T0_GRID)],
    t_star = if (any(coverage >= COVERAGE_THRESHOLD)) min(t0[coverage >= COVERAGE_THRESHOLD]) else NA_real_,
    .groups = "drop"
  )
cat("\nt* (smallest t0 with coverage >=", COVERAGE_THRESHOLD, "; NA = never within",
  max(T0_GRID), "ka):\n")
print(t_star |> mutate(across(starts_with("coverage"), \(x) round(x, 3))), n = Inf)
write_csv(t_star, "output/pl2/climatehistory_tstar.csv", append = FALSE)

cat("\nHeadline (delta, bio10 x bio12):\n")
curves |>
  filter(version == "delta", axes == "primary", t0 %in% REPORT_AGES) |>
  select(scenario, definition, t0, coverage, lo, hi) |>
  mutate(across(c(coverage, lo, hi), \(x) round(x, 3))) |>
  pivot_wider(names_from = t0, values_from = c(coverage), names_prefix = "t0=", id_cols = c(scenario, definition)) |>
  print(n = Inf)

## Per-pixel headroom ####

pixel_out <- expand_grid(scenario = FUTURE_SCENARIOS, t0 = REPORT_AGES) |>
  mutate(res = map2(scenario, t0, \(ssp, t0) {
    win <- window_mask(t0)
    tibble(
      pixel = pixels$pixel,
      headroom_bio10 = future[[ssp]]$bio10 - row_max_in(history$delta$bio10, win),
      headroom_bio10_raw = future[[ssp]]$bio10 - row_max_in(history$raw$bio10, win)
    )
  })) |>
  unnest(res) |>
  left_join(pixels, by = "pixel")
write_csv(pixel_out, "output/pl2/climatehistory_pixels.csv", append = FALSE)

## Figures ####

scen_cols <- setNames(c("#c1462c", "#e39a3b"), FUTURE_SCENARIOS)

### Coverage curve ####

p_cov <- curves |>
  filter(axes == "primary") |>
  ggplot(aes(t0, coverage, colour = definition, linetype = version)) +
  geom_ribbon(
    aes(ymin = lo, ymax = hi, fill = definition),
    alpha = 0.15, colour = NA, na.rm = TRUE
  ) +
  geom_line() +
  geom_hline(yintercept = COVERAGE_THRESHOLD, colour = "grey50", linewidth = 0.3) +
  geom_hline(
    data = reference, aes(yintercept = coverage, colour = definition),
    linetype = "dotted", inherit.aes = FALSE
  ) +
  facet_wrap(~scenario) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Share of future bog pixels inside the climate experienced since t0",
    subtitle = paste(
      "bio10 x bio12. Bands: 90% GCM-cell bootstrap (delta).",
      "Dotted: today's pooled NO+EU bog climate. Grey: threshold."
    ),
    x = "Assumed fen-to-bog transition age t0 (ka b2k)", y = "Coverage"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical")
print(p_cov)
ggsave("output/pl2/climatehistory_coverage.png", p_cov, width = 9, height = 5.5, dpi = 150)

### Trajectories: bio10 against age ####

traj <- as_tibble(history$delta$bio10, .name_repair = "minimal") |>
  set_names(ages_ka) |>
  mutate(pixel = pixels$pixel, icefree_ka = pixels$icefree_ka) |>
  pivot_longer(-c(pixel, icefree_ka), names_to = "age_ka", values_to = "bio10") |>
  mutate(age_ka = as.numeric(age_ka)) |>
  filter(age_ka <= icefree_ka)
traj_band <- traj |>
  group_by(age_ka) |>
  summarise(q05 = quantile(bio10, 0.05), med = median(bio10), q95 = quantile(bio10, 0.95))
fut_band <- imap_dfr(future, \(f, ssp) tibble(
  scenario = ssp,
  q05 = quantile(f$bio10, 0.05), med = median(f$bio10), q95 = quantile(f$bio10, 0.95)
)) |>
  mutate(age_ka = -0.4 - 0.4 * (row_number() - 1))

p_traj <- ggplot(traj, aes(age_ka, bio10)) +
  geom_line(aes(group = pixel), alpha = 0.04, linewidth = 0.2) +
  geom_ribbon(data = traj_band, aes(y = med, ymin = q05, ymax = q95), fill = "#2b6a3f", alpha = 0.25) +
  geom_line(data = traj_band, aes(y = med), colour = "#2b6a3f", linewidth = 0.8) +
  geom_linerange(data = fut_band, aes(y = med, ymin = q05, ymax = q95, colour = scenario), linewidth = 2) +
  geom_point(data = fut_band, aes(y = med, colour = scenario), size = 2) +
  scale_x_reverse(breaks = 0:10) +
  scale_colour_manual(values = scen_cols) +
  labs(
    title = "Warmest-quarter temperature at today's Norwegian bog pixels, 10 ka to 2071-2100",
    subtitle = "TraCE21k anomaly on CHELSA 1981-2010 (delta). Lines: pixels within their ice-free window; band: 5-95%.",
    x = "Age (ka b2k)", y = "bio10 (degC)", colour = "2071-2100"
  ) +
  theme_minimal()
print(p_traj)
ggsave("output/pl2/climatehistory_trajectories.png", p_traj, width = 9, height = 5, dpi = 150)

### Slices in bio10 x bio12 ####

slice_pts <- map_dfr(REPORT_AGES[REPORT_AGES <= 8], \(a) {
  k <- which.min(abs(ages_ka - a))
  tibble(
    age = paste(a, "ka"),
    bio10 = history$delta$bio10[, k], bio12 = history$delta$bio12[, k],
    alive = pixels$icefree_ka >= a
  )
}) |>
  filter(alive) |>
  mutate(age = fct_inorder(age))
fut_pts <- imap_dfr(future, \(f, ssp) tibble(scenario = ssp, bio10 = f$bio10, bio12 = f$bio12))

p_slices <- ggplot(slice_pts, aes(bio10, bio12)) +
  geom_point(data = fut_pts, aes(colour = scenario), size = 0.4, alpha = 0.4) +
  geom_point(size = 0.4, colour = "#2b6a3f", alpha = 0.6) +
  facet_wrap(~age, nrow = 1) +
  scale_colour_manual(values = scen_cols) +
  scale_y_log10() +
  labs(
    title = "Bog pixels' climate in single TraCE slices (green) against 2071-2100",
    x = "bio10 (degC)", y = "bio12 (mm, log)", colour = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_slices)
ggsave("output/pl2/climatehistory_slices.png", p_slices, width = 11, height = 3.8, dpi = 150)

### Headroom map ####

p_map <- pixel_out |>
  filter(scenario == FUTURE_SCENARIO, t0 %in% c(0, 6, 10)) |>
  mutate(t0 = paste("t0 =", t0, "ka")) |>
  ggplot(aes(lon, lat, colour = headroom_bio10)) +
  geom_point(size = 0.8) +
  facet_wrap(~t0) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b", midpoint = 0) +
  coord_quickmap() +
  labs(
    title = paste("Headroom:", FUTURE_SCENARIO, "bio10 minus the warmest slice since t0 (delta)"),
    subtitle = "> 0: warmer than anything this pixel has experienced since t0",
    colour = "degC", x = NULL, y = NULL
  ) +
  theme_minimal()
print(p_map)
ggsave("output/pl2/climatehistory_headroom_map.png", p_map, width = 10, height = 5.5, dpi = 150)

# sessionInfo ####

sessioninfo::session_info()
