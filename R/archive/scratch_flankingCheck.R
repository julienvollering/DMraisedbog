# Scratch: is the "flanked on both sides" argument a real modelling problem? ####

# Motivated by the 2026-08-17 notebook critique of the 2026-06-05 entry. The entry
# argues that because raised bog is a climatic intermediate flanked by two ecologically
# opposite non-bog classes, a binary model "lumps two opposite flanks into one 0" and
# therefore cannot work -- hence multiple 2-class models. This script tests that.
#
# Part 1 (simulation, seconds): can a forest recover an intermediate band at all, with
#        both flanks lumped into 0? Establishes the mechanism claim.
# Part 2 (real data, minutes): restore the warm flank to the modelling population that
#        the artype_60 mask deleted, and read the response curves along the top climate
#        features. If a unimodal band appears, the flanking argument is empirically dead.
#
# Scratch only -- not part of any pipeline. Nothing is written to output/.

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(terra)
library(randomForestSRC)

set.seed(7)

# Part 1: simulation ####

## A bog band flanked cold and warm, at realistic prevalence, buried in noise.

n_sim <- 60000
gdd_sim <- runif(n_sim, 0, 100)
noise_sim <- matrix(rnorm(n_sim * 20), n_sim, 20)

# Truth: bog occupies gdd in [40, 55]; other-peat below, non-peat above.
class_sim <- ifelse(
  gdd_sim < 40,
  "otherpeat",
  ifelse(gdd_sim > 55, "nonpeat", "bog")
)

# Thin the bog class to ~1.5% prevalence, matching the regional problem.
keep_sim <- rep(TRUE, n_sim)
keep_sim[class_sim == "bog"] <- runif(sum(class_sim == "bog")) < 0.09

sim <- data.frame(
  y = factor(class_sim[keep_sim]),
  gdd = gdd_sim[keep_sim],
  wet = runif(sum(keep_sim), 0, 100),
  noise_sim[keep_sim, ]
)

cat("Simulated class counts:\n")
print(table(sim$y))
cat("Bog prevalence:", round(mean(sim$y == "bog"), 4), "\n\n")

## Grid holding every predictor but gdd at a constant, for response curves.
grid_sim <- data.frame(gdd = seq(2, 98, 4), wet = 50, matrix(0, 25, 20))
names(grid_sim) <- names(sim)[-1]

## (a) Binary, both flanks lumped into 0 -- the design the 2026-06-05 entry rejects.
sim_binary <- sim
sim_binary$y <- factor(ifelse(sim$y == "bog", "1", "0"), levels = c("0", "1"))

fit_sim_binary <- imbalanced(y ~ ., sim_binary, ntree = 300, method = "brf")
p_sim_binary <- predict(fit_sim_binary, grid_sim)$predicted[, "1"]

cat("Binary (both flanks = 0), P(bog) along gdd -- truth is a band at 40-55:\n")
print(round(setNames(p_sim_binary, grid_sim$gdd), 2))
cat("argmax at gdd =", grid_sim$gdd[which.max(p_sim_binary)], "\n\n")

## (b) One 3-class balanced RF. `imbalanced()` is two-class only, but that is a limit
## of RFQ's q* classifier, not of the package: rfsrc() is natively J-class and the
## BRF recipe below (documented under ?rfsrc) is written generically over the levels.
fit_sim_3class <- rfsrc(
  y ~ .,
  sim,
  ntree = 300,
  case.wt = randomForestSRC:::make.wt(sim$y),
  sampsize = randomForestSRC:::make.size(sim$y)
)
p_sim_3class <- predict(fit_sim_3class, grid_sim)$predicted

cat("3-class single model, probabilities along gdd:\n")
print(round(cbind(gdd = grid_sim$gdd, p_sim_3class), 2))
cat("rows sum to 1:", all(abs(rowSums(p_sim_3class) - 1) < 1e-8), "\n\n")

# Part 2: real data, warm flank restored ####

## Rebuild the modelling frame WITHOUT the artype_60 mask ####

# `pl2_createModelingFrame.R` masks the population to artype_60 >= 0.5 (plus presence
# cells) before extracting. That mask is what deleted the warm flank. Here we skip it
# and instead use artype_60 to *label* the two absence flanks -- the 3-class response
# of the 2026-06-05 entry, obtained from data already on disk.

rf_global_current <- rast("output/rf_global_pred_regional_current.tif")
names(rf_global_current) <- "rf_global"
preds_current <- rast(
  "output/predictors_regional_250m_Norway_current_EPSG3035.tif"
)
raster_current <- c(rf_global_current, preds_current)

presence <- read_csv(
  "output/presence_coords_regional.csv",
  show_col_types = FALSE
)
absence <- read_csv(
  "output/absence_coords_regional.csv",
  show_col_types = FALSE
)

## Label the two absence flanks ####

# artype_60 is fractional peatland cover: absences on peatland are the COLD flank
# (other mire), absences off peatland are the WARM flank (forest/heath/agriculture).
# Extract that one layer over all absences first (~15 s) so the flanks can be sampled
# deliberately; extracting the full stack over 1.76e6 points would be wasteful.
absence <- absence |>
  mutate(
    artype_60 = terra::extract(
      raster_current[["artype_60"]],
      absence[c("x", "y")],
      ID = FALSE
    )[[1]]
  ) |>
  drop_na(artype_60) |>
  mutate(flank = if_else(artype_60 >= 0.5, "otherpeat", "nonpeat"))

cat("Absences available unmasked:", nrow(absence), "\n")
print(count(absence, flank))
cat(
  "Warm flank share:",
  round(mean(absence$flank == "nonpeat"), 3),
  "-- this is the pool the artype_60 mask deleted.\n",
  "Cold-flank count should match the 61850 absences pl2's masked frame retains.\n\n"
)

## Compose the absence sample ####

# "stratified" samples the two flanks equally. Under the "natural" (uniform-background)
# ratio the warm flank swamps the cold one, so the loss is dominated by the easy
# contrast -- the stratification failure diagnosed in the 2026-08-17 entry. Flip this to
# "natural" to reproduce that failure; it is a sampling knob, not an architecture choice.
flank_design <- "stratified"
n_absence_sample <- 60000

absence_sample <- switch(
  flank_design,
  stratified = absence |>
    group_by(flank) |>
    slice_sample(n = n_absence_sample %/% 2) |>
    ungroup(),
  natural = absence |> slice_sample(n = n_absence_sample)
)

extract_at <- function(coords) {
  terra::extract(raster_current, coords[c("x", "y")], ID = FALSE) |>
    tibble()
}

df_presence <- extract_at(presence) |>
  mutate(class = "bog", .before = 1)

df_absence <- extract_at(absence_sample) |>
  mutate(class = absence_sample$flank, .before = 1)

rm(raster_current, rf_global_current, preds_current, absence)
gc()

frame_unmasked <- bind_rows(df_presence, df_absence) |>
  select(-starts_with("artype")) |> # artype must not be a predictor (see 2025-09-15)
  drop_na() |>
  mutate(class = factor(class, levels = c("otherpeat", "bog", "nonpeat")))

cat("Unmasked 3-class frame (", flank_design, "flanks ):\n")
print(count(frame_unmasked, class))
cat("\n")

predictors <- setdiff(names(frame_unmasked), "class")

## Which features to inspect ####

# BRF variable importance from pl2, which the notebook (2025-12-06) prefers over RFQ VI.
top_features <- read_csv(
  "output/pl2/weights_feature_data_partitioning.csv",
  show_col_types = FALSE
) |>
  filter(method == "Balanced Random Forest", feature %in% predictors) |>
  slice_max(median, n = 6) |>
  pull(feature)

cat("Inspecting:", paste(top_features, collapse = ", "), "\n\n")

## FOP curves: the model-free view ####

# plotFOP shows the empirical frequency of presence along each feature. If the raw data
# carry a unimodal signal once the warm flank is present, it shows up here before any
# model is fitted.
frame_fop <- frame_unmasked |>
  mutate(RV = as.numeric(class == "bog"), .before = 1) |>
  select(-class) |>
  as.data.frame()

fop <- lapply(top_features, \(f) {
  out <- MIAmaxent::plotFOP(frame_fop, EV = f)
  out$FOPdata |>
    mutate(feature = f, EVoptimum = out$EVoptimum)
})

plot_fop <- bind_rows(fop) |>
  ggplot(aes(intEV, intRV)) +
  geom_line(aes(y = loess), colour = "steelblue", linewidth = 0.8) +
  geom_point(size = 0.8) +
  geom_vline(aes(xintercept = EVoptimum), linetype = 2, colour = "firebrick") +
  facet_wrap(~feature, scales = "free", ncol = 3) +
  labs(
    title = "FOP, warm flank restored",
    subtitle = "Interior peak = unimodal signal present in the raw data",
    x = NULL,
    y = "Frequency of observed presence"
  ) +
  theme_bw()

print(plot_fop)

## Fit the two candidate architectures ####

# Bog:absence imbalance is left as it is -- BRF and case.wt are what handle it. The
# flank composition was already set above by `flank_design`.
frame_fit <- frame_unmasked

## (a) Binary: both flanks lumped into 0, exactly the design under dispute.
frame_binary <- frame_fit |>
  mutate(
    response = factor(if_else(class == "bog", "1", "0"), levels = c("0", "1"))
  ) |>
  select(-class) |>
  as.data.frame()

fit_binary <- imbalanced(
  response ~ .,
  frame_binary,
  ntree = 1000,
  method = "brf",
  do.trace = FALSE # see CLAUDE.md: do.trace crashes Positron in nested contexts
)

## (b) One 3-class balanced RF.
frame_3class <- as.data.frame(frame_fit)

fit_3class <- rfsrc(
  class ~ .,
  frame_3class,
  ntree = 1000,
  case.wt = randomForestSRC:::make.wt(frame_3class$class),
  sampsize = randomForestSRC:::make.size(frame_3class$class),
  do.trace = FALSE
)

cat("\n3-class OOB confusion matrix -- the cold contrast cannot hide here:\n")
print(get.confusion(fit_3class$yvar, fit_3class$predicted.oob))

## Response curves along each top feature ####

# Vary one feature over its observed range, hold the rest at their median. Crude next to
# a full partial dependence, but enough to answer the shape question.
response_curve <- function(feature, data, n_grid = 40) {
  medians <- data |>
    select(all_of(predictors)) |>
    summarise(across(everything(), \(x) median(x, na.rm = TRUE)))

  grid <- medians[rep(1, n_grid), ]
  grid[[feature]] <- seq(
    quantile(data[[feature]], 0.01, na.rm = TRUE),
    quantile(data[[feature]], 0.99, na.rm = TRUE),
    length.out = n_grid
  )
  grid <- as.data.frame(grid)

  binary <- predict(fit_binary, grid)$predicted[, "1"]
  three <- predict(fit_3class, grid)$predicted

  tibble(
    feature = feature,
    value = grid[[feature]],
    `binary: P(bog)` = binary,
    `3-class: P(bog)` = three[, "bog"],
    `3-class: P(other-peat)` = three[, "otherpeat"],
    `3-class: P(non-peat)` = three[, "nonpeat"]
  )
}

curves <- lapply(top_features, response_curve, data = frame_fit) |>
  bind_rows() |>
  pivot_longer(-c(feature, value), names_to = "quantity", values_to = "p")

plot_bog <- curves |>
  filter(quantity %in% c("binary: P(bog)", "3-class: P(bog)")) |>
  ggplot(aes(value, p, colour = quantity)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~feature, scales = "free_x", ncol = 3) +
  labs(
    title = "P(bog): binary with both flanks lumped vs one 3-class model",
    subtitle = "If both are unimodal, the flanking argument is empirically dead",
    x = NULL,
    y = "P(bog)",
    colour = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

plot_flanks <- curves |>
  filter(quantity != "binary: P(bog)") |>
  ggplot(aes(value, p, fill = quantity)) +
  geom_area(position = "stack") +
  facet_wrap(~feature, scales = "free_x", ncol = 3) +
  labs(
    title = "3-class composition: which flank the residual mass goes to",
    x = NULL,
    y = "Probability",
    fill = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(plot_bog / plot_flanks)

## What to conclude ####

# - Both curves unimodal  -> the flanking argument does not motivate multiple models.
#   Keep the 3-class response for the per-flank weighting knob, the 3x3 confusion matrix
#   and the endpoint composition; fit it with ONE rfsrc() call.
# - Binary flat/monotone but 3-class unimodal -> the labels carry information the lumped
#   response destroys. Still one model; still no pairwise decomposition.
# - Neither unimodal -> the problem is upstream (features, scale, absence quality), and
#   no architecture rescues it.

## Conclusion ####

# - Neither of the response curves are unimodel, despite the FOP showing a clear unimodal signal in the raw data.
# --> the problem may be upstream in the form of features (more attention to terrain needed) or scale (EU extent needed to lengthen the warm flank axes).

# sessionInfo ####

sessioninfo::session_info()
