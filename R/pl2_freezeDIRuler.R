# Freeze the dissimilarity-index ruler ####

# PURPOSE: Derives the one DI ruler -- feature weights, feature scaling, and the normalisation constant -- from the pooled training frame, so every distance measured downstream sits on a single axis.

# Plan section 1.4. DI is a VI-weighted, scaled, normalised distance, which means it has
# three tunable parts and therefore three ways to stop being comparable between the places
# it is measured. This script fixes all three, once, and every consumer reads the result
# rather than deriving its own.
#
# WHY THIS IS ITS OWN STEP. The ruler used to be built inside pl2_predict.R, which runs
# LAST. That was fine while DI was only a per-fold diagnostic, but it makes the ruler
# unavailable to pl2_evaluate.R, which runs earlier -- so the cross-validation normalised
# each fold by its own constant instead. Those per-fold constants span 0.0362 to 0.0634, a
# factor of 1.75, so CV DI was not comparable across folds, let alone with the DI of the
# future projection. Hoisting the ruler here removes the circularity: it depends only on
# the modelling frame and the feature weights, both of which exist before any distance is
# measured.
#
# WHAT THAT BUYS. Once CV DI and future DI are on one axis, the cross-validation stops
# being only a verdict on the model and becomes the calibration dataset for reliability:
# skill can be read AT the novelty the projection actually sits at, instead of averaged
# over a CV whose median novelty sits below it. Measured on this ruler the CV spans the
# future -- its 90th percentile of DI (0.799) covers 96.2% of future cells and 100% of
# future bog cells -- which is why the separate pl3 cross-validation was dropped rather
# than written.
#
# THE THREE PARTS, AND WHY EACH IS FROZEN:
#
#  - WEIGHTS. Restricted to projection-dynamic predictors. A predictor identical in the
#    current and future frames contributes exactly zero to a current-to-future distance
#    while still inflating the metric's denominator, so weight on one is dilution, not
#    neutrality. slope + elevation + the nine paleo_* variables hold 55.5% of the
#    unrestricted weight mass and are bit-identical between scenarios. They stay in the
#    MODEL; they are excluded from the METRIC.
#  - SCALING. A property of the feature space, not of a fold. Derived per fold it explodes:
#    `gsp` is near-constant inside the narrow partitions, which put partition-1 rows 1.5
#    million fold-SDs away on that one axis and drove DI to 1.9e5.
#  - NORMALISATION CONSTANT. The grand mean of pairwise training distances. Per-fold it
#    rescales the axis between measurements that are supposed to lie on one curve.

library(readr)
library(dplyr)
library(ggplot2)

source("R/functions.R")

## Configuration ####

# Rows sampled to estimate the normalisation constant. It is the grand mean of the
# pairwise distances, so a subsample of rows estimates it without bias; computing all n^2
# pairs costs hours and buys a decimal place nothing downstream can see.
NORM_SAMPLE <- 2000
SEED <- 42

## Inputs ####

mf <- read_csv("output/pl2/modeling_frame_regional.csv", show_col_types = FALSE)

train <- mf |> filter(scenario == "current")

# The production training set defines the ruler, so the basis is the whole pooled current
# frame -- not the partitioned subset, which drops the rows that fall in partition gaps.
predictor_names <- setdiff(
  names(train),
  c("scenario", "response", "dataset", "x", "y")
)

weights_features <- read_csv(
  "output/pl2/weights_feature_data_partitioning.csv",
  show_col_types = FALSE
)

weights <- weights_features |>
  filter(method == "Balanced Random Forest", dynamic) |>
  select(feature, median) |>
  arrange(desc(median)) |>
  tibble::deframe()

weights <- weights[names(weights) %in% predictor_names]
metric_names <- names(weights)

cat("Training rows:", nrow(train), "\n")
cat(
  "Metric variables:", length(metric_names), "of", length(predictor_names),
  "predictors (projection-dynamic only)\n"
)
cat("Excluded from the metric:", paste(setdiff(predictor_names, metric_names), collapse = ", "), "\n\n")

## Derive the ruler ####

scaling <- compute_feature_scaling(train, metric_names)

# calculate_weighted_di() needs a test set, but nothing is being measured here: the
# quantity being retained is the normalisation constant it derives from `train_data`.
# One row is passed so the nearest-neighbour step is trivial.
cat("Computing the normalisation constant over the pooled training frame...\n")
ruler_fit <- calculate_weighted_di(
  train_data = train[, metric_names],
  test_data = train[1, metric_names],
  weights = weights,
  scaling = scaling,
  norm_sample = NORM_SAMPLE,
  seed = SEED,
  verbose = TRUE
)

## Note on the AOA threshold ####

# The area-of-applicability threshold is NOT derived here, and the reason is worth
# recording. The obvious construction -- each training row's distance to its nearest other
# training row, then the Tukey fence (q75 + 1.5 IQR) -- gives 0.0465 on this training set,
# against a future domain whose median DI is 0.261. It would place essentially the whole
# projection outside the area of applicability, which is not a finding but an artefact:
# with 82k densely packed rows, "distance to nearest training point" measures how DENSE
# the training set is, not how novel a NEW point can be. CAST avoids this by excluding
# same-fold neighbours.
#
# The honest reference distribution is therefore the DI of HELD-OUT rows, measured against
# the training set they were held out of. pl2_evaluate.R produces exactly that, so the
# threshold is computed in pl2_fitErrorProfiles.R where that column lives.

di_ruler <- list(
  weights = weights,
  scaling = scaling,
  train_avg_dist = ruler_fit$train_avg_dist,
  # The metric's variables and the model's, kept apart on purpose: consumers must apply
  # the ruler over `metric_names`, never over the full predictor set.
  metric_names = metric_names,
  predictor_names = predictor_names,
  n_train = nrow(train),
  norm_sample = NORM_SAMPLE,
  source = "pl2_freezeDIRuler.R (plan_EUintegration.md section 1.4)"
)

saveRDS(di_ruler, "output/pl2/di_ruler_production.rds")

## Report ####

cat("\nFrozen DI ruler written to output/pl2/di_ruler_production.rds\n")
cat("  normalisation constant:", round(di_ruler$train_avg_dist, 5), "\n")
cat("  top metric weights:\n")
print(round(head(weights, 8), 5))

tibble(
  quantity = c(
    "training rows",
    "metric variables",
    "normalisation constant",
    "normalisation subsample"
  ),
  value = c(
    format(nrow(train)),
    format(length(metric_names)),
    sprintf("%.5f", di_ruler$train_avg_dist),
    format(NORM_SAMPLE)
  )
) |>
  write_csv("output/pl2/di_ruler_summary.csv", append = FALSE)

## Diagnostics ####

# This script encodes the most contested design decision in the pipeline -- that the
# distance metric runs over projection-dynamic predictors only -- so the evidence for it
# is drawn here rather than left in a log.

diag <- weights_features |>
  mutate(
    status = case_when(
      !dynamic ~ "static (excluded from metric)",
      material ~ "dynamic + material (partition axis pool)",
      TRUE ~ "dynamic, not material"
    ),
    abs_shift = abs(shift_sd)
  )

# 1. Importance against projected shift. The whole dynamic/material scheme is one picture:
#    a predictor is only useful to a projection metric if it BOTH carries signal (y) and
#    moves under the scenario (x). Static predictors sit on the x = 0 spine no matter how
#    important they are, which is exactly the trap slope fell into.
p_quad <- diag |>
  ggplot(aes(x = abs_shift, y = median, colour = status)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey40") +
  geom_point(size = 2.4, alpha = 0.85) +
  ggrepel::geom_text_repel(
    aes(label = feature), size = 2.8, max.overlaps = 12, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "static (excluded from metric)" = "#c1462c",
    "dynamic, not material" = "#c9a227",
    "dynamic + material (partition axis pool)" = "#2b6a3f"
  )) +
  labs(
    title = "Importance is not enough: a metric predictor must also move",
    subtitle = paste(
      "Dashed line = the 0.5 SD materiality bar used to pick the partition axis.",
      "\nStatic predictors (terrain, paleo) sit at zero shift however important they are."
    ),
    x = "|mean projected shift| (SD of the predictor's own spread)",
    y = "variable importance (median over seeds)", colour = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_quad)
ggsave("output/pl2/ruler_importance_vs_shift.png", p_quad, width = 9, height = 6.5, dpi = 150)

# 2. Where the weight mass actually sits. The headline number the restriction rests on.
mass <- diag |>
  group_by(status) |>
  summarise(weight_mass = sum(median), n = n(), .groups = "drop") |>
  mutate(share = weight_mass / sum(weight_mass))

cat("\nWeight mass by status (the restriction's justification):\n")
mass |>
  mutate(share = paste0(round(100 * share, 1), "%"), weight_mass = round(weight_mass, 4)) |>
  as.data.frame() |>
  print(row.names = FALSE)

p_mass <- mass |>
  ggplot(aes(x = "", y = share, fill = status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%s\n%.1f%% of weight\n(%d predictors)", status, 100 * share, n)),
            position = position_stack(vjust = 0.5), size = 3.1, lineheight = 0.95) +
  scale_fill_manual(values = c(
    "static (excluded from metric)" = "#c1462c",
    "dynamic, not material" = "#c9a227",
    "dynamic + material (partition axis pool)" = "#2b6a3f"
  )) +
  coord_flip() +
  labs(
    title = "Share of the unrestricted ruler's weight that cannot register any projected change",
    x = NULL, y = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank())
print(p_mass)
ggsave("output/pl2/ruler_weight_mass.png", p_mass, width = 10, height = 3, dpi = 150)

write_csv(mass, "output/pl2/ruler_weight_mass.csv", append = FALSE)

# sessionInfo ####

sessioninfo::session_info()
