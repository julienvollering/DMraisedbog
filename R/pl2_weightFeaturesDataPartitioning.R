# Feature importance assessment at regional scale for data partitioning ####

# PURPOSE: Fits the production balanced 3-class forest to derive the feature weights used by every weighted distance downstream, flagged dynamic/material.

# Produces the feature weights used by every weighted feature-space distance
# downstream: the partitioning feature (pl2_partitionDataByTopFeature.R) and the DI ruler
# (pl2_freezeDIRuler.R).
#
# Three things changed with the EU integration (plan_EUintegration.md):
#
#  1. The response is a 3-level factor (nonpeat / otherpeat / bog), not binary. RFQ is
#     therefore unavailable -- imbalanced()'s q* classifier is two-class only -- but that
#     is a limit of the q* threshold, not of the package: rfsrc() is natively J-class and
#     the documented balanced-RF recipe (case.wt = make.wt(y), sampsize = make.size(y),
#     under ?rfsrc "imbalanced classification data") is written generically over the
#     number of classes. One model, one training set, one prediction pass
#     (notebook 2026-06-18).
#  2. The frame is pooled over two blocks. `dataset` is a bookkeeping column, never a
#     predictor -- it is the block label the dataset-blocked comparison of section 1.5
#     is built on, and letting the forest split on it would let the model learn "EU" as
#     a shortcut for the warm flank.
#  3. Section 1.4 requires the DI ruler to be frozen across the pl3 sweep, taken once
#     from the *production* fit. So the primary weights here are fitted on the whole
#     pooled frame with the production learner, without absence subsampling: subsampling
#     would make these weights something other than the production fit's. The reps vary
#     only the seed, so `sd` measures the forest's own bootstrap/mtry noise.
#
# Two cross-checks from other architectures (randomForest's balanced forest, a multinomial
# elastic net) are fitted alongside and written under their own `method` labels. They set
# no weight downstream -- every consumer filters on the production label -- but the rank
# agreement between them and the production ruler is the evidence that the ruler is a
# property of the data rather than of one implementation.

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(randomForestSRC)
library(randomForest)
library(glmnet)

source("R/functions.R")

## Configuration ####

RESPONSE_LEVELS <- c("nonpeat", "otherpeat", "bog")
NTREE <- 500

# Reps vary the seed only, so they measure the forest's own bootstrap/mtry noise -- not
# sampling variability in the training set, which is fixed by design here. A handful is
# enough for that, and permutation VI over 85k rows is the expensive part: each rep costs
# 43 variables x NTREE OOB prediction passes over the whole frame.
N_REPS <- 5

# Fewer reps for the two cross-checks. Neither sets the ruler, both are markedly more
# expensive per rep at this row count (single-threaded permutation VI for randomForest,
# a full lambda path x nfolds for cv.glmnet), and cv.glmnet's only between-rep variation
# is its fold assignment.
N_REPS_RF <- 3
N_REPS_GLMNET <- 3
GLMNET_NFOLDS <- 5
GLMNET_NLAMBDA <- 50

# The production learner's label. Downstream scripts default to this string, so it must
# stay attached to the fit whose VI is the frozen ruler.
METHOD_PRODUCTION <- "Balanced Random Forest"

## Read modeling frame ####

mf <- read_csv("output/pl2/modeling_frame_regional.csv")

mf_current <- mf |>
  filter(scenario == "current") |>
  select(-scenario, -x, -y) |>
  mutate(response = factor(response, levels = RESPONSE_LEVELS))

stopifnot(!anyNA(mf_current$response), !anyNA(mf_current$dataset))

# Predictors are everything that is neither the response nor the block label.
predictor_names <- setdiff(names(mf_current), c("response", "dataset"))

# Rows with an incomplete predictor vector cannot contribute to a distance metric and
# would be dropped silently inside the fitting functions; drop them here and say so.
n_incomplete <- sum(!complete.cases(mf_current[, predictor_names]))
if (n_incomplete > 0) {
  cat("Dropping", n_incomplete, "rows with incomplete predictors\n")
  mf_current <- mf_current |>
    filter(complete.cases(mf_current[, predictor_names]))
}

cat(
  "Training rows:",
  nrow(mf_current),
  "over",
  length(predictor_names),
  "predictors\n"
)
mf_current |>
  count(dataset, response) |>
  pivot_wider(names_from = response, values_from = n, values_fill = 0) |>
  print()

## Functions ####

# Balanced 3-class rfsrc -- the production learner (notebook 2026-06-18).
# make.wt()/make.size() are internal but are the package's own documented recipe; they
# generalize over the number of classes, giving a per-tree sample of n_classes x n_rarest
# drawn with class-equalizing probabilities.
get_brf_rfsrc_importance <- function(outcome, predictors, n_reps = N_REPS) {
  train_data <- cbind(outcome = droplevels(outcome), predictors) |>
    as.data.frame()
  wt <- randomForestSRC:::make.wt(train_data$outcome)
  size <- randomForestSRC:::make.size(train_data$outcome)

  map_dfr(
    1:n_reps,
    ~ {
      set.seed(1000 + .x)
      cat("  rfsrc rep", .x, "of", n_reps, "...\n")

      rf <- rfsrc(
        formula = outcome ~ .,
        data = train_data,
        ntree = NTREE,
        case.wt = wt,
        sampsize = size,
        importance = TRUE
      )

      # Keep the per-class columns alongside "all": which features carry the *bog*
      # contrast specifically is the quantity section 1.4 expects to shift along the
      # sweep, and it is free here.
      importance_values <- rf$importance |>
        as.data.frame() |>
        tibble::rownames_to_column("feature") |>
        tibble::as_tibble() |>
        mutate(rep = .x)

      rm(rf)
      gc()

      importance_values
    }
  )
}

summarise_importance <- function(importance_reps) {
  importance_reps |>
    group_by(feature) |>
    summarise(
      median = median(all),
      sd = sd(all),
      .groups = "drop"
    )
}

# Cross-implementation check on the ruler: a different package, a different balancing
# mechanism (per-class sampsize rather than case weights), the same permutation-VI idea.
get_brf_rf_importance <- function(outcome, predictors, n_reps = N_REPS_RF) {
  train_data <- cbind(outcome = droplevels(outcome), predictors) |>
    as.data.frame()
  n_per_class <- min(table(train_data$outcome))
  sampsize <- rep(n_per_class, nlevels(train_data$outcome)) |>
    setNames(levels(train_data$outcome))

  map_dfr(
    1:n_reps,
    ~ {
      set.seed(1000 + .x)
      cat("  randomForest rep", .x, "of", n_reps, "...\n")

      rf <- randomForest(
        formula = outcome ~ .,
        data = train_data,
        ntree = NTREE,
        sampsize = sampsize,
        importance = TRUE,
        replace = TRUE
      )

      importance_values <- tibble(
        rep = .x,
        feature = rownames(rf$importance),
        all = rf$importance[, "MeanDecreaseAccuracy"]
      )

      rm(rf)
      gc()

      importance_values
    }
  )
}

# Multinomial elastic net. Predictors are scaled manually so coefficients are on a
# common standardized scale: glmnet's internal standardization returns coefficients on
# the original scale, where a 1-unit change in bio03 and a 1-unit change in gdd10 are
# not comparable quantities.
get_glmnet_importance <- function(
  outcome,
  predictors,
  n_reps = N_REPS_GLMNET,
  alpha = 0.5 # Elastic Net, between Ridge and Lasso
) {
  y <- droplevels(outcome)
  X_scaled <- predictors |>
    as.matrix() |>
    scale()

  # Class weights so the rare bog class is not simply ignored by the likelihood.
  class_freq <- table(y) / length(y)
  obs_weights <- 1 / as.numeric(class_freq[as.character(y)])

  map_dfr(
    1:n_reps,
    ~ {
      set.seed(1000 + .x)
      cat("  glmnet rep", .x, "of", n_reps, "...\n")

      fit <- cv.glmnet(
        X_scaled,
        y,
        family = "multinomial",
        weights = obs_weights,
        alpha = alpha,
        nfolds = GLMNET_NFOLDS,
        nlambda = GLMNET_NLAMBDA,
        standardize = FALSE # Already scaled manually
      )

      # One coefficient vector per class; a feature's weight is its mean absolute
      # coefficient across classes, so a feature that separates only one class still
      # registers rather than being averaged away by a signed sum.
      coef_mat <- coef(fit, s = "lambda.1se") |>
        map(~ as.matrix(.x)[-1, 1, drop = FALSE]) |> # Drop intercept
        (\(m) do.call(cbind, m))()

      importance_values <- tibble(
        rep = .x,
        feature = rownames(coef_mat),
        all = rowMeans(abs(coef_mat))
      )

      rm(fit, coef_mat)
      gc()

      importance_values
    }
  )
}

## Weights on the pooled frame ####

### Balanced Random Forest (rfsrc) -- production learner, frozen ruler ####

tictoc::tic("BRF (rfsrc), pooled")
importance_brf_rfsrc <- get_brf_rfsrc_importance(
  mf_current$response,
  mf_current[, predictor_names]
)
tictoc::toc()

weights_brf_rfsrc <- summarise_importance(importance_brf_rfsrc)

weights_brf_rfsrc |>
  arrange(desc(median)) |>
  print(n = Inf)

### Balanced Random Forest (randomForest) ####

tictoc::tic("BRF (randomForest), pooled")
importance_brf_rf <- get_brf_rf_importance(
  mf_current$response,
  mf_current[, predictor_names]
)
tictoc::toc()

weights_brf_rf <- summarise_importance(importance_brf_rf)

### Penalized Regression (glmnet) ####

tictoc::tic("Multinomial elastic net, pooled")
importance_glmnet <- get_glmnet_importance(
  mf_current$response,
  mf_current[, predictor_names]
)
tictoc::toc()

weights_glmnet <- summarise_importance(importance_glmnet)

## Which predictors can carry projection novelty at all ####

# The distance metric these weights feed is used to measure how far future Norway sits
# from the training data. A predictor that is identical in the current and future frames
# contributes exactly zero to that distance no matter how much weight it carries -- so
# weight placed on one is not neutral, it is dilution: it inflates the denominator of the
# metric and compresses measured novelty toward zero, biasing the reliability reading in
# the direction that makes the projection look safer than it is.
#
# This is not a hypothetical here. `slope` alone takes 39.9% of the weight mass, and it is
# bit-identical between scenarios; so are `elevation` and all nine `paleo_*` variables,
# which describe deglaciation history. Together that is 55.5% of the ruler, inert.
#
# The per-class importance says the same thing from the other side: slope scores 0.555 on
# the OTHER-PEAT contrast and only 0.033 on the bog contrast, level with bio10. The "all"
# column was being set by the easy contrast, not by the one the research question rests on.
#
# The restriction applies to the METRIC only. All 43 predictors stay in the model.

dynamic_info <- identify_dynamic_predictors(mf, predictor_names)

cat("\nStatic predictors, excluded from the distance metric:\n")
print(dynamic_info$static)
cat("Dynamic predictors retained for the metric:", length(dynamic_info$dynamic), "\n")
cat(
  "Of those, material (mean projected shift >= 0.5 SD), the pool the partition axis is",
  "drawn from:", length(dynamic_info$material), "\n"
)

static_share <- weights_brf_rfsrc |>
  summarise(share = sum(median[feature %in% dynamic_info$static]) / sum(median)) |>
  pull(share)
cat(
  "Share of the production ruler's weight mass on static predictors:",
  paste0(round(100 * static_share, 1), "%\n\n")
)

## Combine results ####

weights_df <- bind_rows(
  weights_brf_rfsrc |> mutate(method = METHOD_PRODUCTION),
  weights_brf_rf |> mutate(method = "Balanced Random Forest (randomForest)"),
  weights_glmnet |> mutate(method = "Penalized Regression (glmnet)")
)

## Ruler agreement ####

# The ruler only has to be stable in *rank* for the partitioning feature and the distance
# metric to mean the same thing; absolute VI scales are not comparable across methods.
rank_agreement <- weights_df |>
  group_by(method) |>
  mutate(rank = rank(-median)) |>
  ungroup() |>
  select(feature, method, rank) |>
  pivot_wider(names_from = method, values_from = rank)

rank_correlations <- rank_agreement |>
  select(-feature) |>
  map_dbl(~ cor(.x, rank_agreement[[METHOD_PRODUCTION]], method = "spearman"))

cat("\nSpearman rank correlation of feature weights, vs the production ruler:\n")
print(round(rank_correlations, 3))

cat("\nTop 10 features by the production ruler, ranked under each method:\n")
rank_agreement |>
  arrange(.data[[METHOD_PRODUCTION]]) |>
  head(10) |>
  print(width = Inf)

## Per-class importance under the production ruler ####

# Reported, not used for weighting. Section 1.4 expects warm-flank variables to gain
# importance as the pl3 cut moves warmer; this is the baseline that shift is read against.
importance_by_class <- importance_brf_rfsrc |>
  select(-rep) |>
  pivot_longer(-feature, names_to = "class", values_to = "importance") |>
  group_by(feature, class) |>
  summarise(median = median(importance), .groups = "drop")

importance_by_class |>
  pivot_wider(names_from = class, values_from = median) |>
  arrange(desc(all)) |>
  head(15) |>
  print(width = Inf)

## Plotting ####

# Features are ordered once, by the production ruler, and that order is held across all
# panels: the question the figure has to answer is whether the other methods agree with
# the ordering the distance metric actually uses, which a per-panel reordering hides.
feature_order <- weights_brf_rfsrc |>
  arrange(median) |>
  pull(feature)

weights_df |>
  mutate(feature = factor(feature, levels = feature_order)) |>
  ggplot(aes(x = feature, y = median)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(
    aes(ymin = pmax(0, median - sd), ymax = median + sd),
    width = 0.2
  ) +
  coord_flip() +
  facet_wrap(~method, scales = "free_x", ncol = 2) +
  labs(
    title = "Feature importance for the weighted distance metric",
    subtitle = "Pooled EU + Norway frame, 3-class response; the production ruler and two cross-implementation checks",
    x = "Feature",
    y = "Importance"
  ) +
  theme_minimal()

# Per-class importance for the production ruler: which contrast each feature carries.
importance_by_class |>
  filter(class != "all") |>
  semi_join(
    weights_brf_rfsrc |> slice_max(median, n = 15) |> select(feature),
    by = "feature"
  ) |>
  mutate(feature = forcats::fct_reorder(feature, median, .fun = max)) |>
  ggplot(aes(x = feature, y = median, fill = class)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Per-class importance, production ruler (top 15 features)",
    x = "Feature",
    y = "Importance",
    fill = "Class"
  ) +
  theme_minimal()

## Diagnostic: is the ruler stable enough to freeze? ####

# The weights are frozen once and reused everywhere downstream, so the question a reviewer
# will ask is whether they are stable enough to bear that. Two things shown together: the
# spread across seeds for each predictor, and whether the predictor can contribute to a
# projection metric at all.

stability <- importance_brf_rfsrc |>
  select(feature, rep, all) |>
  left_join(
    tibble(
      feature = predictor_names,
      status = if_else(
        predictor_names %in% dynamic_info$dynamic,
        if_else(predictor_names %in% dynamic_info$material,
                "dynamic + material", "dynamic, not material"),
        "static"
      )
    ),
    by = "feature"
  ) |>
  mutate(feature = factor(feature, levels = rev(feature_order)))

p_stab <- stability |>
  ggplot(aes(x = feature, y = all, colour = status)) +
  geom_point(alpha = 0.6, size = 1.4) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6) +
  coord_flip() +
  scale_colour_manual(values = c(
    "static" = "#c1462c",
    "dynamic, not material" = "#c9a227",
    "dynamic + material" = "#2b6a3f"
  )) +
  labs(
    title = "Variable importance entering the frozen ruler",
    subtitle = paste(
      "One point per seed; bar = median. Colour marks whether the predictor can",
      "\ncontribute to a projection distance at all (see pl2_freezeDIRuler.R)."
    ),
    x = NULL, y = "permutation importance (all-class)", colour = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_stab)
ggsave(
  "output/pl2/weights_stability.png", p_stab,
  width = 8, height = 9, dpi = 150
)

cat("\nRep-to-rep spread relative to the estimate (top 10 by importance):\n")
weights_brf_rfsrc |>
  arrange(desc(median)) |>
  head(10) |>
  transmute(feature, median = round(median, 5), sd = round(sd, 5),
            cv_pct = round(100 * sd / median, 1)) |>
  as.data.frame() |>
  print(row.names = FALSE)

## Diagnostic: frequency of observed presence along the top features ####

# FOP is a binary-response diagnostic, so it is read here on the bog-vs-rest contrast
# only. It says nothing about the other-peat / non-peat contrast and sets no weight --
# it is here to catch a top-ranked feature whose relationship with bog occurrence is
# flat or an artefact.
training_data_fop <- mf_current |>
  mutate(response = as.numeric(response == "bog")) |>
  select(response, all_of(predictor_names)) |>
  as.data.frame()

top_features <- weights_brf_rfsrc |>
  slice_max(median, n = 4) |>
  arrange(desc(median)) |>
  pull(feature)

cat("\nFOP diagnostic (bog vs rest) for:", paste(top_features, collapse = ", "), "\n")

pairs(training_data_fop[, top_features])
for (feature in top_features) {
  MIAmaxent::plotFOP(training_data_fop, EV = feature)
}

## Saving to file ####

# `dynamic` is written alongside the weights rather than applied to them: this file is
# the record of what the fits found, and the restriction is a property of what the metric
# is FOR. Consumers that build a distance filter on it; consumers that just want VI do not.
weights_df |>
  mutate(
    dynamic = feature %in% dynamic_info$dynamic,
    material = feature %in% dynamic_info$material,
    shift_sd = unname(dynamic_info$shift_sd[feature])
  ) |>
  select(feature, method, median, sd, dynamic, material, shift_sd) |>
  arrange(method, desc(median)) |>
  write_csv("output/pl2/weights_feature_data_partitioning.csv", append = FALSE)

cat("Top 10 of the restricted (projection-dynamic) production ruler:\n")
weights_brf_rfsrc |>
  filter(feature %in% dynamic_info$dynamic) |>
  arrange(desc(median)) |>
  head(10) |>
  print()

importance_by_class |>
  pivot_wider(names_from = class, values_from = median) |>
  arrange(desc(all)) |>
  write_csv("output/pl2/weights_feature_by_class.csv", append = FALSE)

# sessionInfo ####

sessioninfo::session_info()
