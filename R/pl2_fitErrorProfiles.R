# Fit skill-versus-novelty error profiles ####

# PURPOSE: Turns the cross-validation into reliability curves -- skill as a function of DI, with cluster-bootstrap intervals -- plus the area-of-applicability threshold the map is masked with.

# The step that converts a verdict into a calibration. pl2_evaluate.R reports how well the
# model does on average; this asks how that skill DECAYS as a prediction gets further from
# the training data, so a number can be attached to each future pixel rather than one
# accuracy quoted for the whole map.
#
# WHY THIS IS POSSIBLE FROM THE ORDINARY CV. It was not, until two things changed. The EU
# block was pooled into training and the partitions were cut along bio10, which made the
# folds genuinely extrapolative; and pl2_freezeDIRuler.R put every fold's DI on one axis.
# Measured that way the CV spans the future -- its 90th percentile of DI (0.799) covers
# 96.2% of future cells and 100% of future bog cells, and 107k held-out rows (478 bog) sit
# above the future's median bog novelty of 0.398. That is why the separately designed pl3
# cross-validation was dropped: `predictions_cv_topfeature.csv` already IS the calibration
# dataset.
#
# WHAT IS MEASURED, AND WHAT IS NOT (notebook 2026-08-17, amendment to 2026-06-03). The
# targets are deliberately rank- and threshold-based, not mean-based:
#
#  - P-hat is NOT a probability of presence. `case.wt`/`make.size` mean the forest's
#    terminal-node frequencies estimate a posterior under a re-weighted prior, so there is
#    no value it "should" approach at a true bog. Curves of E[P-hat | DI] would measure the
#    one property the learner was indifferent to.
#  - These are class-conditional LIKELIHOODS, not posteriors. Inverting them to
#    P(bog | P-hat, DI) needs the prevalence of bog in the future domain, which is the
#    unknown being estimated. Nothing here may be read as a per-pixel probability of bog.
#
# So three layers per DI bin: discrimination (one-vs-rest AUC), separation (the full
# class-conditional score distributions, kept as quantiles rather than collapsed to means,
# so that a shrinking gap and a compressing scale stay distinguishable), and operating
# characteristics (per-class recall and false-positive rate at the argmax operating point,
# which is the layer pl2_mapReliability.R actually applies).
#
# BINS ARE SIZED BY BOG COUNT, NOT BY DI WIDTH. Equal-width DI bins put almost every bog in
# the first bin and leave the high-novelty end -- the end the projection sits at -- resting
# on a handful of presences. Bin edges are therefore quantiles of DI *among bog rows*.
#
# INTERVALS ARE CLUSTER BOOTSTRAP OVER FOLD PAIRS. Rows within one (train, test) partition
# pair share a model and a training set, so they are not independent; resampling rows would
# understate the uncertainty badly. The resampling unit is the fold pair.

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

source("R/functions.R")

## Configuration ####

RESPONSE_LEVELS <- c("nonpeat", "otherpeat", "bog")

# Bin count is a compromise: enough bins to see a trend, few enough that the rarest class
# still has a usable count in each. With ~1.6k bog rows per CV pass, 8 bins leaves ~200
# bog rows per bin.
N_BINS <- 8

N_BOOT <- 500
SEED <- 42

## Inputs ####

preds <- read_csv(
  "output/pl2/predictions_cv_topfeature.csv",
  show_col_types = FALSE
) |>
  mutate(
    response = factor(response, levels = RESPONSE_LEVELS),
    pred_class = factor(pred_class, levels = RESPONSE_LEVELS),
    fold = paste(train_partition, test_partition, sep = "-")
  )

ruler <- readRDS("output/pl2/di_ruler_production.rds")

cat("CV rows:", nrow(preds), "| fold pairs:", n_distinct(preds$fold), "\n")
cat(
  "DI (frozen ruler) median/q90/q99/max:",
  paste(round(quantile(preds$DI, c(.5, .9, .99, 1)), 3), collapse = " / "), "\n\n"
)

## Area-of-applicability threshold ####

# Derived from HELD-OUT rows, not from within-training distances. Each CV row's DI is its
# distance to the training set it was excluded from, which is the same quantity a
# prediction cell has; the Tukey fence over that distribution is therefore a threshold in
# the units the map is measured in, and lands at ~1.31, leaving 0.3% of the future domain
# outside. (Computing it from training-to-training distances instead gives 0.0465 -- below
# the future's 5th percentile -- because with 82k dense rows that measures training
# density, not attainable novelty.)
di_q <- quantile(preds$DI, c(0.25, 0.75))
aoa_threshold <- unname(di_q[2] + 1.5 * (di_q[2] - di_q[1]))

cat("AOA threshold (Tukey fence of held-out DI):", round(aoa_threshold, 4), "\n")
cat(
  "  CV rows beyond it:", sum(preds$DI > aoa_threshold),
  sprintf("(%.1f%%)\n\n", 100 * mean(preds$DI > aoa_threshold))
)

## Bin by DI, sized on bog count ####

bog_di <- preds$DI[preds$response == "bog"]
breaks <- unique(quantile(bog_di, probs = seq(0, 1, length.out = N_BINS + 1)))
# Open the ends so no row falls outside a bin.
breaks[1] <- -Inf
breaks[length(breaks)] <- Inf

preds <- preds |>
  mutate(di_bin = cut(DI, breaks = breaks, labels = FALSE, include.lowest = TRUE))

## Per-bin metrics ####

prob_cols <- c("prob_nonpeat", "prob_otherpeat", "prob_bog")

bin_metrics <- function(d) {
  prob <- as.matrix(d[, prob_cols])
  colnames(prob) <- RESPONSE_LEVELS
  m <- calculate_multiclass_metrics(
    predicted_class = d$pred_class,
    true_class = d$response,
    predicted_prob = prob,
    levels = RESPONSE_LEVELS
  )
  # `n_true` is dropped before the pivot on purpose: it differs between classes, so leaving
  # it in makes pivot_wider treat it as an identifier and emit one sparse row per class
  # instead of one row per bin. Per-class counts are recoverable from `separation` below.
  by_class <- m$by_class |>
    mutate(fpr = 1 - specificity) |>
    select(class, recall, fpr, auc)

  bind_cols(
    tibble(
      n = nrow(d),
      n_bog = sum(d$response == "bog"),
      DI_median = median(d$DI),
      accuracy = m$summary$accuracy,
      Gmean_macro = m$summary$Gmean_macro,
      macro_auc = m$summary$macro_auc
    ),
    by_class |>
      pivot_wider(
        names_from = class,
        values_from = c(recall, fpr, auc),
        names_glue = "{.value}_{class}"
      )
  )
}

profiles <- preds |>
  group_by(di_bin) |>
  group_modify(~ bin_metrics(.x)) |>
  ungroup()

cat("Error profiles by DI bin:\n")
profiles |>
  select(di_bin, n, n_bog, DI_median, macro_auc, auc_bog, recall_bog, fpr_bog) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

## Cluster bootstrap over fold pairs ####

# Resample the 20 (train, test) pairs with replacement, recompute every bin. A bin can
# lose a class entirely in a resample, in which case its metrics are NA and drop out of
# the quantiles rather than being scored as zero.
set.seed(SEED)
folds <- unique(preds$fold)
by_fold <- split(preds, preds$fold)

boot <- map_dfr(seq_len(N_BOOT), function(b) {
  draw <- bind_rows(by_fold[sample(folds, length(folds), replace = TRUE)])
  draw |>
    group_by(di_bin) |>
    group_modify(~ bin_metrics(.x)) |>
    ungroup() |>
    mutate(rep = b)
})

ci <- boot |>
  select(di_bin, macro_auc, auc_bog, recall_bog, fpr_bog, Gmean_macro) |>
  pivot_longer(-di_bin, names_to = "metric", values_to = "value") |>
  group_by(di_bin, metric) |>
  summarise(
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

profiles_ci <- profiles |>
  select(di_bin, n, n_bog, DI_median, macro_auc, auc_bog, recall_bog, fpr_bog, Gmean_macro) |>
  pivot_longer(
    c(macro_auc, auc_bog, recall_bog, fpr_bog, Gmean_macro),
    names_to = "metric", values_to = "estimate"
  ) |>
  left_join(ci, by = c("di_bin", "metric"))

cat("\nWith cluster-bootstrap 95% intervals (", N_BOOT, "resamples of fold pairs):\n")
profiles_ci |>
  filter(metric %in% c("macro_auc", "recall_bog")) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

## Separation: class-conditional score distributions ####

# Quantiles rather than means, because under extrapolation both conditional distributions
# compress toward the training marginal; a pair of mean lines cannot distinguish "signal
# lost" from "scale squashed", which is precisely the failure that retired the old target.
separation <- preds |>
  group_by(di_bin, response) |>
  summarise(
    n = n(),
    q10 = quantile(prob_bog, 0.10),
    q25 = quantile(prob_bog, 0.25),
    median = median(prob_bog),
    q75 = quantile(prob_bog, 0.75),
    q90 = quantile(prob_bog, 0.90),
    .groups = "drop"
  )

## Save ####

error_profiles <- list(
  profiles = profiles,
  profiles_ci = profiles_ci,
  separation = separation,
  breaks = breaks,
  aoa_threshold = aoa_threshold,
  n_boot = N_BOOT,
  ruler_train_avg_dist = ruler$train_avg_dist,
  source = "pl2_fitErrorProfiles.R"
)

saveRDS(error_profiles, "output/pl2/error_profiles.rds")
write_csv(profiles, "output/pl2/error_profiles_by_bin.csv", append = FALSE)
write_csv(profiles_ci, "output/pl2/error_profiles_ci.csv", append = FALSE)

cat("\nWrote output/pl2/error_profiles.rds (+ two CSVs)\n")

## Figures ####

p_auc <- profiles_ci |>
  filter(metric %in% c("macro_auc", "auc_bog")) |>
  ggplot(aes(x = DI_median, y = estimate, colour = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5, linetype = 3) +
  geom_vline(xintercept = aoa_threshold, linetype = 2) +
  labs(
    title = "Discrimination decays with novelty",
    subtitle = "Dashed line = area-of-applicability threshold; dotted = chance",
    x = "Dissimilarity index (frozen ruler)", y = "AUC", colour = NULL, fill = NULL
  ) +
  theme_minimal()
print(p_auc)

p_oc <- profiles_ci |>
  filter(metric %in% c("recall_bog", "fpr_bog")) |>
  ggplot(aes(x = DI_median, y = estimate, colour = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = aoa_threshold, linetype = 2) +
  labs(
    title = "Operating characteristics for bog, at the argmax operating point",
    subtitle = "This is the layer projected onto the map",
    x = "Dissimilarity index (frozen ruler)", y = "Rate", colour = NULL, fill = NULL
  ) +
  theme_minimal()
print(p_oc)

p_sep <- separation |>
  ggplot(aes(x = factor(di_bin), fill = response)) +
  geom_boxplot(
    aes(ymin = q10, lower = q25, middle = median, upper = q75, ymax = q90),
    stat = "identity", alpha = 0.7
  ) +
  labs(
    title = "Class-conditional score distributions across novelty",
    subtitle = "Gap AND compression, which a pair of mean lines would hide",
    x = "DI bin", y = "P(bog) score", fill = NULL
  ) +
  theme_minimal()
print(p_sep)

ggsave("output/pl2/reliability_curve_auc.png", p_auc, width = 7, height = 4.5, dpi = 150)
ggsave("output/pl2/reliability_curve_oc.png", p_oc, width = 7, height = 4.5, dpi = 150)
ggsave("output/pl2/reliability_separation.png", p_sep, width = 7, height = 4.5, dpi = 150)

# sessionInfo ####

sessioninfo::session_info()
