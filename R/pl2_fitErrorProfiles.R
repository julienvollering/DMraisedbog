# Fit skill-versus-novelty error profiles ####

# PURPOSE: Turns the cross-validation into reliability curves -- skill as a function of signed bio10 offset from the training bogs, with cluster-bootstrap intervals -- plus the area-of-applicability threshold the map is masked with.

# The step that converts a verdict into a calibration. pl2_evaluate.R reports how well the
# model does on average; this asks how that skill DECAYS as a prediction moves away from the
# climate the training bogs occupy, so a number can be attached to each future pixel rather
# than one accuracy quoted for the whole map.
#
# WHY THIS IS POSSIBLE FROM THE ORDINARY CV. It was not, until two things changed. The EU
# block was pooled into training and the partitions were cut along bio10, which made the
# folds genuinely extrapolative; and the partition tiling was completed so that rows beyond
# the warmest training bog stay in the folds instead of being dropped. That is why the
# separately designed pl3 cross-validation was retired:
# `predictions_cv_topfeature.csv` already IS the calibration dataset.
#
# WHICH ARM CALIBRATES (CV_ARM below; notebook 2026-09-14). The predictions file holds two
# arms and they do NOT share a novelty scale -- pairwise trains on one partition, LOPO on k-1
# -- so each is profiled separately, never pooled. PAIRWISE calibrates, because it is the arm
# that simulates the extrapolation the projection makes. A model trained on one partition and
# tested on a distant one puts held-out bogs far from their training bogs: around the Lyngstad
# bogs' future offset (+2.5 to +4 C) pairwise holds 70 bog rows from 4 models, LOPO 9 from 1.
# Pairwise also replicates: every offset bin is scored by 4-6 models, where LOPO's are scored
# by 1-2 (see the next block), so only pairwise can put an interval on the curve.
#
# THE PRICE, stated rather than hidden. Pairwise measures the skill of a model trained on a
# median 5,042 rows, not the production fit's 82,359, and that costs skill: macro AUC 0.736
# vs 0.788, macro G-mean 0.371 vs 0.492 against LOPO on the inside-envelope slice. The map
# therefore reads CONSERVATIVE -- expected skill is a floor on production's, not an estimate
# of it. Borrowing only the interval width from pairwise and centring it on LOPO was
# considered and rejected: the two arms disagree on the SHAPE of the curve, not just its
# level (interior bog recall ~0.45 vs ~0.95; coldest bin 0.43 vs 0.12), so a pairwise band
# wrapped around a LOPO estimate would describe neither. LOPO is kept as COMPARE_ARM and
# drawn on the figures, so the gap between a small and a production-sized model stays visible.
#
# WHY THE AXIS IS SIGNED bio10 OFFSET AND NOT DI (notebook 2026-09-12). DI was the axis
# until it was measured properly, and then it failed in two ways at once. (i) NO RESOLUTION
# WHERE THE MAP IS. LOPO bog rows pile up at low DI -- with k-1 partitions in training
# nearly every held-out bog has a close analogue -- so the top bog-count bin began at DI
# 0.187 while the whole future domain sat above it. Every projected cell read one bin, whose
# bootstrap interval on recall_bog was 0.000-0.961. (ii) IT MIXES OPPOSITE CASES. DI is an
# unsigned scalar over 32 dimensions, so a cell that is oddly COLD and a cell that is oddly
# WARM score alike, and this model treats them oppositely: of the 110 LOPO bog rows inside
# the future's DI window, the 81 from the cold-end fold scored recall 0.099 and the 29 from
# the warm-end fold scored 1.000. Pooling them produced 0.336, a number describing neither.
#
# Signed offset fixes both because it is one-dimensional, SIGNED, and bog-referenced:
#
#   offset = row bio10 - median bio10 of the BOG rows in that fold's TRAINING set
#
# Read directly: "this cell is +3.5 C warmer than the typical bog this model learned from."
# Each fold is measured against what that fold actually saw, exactly as DI was, so the folds
# pool onto one axis. The production reference (15.528 C) is saved here and reused by
# pl2_mapReliability.R, so CV offsets and map offsets share an origin. The training bogs span
# -6.03 to +3.53 on it; the future Norwegian domain spans about -5.4 to +3.6, and future bog
# cells sit at roughly +3.5 -- i.e. AT the warm envelope edge, which is the fact the DI axis
# was hiding.
#
# DI IS STILL COMPUTED AND STILL PROFILED, as `profiles_di`, for two reasons: the
# area-of-applicability threshold is a multivariate question and stays on DI, and the flat DI
# curve is itself evidence (discrimination does not decay with generic novelty). It is no
# longer what the map reads.
#
# WHAT IS MEASURED, AND WHAT IS NOT (notebook 2026-08-17, amendment to 2026-06-03). The
# targets are deliberately rank- and threshold-based, not mean-based:
#
#  - P-hat is NOT a probability of presence. `case.wt`/`make.size` mean the forest's
#    terminal-node frequencies estimate a posterior under a re-weighted prior, so there is
#    no value it "should" approach at a true bog. Curves of E[P-hat | axis] would measure
#    the one property the learner was indifferent to.
#  - These are class-conditional LIKELIHOODS, not posteriors. Inverting them to
#    P(bog | P-hat, axis) needs the prevalence of bog in the future domain, which is the
#    unknown being estimated. Nothing here may be read as a per-pixel probability of bog.
#
# So three layers per bin: discrimination (one-vs-rest AUC), separation (the full
# class-conditional score distributions, kept as quantiles rather than collapsed to means, so
# that a shrinking gap and a compressing scale stay distinguishable), and operating
# characteristics (per-class recall and false-positive rate at the argmax operating point,
# which is the layer pl2_mapReliability.R actually applies).
#
# BINS ARE SIZED BY BOG COUNT, NOT BY AXIS WIDTH. Equal-width bins leave the end the
# projection sits at resting on a handful of presences. Bin edges are therefore quantiles of
# the axis *among bog rows*. On a signed axis that costs interpretability -- the bins are no
# longer round numbers of degrees -- so each bin reports its own offset range alongside the
# median the map interpolates on.
#
# INTERVALS ARE CLUSTER BOOTSTRAP OVER FOLDS. Rows within one fold share a model and a
# training set, so they are not independent; resampling rows would understate the uncertainty
# badly. The resampling unit is the fold. Note LOPO offers only k = 5 units against
# pairwise's 20, so its intervals are both wider and coarser -- the honest price of folds
# whose training sets resemble production, not a defect to tune away.
#
# WHY LOPO CANNOT CARRY INTERVALS ON THIS AXIS (`n_folds` below). An interval asks how much a
# bin's skill would change under a different set of models, which can only be estimated by
# comparing several models scored in that bin. LOPO's partitions were cut on bio10 and so are
# the offset bins, so a bin and a fold are nearly the same object: under LOPO, bins 1, 3, 6
# and 8 were scored by exactly one model and none by more than two. A single-model bin has
# nothing to resample -- every bootstrap draw is the same rows duplicated, and recall and AUC
# are invariant to duplication -- so the bootstrap returns a ZERO-WIDTH interval, which reads
# as certainty when it is the opposite. Pairwise does not have this problem, because each test
# partition is scored by k-1 models trained on different partitions, each of which puts the
# same rows at a different offset. The guard below stays regardless: any bin with fewer than
# MIN_FOLDS_FOR_CI models reports NA rather than a fabricated band.

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

source("R/functions.R")

## Configuration ####

RESPONSE_LEVELS <- c("nonpeat", "otherpeat", "bog")

# The CV arm the curves are fitted on, and the arm drawn alongside it for reference. See
# WHICH ARM CALIBRATES above.
CV_ARM <- "pairwise"
COMPARE_ARM <- "lopo"

# The axis the map reads. See WHY THE AXIS IS SIGNED bio10 OFFSET above.
PROFILE_AXIS <- "offset"

# The predictor the partitions were cut on, and therefore the axis the projection travels.
AXIS_FEATURE <- "bio10"

# Bin count is a compromise: enough bins to see a trend, few enough that the rarest class
# still has a usable count in each. Pairwise scores each of ~1.6k bog cells four times (once
# per model that did not train on it), so 8 bins leaves ~800 bog rows per bin; the reference
# LOPO arm, one score per cell, gets ~200.
N_BINS <- 8

N_BOOT <- 500
SEED <- 42

# A fold-level bootstrap needs at least two folds in a bin to say anything. Below this the
# interval is suppressed rather than reported as zero width. See the header.
MIN_FOLDS_FOR_CI <- 2

## Inputs ####

# The predictions must postdate the frame they were scored on, and the ruler the weights it
# was frozen from; otherwise the curves below describe a run that no longer exists (see
# assert_fresher() in R/functions.R for the incident behind this).
assert_fresher(
  "output/pl2/predictions_cv_topfeature.csv",
  "output/pl2/modeling_frame_regional_partitioned_topfeature.csv"
)
assert_fresher(
  "output/pl2/predictions_cv_topfeature.csv",
  "output/pl2/di_ruler_production.rds"
)
assert_fresher(
  "output/pl2/di_ruler_production.rds",
  "output/pl2/weights_feature_data_partitioning.csv"
)

# Both arms are read, because the reference arm is profiled too. They are split apart before
# anything is pooled.
preds <- read_csv(
  "output/pl2/predictions_cv_topfeature.csv",
  show_col_types = FALSE
) |>
  filter(arm %in% c(CV_ARM, COMPARE_ARM)) |>
  mutate(
    response = factor(response, levels = RESPONSE_LEVELS),
    pred_class = factor(pred_class, levels = RESPONSE_LEVELS),
    fold = paste(train_partition, test_partition, sep = "-")
  )

ruler <- readRDS("output/pl2/di_ruler_production.rds")

stopifnot(nrow(preds) > 0)

# The axis feature and the partition labels are not carried in the predictions file, so they
# come back from the partitioned frame. Joined on the full row key; the row count must not
# move, which asserts the key is unique on both sides.
frame <- read_csv(
  "output/pl2/modeling_frame_regional_partitioned_topfeature.csv",
  col_select = all_of(c("scenario", "response", "dataset", "x", "y", "partition", AXIS_FEATURE)),
  show_col_types = FALSE
) |>
  filter(scenario == "current") |>
  select(-scenario)

n_before <- nrow(preds)
preds <- preds |>
  left_join(frame, by = c("x", "y", "dataset", "response"))
stopifnot(nrow(preds) == n_before, !anyNA(preds[[AXIS_FEATURE]]))

## Signed offset from each fold's training bogs ####

# Per fold, the reference is the median axis value of the BOG rows the fold trained on --
# not of all its training rows, and not of the whole dataset. `train_partition` is an integer
# for pairwise folds and "all-but-k" for LOPO folds, so the training partitions are recovered
# from the label rather than assumed.
bog_axis <- frame |>
  filter(response == "bog") |>
  select(partition, value = all_of(AXIS_FEATURE))

train_partitions <- function(label) {
  if (grepl("^all-but-", label)) {
    setdiff(sort(unique(bog_axis$partition)), as.integer(sub("^all-but-", "", label)))
  } else {
    as.integer(label)
  }
}

fold_ref <- preds |>
  distinct(arm, fold, train_partition) |>
  mutate(ref = map_dbl(
    train_partition,
    ~ median(bog_axis$value[bog_axis$partition %in% train_partitions(.x)])
  ))

# A fold whose training set held no bog rows would get an NA reference and every offset
# it scored would vanish from the bins without a word. Every fold has one by design.
stopifnot(all(is.finite(fold_ref$ref)))

# The production reference, for pl2_mapReliability.R: the median over ALL training bogs,
# which is what the production fit saw. Frozen here so CV and map offsets share an origin.
axis_ref_production <- median(bog_axis$value)

preds <- preds |>
  left_join(fold_ref |> select(arm, fold, ref), by = c("arm", "fold")) |>
  mutate(offset = .data[[AXIS_FEATURE]] - ref)

compare <- preds |> filter(arm == COMPARE_ARM)
preds <- preds |> filter(arm == CV_ARM)

cat("CV arm:", CV_ARM, "| rows:", nrow(preds),
    "| folds:", n_distinct(preds$fold),
    "|| reference arm:", COMPARE_ARM, "| rows:", nrow(compare), "\n")
cat("Axis:", PROFILE_AXIS, "on", AXIS_FEATURE,
    "| production reference:", round(axis_ref_production, 3), "\n")
cat("Per-fold reference (", CV_ARM, "):",
    paste(sprintf("%s=%.2f", fold_ref$fold[fold_ref$arm == CV_ARM],
                  fold_ref$ref[fold_ref$arm == CV_ARM]), collapse = " "), "\n")
cat("Held-out offset (q05/median/q95/max):",
    paste(round(quantile(preds$offset, c(.05, .5, .95, 1)), 2), collapse = " / "), "\n")
cat("  among bog rows:",
    paste(round(quantile(preds$offset[preds$response == "bog"], c(.05, .5, .95, 1)), 2),
          collapse = " / "), "\n")
cat("DI (frozen ruler) median/q90/q99/max:",
    paste(round(quantile(preds$DI, c(.5, .9, .99, 1)), 3), collapse = " / "), "\n\n")

## Area-of-applicability threshold ####

# Stays on DI, because applicability is a multivariate question: a cell can be unremarkable
# on bio10 and still sit outside the training data on seasonality or paleo-climate.
#
# Derived from HELD-OUT rows, not from within-training distances. Each CV row's DI is its
# distance to the training set it was excluded from, which is the same quantity a prediction
# cell has; the Tukey fence over that distribution is therefore a threshold in the units the
# map is measured in. (Computing it from training-to-training distances instead gives 0.0465
# -- below the future's 5th percentile -- because with 82k dense rows that measures training
# density, not attainable novelty.)
di_q <- quantile(preds$DI, c(0.25, 0.75))
aoa_threshold <- unname(di_q[2] + 1.5 * (di_q[2] - di_q[1]))

cat("AOA threshold (Tukey fence of held-out DI):", round(aoa_threshold, 4), "\n")
cat(
  "  CV rows beyond it:", sum(preds$DI > aoa_threshold),
  sprintf("(%.1f%%)\n\n", 100 * mean(preds$DI > aoa_threshold))
)

## Per-bin metrics ####

prob_cols <- c("prob_nonpeat", "prob_otherpeat", "prob_bog")

bin_metrics <- function(d, axis) {
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

  v <- d[[axis]]
  bind_cols(
    tibble(
      n = nrow(d),
      n_bog = sum(d$response == "bog"),
      # How many independently fitted models this bin's rows came from. On the offset axis
      # this is usually 1 or 2, which is what makes the interval unestimable -- see header.
      n_folds = n_distinct(d$fold),
      axis_median = median(v),
      axis_lo = min(v),
      axis_hi = max(v),
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

# Bin edges are quantiles of the axis among bog rows, ends opened so no row falls outside.
axis_breaks <- function(d, axis) {
  b <- unique(quantile(
    d[[axis]][d$response == "bog"],
    probs = seq(0, 1, length.out = N_BINS + 1)
  ))
  b[1] <- -Inf
  b[length(b)] <- Inf
  b
}

BOOT_METRICS <- c("macro_auc", "auc_bog", "recall_bog", "fpr_bog", "Gmean_macro")

# One arm and one axis in, the three layers out. Called for the calibrating arm on the offset
# axis (what the map reads), on DI (retained comparison), and for the reference arm on the
# offset axis without a bootstrap, since its bins are mostly single-model anyway.
fit_profiles <- function(data, axis, boot = TRUE) {
  breaks <- axis_breaks(data, axis)
  d <- data |>
    mutate(bin = cut(.data[[axis]], breaks = breaks, labels = FALSE, include.lowest = TRUE))

  profiles <- d |>
    group_by(bin) |>
    group_modify(~ bin_metrics(.x, axis)) |>
    ungroup()

  if (boot) {
    # Resample the arm's folds with replacement, recompute every bin. A bin can lose a class
    # entirely in a resample, in which case its metrics are NA and drop out of the quantiles
    # rather than being scored as zero.
    set.seed(SEED)
    folds <- unique(d$fold)
    by_fold <- split(d, d$fold)

    reps <- map_dfr(seq_len(N_BOOT), function(b) {
      draw <- bind_rows(by_fold[sample(folds, length(folds), replace = TRUE)])
      draw |>
        group_by(bin) |>
        group_modify(~ bin_metrics(.x, axis)) |>
        ungroup() |>
        mutate(rep = b)
    })

    ci <- reps |>
      select(bin, all_of(BOOT_METRICS)) |>
      pivot_longer(-bin, names_to = "metric", values_to = "value") |>
      group_by(bin, metric) |>
      summarise(
        lower = quantile(value, 0.025, na.rm = TRUE),
        upper = quantile(value, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    ci <- tidyr::expand_grid(bin = profiles$bin, metric = BOOT_METRICS) |>
      mutate(lower = NA_real_, upper = NA_real_)
  }

  # Suppress the interval where the bin had too few independent models for a fold-level
  # resample to mean anything; a zero-width band would read as certainty.
  profiles_ci <- profiles |>
    select(bin, n, n_bog, n_folds, axis_median, axis_lo, axis_hi, all_of(BOOT_METRICS)) |>
    pivot_longer(all_of(BOOT_METRICS), names_to = "metric", values_to = "estimate") |>
    left_join(ci, by = c("bin", "metric")) |>
    mutate(across(c(lower, upper), ~ if_else(n_folds < MIN_FOLDS_FOR_CI, NA_real_, .x)))

  # Quantiles rather than means, because under extrapolation both conditional distributions
  # compress toward the training marginal; a pair of mean lines cannot distinguish "signal
  # lost" from "scale squashed", which is precisely the failure that retired the old target.
  separation <- d |>
    group_by(bin, response) |>
    summarise(
      n = n(),
      q10 = quantile(prob_bog, 0.10),
      q25 = quantile(prob_bog, 0.25),
      median = median(prob_bog),
      q75 = quantile(prob_bog, 0.75),
      q90 = quantile(prob_bog, 0.90),
      .groups = "drop"
    )

  list(
    profiles = profiles, profiles_ci = profiles_ci, separation = separation,
    breaks = breaks, axis = axis
  )
}

fit_offset <- fit_profiles(preds, PROFILE_AXIS)
fit_di <- fit_profiles(preds, "DI")
fit_compare <- fit_profiles(compare, PROFILE_AXIS, boot = FALSE)

cat("Error profiles by", PROFILE_AXIS, "bin:\n")
fit_offset$profiles |>
  select(bin, n, n_bog, n_folds, axis_lo, axis_hi, axis_median,
         macro_auc, auc_bog, recall_bog, fpr_bog) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

cat("\nWith cluster-bootstrap 95% intervals (", N_BOOT, "resamples of folds;",
    "NA where n_folds <", MIN_FOLDS_FOR_CI, "):\n")
fit_offset$profiles_ci |>
  filter(metric %in% c("macro_auc", "recall_bog")) |>
  select(bin, n_bog, n_folds, axis_median, metric, estimate, lower, upper) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

cat("\nBins with an estimable interval:",
    sum(fit_offset$profiles$n_folds >= MIN_FOLDS_FOR_CI), "of", nrow(fit_offset$profiles), "\n")

cat("\nRetained DI profiles, for comparison:\n")
fit_di$profiles |>
  select(bin, n, n_bog, axis_median, macro_auc, auc_bog, recall_bog, fpr_bog) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

cat("\nReference arm (", COMPARE_ARM, ") on the offset axis, point estimates only:\n")
fit_compare$profiles |>
  select(bin, n_bog, n_folds, axis_median, macro_auc, auc_bog, recall_bog, fpr_bog) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

## Save ####

error_profiles <- list(
  profiles = fit_offset$profiles,
  profiles_ci = fit_offset$profiles_ci,
  separation = fit_offset$separation,
  breaks = fit_offset$breaks,
  profiles_di = fit_di$profiles,
  profiles_ci_di = fit_di$profiles_ci,
  separation_di = fit_di$separation,
  breaks_di = fit_di$breaks,
  profiles_compare = fit_compare$profiles,
  compare_arm = COMPARE_ARM,
  axis = PROFILE_AXIS,
  axis_feature = AXIS_FEATURE,
  axis_ref_production = axis_ref_production,
  fold_ref = fold_ref,
  aoa_threshold = aoa_threshold,
  cv_arm = CV_ARM,
  n_boot = N_BOOT,
  ruler_train_avg_dist = ruler$train_avg_dist,
  source = "pl2_fitErrorProfiles.R"
)

saveRDS(error_profiles, "output/pl2/error_profiles.rds")
write_csv(fit_offset$profiles, "output/pl2/error_profiles_by_bin.csv", append = FALSE)
write_csv(fit_offset$profiles_ci, "output/pl2/error_profiles_ci.csv", append = FALSE)
write_csv(fit_di$profiles, "output/pl2/error_profiles_by_bin_di.csv", append = FALSE)
write_csv(fit_compare$profiles, "output/pl2/error_profiles_by_bin_compare.csv", append = FALSE)

cat("\nWrote output/pl2/error_profiles.rds (+ four CSVs)\n")

## Figures ####

# The warm edge of the training bogs, in offset units: the point beyond which no fold had a
# labelled bog to learn from. Worth a line on every panel.
bog_edge <- max(bog_axis$value) - axis_ref_production

axis_lab <- sprintf("%s offset from training bogs (C)", AXIS_FEATURE)

# The reference arm as thin dashed lines in the same colours: where a production-sized model
# sits against the small-model curve the map reads.
compare_long <- function(metrics) {
  fit_compare$profiles |>
    select(axis_median, all_of(metrics)) |>
    pivot_longer(all_of(metrics), names_to = "metric", values_to = "estimate")
}
arm_caption <- sprintf("Bands and solid lines: %s (read by the map). Thin dashed: %s, for reference.",
                       CV_ARM, COMPARE_ARM)

p_auc <- fit_offset$profiles_ci |>
  filter(metric %in% c("macro_auc", "auc_bog")) |>
  ggplot(aes(x = axis_median, y = estimate, colour = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_line(data = compare_long(c("macro_auc", "auc_bog")), linetype = 2, linewidth = 0.5) +
  geom_hline(yintercept = 0.5, linetype = 3) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_vline(xintercept = bog_edge, linetype = 2, colour = "grey40") +
  labs(
    title = "Discrimination against warming beyond the training bogs",
    subtitle = "Grey solid = the typical training bog; grey dashed = the warmest; dotted = chance",
    caption = arm_caption,
    x = axis_lab, y = "AUC", colour = NULL, fill = NULL
  ) +
  theme_minimal()
print(p_auc)

p_oc <- fit_offset$profiles_ci |>
  filter(metric %in% c("recall_bog", "fpr_bog")) |>
  ggplot(aes(x = axis_median, y = estimate, colour = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_line(data = compare_long(c("recall_bog", "fpr_bog")), linetype = 2, linewidth = 0.5) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_vline(xintercept = bog_edge, linetype = 2, colour = "grey40") +
  labs(
    title = "Operating characteristics for bog, at the argmax operating point",
    subtitle = "This is the layer projected onto the map",
    caption = arm_caption,
    x = axis_lab, y = "Rate", colour = NULL, fill = NULL
  ) +
  theme_minimal()
print(p_oc)

p_sep <- fit_offset$separation |>
  left_join(fit_offset$profiles |> select(bin, axis_median), by = "bin") |>
  mutate(bin_lab = sprintf("%+.1f", axis_median)) |>
  ggplot(aes(x = reorder(bin_lab, axis_median), fill = response)) +
  geom_boxplot(
    aes(ymin = q10, lower = q25, middle = median, upper = q75, ymax = q90),
    stat = "identity", alpha = 0.7
  ) +
  labs(
    title = "Class-conditional score distributions across the offset axis",
    subtitle = "Gap AND compression, which a pair of mean lines would hide",
    x = axis_lab, y = "P(bog) score", fill = NULL
  ) +
  theme_minimal()
print(p_sep)

p_di <- fit_di$profiles_ci |>
  filter(metric %in% c("macro_auc", "recall_bog")) |>
  ggplot(aes(x = axis_median, y = estimate, colour = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5, linetype = 3) +
  geom_vline(xintercept = aoa_threshold, linetype = 2) +
  labs(
    title = "The retired axis, kept as evidence",
    subtitle = "Discrimination is flat in generic novelty; the operating point is not",
    x = "Dissimilarity index (frozen ruler)", y = "Rate", colour = NULL, fill = NULL
  ) +
  theme_minimal()
print(p_di)

ggsave("output/pl2/reliability_curve_auc.png", p_auc, width = 7, height = 4.5, dpi = 150)
ggsave("output/pl2/reliability_curve_oc.png", p_oc, width = 7, height = 4.5, dpi = 150)
ggsave("output/pl2/reliability_separation.png", p_sep, width = 7, height = 4.5, dpi = 150)
ggsave("output/pl2/reliability_curve_di.png", p_di, width = 7, height = 4.5, dpi = 150)

# sessionInfo ####

sessioninfo::session_info()
