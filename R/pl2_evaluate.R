# Model Evaluation on Test Partition ####

# PURPOSE: Pairwise cross-validation over the partitions: per-class skill, macro G-mean, 3x3 confusion, and per-row DI, reported overall and per dataset block.

# Evaluates model prediction through pairwise CV of partitions, using default
# hyperparameters. All pairwise combinations (20 train-test pairs for 5 partitions).
#
# Two things changed with the EU integration (plan_EUintegration.md):
#
#  1. The BRF-vs-RFQ comparison is gone. RFQ's q* classifier is two-class only, and the
#     response is now 3-class, so there is one learner: a balanced `rfsrc()` using the
#     package's own recipe (case.wt = make.wt(y), sampsize = make.size(y)). One model,
#     one training set, one prediction pass (notebook 2026-06-18). With three classes and
#     a balanced forest there is also no threshold to optimise -- the prediction is the
#     argmax of the probability simplex, which is what the G-mean-optimised threshold was
#     approximating in the binary case.
#  2. Metrics are reported overall AND split by the test row's `dataset`. That split is
#     the dataset-blocked comparison of section 1.5 / section 7.1 arm 1: it says where
#     the pooled model's skill actually lands, rather than letting a good EU score cover
#     a poor Norwegian one or the reverse.
#
# Main outputs (for each train-test pair):
# - Per-row predictions, class probabilities and dissimilarity index for every test cell.
#   DI is per row, never per fold -- fold-level DI summaries in the metrics table are
#   summaries OF this column, and skill can be binned over it at any resolution
#   (plan_EUintegration.md section 1.3).
# - Aggregate per-class recall/precision/specificity/AUC, macro G-mean, and the full
#   3x3 confusion matrix, across the test partition and within each dataset block.
library(readr)
library(dplyr)
library(randomForestSRC)
library(tidyr)
library(purrr)
library(ggplot2)

source("R/functions.R")

## Configuration ####

RESPONSE_LEVELS <- c("nonpeat", "otherpeat", "bog")
NTREE <- 1000

## Load data ####

input_file <- paste0(
    "output/pl2/modeling_frame_regional_partitioned_",
    "topfeature",
    ".csv"
)

cat("Loading data from:", input_file, "\n")
mf <- read_csv(input_file, show_col_types = FALSE)

# Extract training data (current scenario only)
train_data <- mf |>
    filter(scenario == "current") |>
    dplyr::select(-scenario) |> # Keep partition, response, dataset, predictors, x, y
    mutate(response = factor(response, levels = RESPONSE_LEVELS))

cat("Total observations:", nrow(train_data), "\n")
cat(
    "Partitions:",
    paste(sort(unique(train_data$partition)), collapse = ", "),
    "\n\n"
)
print(table(train_data$dataset, train_data$response))
cat("\n")

# Get predictor names. `dataset` is excluded alongside the coordinates and the partition
# label: it is bookkeeping used to split the *reporting*, and a forest allowed to split
# on it could learn "EU" as a shortcut for the warm flank instead of learning climate.
# `envelope_side` joins that list: it is a label describing where a row sits on the
# stratifying axis, so a forest allowed to split on it would be reading the fold design
# rather than the climate.
predictor_names <- setdiff(
    names(train_data),
    c("response", "partition", "envelope_side", "dataset", "x", "y")
)
cat("Number of predictors:", length(predictor_names), "\n\n")

## Load the frozen DI ruler ####

# All three parts of the metric -- weights, scaling, normalisation constant -- come from
# pl2_freezeDIRuler.R, which derives them once from the pooled training frame. Nothing is
# derived per fold here.
#
# This is what makes the cross-validation double as the calibration dataset for
# reliability. Normalising by each fold's own constant (the previous behaviour) left DI
# incomparable between folds -- those constants span a factor of 1.75 -- and incomparable
# with the DI of the future projection, so a skill-vs-novelty curve could not be read off
# it. On one axis the CV covers the future: its 90th percentile of DI (0.799) covers 96.2%
# of future cells and 100% of future bog cells.
ruler <- readRDS("output/pl2/di_ruler_production.rds")
weights <- ruler$weights
stopifnot(all(names(weights) %in% predictor_names))

cat("Frozen DI ruler:", length(weights), "metric variables |",
    "normalisation constant", round(ruler$train_avg_dist, 5), "\n")
cat("Top 5 metric weights:\n")
print(round(head(weights, 5), 5))
cat("\n")

## Cross-validation: two arms ####

# The pairwise arm trains on ONE partition. That is a genuine extrapolation test, but it
# also trains on between 1.8k and 59k rows against a production model that trains on all
# 82,359 -- so its skill confounds "the model is extrapolating" with "the model is
# starved". The two cannot be separated from a single arm, and the gap between pairwise
# macro AUC (~0.75) and production OOB bog AUC (0.967) is exactly the quantity a reader
# will ask about.
#
# The LOPO arm trains on the other k-1 partitions, so it is the closest analogue to the
# production fit that still holds out a whole climate band. Reading the two together:
#
#   LOPO ~ pairwise  -> training size was not the constraint; the pairwise number IS the
#                       extrapolation cost, and is reportable as-is.
#   LOPO >> pairwise -> the pairwise number understates production skill, and the metrics
#                       section should quote LOPO with the gap explained.
#
# LOPO DI is also systematically LOWER, because a test row has k-1 partitions of training
# rows to be near rather than one. That is not a defect: the projection's DI is measured
# against the whole training frame, so LOPO's DI axis is the one that lines up with the
# map's. It is also why pooling the two arms' rows into a single reliability curve is not
# safe -- see pl2_fitErrorProfiles.R.

unique_partitions <- sort(unique(train_data$partition[
    !is.na(train_data$partition)
]))
n_partitions <- length(unique_partitions)

# One fold, whichever arm it belongs to. Everything that defines the learner and the
# metric lives here, so the two arms cannot drift apart: same NTREE, same balanced recipe,
# same frozen ruler, same argmax rule.
run_fold <- function(train_fold, test_fold, arm, train_label, test_label, seed) {
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("[", arm, "] train:", train_label, "| test:", test_label, "\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")

    cat("Training observations:", nrow(train_fold), "\n")
    print(table(train_fold$dataset, train_fold$response))
    cat("\nTest observations:", nrow(test_fold), "\n")
    print(table(test_fold$dataset, test_fold$response))
    cat("\n")

    ## Weighted Euclidean DI, on the frozen ruler ####

    # Per test row. Scaling AND normalisation constant both come from the production
    # training frame, so this DI is comparable across folds and with the DI of the future
    # projection. The nearest-neighbour distance is still measured against THIS fold's
    # training rows -- that part is what makes it a cross-validation rather than an
    # in-sample distance, and it is also why the arms sit on different DI scales.
    DI_result <- calculate_weighted_di(
        train_data = train_fold[, names(weights)],
        test_data = test_fold[, names(weights)],
        weights = weights,
        scaling = ruler$scaling,
        train_avg_dist = ruler$train_avg_dist,
        verbose = TRUE
    )

    cat(
        "  Test DI - mean:", round(mean(DI_result$DI), 4),
        ", median:", round(median(DI_result$DI), 4),
        ", max:", round(max(DI_result$DI), 4), "\n"
    )

    ## Train and evaluate the balanced 3-class forest ####

    cat("Training balanced 3-class RF (ntree =", NTREE, ")...\n")

    fold_train_df <- train_fold |>
        dplyr::select(response, all_of(predictor_names)) |>
        as.data.frame()
    fold_train_df$response <- droplevels(fold_train_df$response)

    # Seeded per fold (rfsrc takes a negative integer), so two runs of unchanged code give
    # identical CV numbers and the run-to-run diff in RUNALL separates code changes from
    # forest noise.
    model_brf <- rfsrc(
        formula = response ~ .,
        data = fold_train_df,
        ntree = NTREE,
        case.wt = randomForestSRC:::make.wt(fold_train_df$response),
        sampsize = randomForestSRC:::make.size(fold_train_df$response),
        importance = FALSE,
        seed = seed
    )

    cat("  Per-tree balanced sample size:", model_brf$sampsize, "\n")

    test_pred <- predict(
        model_brf,
        newdata = test_fold |>
            dplyr::select(all_of(predictor_names)) |>
            as.data.frame()
    )

    # Probability simplex over the classes the fold could learn; a class absent from the
    # train fold gets probability 0 rather than a missing column.
    prob_matrix <- matrix(
        0,
        nrow = nrow(test_fold),
        ncol = length(RESPONSE_LEVELS),
        dimnames = list(NULL, RESPONSE_LEVELS)
    )
    prob_matrix[, colnames(test_pred$predicted)] <- test_pred$predicted

    # No threshold: the balanced forest has already corrected the prior, so the label is
    # the argmax of the simplex.
    test_pred_class <- RESPONSE_LEVELS[max.col(prob_matrix, ties.method = "first")]

    ## Metrics: overall, by dataset block, and by position on the stratifying axis ####

    metrics_overall <- calculate_multiclass_metrics(
        predicted_class = test_pred_class,
        true_class = test_fold$response,
        predicted_prob = prob_matrix,
        levels = RESPONSE_LEVELS
    )

    cat("\n  Confusion matrix (rows = actual):\n")
    print(metrics_overall$confusion)
    cat("\n  Macro G-mean:", round(metrics_overall$summary$Gmean_macro, 4), "\n\n")

    # Same metrics restricted to a subset of the test rows. Used for both the dataset
    # blocks and the envelope sides so the two are computed identically. Returns a `slice`
    # column rather than a grouping-specific one; `slice_type` records which grouping it
    # came from.
    slice_metrics <- function(keep, slice_type, slice_val) {
        if (!any(keep)) {
            return(NULL)
        }
        m <- calculate_multiclass_metrics(
            predicted_class = test_pred_class[keep],
            true_class = test_fold$response[keep],
            predicted_prob = prob_matrix[keep, , drop = FALSE],
            levels = RESPONSE_LEVELS
        )
        recalls <- m$by_class |>
            dplyr::select(class, recall) |>
            pivot_wider(
                names_from = class, values_from = recall, names_prefix = "recall_"
            )
        fprs <- m$by_class |>
            mutate(fpr = 1 - specificity) |>
            dplyr::select(class, fpr) |>
            pivot_wider(names_from = class, values_from = fpr, names_prefix = "fpr_")
        bind_cols(
            tibble(slice_type = slice_type, slice = slice_val),
            m$summary, recalls, fprs
        )
    }

    metrics_slices <- bind_rows(
        slice_metrics(rep(TRUE, nrow(test_fold)), "all", "all"),
        sort(unique(test_fold$dataset)) |>
            map_dfr(function(ds) {
                slice_metrics(test_fold$dataset == ds, "dataset", ds)
            }),
        # Position relative to the presence envelope. `above` holds no presences at all --
        # it is the region beyond the warmest bog, which is where the projection goes --
        # so recall_bog is NA there by construction and fpr_bog is the number that
        # matters: it says whether the model correctly withholds the bog label outside the
        # envelope. `below` is the opposite flank and must not be pooled with it.
        sort(unique(test_fold$envelope_side)) |>
            map_dfr(function(es) {
                slice_metrics(test_fold$envelope_side == es, "envelope_side", es)
            })
    )

    cat("  By slice:\n")
    print(as.data.frame(metrics_slices), row.names = FALSE)
    cat("\n")

    ## Assemble ####

    predictions_test <- tibble(
        arm = arm,
        train_partition = train_label,
        test_partition = test_label,
        x = test_fold$x,
        y = test_fold$y,
        dataset = test_fold$dataset,
        envelope_side = test_fold$envelope_side,
        response = as.character(test_fold$response),
        pred_class = test_pred_class,
        prob_nonpeat = prob_matrix[, "nonpeat"],
        prob_otherpeat = prob_matrix[, "otherpeat"],
        prob_bog = prob_matrix[, "bog"],
        DI = DI_result$DI
    )

    metrics_test <- metrics_slices |>
        mutate(
            arm = arm,
            train_partition = train_label,
            test_partition = test_label,
            .before = 1
        ) |>
        mutate(
            DI_mean = mean(DI_result$DI),
            DI_median = median(DI_result$DI),
            DI_max = max(DI_result$DI),
            n_train = nrow(train_fold)
        )

    confusion_test <- metrics_overall$confusion |>
        as.data.frame() |>
        mutate(
            arm = arm,
            train_partition = train_label,
            test_partition = test_label,
            .before = 1
        )

    # What each fold trained and tested on, by block and class. The tables printed at the
    # top of every fold end up only in the knit; the size gap between the arms (median
    # 5,042 vs 77,317 training rows) is the argument behind which arm calibrates, so it
    # must exist as a file.
    composition_test <- bind_rows(
        train_fold |> count(dataset, response, name = "n") |> mutate(role = "train"),
        test_fold |> count(dataset, response, name = "n") |> mutate(role = "test")
    ) |>
        mutate(
            arm = arm,
            train_partition = train_label,
            test_partition = test_label,
            .before = 1
        )

    list(
        predictions = predictions_test,
        metrics = metrics_test,
        confusion = confusion_test,
        composition = composition_test
    )
}

# Initialize results storage
predictions_all <- list()
metrics_all <- list()
confusion_all <- list()
composition_all <- list()

collect <- function(res) {
    predictions_all[[length(predictions_all) + 1]] <<- res$predictions
    metrics_all[[length(metrics_all) + 1]] <<- res$metrics
    confusion_all[[length(confusion_all) + 1]] <<- res$confusion
    composition_all[[length(composition_all) + 1]] <<- res$composition
}

## Arm 1: pairwise ####

pairwise_combinations <- expand.grid(
    test_partition = unique_partitions,
    train_partition = unique_partitions
) |>
    filter(test_partition != train_partition) # Exclude same partition

cat("Arm 1 -- pairwise\n")
cat("  Partitions:", n_partitions, "\n")
cat("  Train-test pairs:", nrow(pairwise_combinations), "\n\n")

for (i in seq_len(nrow(pairwise_combinations))) {
    test_partition <- pairwise_combinations$test_partition[i]
    train_partition <- pairwise_combinations$train_partition[i]

    collect(run_fold(
        train_fold = train_data |> filter(partition == train_partition),
        test_fold = train_data |> filter(partition == test_partition),
        arm = "pairwise",
        train_label = as.character(train_partition),
        test_label = as.character(test_partition),
        seed = -(1000L + i)
    ))
}

## Arm 2: leave-one-partition-out ####

cat("\nArm 2 -- leave-one-partition-out\n")
cat("  Folds:", n_partitions, "\n\n")

for (p in unique_partitions) {
    collect(run_fold(
        train_fold = train_data |> filter(!is.na(partition), partition != p),
        test_fold = train_data |> filter(partition == p),
        arm = "lopo",
        train_label = paste0("all-but-", p),
        test_label = as.character(p),
        seed = -(2000L + as.integer(p))
    ))
}

## Combine and save results ####

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SAVING RESULTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

predictions_combined <- bind_rows(predictions_all)
metrics_combined <- bind_rows(metrics_all)
confusion_combined <- bind_rows(confusion_all)

# Row accounting. LOPO scores every training row exactly once; pairwise scores every row
# once per training partition other than its own, i.e. k - 1 times. A row missing from
# either arm means a partition was silently dropped somewhere upstream.
scored <- predictions_combined |>
    count(arm, x, y, dataset, response, name = "times")
stopifnot(
    sum(scored$arm == "lopo") == nrow(train_data),
    all(scored$times[scored$arm == "lopo"] == 1L),
    sum(scored$arm == "pairwise") == nrow(train_data),
    all(scored$times[scored$arm == "pairwise"] == n_partitions - 1L)
)
cat("Row accounting: every training row scored once under LOPO and",
    n_partitions - 1L, "times under pairwise\n\n")

# Save predictions
predictions_file <- "output/pl2/predictions_cv_topfeature.csv"
write_csv(predictions_combined, predictions_file)
cat("Saved predictions to:", predictions_file, "\n")
cat("  Rows:", nrow(predictions_combined), "\n\n")

# Save metrics
metrics_file <- "output/pl2/metrics_cv_topfeature.csv"
write_csv(metrics_combined, metrics_file)
cat("Saved metrics to:", metrics_file, "\n")
cat("  Rows:", nrow(metrics_combined), "\n\n")

# Save confusion matrices
confusion_file <- "output/pl2/confusion_cv_topfeature.csv"
write_csv(confusion_combined, confusion_file)
cat("Saved confusion matrices to:", confusion_file, "\n\n")

# Save fold composition
composition_file <- "output/pl2/cv_fold_composition.csv"
write_csv(bind_rows(composition_all), composition_file)
cat("Saved fold composition to:", composition_file, "\n\n")

## Does this cross-validation actually reach the projection? ####

# The load-bearing check for the whole architecture. A separate pl3 cross-validation was
# dropped on the grounds that these folds already span the novelty the projection sits at;
# if that is false, the reliability curves fitted downstream are extrapolating and the
# claim collapses. Both distributions are on the frozen ruler, so they are comparable by
# construction.

future_rows <- mf |> filter(scenario == "future")
train_rows <- train_data

ruled <- function(d) sweep(
    scale(
        as.matrix(d[, names(weights)]),
        center = ruler$scaling[["center"]][names(weights)],
        scale = ruler$scaling[["scale"]][names(weights)]
    ),
    2, pmax(weights, 0), "*"
)

di_future <- FNN::get.knnx(
    ruled(train_rows), ruled(future_rows), k = 1, algorithm = "kd_tree"
)$nn.dist[, 1] / ruler$train_avg_dist

bog_future <- train_rows |>
    filter(response == "bog", dataset == "NO") |>
    select(x, y) |>
    inner_join(future_rows, by = c("x", "y"))
di_future_bog <- FNN::get.knnx(
    ruled(train_rows), ruled(bog_future), k = 1, algorithm = "kd_tree"
)$nn.dist[, 1] / ruler$train_avg_dist

# Per arm, because the arms do NOT share a DI scale: a pairwise test row is measured
# against one partition of training rows, a LOPO row against k-1 of them. The projection
# DI is measured against the whole training frame, so it is the LOPO row of this table
# that has to cover the future -- the pairwise row covering it is not sufficient, and a
# pooled quantile here would hide exactly that distinction.
di_probs <- c(.5, .9, .99, 1)
di_labels <- c("median", "q90", "q99", "max")

# Built per arm with map_dfr rather than a grouped summarise: a column named `quantile`
# would mask stats::quantile() under data masking, and the multi-row-per-group idiom is
# deprecated in current dplyr anyway.
coverage <- sort(unique(predictions_combined$arm)) |>
    map_dfr(function(a) {
        di_a <- predictions_combined$DI[predictions_combined$arm == a]
        tibble(
            arm = a,
            quantile = di_labels,
            cv_heldout = stats::quantile(di_a, di_probs),
            future_all = stats::quantile(di_future, di_probs),
            future_bog = stats::quantile(di_future_bog, di_probs)
        )
    })
cat("\nDI on the frozen ruler -- does the CV span the projection?\n")
coverage |> mutate(across(where(is.numeric), ~ round(.x, 3))) |> as.data.frame() |>
    print(row.names = FALSE)

cat("\nShare of the projection the CV's novelty range covers:\n")
for (a in sort(unique(predictions_combined$arm))) {
    di_a <- predictions_combined$DI[predictions_combined$arm == a]
    for (q in c(.5, .9, .99)) {
        hi <- quantile(di_a, q)
        cat(sprintf(
            "  [%-8s] CV %2.0f%% quantile = %.3f -> %5.1f%% future cells, %5.1f%% future bog cells\n",
            a, 100 * q, hi, 100 * mean(di_future <= hi), 100 * mean(di_future_bog <= hi)
        ))
    }
}

di_long <- bind_rows(
    predictions_combined |>
        transmute(DI, set = paste0("CV held-out rows (", arm, ")")),
    tibble(DI = di_future, set = "future projection domain"),
    tibble(DI = di_future_bog, set = "future, at current bog cells")
)

p_cov <- di_long |>
    ggplot(aes(x = DI, fill = set, colour = set)) +
    geom_density(alpha = 0.28, linewidth = 0.7) +
    coord_cartesian(xlim = c(0, quantile(di_long$DI, 0.999))) +
    scale_fill_manual(values = c(
        "CV held-out rows (pairwise)" = "#4c72b0",
        "CV held-out rows (lopo)" = "#7fa8d6",
        "future projection domain" = "#c1462c",
        "future, at current bog cells" = "#2b6a3f"
    )) +
    scale_colour_manual(values = c(
        "CV held-out rows (pairwise)" = "#4c72b0",
        "CV held-out rows (lopo)" = "#7fa8d6",
        "future projection domain" = "#c1462c",
        "future, at current bog cells" = "#2b6a3f"
    )) +
    labs(
        title = "Does the cross-validation reach the novelty the projection sits at?",
        subtitle = paste(
            "All three on the frozen ruler. If the blue distribution did not cover the",
            "\nred and green ones, the reliability curves would be extrapolating."
        ),
        x = "dissimilarity index", y = "density", fill = NULL, colour = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
print(p_cov)
ggsave("output/pl2/cv_novelty_coverage.png", p_cov, width = 9, height = 5, dpi = 150)
write_csv(coverage, "output/pl2/cv_novelty_coverage.csv", append = FALSE)

# A caution that belongs next to the figure: covering the projection's DI RANGE is not the
# same as probing the axis it travels. The CV's most novel bog cells are cooler, not
# warmer, so the skill-vs-DI curve says little about the warm edge specifically.
cat("\nWhere the CV's most novel BOG cells actually sit (bio10):\n")
predictions_combined |>
    filter(response == "bog", arm == "lopo") |>
    left_join(train_rows |> select(x, y, bio10), by = c("x", "y")) |>
    mutate(di_quartile = ntile(DI, 4)) |>
    group_by(di_quartile) |>
    summarise(
        n = n(), DI_median = round(median(DI), 3),
        bio10_median = round(median(bio10), 2),
        pct_above_17 = round(100 * mean(bio10 > 17), 1), .groups = "drop"
    ) |>
    as.data.frame() |>
    print(row.names = FALSE)

# Skill as a function of novelty is NOT read here. pl2_fitErrorProfiles.R fits it from the
# per-row predictions saved above, on the signed bio10 offset and on DI, with fold-level
# intervals; the summaries below are per-fold means.

## Summary statistics ####

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SUMMARY STATISTICS ACROSS ALL CV FOLDS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

metrics_summary <- metrics_combined |>
    group_by(arm, slice_type, slice) |>
    summarise(
        n_folds = n(),
        mean_Gmean = mean(Gmean_macro, na.rm = TRUE),
        sd_Gmean = sd(Gmean_macro, na.rm = TRUE),
        mean_balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
        mean_macro_auc = mean(macro_auc, na.rm = TRUE),
        mean_recall_bog = mean(recall_bog, na.rm = TRUE),
        mean_recall_otherpeat = mean(recall_otherpeat, na.rm = TRUE),
        mean_recall_nonpeat = mean(recall_nonpeat, na.rm = TRUE),
        .groups = "drop"
    )

cat("Performance by arm and slice (mean across folds):\n\n")
print(as.data.frame(metrics_summary), row.names = FALSE)

# The comparison the second arm exists for. If these rows are close, training-set size was
# not what held the pairwise number down and it is reportable as the extrapolation cost.
# If LOPO is well clear of pairwise, the pairwise number understates the model that is
# actually projected, and the metrics section should quote LOPO.
cat("\n\nArm comparison, overall slice -- extrapolation cost vs training size:\n\n")
metrics_combined |>
    filter(slice == "all") |>
    group_by(arm) |>
    summarise(
        n_folds = n(),
        median_n_train = median(n_train),
        Gmean_macro = round(mean(Gmean_macro, na.rm = TRUE), 4),
        macro_auc = round(mean(macro_auc, na.rm = TRUE), 4),
        recall_bog = round(mean(recall_bog, na.rm = TRUE), 4),
        DI_mean = round(mean(DI_mean, na.rm = TRUE), 4),
        .groups = "drop"
    ) |>
    as.data.frame() |>
    print(row.names = FALSE)

# The tails, which no fold could previously reach. Neither holds a presence, so recall_bog
# is NA by construction and fpr_bog carries the information: whether the model withholds
# the bog label beyond the envelope it was trained in.
cat("\n\nSkill beyond the presence envelope (fpr_bog is the readable column):\n\n")
metrics_combined |>
    filter(slice_type == "envelope_side") |>
    group_by(arm, slice) |>
    summarise(
        n_folds = n(),
        n_rows = sum(n),
        recall_bog = round(mean(recall_bog, na.rm = TRUE), 4),
        fpr_bog = round(mean(fpr_bog, na.rm = TRUE), 4),
        recall_nonpeat = round(mean(recall_nonpeat, na.rm = TRUE), 4),
        .groups = "drop"
    ) |>
    as.data.frame() |>
    print(row.names = FALSE)

cat("\n\nDissimilarity Index summary:\n\n")
DI_summary <- metrics_combined |>
    filter(slice == "all") |>
    group_by(arm) |>
    summarise(
        mean_DI_mean = mean(DI_mean, na.rm = TRUE),
        sd_DI_mean = sd(DI_mean, na.rm = TRUE),
        mean_DI_median = mean(DI_median, na.rm = TRUE),
        mean_DI_max = mean(DI_max, na.rm = TRUE),
        .groups = "drop"
    )

print(as.data.frame(DI_summary), row.names = FALSE)

cat("\n\nPooled confusion matrix across all folds (rows = actual):\n")
confusion_pooled <- confusion_combined |>
    group_by(arm, Actual, Predicted) |>
    summarise(n = sum(Freq), .groups = "drop") |>
    pivot_wider(names_from = Predicted, values_from = n, values_fill = 0)

print(as.data.frame(confusion_pooled), row.names = FALSE)

# sessionInfo ####

sessioninfo::session_info()
