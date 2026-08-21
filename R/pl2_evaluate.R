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
predictor_names <- setdiff(
    names(train_data),
    c("response", "partition", "dataset", "x", "y")
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

## Pairwise cross-validation ####

unique_partitions <- sort(unique(train_data$partition[
    !is.na(train_data$partition)
]))
n_partitions <- length(unique_partitions)

# Create all pairwise combinations (test_partition, train_partition)
pairwise_combinations <- expand.grid(
    test_partition = unique_partitions,
    train_partition = unique_partitions
) |>
    filter(test_partition != train_partition) # Exclude same partition

n_combinations <- nrow(pairwise_combinations)

cat("Running pairwise cross-validation with all combinations\n")
cat("  Partitions:", n_partitions, "\n")
cat("  Train-test pairs:", n_combinations, "\n\n")

# Initialize results storage
predictions_all <- list()
metrics_all <- list()
confusion_all <- list()

for (i in 1:n_combinations) {
    test_partition <- pairwise_combinations$test_partition[i]
    train_partition <- pairwise_combinations$train_partition[i]

    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat(
        "CV Iteration",
        i,
        "/",
        n_combinations,
        "- Train partition:",
        train_partition,
        "| Test partition:",
        test_partition,
        "\n"
    )
    cat(paste(rep("=", 80), collapse = ""), "\n\n")

    # Split data - train on single partition, test on another single partition
    train_fold <- train_data |> filter(partition == train_partition)
    test_fold <- train_data |> filter(partition == test_partition)

    cat("Training observations:", nrow(train_fold), "\n")
    print(table(train_fold$dataset, train_fold$response))
    cat("\nTest observations:", nrow(test_fold), "\n")
    print(table(test_fold$dataset, test_fold$response))
    cat("\n")

    ## Calculate weighted Euclidean DI ####

    cat(
        "Calculating weighted Euclidean DI from train partition",
        train_partition,
        "to test partition",
        test_partition,
        "...\n"
    )

    # Per test row, on the frozen ruler: scaling AND normalisation constant both come from
    # the production training frame, so this DI is directly comparable across folds and
    # with the DI of the future projection. The nearest-neighbour distance is still
    # measured against THIS fold's training rows -- that part is what makes it a
    # cross-validation rather than an in-sample distance.
    DI_result <- calculate_weighted_di(
        train_data = train_fold[, names(weights)],
        test_data = test_fold[, names(weights)],
        weights = weights,
        scaling = ruler$scaling,
        train_avg_dist = ruler$train_avg_dist,
        verbose = TRUE
    )

    cat(
        "  Test DI - mean:",
        round(mean(DI_result$DI), 4),
        ", median:",
        round(median(DI_result$DI), 4),
        ", max:",
        round(max(DI_result$DI), 4),
        "\n"
    )
    cat("  Training avg distance:", round(DI_result$train_avg_dist, 4), "\n\n")

    ## Train and evaluate the balanced 3-class forest ####

    cat("Training balanced 3-class RF (ntree =", NTREE, ")...\n")

    fold_train_df <- train_fold |>
        select(response, all_of(predictor_names)) |>
        as.data.frame()
    fold_train_df$response <- droplevels(fold_train_df$response)

    model_brf <- rfsrc(
        formula = response ~ .,
        data = fold_train_df,
        ntree = NTREE,
        case.wt = randomForestSRC:::make.wt(fold_train_df$response),
        sampsize = randomForestSRC:::make.size(fold_train_df$response),
        importance = FALSE
    )

    cat("  Per-tree balanced sample size:", model_brf$sampsize, "\n")

    # Predict on test set
    test_pred <- predict(
        model_brf,
        newdata = test_fold |>
            select(all_of(predictor_names)) |>
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

    ## Metrics, overall and by dataset block ####

    metrics_overall <- calculate_multiclass_metrics(
        predicted_class = test_pred_class,
        true_class = test_fold$response,
        predicted_prob = prob_matrix,
        levels = RESPONSE_LEVELS
    )

    cat("\n  Confusion matrix (rows = actual):\n")
    print(metrics_overall$confusion)
    cat("\n  Per-class metrics:\n")
    print(as.data.frame(metrics_overall$by_class), row.names = FALSE)
    cat("\n  Macro G-mean:", round(metrics_overall$summary$Gmean_macro, 4), "\n\n")

    metrics_by_dataset <- unique(test_fold$dataset) |>
        sort() |>
        map_dfr(function(ds) {
            keep <- test_fold$dataset == ds
            m <- calculate_multiclass_metrics(
                predicted_class = test_pred_class[keep],
                true_class = test_fold$response[keep],
                predicted_prob = prob_matrix[keep, , drop = FALSE],
                levels = RESPONSE_LEVELS
            )
            bind_cols(tibble(dataset = ds), m$summary) |>
                bind_cols(
                    m$by_class |>
                        select(class, recall) |>
                        pivot_wider(
                            names_from = class,
                            values_from = recall,
                            names_prefix = "recall_"
                        )
                )
        })

    cat("  By dataset block:\n")
    print(as.data.frame(metrics_by_dataset), row.names = FALSE)
    cat("\n")

    ## Store predictions ####

    predictions_test <- tibble(
        train_partition = train_partition,
        test_partition = test_partition,
        x = test_fold$x,
        y = test_fold$y,
        dataset = test_fold$dataset,
        response = as.character(test_fold$response),
        pred_class = test_pred_class,
        prob_nonpeat = prob_matrix[, "nonpeat"],
        prob_otherpeat = prob_matrix[, "otherpeat"],
        prob_bog = prob_matrix[, "bog"],
        DI = DI_result$DI
    )

    predictions_all[[i]] <- predictions_test

    ## Store aggregate metrics ####

    metrics_test <- bind_cols(
        tibble(
            train_partition = train_partition,
            test_partition = test_partition,
            dataset = "all"
        ),
        metrics_overall$summary
    ) |>
        bind_cols(
            metrics_overall$by_class |>
                select(class, recall) |>
                pivot_wider(
                    names_from = class,
                    values_from = recall,
                    names_prefix = "recall_"
                )
        ) |>
        bind_rows(
            metrics_by_dataset |>
                mutate(
                    train_partition = train_partition,
                    test_partition = test_partition,
                    .before = 1
                )
        ) |>
        mutate(
            DI_mean = mean(DI_result$DI),
            DI_median = median(DI_result$DI),
            DI_max = max(DI_result$DI),
            train_avg_dist = DI_result$train_avg_dist
        )

    metrics_all[[i]] <- metrics_test

    confusion_all[[i]] <- metrics_overall$confusion |>
        as.data.frame() |>
        mutate(
            train_partition = train_partition,
            test_partition = test_partition,
            .before = 1
        )
}

## Combine and save results ####

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SAVING RESULTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

predictions_combined <- bind_rows(predictions_all)
metrics_combined <- bind_rows(metrics_all)
confusion_combined <- bind_rows(confusion_all)

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

## Does this cross-validation actually reach the projection? ####

# The load-bearing check for the whole architecture. A separate pl3 cross-validation was
# dropped on the grounds that these folds already span the novelty the projection sits at;
# if that is false, the reliability curves fitted downstream are extrapolating and the
# claim collapses. Both distributions are on the frozen ruler, so they are comparable by
# construction.

future_rows <- read_csv(input_file, show_col_types = FALSE) |>
    filter(scenario == "future")
train_rows <- read_csv(input_file, show_col_types = FALSE) |>
    filter(scenario == "current")

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

coverage <- tibble(
    quantile = c("median", "q90", "q99", "max"),
    cv_heldout = quantile(predictions_combined$DI, c(.5, .9, .99, 1)),
    future_all = quantile(di_future, c(.5, .9, .99, 1)),
    future_bog = quantile(di_future_bog, c(.5, .9, .99, 1))
)
cat("\nDI on the frozen ruler -- does the CV span the projection?\n")
coverage |> mutate(across(-quantile, ~ round(.x, 3))) |> as.data.frame() |>
    print(row.names = FALSE)

cat("\nShare of the projection the CV's novelty range covers:\n")
for (q in c(.5, .9, .99)) {
    hi <- quantile(predictions_combined$DI, q)
    cat(sprintf(
        "  CV %2.0f%% quantile = %.3f -> %5.1f%% of future cells, %5.1f%% of future bog cells\n",
        100 * q, hi, 100 * mean(di_future <= hi), 100 * mean(di_future_bog <= hi)
    ))
}

di_long <- bind_rows(
    tibble(DI = predictions_combined$DI, set = "CV held-out rows"),
    tibble(DI = di_future, set = "future projection domain"),
    tibble(DI = di_future_bog, set = "future, at current bog cells")
)

p_cov <- di_long |>
    ggplot(aes(x = DI, fill = set, colour = set)) +
    geom_density(alpha = 0.28, linewidth = 0.7) +
    coord_cartesian(xlim = c(0, quantile(di_long$DI, 0.999))) +
    scale_fill_manual(values = c(
        "CV held-out rows" = "#4c72b0",
        "future projection domain" = "#c1462c",
        "future, at current bog cells" = "#2b6a3f"
    )) +
    scale_colour_manual(values = c(
        "CV held-out rows" = "#4c72b0",
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
    filter(response == "bog") |>
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

## Skill vs novelty ####

# The fold-level view: one point per train-test pair.
plot(
    Gmean_macro ~ DI_mean,
    data = filter(metrics_combined, dataset == "all"),
    xlab = "Mean Dissimilarity Index",
    ylab = "Macro G-mean",
    main = "Skill vs novelty, by fold"
)

# The per-row view, which is the one section 1.3 actually calls for: DI is a per-pixel
# quantity, so skill can be binned over it as finely as the sample supports rather than
# read off one point per fold. Counts are printed with it -- a bin's skill is not
# interpretable without them, and the high-DI bins are the thin ones.
skill_by_di_bin <- predictions_combined |>
    mutate(
        DI_bin = cut(DI, breaks = quantile(DI, probs = seq(0, 1, 0.05)),
                     include.lowest = TRUE)
    ) |>
    group_by(DI_bin) |>
    summarise(
        n = n(),
        DI_mid = median(DI),
        n_bog = sum(response == "bog"),
        accuracy = mean(pred_class == response),
        recall_bog = {
            is_bog <- response == "bog"
            if (any(is_bog)) mean(pred_class[is_bog] == "bog") else NA_real_
        },
        .groups = "drop"
    )

cat("Skill by DI ventile (per-row binning):\n")
print(as.data.frame(skill_by_di_bin), row.names = FALSE)

plot(
    accuracy ~ DI_mid,
    data = skill_by_di_bin,
    type = "b",
    xlab = "Dissimilarity Index (bin median)",
    ylab = "Accuracy",
    main = "Skill vs novelty, per-row binning"
)

## Summary statistics ####

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SUMMARY STATISTICS ACROSS ALL CV FOLDS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

metrics_summary <- metrics_combined |>
    group_by(dataset) |>
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

cat("Performance by dataset block (mean across folds):\n\n")
print(as.data.frame(metrics_summary), row.names = FALSE)

cat("\n\nDissimilarity Index summary:\n\n")
DI_summary <- metrics_combined |>
    filter(dataset == "all") |>
    summarise(
        mean_DI_mean = mean(DI_mean, na.rm = TRUE),
        sd_DI_mean = sd(DI_mean, na.rm = TRUE),
        mean_DI_median = mean(DI_median, na.rm = TRUE),
        mean_DI_max = mean(DI_max, na.rm = TRUE),
        mean_train_avg_dist = mean(train_avg_dist, na.rm = TRUE)
    )

print(as.data.frame(DI_summary), row.names = FALSE)

cat("\n\nPooled confusion matrix across all folds (rows = actual):\n")
confusion_pooled <- confusion_combined |>
    group_by(Actual, Predicted) |>
    summarise(n = sum(Freq), .groups = "drop") |>
    pivot_wider(names_from = Predicted, values_from = n, values_fill = 0)

print(as.data.frame(confusion_pooled), row.names = FALSE)

# sessionInfo ####

sessioninfo::session_info()
