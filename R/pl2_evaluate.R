# Model Evaluation on Test Partition ####

# This script evaluates model prediction through pairwise CV of partitions, using default
# hyperparameters. Implements all pairwise combinations (20 train-test pairs for 5 partitions).
#
# It compares two algorithms for imbalanced classification:
# Balanced Random Forest (BRF) and Quantile Random Forest (RFQ).
# For BRF the classification threshold is optimized based on G-mean. For RFQ the implicit
# threshold is used.
#
# Main outputs (for each train-test pair):
# - (Individual) predictions and dissimilarity index (weighted Euclidean distance) for all cells in test partition
# - (Aggregate) TPR, FPR, TNR, FNR, G-mean, across the test partition
library(readr)
library(dplyr)
library(randomForest)
library(randomForestSRC)
library(ROCR)
library(tidyr)
library(purrr)

source("R/functions.R")

## Configuration ####

# Fixed hyperparameters (default values)

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
    dplyr::select(-scenario) # Keep partition, response, predictors, x, y

cat("Total observations:", nrow(train_data), "\n")
cat(
    "Partitions:",
    paste(sort(unique(train_data$partition)), collapse = ", "),
    "\n\n"
)

# Get predictor names (exclude response, partition, x, y)
predictor_names <- setdiff(
    names(train_data),
    c("response", "partition", "x", "y")
)
cat("Number of predictors:", length(predictor_names), "\n\n")

## Load feature weights for trainDI ####

weights_file <- "output/pl2/weights_feature_data_partitioning.csv"
if (file.exists(weights_file)) {
    weights_features <- read_csv(weights_file, show_col_types = FALSE)
    weighting_method <- "Balanced Random Forest"

    weights <- weights_features |>
        filter(method == weighting_method) |>
        select(feature, median) |>
        arrange(desc(median)) |>
        tibble::deframe()

    # Filter to only predictors in modeling frame
    weights <- weights[names(weights) %in% predictor_names]

    cat("Loaded feature weights for", length(weights), "predictors\n")
    cat("Top 5 features by importance:\n")
    print(head(weights, 5))
} else {
    # Equal weights if no weight file exists
    weights <- setNames(rep(1, length(predictor_names)), predictor_names)
    cat("Using equal weights for all predictors (weights file not found)\n")
}
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
    train_fold <- train_data |>
        filter(partition == train_partition) |>
        select(-x, -y)

    test_fold <- train_data |>
        filter(partition == test_partition) |>
        select(-x, -y)

    test_fold_with_coords <- train_data |>
        filter(partition == test_partition)

    # Convert response to factor with explicit levels
    train_fold$response <- factor(
        as.character(train_fold$response),
        levels = c("0", "1")
    )
    test_fold$response <- factor(
        as.character(test_fold$response),
        levels = c("0", "1")
    )

    cat(
        "Training observations:",
        nrow(train_fold),
        "(prevalence =",
        round(mean(as.numeric(as.character(train_fold$response))), 4),
        ")\n"
    )
    cat(
        "Test observations:",
        nrow(test_fold),
        "(prevalence =",
        round(mean(as.numeric(as.character(test_fold$response))), 4),
        ")\n\n"
    )

    ## Calculate weighted Euclidean DI ####

    cat(
        "Calculating weighted Euclidean DI from train partition",
        train_partition,
        "to test partition",
        test_partition,
        "...\n"
    )

    # Calculate DI using custom weighted Euclidean distance function
    DI_result <- calculate_weighted_di(
        train_data = train_fold |> select(-partition, -response),
        test_data = test_fold |> select(-partition, -response),
        weights = weights,
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

    ## Train and evaluate BRF ####

    cat("Training BRF model (ntree = 1000)...\n")

    # Calculate balanced sample size
    n_presence <- sum(train_fold$response == "1")
    n_absence <- sum(train_fold$response == "0")
    balanced_size <- min(n_presence, n_absence)

    cat("  Balanced sampling size:", balanced_size, "\n")

    model_brf <- randomForest::randomForest(
        formula = response ~ .,
        data = train_fold |> select(-partition),
        ntree = 1000,
        sampsize = c(balanced_size, balanced_size),
        replace = TRUE,
        importance = FALSE,
        do.trace = FALSE
    )

    # Get predictions on training data for threshold optimization
    train_pred_prob <- predict(
        model_brf,
        newdata = train_fold |> select(-partition),
        type = "prob"
    )[, "1"]
    train_true <- as.numeric(as.character(train_fold$response))

    # Optimize threshold using G-mean
    thresholds_brf <- calculate_classification_thresholds(
        predicted_prob = train_pred_prob,
        true_labels = train_true,
        method = "gmean"
    )

    cat(
        "  BRF optimal threshold (G-mean):",
        round(thresholds_brf$threshold_default, 4),
        "\n"
    )

    # Predict on test set
    test_pred_prob_brf <- predict(
        model_brf,
        newdata = test_fold |> select(-partition),
        type = "prob"
    )[, "1"]
    test_pred_class_brf <- ifelse(
        test_pred_prob_brf >= thresholds_brf$threshold_default,
        "1",
        "0"
    )

    # Calculate metrics
    metrics_brf <- get_metrics_at_threshold(
        predicted_prob = test_pred_prob_brf,
        true_class = as.character(test_fold$response),
        threshold = thresholds_brf$threshold_default
    )

    cat("  BRF test metrics:\n")
    cat(
        "    TPR:",
        round(metrics_brf$TPR, 4),
        "| FPR:",
        round(metrics_brf$FPR, 4),
        "\n"
    )
    cat(
        "    TNR:",
        round(metrics_brf$TNR, 4),
        "| FNR:",
        round(metrics_brf$FNR, 4),
        "\n"
    )
    cat("    G-mean:", round(metrics_brf$Gmean, 4), "\n\n")

    ## Train and evaluate RFQ ####

    cat("Training RFQ model (ntree = 3000)...\n")

    model_rfq <- randomForestSRC::imbalanced(
        formula = response ~ .,
        data = train_fold |> select(-partition) |> as.data.frame(),
        ntree = 3000,
        importance = FALSE,
        do.trace = FALSE
    )

    # Predict on test set
    test_pred_rfq <- predict(
        model_rfq,
        newdata = test_fold |> select(-partition) |> as.data.frame()
    )
    test_pred_prob_rfq <- test_pred_rfq$predicted[, "1"]
    test_pred_class_rfq <- as.character(test_pred_rfq$class)

    # RFQ uses training prevalence as implicit threshold
    threshold_rfq_implicit <- n_presence / (n_presence + n_absence)
    cat(
        "  RFQ implicit threshold (train prevalence):",
        round(threshold_rfq_implicit, 4),
        "\n"
    )

    # Calculate metrics using model's class predictions
    gmean_result <- calculate_gmean(
        test_pred_class_rfq,
        as.character(test_fold$response)
    )

    # Build full metrics (FPR = 1 - TNR, FNR = 1 - TPR)
    metrics_rfq <- data.frame(
        TPR = gmean_result$tpr,
        FPR = 1 - gmean_result$tnr,
        TNR = gmean_result$tnr,
        FNR = 1 - gmean_result$tpr,
        Gmean = gmean_result$gmean
    )

    cat("  RFQ test metrics:\n")
    cat(
        "    TPR:",
        round(metrics_rfq$TPR, 4),
        "| FPR:",
        round(metrics_rfq$FPR, 4),
        "\n"
    )
    cat(
        "    TNR:",
        round(metrics_rfq$TNR, 4),
        "| FNR:",
        round(metrics_rfq$FNR, 4),
        "\n"
    )
    cat("    G-mean:", round(metrics_rfq$Gmean, 4), "\n\n")

    ## Store predictions ####

    predictions_test <- tibble(
        train_partition = train_partition,
        test_partition = test_partition,
        x = test_fold_with_coords$x,
        y = test_fold_with_coords$y,
        response = as.numeric(as.character(test_fold$response)),
        pred_prob_brf = test_pred_prob_brf,
        pred_class_brf = as.numeric(test_pred_class_brf),
        threshold_brf = thresholds_brf$threshold_default,
        pred_prob_rfq = test_pred_prob_rfq,
        pred_class_rfq = as.numeric(test_pred_class_rfq),
        threshold_rfq = threshold_rfq_implicit,
        DI = DI_result$DI
    )

    predictions_all[[i]] <- predictions_test

    ## Store aggregate metrics ####

    metrics_test <- tibble(
        train_partition = train_partition,
        test_partition = test_partition,
        model = c("BRF", "RFQ"),
        n_test = nrow(test_fold),
        n_presence = sum(test_fold$response == "1"),
        n_absence = sum(test_fold$response == "0"),
        prevalence = sum(test_fold$response == "1") / nrow(test_fold),
        threshold = c(thresholds_brf$threshold_default, threshold_rfq_implicit),
        TPR = c(metrics_brf$TPR, metrics_rfq$TPR),
        FPR = c(metrics_brf$FPR, metrics_rfq$FPR),
        TNR = c(metrics_brf$TNR, metrics_rfq$TNR),
        FNR = c(metrics_brf$FNR, metrics_rfq$FNR),
        Gmean = c(metrics_brf$Gmean, metrics_rfq$Gmean),
        DI_mean = mean(DI_result$DI),
        DI_median = median(DI_result$DI),
        DI_max = max(DI_result$DI),
        train_avg_dist = DI_result$train_avg_dist
    )

    metrics_all[[i]] <- metrics_test
}

## Combine and save results ####

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SAVING RESULTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

predictions_combined <- bind_rows(predictions_all)
metrics_combined <- bind_rows(metrics_all)

plot(
    Gmean ~ DI_mean,
    data = filter(metrics_combined, model == "RFQ"),
    xlab = "Mean Dissimilarity Index",
    ylab = "G-mean",
    main = "G-mean vs Dissimilarity Index (RFQ Model)"
)

plot(
    Gmean ~ DI_mean,
    data = filter(metrics_combined, model == "BRF"),
    xlab = "Mean Dissimilarity Index",
    ylab = "G-mean",
    main = "G-mean vs Dissimilarity Index (BRF Model)"
)

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

## Summary statistics ####

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SUMMARY STATISTICS ACROSS ALL CV FOLDS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

metrics_summary <- metrics_combined |>
    group_by(model) |>
    summarise(
        n_folds = n(),
        mean_Gmean = mean(Gmean, na.rm = TRUE),
        sd_Gmean = sd(Gmean, na.rm = TRUE),
        mean_TPR = mean(TPR, na.rm = TRUE),
        sd_TPR = sd(TPR, na.rm = TRUE),
        mean_FPR = mean(FPR, na.rm = TRUE),
        sd_FPR = sd(FPR, na.rm = TRUE),
        mean_TNR = mean(TNR, na.rm = TRUE),
        sd_TNR = sd(TNR, na.rm = TRUE),
        mean_FNR = mean(FNR, na.rm = TRUE),
        sd_FNR = sd(FNR, na.rm = TRUE),
        .groups = "drop"
    )

cat("Performance metrics summary (mean ± SD):\n\n")
print(metrics_summary, n = Inf)

cat("\n\nDissimilarity Index summary:\n\n")
DI_summary <- metrics_combined |>
    group_by(model) |>
    summarise(
        mean_DI_mean = mean(DI_mean, na.rm = TRUE),
        sd_DI_mean = sd(DI_mean, na.rm = TRUE),
        mean_DI_median = mean(DI_median, na.rm = TRUE),
        mean_DI_max = mean(DI_max, na.rm = TRUE),
        mean_train_avg_dist = mean(train_avg_dist, na.rm = TRUE),
        .groups = "drop"
    ) |>
    distinct()

print(DI_summary)

cat("\n\nModel comparison (BRF vs RFQ):\n")
cat(
    "  BRF mean G-mean:",
    round(metrics_summary$mean_Gmean[metrics_summary$model == "BRF"], 4),
    "±",
    round(metrics_summary$sd_Gmean[metrics_summary$model == "BRF"], 4),
    "\n"
)
cat(
    "  RFQ mean G-mean:",
    round(metrics_summary$mean_Gmean[metrics_summary$model == "RFQ"], 4),
    "±",
    round(metrics_summary$sd_Gmean[metrics_summary$model == "RFQ"], 4),
    "\n\n"
)

# sessionInfo ####

sessioninfo::session_info()
