# Model Evaluation on Test Partition ####

# This script trains final models on pooled CV partitions (2-5) using optimal
# hyperparameters from tuning, then evaluates performance on held-out test
# partition (1) across three partitioning schemes (importance, pca, maxshift).
#
# Evaluation includes metrics at three different classification thresholds:
# 1. Default/Optimal: G-mean optimal (BRF) or implicit threshold (RFQ)
# 2. High sensitivity: Threshold set to achieve 95% sensitivity on training data
# 3. High specificity: Threshold set to achieve 95% specificity on training data
#
# Metrics calculated: TPR (sensitivity), FPR, TNR (specificity), FNR, G-mean

library(readr)
library(dplyr)
library(randomForest)
library(randomForestSRC)
library(ROCR)

source("R/functions.R")

## Configuration ####

# Select which partitioning schemes to evaluate
# Options: "importance", "pca", "maxshift"
partitioning_schemes <- c("importance")

cat("Partitioning schemes to evaluate:", paste(partitioning_schemes, collapse = ", "), "\n\n")

## Main loop over partitioning schemes ####

for (scheme in partitioning_schemes) {

  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("Evaluating scheme:", toupper(scheme), "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  ## Load data ####

  input_file <- paste0("output/pl2/modeling_frame_regional_partitioned_", scheme, ".csv")
  cat("Loading data from:", input_file, "\n")

  mf <- read_csv(input_file, show_col_types = FALSE)

  # Extract all current scenario data
  all_data <- mf |>
    filter(scenario == "current") |>
    select(-scenario, -x, -y)

  cat("Total observations:", nrow(all_data), "\n")
  cat("Partitions:", paste(sort(unique(all_data$partition)), collapse = ", "), "\n\n")

  ### Select test partition based on prevalence ####

  # Calculate overall prevalence
  overall_prevalence <- sum(all_data$response == 1) / nrow(all_data)

  # Calculate prevalence for each partition
  all_partition_stats <- all_data |>
    group_by(partition) |>
    summarise(
      n_obs = n(),
      n_presence = sum(response == 1),
      n_absence = sum(response == 0),
      prevalence = n_presence / n_obs,
      prev_diff = abs(prevalence - overall_prevalence),
      .groups = "drop"
    )

  cat("Overall prevalence:", round(overall_prevalence, 4), "\n\n")
  cat("Partition prevalences:\n")
  print(all_partition_stats |> select(partition, n_obs, prevalence, prev_diff))

  # Select partition with prevalence closest to overall prevalence as test partition
  test_partition <- all_partition_stats |>
    slice_min(prev_diff, n = 1) |>
    pull(partition)

  cat("\nTest partition selected:", test_partition,
      "(prevalence =", round(all_partition_stats$prevalence[test_partition], 4),
      ", closest to overall)\n")

  # CV folds (remaining partitions used for training)
  cv_partitions <- setdiff(1:5, test_partition)

  cat("Training partitions (pooled):", paste(cv_partitions, collapse = ", "), "\n\n")

  # Split into training (pooled CV partitions) and test
  train_data <- all_data |>
    filter(partition %in% cv_partitions) |>
    select(-partition) |>
    as.data.frame()

  test_data <- all_data |>
    filter(partition == test_partition) |>
    select(-partition) |>
    as.data.frame()

  cat("Training data (pooled CV partitions):", nrow(train_data), "obs\n")
  cat("  Presences:", sum(train_data$response == 1), "\n")
  cat("  Absences:", sum(train_data$response == 0), "\n")
  cat("Test data (partition", test_partition, "):", nrow(test_data), "obs\n")
  cat("  Presences:", sum(test_data$response == 1), "\n")
  cat("  Absences:", sum(test_data$response == 0), "\n\n")

  ## Load optimal hyperparameters ####

  hyperparams_file <- paste0("output/pl2/hyperparameters_", scheme, ".csv")
  hyperparams <- read_csv(hyperparams_file, show_col_types = FALSE)

  model_type <- hyperparams$model_type
  mtry <- hyperparams$mtry
  nodesize <- hyperparams$nodesize

  cat("Optimal hyperparameters:\n")
  cat("  Model type:", model_type, "\n")
  cat("  mtry:", mtry, "\n")
  cat("  nodesize:", nodesize, "\n\n")

  ## Prepare data ####

  # Convert response to factor with explicit levels
  train_data$response <- factor(as.character(train_data$response), levels = c("0", "1"))
  test_data$response <- factor(as.character(test_data$response), levels = c("0", "1"))

  # Ensure ar50 is properly handled as factor if it exists
  if("ar50" %in% names(train_data)) {
    train_data$ar50 <- as.factor(as.character(train_data$ar50))
    test_data$ar50 <- as.factor(as.character(test_data$ar50))
  }

  ## Train final model ####

  cat("Training final", model_type, "model on pooled CV partitions...\n")

  if (model_type == "BRF") {
    # Calculate class sizes for balanced sampling
    n_presence <- sum(train_data$response == "1")
    n_absence <- sum(train_data$response == "0")
    balanced_size <- min(n_presence, n_absence)

    # Train BRF model with balanced sampling and importance
    model <- randomForest::randomForest(
      formula = response ~ .,
      data = train_data,
      ntree = 1000,
      mtry = mtry,
      nodesize = nodesize,
      sampsize = c(balanced_size, balanced_size),
      replace = TRUE,
      importance = TRUE,
      do.trace = FALSE
    )

  } else if (model_type == "RFQ") {
    # Train RFQ model with importance
    model <- randomForestSRC::imbalanced(
      formula = response ~ .,
      data = train_data,
      ntree = 3000,
      mtry = mtry,
      nodesize = nodesize,
      importance = TRUE,
      do.trace = FALSE
    )
  }

  cat("Model training complete.\n\n")

  ## Calculate thresholds on training data ####

  cat("Calculating classification thresholds on training data...\n")

  if (model_type == "BRF") {
    # Get probabilities on training data
    train_pred_prob <- predict(model, newdata = train_data, type = "prob")[, "1"]

    # Calculate thresholds using ROCR
    thresholds <- calculate_classification_thresholds(
      predicted_prob = train_pred_prob,
      true_labels = train_data$response,
      method = "gmean",
      target_sensitivity = 0.95,
      target_specificity = 0.95
    )

  } else if (model_type == "RFQ") {
    # Get probabilities and class predictions on training data
    train_pred <- predict(model, newdata = train_data)
    train_pred_prob <- train_pred$predicted[, "1"]
    train_pred_class <- train_pred$class

    # Calculate thresholds using ROCR
    thresholds <- calculate_classification_thresholds(
      predicted_prob = train_pred_prob,
      true_labels = train_data$response,
      method = "implicit",
      implicit_class = train_pred_class,
      target_sensitivity = 0.95,
      target_specificity = 0.95
    )
  }

  # Extract thresholds from result
  threshold_default <- thresholds$threshold_default
  threshold_sens95 <- thresholds$threshold_sens
  threshold_spec95 <- thresholds$threshold_spec

  cat("Thresholds calculated:\n")
  cat("  Default/Optimal:", round(threshold_default, 4), "\n")
  cat("  High sensitivity (95%):", round(threshold_sens95, 4), "\n")
  cat("  High specificity (95%):", round(threshold_spec95, 4), "\n\n")

  ## Save thresholds with training performance ####

  thresholds_summary <- data.frame(
    threshold_type = c("default", "sens95", "spec95"),
    threshold_value = c(threshold_default, threshold_sens95, threshold_spec95),
    stringsAsFactors = FALSE
  )

  # Add training set performance for each threshold
  train_metrics_list <- list(
    get_metrics_at_threshold(train_pred_prob, train_data$response, threshold_default),
    get_metrics_at_threshold(train_pred_prob, train_data$response, threshold_sens95),
    get_metrics_at_threshold(train_pred_prob, train_data$response, threshold_spec95)
  )

  thresholds_summary <- cbind(thresholds_summary, bind_rows(train_metrics_list))
  names(thresholds_summary)[3:7] <- paste0("train_", names(thresholds_summary)[3:7])

  thresholds_file <- paste0("output/pl2/thresholds_", scheme, ".csv")
  write_csv(thresholds_summary, thresholds_file)
  cat("Saved threshold summary to:", thresholds_file, "\n\n")

  ## Evaluate on test partition ####

  cat("Evaluating on test partition...\n")

  # Get predictions on test data
  if (model_type == "BRF") {
    test_pred_prob <- predict(model, newdata = test_data, type = "prob")[, "1"]
  } else if (model_type == "RFQ") {
    test_pred <- predict(model, newdata = test_data)
    test_pred_prob <- test_pred$predicted[, "1"]
  }

  # Calculate metrics at each threshold
  test_metrics_list <- list(
    get_metrics_at_threshold(test_pred_prob, test_data$response, threshold_default),
    get_metrics_at_threshold(test_pred_prob, test_data$response, threshold_sens95),
    get_metrics_at_threshold(test_pred_prob, test_data$response, threshold_spec95)
  )

  evaluation_results <- data.frame(
    scheme = scheme,
    model_type = model_type,
    threshold_type = c("default", "sens95", "spec95"),
    threshold_value = c(threshold_default, threshold_sens95, threshold_spec95),
    stringsAsFactors = FALSE
  )

  evaluation_results <- cbind(evaluation_results, bind_rows(test_metrics_list))

  cat("\nTest partition evaluation results:\n")
  print(evaluation_results)

  ## Save results ####

  # Save model
  model_file <- paste0("output/pl2/model_", scheme, ".rds")
  saveRDS(model, model_file)
  cat("\nSaved model to:", model_file, "\n")

  # Save evaluation results
  eval_file <- paste0("output/pl2/evaluation_", scheme, ".csv")
  write_csv(evaluation_results, eval_file)
  cat("Saved evaluation results to:", eval_file, "\n")

  cat("\nCompleted evaluation for scheme:", toupper(scheme), "\n")
}


# sessionInfo ####

sessioninfo::session_info()
