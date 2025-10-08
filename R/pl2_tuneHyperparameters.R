# BRF vs RFQ Model Comparison and Hyperparameter Tuning ####

# This script implements comprehensive hyperparameter tuning that compares
# Balanced Random Forest (BRF) and Quantile Random Forest (RFQ) approaches
# across three data partitioning schemes (importance, pca, maxshift).
#
# CV Strategy:
# - Partition 1 held out entirely as final test set (never used in tuning)
# - Partitions 2-5 used for 4-fold LOO CV during hyperparameter tuning
# - Each CV iteration: 1 of 4 partitions held out for validation, 3 for training
#
# Optimal hyperparameters selected using weighted average G-mean by partition size
# to account for imbalanced CV fold sizes.

library(readr)
library(dplyr)
library(purrr)
library(tidyr)

source("R/functions.R")

## Configuration ####

# Select which partitioning schemes to tune
# Options: "importance", "pca", "maxshift"
partitioning_schemes <- c("importance", "pca", "maxshift")

cat("Partitioning schemes to process:", paste(partitioning_schemes, collapse = ", "), "\n\n")

## Main loop over partitioning schemes ####

for (scheme in partitioning_schemes) {

  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("Processing partitioning scheme:", toupper(scheme), "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  ### Load partitioned data ####

  input_file <- paste0("output/pl2/modeling_frame_regional_partitioned_", scheme, ".csv")
  cat("Loading data from:", input_file, "\n")

  mf <- read_csv(input_file, show_col_types = FALSE)

  # Extract training data (current scenario only)
  train_data <- mf |>
    filter(scenario == "current") |>
    select(-scenario, -x, -y)  # Keep partition, response, and predictors

  cat("Total observations:", nrow(train_data), "\n")
  cat("Partitions:", paste(sort(unique(train_data$partition)), collapse = ", "), "\n\n")

  ### Select test partition based on prevalence ####

  # Calculate overall prevalence
  overall_prevalence <- sum(train_data$response == 1) / nrow(train_data)

  # Calculate prevalence for each partition
  all_partition_stats <- train_data |>
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

  # CV folds (remaining partitions used for tuning)
  cv_partitions <- setdiff(1:5, test_partition)

  cat("CV partitions (for tuning):", paste(cv_partitions, collapse = ", "), "\n\n")

  # Extract CV and test partition statistics from all_partition_stats
  partition_stats <- all_partition_stats |>
    filter(partition %in% cv_partitions) |>
    select(partition, n_obs, n_presence, n_absence, prevalence)

  test_stats <- all_partition_stats |>
    filter(partition == test_partition) |>
    select(n_obs, n_presence, n_absence, prevalence)

  cat("CV partition statistics:\n")
  print(partition_stats)

  cat("\nTest partition (", test_partition, ") statistics:\n", sep = "")
  print(test_stats)

  ### Create parameter grid ####

  # Define hyperparameter grid for tuning
  n_predictors <- ncol(train_data) - 1  # Exclude partition column (response is included as predictor count)
  n_predictors <- n_predictors - 1  # Now exclude response
  mtry_values <- c(4, 7, 11) # seq(floor(sqrt(n_predictors)), ceiling(n_predictors/2), by = 1)
  nodesize_values <- c(1, 20)
  model_types <- c("BRF", "RFQ")

  # Create complete parameter combination grid
  param_grid <- expand_grid(
    cv_fold = cv_partitions,  # Only CV partitions (exclude test partition)
    model_type = model_types,
    mtry = mtry_values,
    nodesize = nodesize_values
  ) |>
    mutate(
      job_id = row_number()
    )

  cat("Total parameter combinations:", nrow(param_grid), "\n")
  cat("Jobs per CV fold:", nrow(param_grid) / length(cv_partitions), "\n")

  ### Check parameter grid progress ####

  # Initialize progress tracking CSV file
  progress_file <- paste0("output/pl2/hyperparameter_tuning_progress_", scheme, ".csv")

  if(file.exists(progress_file)) {
    cat("\nFound existing progress file.\nChecking for completed jobs...\n")

    # Read existing results with explicit column types to ensure timestamp is character
    completed_runs <- read_csv(progress_file, show_col_types = FALSE,
                             col_types = cols(timestamp = col_character()))

    if(nrow(completed_runs) > 0) {
      cat("Found", nrow(completed_runs), "completed jobs\n")

      # Filter out completed parameter combinations
      param_grid <- param_grid |>
        anti_join(completed_runs, by = c("cv_fold", "model_type", "mtry", "nodesize"))

      cat("Remaining jobs to complete:", nrow(param_grid), "\n")

      if(nrow(param_grid) == 0) {
        cat("All parameter combinations already completed!\n")
        cat("Loading results from progress file...\n")
        all_results <- completed_runs
      }
    } else {
      cat("Progress file exists but is empty. Starting fresh...\n")
    }
  } else {
    cat("\nNo existing progress file. Starting fresh...\n")
    # Create empty file with headers
    empty_result <- data.frame(
      timestamp = character(0),
      model_type = character(0),
      cv_fold = numeric(0),
      mtry = numeric(0),
      nodesize = numeric(0),
      gmean = numeric(0),
      tpr = numeric(0),
      tnr = numeric(0),
      stringsAsFactors = FALSE
    )
    write_csv(empty_result, progress_file)
  }

  ### Tune remaining parameter grid ####
  if(nrow(param_grid) > 0) {

    cat("\nStarting hyperparameter tuning...\n")
    cat("Running", nrow(param_grid), "parameter combinations\n")

    # Use purrr::pmap to run each parameter combination
    tuning_results <- param_grid |>
      mutate(
        results = pmap(
          list(cv_fold, model_type, mtry, nodesize, job_id),
          ~train_single_model(
            cv_fold = ..1,
            model_type = ..2,
            mtry = ..3,
            nodesize = ..4,
            job_id = ..5,
            train_data = train_data,
            test_partition = test_partition
          )
        )
      ) |>
      select(results) |>  # Keep only results column to avoid duplicates
      pull(results) |>    # Extract list of data.frames
      bind_rows()         # Combine into single data.frame

    cat("\nHyperparameter tuning completed!\n")

    # Combine with any existing results
    if(exists("completed_runs") && nrow(completed_runs) > 0) {
      all_results <- bind_rows(completed_runs, tuning_results)
    } else {
      all_results <- tuning_results
    }
  } else {
    # All jobs already completed - use existing results
    all_results <- completed_runs
  }

  ### Aggregate across CV folds with WEIGHTED averaging ####

  # Join partition sizes to results
  all_results_weighted <- all_results |>
    left_join(partition_stats |> select(partition, n_obs),
              by = c("cv_fold" = "partition"))

  # Calculate WEIGHTED average performance across folds
  avg_results <- all_results_weighted |>
    group_by(model_type, mtry, nodesize) |>
    summarize(
      weighted_mean_gmean = weighted.mean(gmean, w = n_obs),
      weighted_mean_tpr = weighted.mean(tpr, w = n_obs),
      weighted_mean_tnr = weighted.mean(tnr, w = n_obs),
      # Also keep unweighted for comparison
      mean_gmean = mean(gmean),
      sd_gmean = sd(gmean),
      min_gmean = min(gmean),
      max_gmean = max(gmean),
      .groups = "drop"
    ) |>
    arrange(desc(weighted_mean_gmean))

  cat("Top 10 parameter combinations by WEIGHTED mean G-mean:\n")
  print(head(avg_results, 10))

  # Find the best performing configuration overall (by weighted mean)
  best_config <- avg_results |>
    slice_max(weighted_mean_gmean, n = 1)

  if(nrow(best_config) > 1) {
    # If there's a tie, select the one with smallest standard deviation
    best_config <- best_config |>
      slice_min(sd_gmean, n = 1)
  }

  final_model_type <- best_config$model_type
  final_mtry <- best_config$mtry
  final_nodesize <- best_config$nodesize
  final_weighted_mean_gmean <- best_config$weighted_mean_gmean
  final_mean_gmean <- best_config$mean_gmean
  final_sd_gmean <- best_config$sd_gmean

  ### Save hyperparameter tuning results ####

  # Save aggregated results
  avg_results_file <- paste0("output/pl2/hyperparameter_tuning_summary_", scheme, ".csv")
  write_csv(avg_results, avg_results_file)

  # Save final hyperparameters
  final_hyperparams <- data.frame(
    scheme = scheme,
    test_partition = test_partition,
    cv_partitions = paste(cv_partitions, collapse = ","),
    model_type = final_model_type,
    mtry = final_mtry,
    nodesize = final_nodesize,
    weighted_mean_gmean = final_weighted_mean_gmean,
    mean_gmean = final_mean_gmean,
    sd_gmean = final_sd_gmean,
    min_gmean = best_config$min_gmean,
    max_gmean = best_config$max_gmean,
    weighted_mean_tpr = best_config$weighted_mean_tpr,
    weighted_mean_tnr = best_config$weighted_mean_tnr,
    stringsAsFactors = FALSE
  )

  final_hyperparams_file <- paste0("output/pl2/hyperparameters_", scheme, ".csv")
  write_csv(final_hyperparams, final_hyperparams_file)
  cat("Saved final hyperparameters to:", final_hyperparams_file, "\n")

  cat("\nCompleted scheme:", toupper(scheme), "\n")
}

## Summary across schemes ####

cat("Summary across all schemes:\n\n")
for (scheme in partitioning_schemes) {
  final_file <- paste0("output/pl2/hyperparameters_", scheme, ".csv")
  if (file.exists(final_file)) {
    results <- read_csv(final_file, show_col_types = FALSE)
    cat(sprintf("%-15s: %s, mtry=%d, nodesize=%d, weighted G-mean=%.4f\n",
                toupper(scheme),
                results$model_type,
                results$mtry,
                results$nodesize,
                results$weighted_mean_gmean))
  }
}

# sessionInfo ####

sessioninfo::session_info()
