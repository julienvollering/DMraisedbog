# BRF vs QRF Model Comparison and Hyperparameter Tuning ####

# This script implements comprehensive hyperparameter tuning that compares
# Balanced Random Forest (BRF) and Quantile Random Forest (QRF) approaches
# as if model type were a hyperparameter. BRF uses balanced sampling with
# optimal threshold selection, while QRF uses the imbalanced() approach.
# The best performing model type and hyperparameters are selected based on
# cross-validated G-mean performance.

library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(randomForestSRC)
library(randomForest)
library(ROCR)

## Load partitioned data ####

# Read main partitioned dataset
mf_partitioned <- read_csv("output/pl2/modeling_frame_regional_partitioned.csv")

mf_current <- mf_partitioned |> 
  filter(scenario == "current") |>
  select(-scenario)  # Keep partition columns and coordinates for spatial CV

# Extract training partition for hyperparameter tuning
train_data <- mf_current |>
  filter(outer == "train") |>
  select(-outer, -x, -y)  # Remove partition column and coordinates

cat("\nTraining data for hyperparameter tuning:\n")
cat("Observations:", nrow(train_data), "\n")
cat("Inner CV folds:", table(train_data$inner), "\n")
cat("Response variable class:", class(train_data$response), "\n")
cat("Response variable values:", paste(unique(train_data$response), collapse = ", "), "\n")
cat("Prevalence:", round(mean(train_data$response == 1) * 100, 2), "%\n")

## Helper functions ####

# Function to calculate G-mean from predictions and true values
# Robust to factor level ordering by explicitly identifying positive/negative classes
calculate_gmean <- function(predicted_class, true_class) {
  # Convert to character to avoid factor level issues
  predicted_class <- as.character(predicted_class)
  true_class <- as.character(true_class)

  # Create confusion matrix
  cm <- table(Actual = true_class, Predicted = predicted_class)

  if(nrow(cm) == 2 && ncol(cm) == 2) {
    # Explicitly identify which row/column corresponds to positive class ("1")
    pos_row <- which(rownames(cm) == "1")
    pos_col <- which(colnames(cm) == "1")
    neg_row <- which(rownames(cm) == "0")
    neg_col <- which(colnames(cm) == "0")

    # Handle case where positive or negative class is missing
    if(length(pos_row) == 0 || length(pos_col) == 0 ||
       length(neg_row) == 0 || length(neg_col) == 0) {
      tpr <- 0
      tnr <- 0
      gmean <- 0
    } else {
      # Extract confusion matrix values using explicit indexing
      tp <- cm[pos_row, pos_col]  # True positives
      fn <- cm[pos_row, neg_col]  # False negatives
      fp <- cm[neg_row, pos_col]  # False positives
      tn <- cm[neg_row, neg_col]  # True negatives

      # Calculate TPR and TNR
      tpr <- tp / (tp + fn)  # Sensitivity
      tnr <- tn / (tn + fp)  # Specificity

      # Calculate G-mean
      gmean <- sqrt(tpr * tnr)
    }
  } else {
    # Handle edge cases (only one class present)
    tpr <- 0
    tnr <- 0
    gmean <- 0
  }

  return(list(gmean = gmean, tpr = tpr, tnr = tnr))
}

## Model training functions ####

# Function to train a BRF model with threshold optimization
train_brf_model <- function(cv_fold, mtry, nodesize, job_id, train_data) {

  cat("Job", job_id, "- BRF Fold", cv_fold, "mtry:", mtry, "nodesize:", nodesize, "...\n")

  # Split data by inner CV fold
  train_fold <- train_data |>
    filter(inner != cv_fold) |>
    select(-inner) |>
    as.data.frame()

  val_fold <- train_data |>
    filter(inner == cv_fold) |>
    select(-inner) |>
    as.data.frame()

  # Convert response to factor with explicit levels (0, 1) for consistency
  train_fold$response <- factor(as.character(train_fold$response), levels = c("0", "1"))
  val_fold$response <- factor(as.character(val_fold$response), levels = c("0", "1"))

  # Ensure ar50 is properly handled as factor if it exists
  if("ar50" %in% names(train_fold)) {
    train_fold$ar50 <- as.factor(as.character(train_fold$ar50))
    val_fold$ar50 <- as.factor(as.character(val_fold$ar50))
  }

  # Calculate class sizes for balanced sampling
  n_presence <- sum(train_fold$response == "1")
  n_absence <- sum(train_fold$response == "0")
  balanced_size <- min(n_presence, n_absence)

  # Train BRF model with balanced sampling
  model <- randomForest(
    formula = response ~ .,
    data = train_fold,
    ntree = 1000,
    mtry = mtry,
    nodesize = nodesize,
    sampsize = c(balanced_size, balanced_size),  # Balanced sampling
    replace = TRUE,
    importance = FALSE,
    do.trace = FALSE
  )

  # Get predictions on training fold for threshold optimization
  train_pred_prob <- predict(model, newdata = train_fold, type = "prob")[, "1"]
  train_true <- as.numeric(as.character(train_fold$response))

  # Optimize threshold using ROCR on training folds
  pred_obj <- prediction(train_pred_prob, train_true)

  # Calculate G-mean for all possible thresholds
  tpr_values <- performance(pred_obj, "tpr")@y.values[[1]]
  tnr_values <- performance(pred_obj, "tnr")@y.values[[1]]
  gmean_values <- sqrt(tpr_values * tnr_values)

  # Find optimal threshold
  optimal_idx <- which.max(gmean_values)
  optimal_threshold <- performance(pred_obj, "tpr")@x.values[[1]][optimal_idx]

  # Get predictions on validation fold
  val_pred_prob <- predict(model, newdata = val_fold, type = "prob")[, "1"]
  val_pred_class <- ifelse(val_pred_prob >= optimal_threshold, "1", "0")
  val_true <- as.character(val_fold$response)

  # Calculate G-mean on validation set with optimal threshold
  metrics <- calculate_gmean(val_pred_class, val_true)
  gmean <- metrics$gmean
  tpr <- metrics$tpr
  tnr <- metrics$tnr

  # Create result row
  result <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    model_type = "BRF",
    cv_fold = cv_fold,
    mtry = mtry,
    nodesize = nodesize,
    gmean = gmean,
    tpr = tpr,
    tnr = tnr,
    stringsAsFactors = FALSE
  )

  # Append to CSV immediately
  write_csv(result, progress_file, append = TRUE)
  cat("  Job", job_id, "completed - G-mean:", round(gmean, 4), "\n")

  return(result)
}

# Function to train a QRF model (existing approach)
train_qrf_model <- function(cv_fold, mtry, nodesize, job_id, train_data) {

  cat("Job", job_id, "- QRF Fold", cv_fold, "mtry:", mtry, "nodesize:", nodesize, "...\n")

  # Split data by inner CV fold and convert to data.frames
  train_fold <- train_data |>
    filter(inner != cv_fold) |>
    select(-inner) |>
    as.data.frame()  # Convert tibble to data.frame for imbalanced()

  val_fold <- train_data |>
    filter(inner == cv_fold) |>
    select(-inner) |>
    as.data.frame()  # Convert tibble to data.frame for imbalanced()

  # Convert response to factor with explicit levels (0, 1) for consistency
  train_fold$response <- factor(as.character(train_fold$response), levels = c("0", "1"))
  val_fold$response <- factor(as.character(val_fold$response), levels = c("0", "1"))

  # Ensure ar50 is properly handled as factor if it exists
  if("ar50" %in% names(train_fold)) {
    train_fold$ar50 <- as.factor(as.character(train_fold$ar50))
    val_fold$ar50 <- as.factor(as.character(val_fold$ar50))
  }

  # Train model with current parameters
  model <- imbalanced(
    formula = response ~ .,
    data = train_fold,
    ntree = 3000,
    mtry = mtry,
    nodesize = nodesize,
    importance = FALSE,  # Skip importance for speed during tuning
    do.trace = FALSE
  )

  # Predict on validation set
  pred <- predict(model, newdata = val_fold)

  # Calculate G-mean on validation set
  pred_class <- pred$class
  true_class <- val_fold$response

  # Calculate G-mean using common function
  metrics <- calculate_gmean(pred_class, true_class)
  gmean <- metrics$gmean
  tpr <- metrics$tpr
  tnr <- metrics$tnr

  # Create result row
  result <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    model_type = "QRF",
    cv_fold = cv_fold,
    mtry = mtry,
    nodesize = nodesize,
    gmean = gmean,
    tpr = tpr,
    tnr = tnr,
    stringsAsFactors = FALSE
  )

  # Append to CSV immediately
  write_csv(result, progress_file, append = TRUE)
  cat("  Job", job_id, "completed - G-mean:", round(gmean, 4), "\n")

  return(result)
}

# Unified function to train either model type
train_single_model <- function(cv_fold, model_type, mtry, nodesize, job_id, train_data) {
  if(model_type == "BRF") {
    return(train_brf_model(cv_fold, mtry, nodesize, job_id, train_data))
  } else if(model_type == "QRF") {
    return(train_qrf_model(cv_fold, mtry, nodesize, job_id, train_data))
  } else {
    stop("Unknown model_type: ", model_type)
  }
}

## Create parameter grid ####

# Define hyperparameter grid for tuning
n_predictors <- ncol(train_data) - 2  # Exclude response and inner columns
mtry_values <- seq(floor(sqrt(n_predictors)), ceiling(n_predictors/2), by = 1)
nodesize_values <- c(1, 5, 20)
model_types <- c("BRF", "QRF")

cat("Number of predictors:", n_predictors, "\n")
cat("Model types to test:", paste(model_types, collapse = ", "), "\n")
cat("mtry values to test:", paste(mtry_values, collapse = ", "), "\n")
cat("nodesize values to test:", paste(nodesize_values, collapse = ", "), "\n")

cat("\nCreating parameter grid for hyperparameter tuning...\n")

# Create complete parameter combination grid
param_grid <- expand_grid(
  cv_fold = 1:3,
  model_type = model_types,
  mtry = mtry_values,
  nodesize = nodesize_values
) |>
  mutate(
    job_id = row_number()
  )

cat("Total parameter combinations:", nrow(param_grid), "\n")
cat("Jobs per CV fold:", nrow(param_grid) / 3, "\n")
cat("Jobs per model type:", nrow(param_grid) / (3 * length(model_types)), "\n")

### Runtime estimation ####

cat("\nEstimating runtime with sample models...\n")

# Prepare small test dataset (first fold) for timing
test_data <- train_data |>
  filter(inner != 1) |>
  select(-inner) |>
  mutate(response = as.factor(as.character(response))) |>
  as.data.frame()

cat("Test data size:", nrow(test_data), "observations\n")

# Time QRF model
cat("Timing QRF model (3000 trees)...\n")
start_time <- Sys.time()
test_model_qrf <- imbalanced(
  formula = response ~ .,
  data = test_data,
  ntree = 3000,
  mtry = mtry_values[1],
  nodesize = nodesize_values[1],
  importance = FALSE,
  do.trace = FALSE
)
end_time <- Sys.time()
runtime_qrf <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Time BRF model
cat("Timing BRF model (1000 trees)...\n")
n_presence <- sum(test_data$response == "1")
n_absence <- sum(test_data$response == "0")
balanced_size <- min(n_presence, n_absence)

start_time <- Sys.time()
test_model_brf <- randomForest(
  formula = response ~ .,
  data = test_data,
  ntree = 1000,
  mtry = mtry_values[1],
  nodesize = nodesize_values[1],
  sampsize = c(balanced_size, balanced_size),
  replace = TRUE,
  importance = FALSE,
  do.trace = FALSE
)
end_time <- Sys.time()
runtime_brf <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("Runtime for QRF (3000 trees):", round(runtime_qrf, 2), "seconds\n")
cat("Runtime for BRF (1000 trees):", round(runtime_brf, 2), "seconds\n")

# Estimate total tuning time
n_qrf_jobs <- sum(param_grid$model_type == "QRF")
n_brf_jobs <- sum(param_grid$model_type == "BRF")
estimated_total_time <- (n_qrf_jobs * runtime_qrf) + (n_brf_jobs * runtime_brf * 1.2)  # 20% overhead for threshold optimization

cat("Estimated total tuning time:", round(estimated_total_time / 60, 1), "minutes\n")

## Check paramater grid progress ####

# Initialize progress tracking CSV file

progress_file <- "output/pl2/hyperparameter_tuning_progress.csv"

if(file.exists(progress_file)) {
  cat("Found existing progress file.\nChecking for completed jobs...\n")
  
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
  cat("No existing progress file. Starting fresh...\n")
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

## Tune remaining paramter grid ####
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
          train_data = train_data
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

## Aggregate across CV folds and select best model ####

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("Cross-validation tuning results: BRF vs QRF comparison\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Calculate average performance across folds for each model type and hyperparameter combination
avg_results <- all_results %>%
  group_by(model_type, mtry, nodesize) %>%
  summarize(
    mean_gmean = mean(gmean),
    sd_gmean = sd(gmean),
    mean_tpr = mean(tpr),
    mean_tnr = mean(tnr),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_gmean))

cat("\nTop 10 model configurations by mean G-mean:\n")
print(head(avg_results, 10))

# Find the best performing configuration overall
best_config <- avg_results %>%
  slice_max(mean_gmean, n = 1)

if(nrow(best_config) > 1) {
  # If there's a tie, select the one with smallest standard deviation
  best_config <- best_config %>%
    slice_min(sd_gmean, n = 1)

  if(nrow(best_config) > 1) {
    # If still tied, prefer BRF for interpretability
    best_config <- best_config %>%
      arrange(model_type) %>%
      slice(1)
  }
}

final_model_type <- best_config$model_type
final_mtry <- best_config$mtry
final_nodesize <- best_config$nodesize
final_mean_gmean <- best_config$mean_gmean
final_sd_gmean <- best_config$sd_gmean

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("FINAL SELECTED MODEL CONFIGURATION:\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Model type:", final_model_type, "\n")
if(final_model_type == "BRF") {
  cat("ntree: 1000 (fixed for BRF)\n")
} else {
  cat("ntree: 3000 (fixed for QRF)\n")
}
cat("mtry:", final_mtry, "\n")
cat("nodesize:", final_nodesize, "\n")
cat("CV G-mean (mean ± sd):", round(final_mean_gmean, 4), "±", round(final_sd_gmean, 4), "\n")

# Show comparison between best BRF and QRF
cat("\n", paste(rep("-", 50), collapse = ""), "\n")
cat("Model type comparison (best of each):\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

best_brf <- avg_results %>% filter(model_type == "BRF") %>% slice_max(mean_gmean, n = 1)
best_qrf <- avg_results %>% filter(model_type == "QRF") %>% slice_max(mean_gmean, n = 1)

if(nrow(best_brf) > 0) {
  cat("Best BRF: G-mean =", round(best_brf$mean_gmean, 4),
      "mtry =", best_brf$mtry, "nodesize =", best_brf$nodesize, "\n")
}

if(nrow(best_qrf) > 0) {
  cat("Best QRF: G-mean =", round(best_qrf$mean_gmean, 4),
      "mtry =", best_qrf$mtry, "nodesize =", best_qrf$nodesize, "\n")
}

## Save hyperparameter tuning results ####

# Save all tuning results and CV results by fold
write_csv(all_results, "output/pl2/hyperparameter_tuning_all_results.csv", 
          append = FALSE)

# Save final hyperparameters and model selection results

# Create comprehensive final hyperparameters file
final_params <- data.frame(
  parameter = c("model_type", "ntree", "mtry", "nodesize", "cv_mean_gmean", "cv_sd_gmean"),
  value = c(
    final_model_type,
    ifelse(final_model_type == "BRF", "1000", "3000"),
    as.character(final_mtry),
    as.character(final_nodesize),
    as.character(round(final_mean_gmean, 6)),
    as.character(round(final_sd_gmean, 6))
  ),
  stringsAsFactors = FALSE
)


write_csv(final_params, "output/pl2/final_hyperparameters.csv", append = FALSE)

# Save aggregated results for analysis
write_csv(avg_results, "output/pl2/hyperparameter_comparison_summary.csv", append = FALSE)

cat("\nSaved results to:\n")
cat("- output/pl2/hyperparameter_tuning_all_results.csv (detailed CV results)\n")
cat("- output/pl2/final_hyperparameters.csv (winning configuration)\n")
cat("- output/pl2/hyperparameter_comparison_summary.csv (model comparison summary)\n")

# sessionInfo ####

sessioninfo::session_info()
