# Random Forest Quantile Classification for Local Scale Modeling ####

# This script implements hyperparameter tuning for random forest quantile 
# classification using randomForestSRC::tune.rfsrc() with imbalanced classification.
# Uses inner CV folds from partitionData.R for hyperparameter optimization.

library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(randomForestSRC)

## Load partitioned data ####

# Read main partitioned dataset
training_data_partitioned <- read_csv("output/data_partitioned.csv")

cat("Loaded partitioned data:\n")
cat("Total observations:", nrow(training_data_partitioned), "\n")
cat("Outer partition summary:\n")
print(table(training_data_partitioned$outer))

# Extract training partition for hyperparameter tuning
train_data <- training_data_partitioned |>
  filter(outer == "train") |>
  select(-outer, -x, -y)  # Remove partition column and coordinates

cat("\nTraining data for hyperparameter tuning:\n")
cat("Observations:", nrow(train_data), "\n")
cat("Inner CV folds:", table(train_data$inner), "\n")
cat("Response variable class:", class(train_data$response), "\n")
cat("Response variable values:", paste(unique(train_data$response), collapse = ", "), "\n")
cat("Prevalence:", round(mean(train_data$response == 1) * 100, 2), "%\n")

## Custom hyperparameter tuning with imbalanced() ####

# Function to train a single model with given parameters
train_single_model <- function(cv_fold, mtry, nodesize, job_id, train_data) {
  
  cat("Job", job_id, "- Fold", cv_fold, "mtry:", mtry, "nodesize:", nodesize, "...\n")
  
  # Split data by inner CV fold and convert to data.frames
  train_fold <- train_data |> 
    filter(inner != cv_fold) |> 
    select(-inner) |>
    as.data.frame()  # Convert tibble to data.frame for imbalanced()
  
  val_fold <- train_data |> 
    filter(inner == cv_fold) |> 
    select(-inner) |>
    as.data.frame()  # Convert tibble to data.frame for imbalanced()
  
  # Convert response to factor explicitly
  train_fold$response <- as.factor(as.character(train_fold$response))
  val_fold$response <- as.factor(as.character(val_fold$response))
  
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
  
  # Create confusion matrix
  cm <- table(Actual = true_class, Predicted = pred_class)
  if(nrow(cm) == 2 && ncol(cm) == 2) {
    tn <- cm[1,1]  # True negatives
    fp <- cm[1,2]  # False positives  
    fn <- cm[2,1]  # False negatives
    tp <- cm[2,2]  # True positives
    
    # Calculate TPR and TNR
    tpr <- tp / (tp + fn)  # Sensitivity
    tnr <- tn / (tn + fp)  # Specificity
    
    # Calculate G-mean
    gmean <- sqrt(tpr * tnr)
  } else {
    tpr <- 0  # Handle edge cases
    tnr <- 0
    gmean <- 0
  }
  
  # Create result row
  result <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
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

## Create parameter grid ####

# Define hyperparameter grid for tuning
n_predictors <- ncol(train_data) - 2  # Exclude response and inner columns
mtry_values <- seq(floor(sqrt(n_predictors)), ceiling(n_predictors/3), by = 1)
nodesize_values <- c(1, 5, 20)

cat("Number of predictors:", n_predictors, "\n")
cat("mtry values to test:", paste(mtry_values, collapse = ", "), "\n")
cat("nodesize values to test:", paste(nodesize_values, collapse = ", "), "\n")

cat("\nCreating parameter grid for hyperparameter tuning...\n")

# Create complete parameter combination grid
param_grid <- expand_grid(
  cv_fold = 1:3,
  mtry = mtry_values, 
  nodesize = nodesize_values
) |>
  mutate(
    job_id = row_number()
  )

cat("Total parameter combinations:", nrow(param_grid), "\n")
cat("Jobs per CV fold:", nrow(param_grid) / 3, "\n")

### Runtime estimation ####

cat("\nEstimating runtime with a single RFQ model (500 trees)...\n")

# Prepare small test dataset (first fold) for timing
test_data <- train_data |>
  filter(inner != 1) |>
  select(-inner) |>
  mutate(response = as.factor(as.character(response))) |>
  as.data.frame()

# Ensure ar50 is properly handled as factor if it exists
if("ar50" %in% names(test_data)) {
  test_data$ar50 <- as.factor(as.character(test_data$ar50))
}

cat("Test data size:", nrow(test_data), "observations\n")

# Time a single model run
start_time <- Sys.time()

test_model <- imbalanced(
  formula = response ~ .,
  data = test_data,
  ntree = 500,  # Reduced trees for timing
  mtry = mtry_values[1],  # Use first mtry value
  nodesize = nodesize_values[1],  # Use first nodesize value
  importance = FALSE,
  do.trace = FALSE
)

end_time <- Sys.time()
runtime_500_trees <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat("Runtime for 500 trees:", round(runtime_500_trees, 2), "minutes\n")

# Estimate total tuning time
n_calls_to_imbalanced <- nrow(param_grid)
scaling_factor <- 3000 / 500
estimated_total_time <- n_calls_to_imbalanced * runtime_500_trees * scaling_factor

cat("Estimated total tuning time:", round(estimated_total_time / 60, 1), "hours\n")

## Check paramater grid progress ####

# Initialize progress tracking CSV file

progress_file <- "output/hyperparameter_tuning_progress.csv"

if(file.exists(progress_file)) {
  cat("Found existing progress file.\nChecking for completed jobs...\n")
  
  # Read existing results
  completed_runs <- read_csv(progress_file, show_col_types = FALSE)
  
  if(nrow(completed_runs) > 0) {
    cat("Found", nrow(completed_runs), "completed jobs\n")
    
    # Filter out completed parameter combinations
    param_grid <- param_grid |>
      anti_join(completed_runs, by = c("cv_fold", "mtry", "nodesize"))
    
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
        list(cv_fold, mtry, nodesize, job_id),
        ~train_single_model(
          cv_fold = ..1, 
          mtry = ..2, 
          nodesize = ..3, 
          job_id = ..4,
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

## Aggregate across CV folds ####

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Cross-validation tuning results:\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Find best parameters by CV fold
best_params_by_fold <- all_results |>
  group_by(cv_fold) |>
  slice_max(gmean, n = 1, with_ties = FALSE) |>
  select(cv_fold, mtry, nodesize, gmean) |>
  rename(best_mtry = mtry, best_nodesize = nodesize, best_gmean = gmean) |>
  ungroup()

print(best_params_by_fold)

# Calculate average performance across folds
avg_results <- best_params_by_fold |>
  summarise(
    mean_mtry = mean(best_mtry),
    mean_nodesize = mean(best_nodesize),
    mean_gmean = mean(best_gmean),
    sd_gmean = sd(best_gmean),
    .groups = 'drop'
  )

cat("\nAverage results across CV folds:\n")
cat("Mean optimal mtry:", round(avg_results$mean_mtry, 1), "\n")
cat("Mean optimal nodesize:", round(avg_results$mean_nodesize, 1), "\n")
cat("Mean G-mean:", round(avg_results$mean_gmean, 4), "\n")
cat("G-mean std dev:", round(avg_results$sd_gmean, 4), "\n")

# Average the best parameters in each fold
final_mtry <- round(avg_results$mean_mtry)
final_nodesize <- round(avg_results$mean_nodesize)

cat("\nFinal selected hyperparameters:\n")
cat("ntree: 3000 (fixed)\n")
cat("mtry:", final_mtry, "\n")
cat("nodesize:", final_nodesize, "\n")

## Save hyperparameter tuning results ####

# Save all tuning results and CV results by fold
write_csv(all_results, "output/hyperparameter_tuning_all_results.csv")
write_csv(best_params_by_fold, "output/hyperparameter_tuning_by_fold.csv")

# Save final hyperparameters
final_params <- data.frame(
  parameter = c("ntree", "mtry", "nodesize", "cv_mean_gmean", "cv_sd_gmean"),
  value = c(3000, final_mtry, final_nodesize, 
            round(avg_results$mean_gmean, 4), round(avg_results$sd_gmean, 4))
)
write_csv(final_params, "output/final_hyperparameters.csv")

cat("\nSaved hyperparameter tuning results:\n")
cat("- Progress tracking: output/hyperparameter_tuning_progress.csv\n")
cat("- All tuning results: output/hyperparameter_tuning_all_results.csv\n")
cat("- CV tuning by fold: output/hyperparameter_tuning_by_fold.csv\n")
cat("- Final hyperparameters: output/final_hyperparameters.csv\n")

cat("\nHyperparameter tuning completed successfully!\n")
cat("Next step: Run 'R/modelLocalScale.R' to train the final model with these optimized parameters.\n")

# sessionInfo ####

sessioninfo::session_info()
