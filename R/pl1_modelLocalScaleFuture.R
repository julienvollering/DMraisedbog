# Local Scale Random Forest Current and Future Predictions ####

# This script trains a random forest quantile classification model using
# train + test partitions (80% of data), refits calibration on the calib
# partition, and generates predictions over Norway using both current
# climate predictors (1981-2010) and future climate predictors under
# SSP585 scenario (2071-2100).

library(readr)
library(dplyr)
library(purrr)
library(randomForestSRC)
library(terra)
library(probably)
library(yardstick)

## Load data and hyperparameters ####

cat("Loading partitioned data and CV-optimized hyperparameters...\n")

# Read main partitioned dataset
training_data_partitioned <- read_csv("output/data_partitioned.csv")

# Extract training + test partitions for final model (80% of data)
train_test_data <- training_data_partitioned |>
  filter(outer %in% c("train", "test")) |>
  select(-outer, -inner, -x, -y)  # Remove partition columns and coordinates

# Extract calibration partition
calib_data <- training_data_partitioned |>
  filter(outer == "calib") |>
  select(-outer, -inner, -x, -y)  # Remove partition columns and coordinates

cat("Training + test data loaded:\n")
cat("Observations:", nrow(train_test_data), "\n")
cat("Response variable class:", class(train_test_data$response), "\n")
cat("Response variable values:", paste(unique(train_test_data$response), collapse = ", "), "\n")
cat("Prevalence:", round(mean(train_test_data$response == 1) * 100, 2), "%\n")

cat("Calibration data:\n")
cat("Observations:", nrow(calib_data), "\n")
cat("Prevalence:", round(mean(calib_data$response == 1) * 100, 2), "%\n")

# Load CV-optimized hyperparameters
if(!file.exists("output/final_hyperparameters.csv")) {
  ntree <- 3000  # Default value if not tuned
  final_mtry <- 10  # Default value if not tuned
  final_nodesize <- 5  # Default value if not tuned
  
  cat("final_hyperparameters.csv not found. Using script hyperparameters.\n")
  cat("ntree:", ntree, "\n")
  cat("mtry:", final_mtry, "\n")
  cat("nodesize:", final_nodesize, "\n")

} else {
  cat("final_hyperparameters.csv found. Loading hyperparameters...\n")
  final_params <- read_csv("output/final_hyperparameters.csv", col_types = "cd")
  
  # Extract hyperparameter values
  ntree <- as.numeric(final_params$value[final_params$parameter == "ntree"])
  final_mtry <- as.numeric(final_params$value[final_params$parameter == "mtry"])
  final_nodesize <- as.numeric(final_params$value[final_params$parameter == "nodesize"])
  cv_mean_gmean <- as.numeric(final_params$value[final_params$parameter == "cv_mean_gmean"])
  cv_sd_gmean <- as.numeric(final_params$value[final_params$parameter == "cv_sd_gmean"])
  
  cat("ntree:", ntree, "\n")
  cat("mtry:", final_mtry, "\n")
  cat("nodesize:", final_nodesize, "\n")
  cat("CV G-mean (mean ± sd):", round(cv_mean_gmean, 4), "±", round(cv_sd_gmean, 4), "\n")
}

## Train final model on 80% dataset ####

# Prepare full training data as data.frame (imbalanced() doesn't handle tibbles)
full_train_data <- train_test_data |>
  mutate(response = factor(response, levels = c("1", "0"))) |> # make "1" the positive class
  as.data.frame()  # Convert tibble to data.frame for imbalanced()

# Ensure ar50 is properly handled as factor if it exists
if("ar50" %in% names(full_train_data)) {
  full_train_data$ar50 <- as.factor(as.character(full_train_data$ar50))
}

n_train <- 20e3  # Use more observations since we have more data
cat("Training final model on", n_train, "observations from train+test partitions (80% of data)...\n")

# Randomly subsample training data for runtime efficiency
set.seed(42)
large_train_data <- full_train_data |> 
  slice_sample(n = min(n_train, nrow(full_train_data)))

# Train final model with selected hyperparameters
start_time <- Sys.time()
final_model_80split <- imbalanced(
  formula = response ~ .,
  data = large_train_data,
  ntree = ntree,
  mtry = final_mtry,
  nodesize = final_nodesize,
  importance = TRUE,  # Calculate variable importance
  do.trace = 60,
  seed = -42
)
end_time <- Sys.time()
training_duration <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat("Final model training completed in", round(training_duration, 0), "minutes.\n")
cat("Trained on", nrow(large_train_data), "of", nrow(full_train_data),
    "(", round(nrow(large_train_data) / nrow(full_train_data) * 100, 1), 
    "% ) observations.\n")

final_model_80split

## Refit calibration on calibration partition ####

cat("\nRefitting calibration on calibration partition...\n")

# Prepare calibration data as data.frame for randomForestSRC
calib_df <- calib_data |>
  mutate(response = factor(response, levels = c("1", "0"))) |>
  as.data.frame()

# Ensure ar50 is properly handled as factor if it exists
if("ar50" %in% names(calib_df)) {
  calib_df$ar50 <- as.factor(as.character(calib_df$ar50))
}

# Get calibration predictions (probability for positive class)
calib_preds <- predict.rfsrc(final_model_80split, calib_df, importance = FALSE)

# Create calibration data frame for probably package
calib_results <- tibble(
  truth = factor(calib_data$response, levels = c("1", "0")),
  .pred_1 = calib_preds$predicted[,1],  # Probability for positive class ("1")
  .pred_0 = calib_preds$predicted[,2],  # Probability for negative class ("0")
  class = calib_preds$class
) 

cat("Calibration predictions generated.\n")
cat("Raw prediction range:", round(min(calib_preds$predicted[,1]), 4), 
    "to", round(max(calib_preds$predicted[,1]), 4), "\n")

# Fit beta calibration method on calibration set
tryCatch({
  final_cal_80split <- cal_estimate_beta(
    calib_results,
    truth = truth,
    estimate = dplyr::starts_with(".pred_"),
    event_level = "first"
  )
  cat("Beta calibration fitted successfully.\n")
  print(final_cal_80split)
}, error = function(e) {
  cat("Beta calibration failed, using raw predictions.\n")
  final_cal_80split <<- NULL
})

## Load current and future predictors ####

cat("\nLoading current and future predictor rasters...\n")

# Load global model predictions for both time periods
rf_global_current <- rast("output/rf_global_pred_regional_current.tif")
names(rf_global_current) <- "rf_global"
cat("Global model current predictions loaded.\n")

rf_global_future_files <- list.files("output", pattern = "rf_global.*future", full.names = TRUE)
if(length(rf_global_future_files) > 0) {
  rf_global_future <- rast(rf_global_future_files[1])
  names(rf_global_future) <- "rf_global"
  cat("Global model future predictions loaded:", rf_global_future_files[1], "\n")
} else {
  # If future global predictions don't exist, use current as placeholder
  cat("Future global predictions not found, using current predictions as placeholder...\n")
  rf_global_future <- rf_global_current
}

# Load regional predictors for both time periods
preds_nor_250m_current <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")
cat("Regional current predictors loaded.\n")

preds_nor_250m_future <- rast("output/predictors_regional_250m_Norway_future_EPSG3035.tif")
cat("Regional future predictors loaded.\n")

# Combine predictors for each time period
predictors_current <- c(rf_global_current, preds_nor_250m_current)
predictors_future <- c(rf_global_future, preds_nor_250m_future)
cat("Combined current predictor stack created with", nlyr(predictors_current), "layers.\n")
cat("Combined future predictor stack created with", nlyr(predictors_future), "layers.\n")

## Generate current and future predictions over full Norway extent ####

cat("\nGenerating current and future predictions over Norway...\n")

# Define prediction function that handles calibration
predfun <- function(model, data) {
  v <- predict.rfsrc(model, data, importance = FALSE)
  raw_prob <- v$predicted[,1]  # Probability for positive class
  raw_class <- v$class
  
  # Apply calibration if available
  if(!is.null(final_cal_80split)) {
    # Create data frame for calibration
    cal_data <- tibble(
      truth = factor("1", levels = c("1", "0")),  # Dummy truth for cal_apply
      .pred_1 = raw_prob,
      .pred_0 = v$predicted[,2]
    )
    # Apply calibration
    cal_result <- cal_apply(cal_data, final_cal_80split)
    calibrated_prob <- cal_result$.pred_1
    calibrated_class <- if_else(calibrated_prob >= 0.5, "1", "0")
    
    return(cbind(
      p_raw = raw_prob, 
      p_cal = calibrated_prob,
      class_raw = as.numeric(raw_class) - 1,  # Convert factor to 0/1
      class_cal = as.numeric(calibrated_class == "1")
    ))
  } else {
    return(cbind(
      p_raw = raw_prob, 
      p_cal = raw_prob,  # Same as raw if no calibration
      class_raw = as.numeric(raw_class) - 1,
      class_cal = as.numeric(raw_class) - 1
    ))
  }
}

# Make spatial predictions for current conditions
cat("Predicting over full Norway extent for current conditions (this may take time)...\n")

start_pred_time_current <- Sys.time()

spatial_preds_current <- terra::predict(
  predictors_current, 
  final_model_80split, 
  fun = predfun, 
  na.rm = TRUE,
  filename = "output/rf_local_pred_current_all_layers.tif",
  overwrite = TRUE
)

end_pred_time_current <- Sys.time()
prediction_duration_current <- as.numeric(difftime(end_pred_time_current, start_pred_time_current, units = "mins"))
cat("Current predictions completed in", round(prediction_duration_current, 1), "minutes.\n")

# Make spatial predictions for future conditions
cat("Predicting over full Norway extent for future conditions (this may take time)...\n")

start_pred_time_future <- Sys.time()

spatial_preds_future <- terra::predict(
  predictors_future, 
  final_model_80split, 
  fun = predfun, 
  na.rm = TRUE,
  filename = "output/rf_local_pred_future_all_layers.tif",
  overwrite = TRUE
)

end_pred_time_future <- Sys.time()
prediction_duration_future <- as.numeric(difftime(end_pred_time_future, start_pred_time_future, units = "mins"))
cat("Future predictions completed in", round(prediction_duration_future, 1), "minutes.\n")

# Summary statistics for both time periods
cat("\n=== CURRENT CONDITIONS SUMMARY ===\n")
cat("Raw positive predictions:", global(spatial_preds_current$class_raw, "sum", na.rm=TRUE)[,1], "pixels\n")
if(!is.null(final_cal_80split)) {
  cat("Calibrated positive predictions:", global(spatial_preds_current$class_cal, "sum", na.rm=TRUE)[,1], "pixels\n")
}

cat("\n=== FUTURE CONDITIONS SUMMARY ===\n")
cat("Raw positive predictions:", global(spatial_preds_future$class_raw, "sum", na.rm=TRUE)[,1], "pixels\n")
if(!is.null(final_cal_80split)) {
  cat("Calibrated positive predictions:", global(spatial_preds_future$class_cal, "sum", na.rm=TRUE)[,1], "pixels\n")
}

# Calculate change between time periods
if(!is.null(final_cal_80split)) {
  current_pixels <- global(spatial_preds_current$class_cal, "sum", na.rm=TRUE)[,1]
  future_pixels <- global(spatial_preds_future$class_cal, "sum", na.rm=TRUE)[,1]
  change_pixels <- future_pixels - current_pixels
  change_percent <- (change_pixels / current_pixels) * 100
  
  cat("\n=== PREDICTED CHANGE ===\n")
  cat("Change in suitable pixels:", change_pixels, "(", round(change_percent, 1), "% change)\n")
}

# View variable importance
cat("\nVariable importance (top 15):\n")
final_model_80split$importance[,1] |> 
  tibble::enframe(name = "variable", value = "importance") |> 
  arrange(desc(importance)) |> 
  slice_head(n = 15) |>
  print(n = 15)

## Save results ####

cat("\nSaving results...\n")

# Save the trained model
model_path <- "output/final_model_local_80split.rds"
saveRDS(final_model_80split, model_path)
cat("Final model saved to:", model_path, "\n")

# Save calibration object if applicable
if(!is.null(final_cal_80split)) {
  cal_path <- "output/calibration80split.rds"
  saveRDS(final_cal_80split, cal_path)
  cat("Calibration object saved to:", cal_path, "\n")
} else {
  cat("No calibration object to save (using raw predictions).\n")
}

cat("\nOutput files created:\n")
cat("- Current predictions: output/rf_local_pred_current_all_layers.tif\n")
cat("- Future predictions: output/rf_local_pred_future_all_layers.tif\n")
cat("- Total prediction time:", round(prediction_duration_current + prediction_duration_future, 1), "minutes\n")

# sessionInfo ####

sessioninfo::session_info()
