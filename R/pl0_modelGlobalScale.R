library(tidyverse)
library(terra)
library(randomForest)
library(randomForestSRC)

# Load and prepare presence data ####

# Load combined presence data
presence_no <- read_csv("output/presence_coords_global_no_thinned.csv")
presence_eu <- read_csv("output/presence_coords_global_eu_thinned.csv")
presence_combined <- bind_rows(presence_no, presence_eu)

n_presence <- nrow(presence_combined)

# Load global predictors
predictors_global <- rast("output/predictors_global_5km_EUNorway_EPSG3035.tif")

# Convert raster to data frame for efficient processing
coords_global <- as.data.frame(predictors_global, xy = TRUE, na.rm = TRUE)

# Create training dataset ####

# Extract predictor values for presence points
presence_values <- terra::extract(
  predictors_global,
  presence_combined[c("x", "y")]
)
presence_data <- bind_cols(presence_combined, presence_values[, -1]) # Remove ID column

# Remove any presence points with NA values
presence_data_complete <- presence_data |>
  drop_na()

# Create presence data frame
presence_df <- presence_data_complete |>
  mutate(response = 1) |>
  select(-x, -y) # Remove coordinates, keep only predictors and response

# Create background data frame (sample for computational efficiency)
set.seed(123) # For reproducibility
background_df <- coords_global |>
  mutate(response = 0) |>
  sample_n(1e4) |> # Sample 10000 background points
  select(-x, -y) # Remove coordinates, keep only predictors and response

# Combine into training dataset
training_data <- bind_rows(presence_df, background_df) |>
  mutate(response = factor(response))
# Ensure explicit levels for response variable
training_data$response <- factor(
  as.character(training_data$response),
  levels = c("0", "1")
)

# Convert to data frame for `randomForestSRC`
training_data <- as.data.frame(training_data)

# Get final counts
n_final_presence <- sum(training_data$response == "1")
n_final_background <- sum(training_data$response == "0")

cat(
  "Final training data - Presences:",
  n_final_presence,
  "Background:",
  n_final_background,
  "\n"
)

# Fit balanced random forest model (BRF) ####

set.seed(456) # For reproducibility
rf_model_brf <- randomForest(
  response ~ .,
  data = training_data,
  ntree = 1000,
  sampsize = c(n_final_presence, n_final_presence), # Equal sampling for both classes
  replace = TRUE, # With replacement for balanced sampling
  importance = TRUE
)

print(rf_model_brf)

# Fit quantile random forest model (QRF) ####

## Estimate run time for QRF model ####
rt <- system.time({
  temp_model <- randomForestSRC::imbalanced(
    response ~ .,
    data = training_data,
    ntree = 1000,
    importance = FALSE
  )
})

cat(
  "\nEstimated time for QRF model fitting with 3000 trees:",
  round(((3000 / 1000) * rt[3]) / 60, 1),
  "minutes\n"
)

## Actual model fitting ####
set.seed(456) # For reproducibility
rf_model_qrf <- randomForestSRC::imbalanced(
  response ~ .,
  data = training_data,
  ntree = 3000, # More trees for quantile regression
  importance = TRUE
)

print(rf_model_qrf)

# Model evaluation using Norwegian regional data ####

library(pROC)

# Function to evaluate a model
evaluate_model <- function(pred_raster, model_name) {
  # Extract predictions for evaluation points
  presence_predictions <- terra::extract(
    pred_raster,
    presence_coords_regional[c("x", "y")],
    ID = FALSE
  )
  absence_predictions <- terra::extract(
    pred_raster,
    absence_coords_regional[c("x", "y")],
    ID = FALSE
  )

  # Create evaluation dataframe
  eval_data <- bind_rows(
    # Presence data (observed = 1)
    bind_cols(presence_coords_regional, presence_predictions) |>
      rename(predicted_prob = 3) |> # Rename the prediction column
      mutate(observed = 1),
    # Absence data (observed = 0)
    bind_cols(absence_coords_regional, absence_predictions) |>
      rename(predicted_prob = 3) |> # Rename the prediction column
      mutate(observed = 0)
  ) |>
    drop_na() # Remove any points with NA prediction values

  # Calculate AUC
  auc_result <- roc(eval_data$observed, eval_data$predicted_prob)
  auc_value <- auc(auc_result)

  # TSS - find optimal threshold using Youden's index
  optimal_threshold <- coords(auc_result, "best", ret = "threshold") |>
    as.numeric()
  eval_data <- eval_data |>
    mutate(predicted_binary = predicted_prob >= optimal_threshold)

  # Calculate confusion matrix components
  tp <- sum(eval_data$observed == 1 & eval_data$predicted_binary == TRUE)
  tn <- sum(eval_data$observed == 0 & eval_data$predicted_binary == FALSE)
  fp <- sum(eval_data$observed == 0 & eval_data$predicted_binary == TRUE)
  fn <- sum(eval_data$observed == 1 & eval_data$predicted_binary == FALSE)

  # Calculate metrics
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  tss_value <- sensitivity + specificity - 1

  # Return results
  data.frame(
    Model = model_name,
    Metric = c("AUC", "TSS"),
    Value = c(auc_value, tss_value)
  )
}

# Load Norwegian presence-absence data
presence_coords_regional <- read_csv("output/presence_coords_regional.csv")
absence_coords_regional <- read_csv("output/absence_coords_regional.csv")

# Load regional predictors for predictions
predictors_current <- rast(
  "output/predictors_regional_250m_Norway_current_EPSG3035.tif"
)
predictors_future <- rast(
  "output/predictors_regional_250m_Norway_future_EPSG3035.tif"
)

# Generate predictions for evaluation
pred_current_brf <- terra::predict(
  predictors_current,
  rf_model_brf,
  type = "prob",
  index = 2
)

# Define custom prediction function for QRF
predict_qrf <- function(model, data) {
  pred <- predict.rfsrc(model, newdata = data, importance = FALSE)
  return(pred$predicted[, "1"]) # Extract probability for class "1"
}
pred_current_qrf <- terra::predict(
  predictors_current,
  rf_model_qrf,
  fun = predict_qrf,
  na.rm = TRUE
)

# Evaluate BRF model
cat("\nEvaluating BRF model...\n")
eval_brf <- evaluate_model(pred_current_brf, "BRF")
print(eval_brf)

# Evaluate QRF model
cat("\nEvaluating QRF model...\n")
eval_qrf <- evaluate_model(pred_current_qrf, "QRF")
print(eval_qrf)

# Combine results
evaluation_results <- bind_rows(eval_brf, eval_qrf)

cat("\n=== Model Evaluation Results ===\n")
print(evaluation_results)

# Save evaluation results
write_csv(evaluation_results, "output/model_global_evaluation_metrics.csv")

# Choose algorithm based on evaluation

chosen <- evaluation_results |>
  pivot_wider(names_from = Metric, values_from = Value) |>
  arrange(desc(AUC), desc(TSS)) |>
  slice(1) |>
  pull(Model)

# Generate predictions ####

if (chosen == "BRF") {
  cat("\nChosen model: Balanced Random Forest (BRF)\n")
  # Generate BRF predictions (probability of class "1")
  pred_global_brf <- terra::predict(
    predictors_global,
    rf_model_brf,
    type = "prob",
    index = 2
  )
  pred_future_brf <- terra::predict(
    predictors_future,
    rf_model_brf,
    type = "prob",
    index = 2
  )
} else if (chosen == "QRF") {
  cat("\nChosen model: Quantile Random Forest (QRF)\n")
  # Generate QRF predictions (probability of class "1")
  # Note: randomForestSRC requires custom prediction function
  pred_global_qrf <- terra::predict(
    predictors_global,
    rf_model_qrf,
    fun = predict_qrf,
    na.rm = TRUE
  )
  pred_future_qrf <- terra::predict(
    predictors_future,
    rf_model_qrf,
    fun = predict_qrf,
    na.rm = TRUE
  )
}

# Save outputs ####

if (chosen == "BRF") {
  # Save BRF predictions
  writeRaster(
    pred_global_brf,
    "output/rf_global_pred_global_current.tif",
    overwrite = TRUE
  )
  writeRaster(
    pred_current_brf,
    "output/rf_global_pred_regional_current.tif",
    overwrite = TRUE
  )
  writeRaster(
    pred_future_brf,
    "output/rf_global_pred_regional_future.tif",
    overwrite = TRUE
  )
  saveRDS(rf_model_brf, "output/rf_global_model.rds")
} else if (chosen == "QRF") {
  # Save QRF predictions
  writeRaster(
    pred_global_qrf,
    "output/rf_global_pred_global_current.tif",
    overwrite = TRUE
  )
  writeRaster(
    pred_current_qrf,
    "output/rf_global_pred_regional_current.tif",
    overwrite = TRUE
  )
  writeRaster(
    pred_future_qrf,
    "output/rf_global_pred_regional_future.tif",
    overwrite = TRUE
  )
  saveRDS(rf_model_qrf, "output/rf_global_model.rds")
}

# sessionInfo ####

sessioninfo::session_info()
