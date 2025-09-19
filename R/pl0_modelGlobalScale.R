library(tidyverse)
library(terra)
library(randomForest)

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
presence_values <- terra::extract(predictors_global, presence_combined[c("x", "y")])
presence_data <- bind_cols(presence_combined, presence_values[, -1])  # Remove ID column

# Remove any presence points with NA values
presence_data_complete <- presence_data |>
  drop_na()

# Create presence data frame
presence_df <- presence_data_complete |>
  mutate(response = 1) |>
  select(-x, -y)  # Remove coordinates, keep only predictors and response

# Create background data frame (all non-presence cells)
background_df <- coords_global |>
  mutate(response = 0) |>
  select(-x, -y)  # Remove coordinates, keep only predictors and response

# Combine into training dataset
training_data <- bind_rows(presence_df, background_df) |>
  mutate(response = factor(response))

# Get final counts
n_final_presence <- sum(training_data$response == "1")
n_final_background <- sum(training_data$response == "0")

cat("Final training data - Presences:", n_final_presence, "Background:", n_final_background, "\n")

# Fit balanced random forest model ####

set.seed(456)  # For reproducibility
rf_model <- randomForest(
  response ~ .,
  data = training_data,
  ntree = 1000,
  sampsize = c(n_final_presence, n_final_presence),  # Equal sampling for both classes
  replace = TRUE,  # With replacement for balanced sampling
  importance = TRUE
)

print(rf_model)

# Generate predictions ####

# Load regional predictors for predictions
predictors_current <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")
predictors_future <- rast("output/predictors_regional_250m_Norway_future_EPSG3035.tif")

# Generate predictions for current and future scenarios
pred_global <- terra::predict(predictors_global, rf_model, type = "prob", index = 2)
pred_current <- terra::predict(predictors_current, rf_model, type = "prob", index = 2)
pred_future <- terra::predict(predictors_future, rf_model, type = "prob", index = 2)

# Save outputs ####

writeRaster(pred_global, "output/rf_global_pred_global_current.tif", overwrite = TRUE)
writeRaster(pred_current, "output/rf_global_pred_regional_current.tif", overwrite = TRUE)
writeRaster(pred_future, "output/rf_global_pred_regional_future.tif", overwrite = TRUE)
saveRDS(rf_model, "output/rf_global_model.rds")

# Model evaluation using Norwegian regional data ####

# Load Norwegian presence-absence data
presence_coords_regional <- read_csv("output/presence_coords_regional.csv")
absence_coords_regional <- read_csv("output/absence_coords_regional.csv")

# Load the regional prediction raster
pred_regional <- rast("output/rf_global_pred_regional_current.tif")

# Extract predictions for evaluation points
presence_predictions <- terra::extract(pred_regional, presence_coords_regional[c("x", "y")])
absence_predictions <- terra::extract(pred_regional, absence_coords_regional[c("x", "y")])

# Create evaluation dataframe
eval_data <- bind_rows(
  # Presence data (observed = 1)
  bind_cols(presence_coords_regional, presence_predictions[, -1]) |>
    rename(predicted_prob = 3) |>  # Rename the prediction column
    mutate(observed = 1),
  # Absence data (observed = 0)  
  bind_cols(absence_coords_regional, absence_predictions[, -1]) |>
    rename(predicted_prob = 3) |>  # Rename the prediction column
    mutate(observed = 0)
) |>
drop_na()  # Remove any points with NA prediction values

# Calculate evaluation metrics
# AUC using pROC package
library(pROC)
auc_result <- roc(eval_data$observed, eval_data$predicted_prob)
auc_value <- auc(auc_result)

# TSS (True Skill Statistic) - find optimal threshold using Youden's index
optimal_threshold <- coords(auc_result, "best", ret = "threshold") |> 
  as.numeric()
eval_data <- eval_data |> 
  mutate(predicted_binary = predicted_prob >= optimal_threshold)

# Calculate confusion matrix components
tp <- sum(eval_data$observed == 1 & eval_data$predicted_binary == TRUE)  # True positives
tn <- sum(eval_data$observed == 0 & eval_data$predicted_binary == FALSE)  # True negatives
fp <- sum(eval_data$observed == 0 & eval_data$predicted_binary == TRUE)  # False positives
fn <- sum(eval_data$observed == 1 & eval_data$predicted_binary == FALSE)  # False negatives

# Calculate sensitivity and specificity
sensitivity <- tp / (tp + fn)  # True positive rate
specificity <- tn / (tn + fp)  # True negative rate

# TSS = Sensitivity + Specificity - 1
tss_value <- sensitivity + specificity - 1

# Create results summary
evaluation_results <- data.frame(
  Metric = c("AUC", "TSS"),
  Value = c(auc_value, tss_value)
)

print("Model Evaluation Results:")
print(evaluation_results)

# Save evaluation results
write_csv(evaluation_results, "output/model_global_evaluation_metrics.csv")

## Session Information ####
sessionInfo()