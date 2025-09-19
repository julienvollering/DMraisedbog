# Area of Applicability Assessment for Local-Scale Future Predictions
# Based on environmental distances between CV folds during training
# Using CAST::aoa() with variable importance weighting

# Load required packages
library(tidyverse)
library(terra)
library(CAST)
library(randomForestSRC)

# Load the trained local-scale model and training data
final_model <- readRDS("output/pl1/final_model_local_80split.rds")
var_importance <- final_model$importance[, 1]
training_data <- read_csv("output/pl1/data_partitioned.csv")

# Remove spatial coordinates and partition columns for environmental space
predictor_cols <- c("rf_global", paste0("bio", 1:19), "elevation", "slope", 
                   paste0("artype_", c(10, 20, 30, 50, 60, 70, 81)))

# Prepare training environmental data - only use outer == "train" partition
available_predictors <- predictor_cols[predictor_cols %in% colnames(training_data)]
training_env <- training_data |>
  filter(outer == "train") |>
  select(all_of(available_predictors), inner) |>
  drop_na()

# Match variable importance to available predictors
importance_weights <- var_importance[available_predictors] |> 
  t() |> 
  as.data.frame()

# Load future scenario predictors
rf_global_future <- rast("output/rf_global_pred_regional_future.tif") 
names(rf_global_future) <- "rf_global"
preds_nor_250m_future <- rast("output/predictors_regional_250m_Norway_future_EPSG3035.tif")
predictors_future <- c(rf_global_future, preds_nor_250m_future)
cat("Combined predictor stack created with", nlyr(predictors_future), "layers.\n")

set.seed(42) # For reproducibility
training_env_smaller <- training_env |> 
  slice_sample(n = min(10e3, nrow(training_env)), replace = FALSE)
cat("Using", nrow(training_env_smaller), "training samples for AOA calculation.\n")
cat("N per fold (inner):\n")
print(table(training_env_smaller$inner))

# Test AOA calculation on small subset (1%) to estimate runtime
cat("Testing AOA on 1% sample to estimate runtime...\n")
predictors_sample <- spatSample(
  predictors_future, 
  size = ncell(predictors_future) * 0.01, 
  method = "regular", as.raster = TRUE)

start_time <- Sys.time()
aoa_test <- aoa(newdata = predictors_sample,
                train = training_env_smaller[, available_predictors], 
                variables = available_predictors,
                weight = importance_weights,
                CVtest = training_env_smaller$inner)
test_time <- Sys.time() - start_time

# Estimate full runtime
full_cells <- ncell(predictors_future)
sample_cells <- ncell(predictors_sample)
estimated_time <- test_time * (full_cells / sample_cells)

cat("Test completed in", round(as.numeric(test_time, units = "secs"), 2), "seconds\n")
cat("Sample cells:", sample_cells, "Full cells:", full_cells, "\n")
cat("Estimated full runtime:", round(as.numeric(estimated_time, units = "mins"), 1), "minutes\n")

# Calculate Area of Applicability for future scenario (full data)
cat("Proceeding with full AOA calculation...\n")
aoa_future <- aoa(newdata = predictors_future,
                  train = training_env_smaller[, available_predictors], 
                  variables = available_predictors,
                  weight = importance_weights,
                  CVtest = training_env_smaller$inner)

plot(aoa_future)

# Save AOA object
saveRDS(aoa_future, "output/pl1/aoa_future_local_scale.rds")

# Save AOA results as rasters
writeRaster(aoa_future$AOA,
           "output/pl1/aoa_future_local_scale.tif", 
           overwrite = TRUE)

writeRaster(aoa_future$DI,
           "output/pl1/dissimilarity_index_future_local_scale.tif",
           overwrite = TRUE)

# Calculate summary statistics for future scenario
future_aoa_area <- sum(values(aoa_future$AOA), na.rm = TRUE) * res(aoa_future$AOA)[1]^2 / 1e6 # km²
future_total_area <- sum(!is.na(values(aoa_future$AOA))) * res(aoa_future$AOA)[1]^2 / 1e6 # km²
future_aoa_percent <- future_aoa_area / future_total_area * 100

future_di_stats <- c(
  min = min(values(aoa_future$DI), na.rm = TRUE),
  mean = mean(values(aoa_future$DI), na.rm = TRUE),
  max = max(values(aoa_future$DI), na.rm = TRUE), 
  sd = sd(values(aoa_future$DI), na.rm = TRUE)
)

# Create summary data frame
aoa_summary <- data.frame(
  scenario = "future",
  aoa_area_km2 = future_aoa_area,
  total_area_km2 = future_total_area,
  aoa_percent = future_aoa_percent,
  di_min = future_di_stats[["min"]],
  di_mean = future_di_stats[["mean"]],
  di_max = future_di_stats[["max"]],
  di_sd = future_di_stats[["sd"]]
)

# Save summary
write_csv(aoa_summary, "output/pl1/aoa_summary_local_scale.csv")

# Print summary
print(aoa_summary)
importance_weights |> 
  t() |> 
  as.data.frame() |> 
  rownames_to_column("variable") |> 
  arrange(desc(V1))

# sessionInfo

sessioninfo::session_info()
