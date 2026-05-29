# Production Model Training and Future Predictions ####

# This script trains production models on ALL training data (no partitions)
# and generates current and future predictions for both BRF and RFQ models.

library(readr)
library(dplyr)
library(randomForest)
library(randomForestSRC)
library(terra)

## Load data ####

mf <- read_csv("output/pl2/modeling_frame_regional.csv", show_col_types = FALSE)

# Prepare training data (current scenario only, all observations)
train_data <- mf |>
  filter(scenario == "current") |>
  select(-scenario, -x, -y) |>
  as.data.frame()

cat("Training data loaded:", nrow(train_data), "observations\n")
cat("  Presences:", sum(train_data$response == 1), "\n")
cat("  Absences:", sum(train_data$response == 0), "\n")

# Load scenario rasters for spatial prediction
current_rasters <- rast("output/pl2/scenario_current.tif")
future_rasters <- rast("output/pl2/scenario_future.tif")
cat("Current scenario rasters loaded:", nlyr(current_rasters), "layers\n")
cat("Future scenario rasters loaded:", nlyr(future_rasters), "layers\n\n")

## Prepare model data ####

train_data_model <- train_data
train_data_model$response <- factor(
  as.character(train_data_model$response),
  levels = c("0", "1")
)

## Train BRF model ####

cat("Training BRF model (ntree = 1000)...\n")

n_presence <- sum(train_data_model$response == "1")
n_absence <- sum(train_data_model$response == "0")
balanced_size <- min(n_presence, n_absence)

cat("  Balanced sampling size:", balanced_size, "\n")

start_time <- Sys.time()

model_brf <- randomForest::randomForest(
  formula = response ~ .,
  data = train_data_model,
  ntree = 1000,
  sampsize = c(balanced_size, balanced_size),
  replace = TRUE,
  importance = TRUE,
  do.trace = FALSE
)

cat(
  "Training completed in",
  round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1),
  "minutes.\n\n"
)

# Extract and save variable importance
var_importance <- model_brf$importance[, 1]
importance_df <- data.frame(
  variable = names(var_importance),
  importance = as.numeric(var_importance),
  stringsAsFactors = FALSE
) |>
  arrange(desc(importance))

write_csv(importance_df, "output/pl2/variable_importance_brf.csv")

cat("Top 10 most important variables (BRF):\n")
print(head(importance_df, 10))
cat("\n")

## Train RFQ model ####

cat("Training RFQ model (ntree = 3000)...\n")

start_time <- Sys.time()

model_rfq <- randomForestSRC::imbalanced(
  formula = response ~ .,
  data = train_data_model,
  ntree = 3000,
  # See "R\pl2_weightFeaturesDataPartitioning.R"; RFQ VI is not reliable.
  importance = FALSE,
  do.trace = FALSE,
  seed = -42
)

cat(
  "Training completed in",
  round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1),
  "minutes.\n\n"
)

## Generate spatial predictions ####

cat("Generating spatial predictions...\n")

predfun_brf <- function(model, data) {
  predict(model, newdata = data, type = "prob")[, "1"]
}

predfun_rfq <- function(model, data) {
  pred <- predict.rfsrc(model, newdata = data, importance = FALSE)
  pred$predicted[, "1"]
}

# BRF predictions
cat("Predicting with BRF...\n")
brf_current <- terra::predict(current_rasters, model_brf, fun = predfun_brf, na.rm = TRUE)
brf_future <- terra::predict(future_rasters, model_brf, fun = predfun_brf, na.rm = TRUE)
brf_combined <- c(brf_current, brf_future)
names(brf_combined) <- c("current", "future")
writeRaster(brf_combined, "output/pl2/rf_local_pred_brf.tif", overwrite = TRUE)
cat("Saved BRF predictions to output/pl2/rf_local_pred_brf.tif\n\n")

# RFQ predictions
cat("Predicting with RFQ...\n")
rfq_current <- terra::predict(current_rasters, model_rfq, fun = predfun_rfq, na.rm = TRUE)
rfq_future <- terra::predict(future_rasters, model_rfq, fun = predfun_rfq, na.rm = TRUE)
rfq_combined <- c(rfq_current, rfq_future)
names(rfq_combined) <- c("current", "future")
writeRaster(rfq_combined, "output/pl2/rf_local_pred_rfq.tif", overwrite = TRUE)
cat("Saved RFQ predictions to output/pl2/rf_local_pred_rfq.tif\n\n")

## Save production models ####

saveRDS(model_brf, "output/pl2/model_production_brf.rds")
saveRDS(model_rfq, "output/pl2/model_production_rfq.rds")
cat("Production models saved.\n")

# sessionInfo ####

sessioninfo::session_info()
