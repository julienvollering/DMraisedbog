# Local Scale Random Forest Quantile Classification Model ####

# This script trains the final random forest quantile classification model
# using CV-optimized hyperparameters from tuneHyperparametersLocal.R.
# Uses the full training partition for final model training.

library(readr)
library(dplyr)
library(purrr)
library(randomForestSRC)
library(terra)

## Load data and hyperparameters ####

cat("Loading partitioned data and CV-optimized hyperparameters...\n")

# Read main partitioned dataset
training_data_partitioned <- read_csv("output/training_data_partitioned.csv")

# Extract training partition for final model
train_data <- training_data_partitioned |>
  filter(outer == "train") |>
  select(-outer, -inner, -x, -y)  # Remove partition columns and coordinates

cat("Training data loaded:\n")
cat("Observations:", nrow(train_data), "\n")
cat("Response variable class:", class(train_data$response), "\n")
cat("Response variable values:", paste(unique(train_data$response), collapse = ", "), "\n")
cat("Prevalence:", round(mean(train_data$response == 1) * 100, 2), "%\n")

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

## Train final model with selected parameters ####

# Prepare full training data as data.frame (imbalanced() doesn't handle tibbles)
full_train_data <- train_data |>
  mutate(response = factor(response, levels = c("1", "0"))) |> # make "1" the positive class
  as.data.frame()  # Convert tibble to data.frame for imbalanced()

# Ensure ar50 is properly handled as factor if it exists
if("ar50" %in% names(full_train_data)) {
  full_train_data$ar50 <- as.factor(as.character(full_train_data$ar50))
}

n_train <- 10e3
cat("Training final model on", n_train, "observations from training partition...\n")

# Randomly subsample training data for runtime efficiency
set.seed(42)
small_train_data <- full_train_data |> 
  slice_sample(n = min(n_train, nrow(full_train_data)))

# Train final model with selected hyperparameters
start_time <- Sys.time()
final_model <- imbalanced(
  formula = response ~ .,
  data = small_train_data,
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
cat("Trained on", nrow(small_train_data), "of", nrow(full_train_data),
    "(", round(nrow(small_train_data) / nrow(full_train_data) * 100, 1), 
    "% ) observations.\n")

final_model

## Geographical visualization of model predictions ####
rf_global_current <- rast("output/rf_global_pred_regional_current.tif") 
names(rf_global_current) <- "rf_global"
preds_nor_250m_cur <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")
predictors <- c(rf_global_current, preds_nor_250m_cur)
predictors_se <- crop(predictors, ext(435e4, 447e4, 410e4, 420e4))

# Make spatial predictions for a cropped extent
predfun <- function(model, data) {
  v <- predict.rfsrc(model, data, importance = FALSE)
  cbind(p = v$predicted[,1], class = v$class) # Taking the first of two classes to be the positive class
}
spatial_preds <- terra::predict(
  predictors_se, 
  final_model, 
  fun = predfun, 
  na.rm = TRUE
)
plot(spatial_preds$p)
plot(spatial_preds$class)
prevalence_training <- mean(small_train_data$response == 1)
plot(spatial_preds$p > prevalence_training)

# View variable importance
final_model$importance[,1] |> 
  tibble::enframe(name = "variable", value = "importance") |> 
  arrange(desc(importance)) |> 
  print(n = 29)

## Save results ####

model_path <- "output/final_model_local.rds"
cat("\nSaving final model to", model_path,"\n")

# Save the trained model object (RDS format for R objects)
saveRDS(final_model, model_path)

# sessionInfo ####

sessioninfo::session_info()
