# Production Model Training and Future Predictions ####

# This script trains production models on ALL training data (no partitions)
# and generates future predictions under climate change scenarios.
# Runs for three partitioning schemes (importance, pca, maxshift) using
# scheme-specific optimal hyperparameters from tuning.
#

library(readr)
library(dplyr)
library(randomForest)
library(randomForestSRC)
library(terra)

## Configuration ####

# Select which partitioning schemes to process
partitioning_schemes <- c("importance", "pca", "maxshift")

## Load shared data (used by all schemes) ####

# Read full modeling frame
mf <- read_csv("output/pl2/modeling_frame_regional.csv", show_col_types = FALSE)

# Prepare training data (current scenario only, all observations)
train_data <- mf |>
  filter(scenario == "current") |>
  select(-scenario, -x, -y) |>
  as.data.frame()

cat("Training data loaded:", nrow(train_data), "observations\n")
cat("  Presences:", sum(train_data$response == 1), "\n")
cat("  Absences:", sum(train_data$response == 0), "\n")

prevalence <- mean(train_data$response == 1)
cat("Training data prevalence:", prevalence, "\n")

# Load scenario rasters for spatial prediction
current_rasters <- rast("output/pl2/scenario_current.tif")
future_rasters <- rast("output/pl2/scenario_future.tif")
cat("Current scenario rasters loaded:", nlyr(current_rasters), "layers\n")
cat("Future scenario rasters loaded:", nlyr(future_rasters), "layers\n\n")

## Load all hyperparameters to detect duplicates ####

cat("Loading hyperparameters for all schemes...\n")

all_hyperparams <- lapply(partitioning_schemes, function(scheme) {
  hyperparams_file <- paste0("output/pl2/hyperparameters_", scheme, ".csv")
  hp <- read_csv(hyperparams_file, show_col_types = FALSE)
  hp$scheme <- scheme
  return(hp)
}) |>
  bind_rows()

# Create unique hyperparameter signature for each scheme
all_hyperparams$hp_signature <- paste(
  all_hyperparams$model_type,
  all_hyperparams$mtry,
  all_hyperparams$nodesize,
  sep = "_"
)

cat("\nHyperparameter summary:\n")
print(all_hyperparams[, c("scheme", "model_type", "mtry", "nodesize", "hp_signature")])
cat("\n")

# Identify unique hyperparameter configurations
unique_signatures <- unique(all_hyperparams$hp_signature)
cat("Found", length(unique_signatures), "unique hyperparameter configuration(s)\n\n")

## Prepare training data (common for all models) ####

# Convert response to factor with explicit levels
train_data_model <- train_data
train_data_model$response <- factor(as.character(train_data_model$response),
  levels = c("0", "1")
)

## Train models for unique hyperparameter configurations ####

trained_models <- list()

for (sig in unique_signatures) {

  # Get first scheme with this signature (for labeling)
  sig_schemes <- all_hyperparams$scheme[all_hyperparams$hp_signature == sig]
  hp_row <- all_hyperparams[all_hyperparams$hp_signature == sig, ][1, ]

  model_type <- hp_row$model_type
  mtry <- hp_row$mtry
  nodesize <- hp_row$nodesize

  cat("Training model for hyperparameters:", sig, "\n")
  cat("  (Used by schemes:", paste(sig_schemes, collapse = ", "), ")\n")
  cat("  Model type:", model_type, "\n")
  cat("  mtry:", mtry, "\n")
  cat("  nodesize:", nodesize, "\n")

  start_time <- Sys.time()

  if (model_type == "BRF") {
    # Calculate class sizes for balanced sampling
    n_presence <- sum(train_data_model$response == "1")
    n_absence <- sum(train_data_model$response == "0")
    balanced_size <- min(n_presence, n_absence)

    # Train BRF model with balanced sampling and importance
    production_model <- randomForest::randomForest(
      formula = response ~ .,
      data = train_data_model,
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
    production_model <- randomForestSRC::imbalanced(
      formula = response ~ .,
      data = train_data_model,
      ntree = 3000,
      mtry = mtry,
      nodesize = nodesize,
      importance = TRUE,
      do.trace = 100,
      seed = -42
    )
  }

  end_time <- Sys.time()
  training_duration <- as.numeric(difftime(end_time, start_time, units = "mins"))

  cat("Training completed in", round(training_duration, 1), "minutes.\n\n")

  # Store trained model
  trained_models[[sig]] <- list(
    model = production_model,
    model_type = model_type,
    training_duration = training_duration
  )
}

## Main loop over partitioning schemes ####

for (scheme in partitioning_schemes) {

  cat("Processing scheme:", toupper(scheme), "\n")

  ### Get hyperparameters and corresponding model ####

  hp_row <- all_hyperparams[all_hyperparams$scheme == scheme, ]
  sig <- hp_row$hp_signature
  model_type <- hp_row$model_type

  cat("Using model with hyperparameters:", sig, "\n")

  # Retrieve trained model
  production_model <- trained_models[[sig]]$model
  training_duration <- trained_models[[sig]]$training_duration

  ### Extract and save variable importance ####

  cat("Extracting variable importance...\n")

  if (model_type == "BRF") {
    # BRF: Extract MeanDecreaseAccuracy (first column)
    var_importance <- production_model$importance[, 1]
    importance_df <- data.frame(
      variable = names(var_importance),
      importance = as.numeric(var_importance),
      stringsAsFactors = FALSE
    ) |>
      arrange(desc(importance))
  } else if (model_type == "RFQ") {
    # RFQ: Extract importance values (single column in randomForestSRC)
    var_importance <- production_model$importance[, 1]
    importance_df <- data.frame(
      variable = names(var_importance),
      importance = as.numeric(var_importance),
      stringsAsFactors = FALSE
    ) |>
      arrange(desc(importance))
  }

  # Save variable importance
  importance_file <- paste0("output/pl2/variable_importance_", scheme, ".csv")
  write_csv(importance_df, importance_file)

  # Display top 10 most important variables
  cat("\nTop 10 most important variables:\n")
  print(head(importance_df, 10))
  cat("\n")

  ### Generate spatial predictions ####

  cat("Generating spatial predictions...\n")

  # Define prediction function
  if (model_type == "BRF") {
    predfun <- function(model, data) {
      # Return probability for positive class ("1")
      pred_prob <- predict(model, newdata = data, type = "prob")[, "1"]
      return(pred_prob)
    }
  } else if (model_type == "RFQ") {
    predfun <- function(model, data) {
      # Return probability for positive class ("1")
      pred <- predict.rfsrc(model, newdata = data, importance = FALSE)
      pred_prob <- pred$predicted[, "1"] # Taking 1's as positive class
      return(pred_prob)
    }
  }

  # Generate current scenario predictions
  current_predictions <- terra::predict(
    current_rasters,
    production_model,
    fun = predfun,
    na.rm = TRUE
  )

  # Generate future scenario predictions
  future_predictions <- terra::predict(
    future_rasters,
    production_model,
    fun = predfun,
    na.rm = TRUE
  )

  # Combine into multi-band raster and save
  predictions_combined <- c(current_predictions, future_predictions)
  names(predictions_combined) <- c("current", "future")

  writeRaster(
    predictions_combined,
    filename = paste0("output/pl2/rf_local_pred_", scheme, ".tif"),
    overwrite = TRUE
  )

  ### Save production model ####

  model_file <- paste0("output/pl2/model_production_", scheme, ".rds")
  saveRDS(production_model, model_file)
  cat("Completed scheme:", toupper(scheme), "\n\n")
}

cat("All schemes completed.\n")

# sessionInfo ####

sessioninfo::session_info()
