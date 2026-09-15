# Production Model Training and Future Predictions ####

# PURPOSE: Fits the production model on the whole pooled frame, writes the frozen DI ruler, and projects the 3-class simplex over Norway for both scenarios.

# Trains the production model on ALL training data (no partitions) and generates current
# and future predictions over the Norwegian projection domain.
#
# Changes with the EU integration (plan_EUintegration.md):
#
#  - One model, not two. RFQ is two-class only and the response is 3-class, so the
#    production learner is a balanced `rfsrc()` using the package's own recipe
#    (notebook 2026-06-18). Predictions are a three-layer probability simplex per
#    scenario rather than a single presence probability.
#  - Training is pooled (EU + Norway); prediction is Norway only. That asymmetry is the
#    architecture: the EU block exists to give the model labelled data on the warm flank
#    Norway is moving into, not to be projected itself.
#  - The frozen DI ruler is NOT built here any more. It moved to pl2_freezeDIRuler.R so
#    that pl2_evaluate.R, which runs earlier, can measure its cross-validation on the same
#    axis the projection is read on (section 1.4).

library(readr)
library(dplyr)
library(randomForestSRC)
library(terra)

source("R/functions.R")
source("R/config.R")

## Configuration ####

NTREE <- 1000

record_settings(
  "R/pl2_predict.R",
  ntree = NTREE,
  seed = -42
)

## Load data ####

mf <- read_csv("output/pl2/modeling_frame_regional.csv", show_col_types = FALSE)

# Prepare training data (current scenario only, all observations)
train_data <- mf |>
  filter(scenario == "current") |>
  select(-scenario, -x, -y) |>
  mutate(response = factor(response, levels = RESPONSE_LEVELS)) |>
  as.data.frame()

cat("Training data loaded:", nrow(train_data), "observations\n")
print(table(train_data$dataset, train_data$response))
cat("\n")

# `dataset` is bookkeeping, never a predictor -- and here there is a second reason
# beyond section 1.5's: the prediction rasters carry no such layer, so a model that
# split on it could not be applied to them at all.
predictor_names <- setdiff(names(train_data), c("response", "dataset"))

# Load scenario rasters for spatial prediction
current_rasters <- rast("output/pl2/scenario_current.tif")
future_rasters <- rast("output/pl2/scenario_future.tif")
cat("Current scenario rasters loaded:", nlyr(current_rasters), "layers\n")
cat("Future scenario rasters loaded:", nlyr(future_rasters), "layers\n\n")

# The model can only be applied to a raster stack that carries every predictor it was
# trained on; fail here rather than inside terra::predict's block loop.
missing_layers <- setdiff(predictor_names, names(current_rasters))
if (length(missing_layers) > 0) {
  stop(
    "Predictors missing from the scenario rasters: ",
    paste(missing_layers, collapse = ", ")
  )
}
stopifnot(identical(names(current_rasters), names(future_rasters)))

## Train the production model ####

cat("Training balanced 3-class RF (ntree =", NTREE, ")...\n")

train_data_model <- train_data[, c("response", predictor_names)]

start_time <- Sys.time()

model_brf <- rfsrc(
  formula = response ~ .,
  data = train_data_model,
  ntree = NTREE,
  case.wt = randomForestSRC:::make.wt(train_data_model$response),
  sampsize = randomForestSRC:::make.size(train_data_model$response),
  importance = TRUE,
  seed = -42
)

cat(
  "Training completed in",
  round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1),
  "minutes.\n\n"
)

# Extract and save variable importance, per class and overall. This is the same
# quantity pl2_weightFeaturesDataPartitioning.R medians over repeated seeds to form the
# frozen ruler; saved here so the production fit carries its own record of it.
importance_df <- model_brf$importance |>
  as.data.frame() |>
  tibble::rownames_to_column("variable") |>
  arrange(desc(all))

write_csv(importance_df, "output/pl2/variable_importance_brf.csv")

cat("Top 10 most important variables:\n")
print(head(importance_df, 10))
cat("\n")

cat("Out-of-bag confusion matrix (rows = actual):\n")
oob_class <- RESPONSE_LEVELS[max.col(
  model_brf$predicted.oob[, RESPONSE_LEVELS],
  ties.method = "first"
)]
oob_metrics <- calculate_multiclass_metrics(
  predicted_class = oob_class,
  true_class = train_data_model$response,
  predicted_prob = model_brf$predicted.oob,
  levels = RESPONSE_LEVELS
)
print(oob_metrics$confusion)
cat("\n")
print(as.data.frame(oob_metrics$by_class), row.names = FALSE)

# Saved, not only printed: the production fit's own skill is quoted in the manuscript and
# compared against the cross-validation, so it has to exist as a file, not as knit output.
oob_metrics$by_class |>
  mutate(
    accuracy = oob_metrics$summary$accuracy,
    Gmean_macro = oob_metrics$summary$Gmean_macro,
    macro_auc = oob_metrics$summary$macro_auc
  ) |>
  write_csv("output/pl2/oob_metrics_production.csv", append = FALSE)
oob_metrics$confusion |>
  as.data.frame() |>
  write_csv("output/pl2/oob_confusion_production.csv", append = FALSE)
cat("\nOOB macro G-mean:", round(oob_metrics$summary$Gmean_macro, 4), "\n\n")

## Generate spatial predictions ####

cat("Generating spatial predictions...\n")

# Returns the full simplex, so each scenario yields one layer per class rather than a
# single presence probability. The three layers sum to 1 and the endpoint classes are
# what make the bog layer interpretable as an intermediate.
predfun <- function(model, data) {
  pred <- predict.rfsrc(model, newdata = data, importance = FALSE)
  pred$predicted[, RESPONSE_LEVELS, drop = FALSE]
}

cat("Predicting current...\n")
pred_current <- terra::predict(
  current_rasters,
  model_brf,
  fun = predfun,
  na.rm = TRUE
)
names(pred_current) <- paste0("current_", RESPONSE_LEVELS)

cat("Predicting future...\n")
pred_future <- terra::predict(
  future_rasters,
  model_brf,
  fun = predfun,
  na.rm = TRUE
)
names(pred_future) <- paste0("future_", RESPONSE_LEVELS)

pred_combined <- c(pred_current, pred_future)
writeRaster(
  pred_combined,
  "output/pl2/rf_local_pred_brf.tif",
  overwrite = TRUE
)
cat("Saved predictions to output/pl2/rf_local_pred_brf.tif\n\n")

plot(pred_combined)

## Save production model ####

saveRDS(model_brf, "output/pl2/model_production_brf.rds")
cat("Production model saved.\n")

# sessionInfo ####

sessioninfo::session_info()
