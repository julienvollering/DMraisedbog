# Calibrate and Evaluate Local Scale Model Predictions ####

# This script calibrates the raw predictions from the RFQ model trained in
# modelLocalScale.R, using cross-validation to select
# between beta calibration and no calibration. Uses 5-fold CV within the
# calibration partition to compare methods based on Brier score.

# Custom thresholds (95% specificity and 95% sensitivity) are calculated from
# the combined training and calibration partitions to ensure robust estimates
# based on the full development dataset.

# After calibrating, the model is evaluated on the test partition, with both
# probability-based and class-based metrics using multiple threshold strategies

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(probably)
library(randomForestSRC)
library(yardstick)
library(rsample)
library(purrr)
library(stringr)
library(patchwork)
library(ROCR)

## Load data and trained model ####

cat("Loading partitioned data and trained model...\n")

# Read main partitioned dataset
training_data_partitioned <- read_csv("output/training_data_partitioned.csv")

# Load the trained model
final_model <- readRDS("output/final_model_local_60split.rds")

cat("Model loaded successfully.\n")
cat("Training data dimensions:", nrow(training_data_partitioned), "x", ncol(training_data_partitioned), "\n")

# Extract partitions
train_data <- training_data_partitioned |>
  filter(outer == "train") |>
  select(-outer, -inner, -x, -y)  # Remove partition columns and coordinates

calib_data <- training_data_partitioned |>
  filter(outer == "calib") |>
  select(-outer, -inner, -x, -y)  # Remove partition columns and coordinates

test_data <- training_data_partitioned |>
  filter(outer == "test") |>
  select(-outer, -inner, -x, -y)  # Remove partition columns and coordinates

cat("Training set:", nrow(train_data), "observations\n")
cat("Calibration set:", nrow(calib_data), "observations\n")
cat("Test set:", nrow(test_data), "observations\n")
cat("Training prevalence:", round(mean(train_data$response == 1) * 100, 2), "%\n")
cat("Calibration prevalence:", round(mean(calib_data$response == 1) * 100, 2), "%\n")
cat("Test prevalence:", round(mean(test_data$response == 1) * 100, 2), "%\n")

## Step 1: Get predictions on all partitions (training, calibration, test) ####

cat("\nGenerating predictions on training partition...\n")

# Prepare training data as data.frame for randomForestSRC
train_df <- train_data |>
  mutate(response = factor(response, levels = c("1", "0"))) |>
  as.data.frame()

# Ensure ar50 is properly handled as factor if it exists
if("ar50" %in% names(train_df)) {
  train_df$ar50 <- as.factor(as.character(train_df$ar50))
}

# Get training predictions (probability for positive class)
train_preds <- predict.rfsrc(final_model, train_df, importance = FALSE)

# Create training data frame for probably package
train_results <- tibble(
  truth = factor(train_data$response, levels = c("1", "0")),
  .pred_1 = train_preds$predicted[,1],  # Probability for positive class ("1")
  .pred_0 = train_preds$predicted[,2],  # Probability for negative class ("0")
  class = train_preds$class
)

cat("Predictions generated on training partition.\n")
cat("Raw prediction range:", round(min(train_preds$predicted[,1]), 4),
    "to", round(max(train_preds$predicted[,1]), 4), "\n")

cat("\nGenerating predictions on calibration partition...\n")

# Prepare calibration data as data.frame for randomForestSRC
calib_df <- calib_data |>
  mutate(response = factor(response, levels = c("1", "0"))) |>
  as.data.frame()

# Ensure ar50 is properly handled as factor if it exists
if("ar50" %in% names(calib_df)) {
  calib_df$ar50 <- as.factor(as.character(calib_df$ar50))
}

# Get calibration predictions (probability for positive class)
calib_preds <- predict.rfsrc(final_model, calib_df, importance = FALSE)

# Create calibration data frame for probably package
calib_results <- tibble(
  truth = factor(calib_data$response, levels = c("1", "0")),
  .pred_1 = calib_preds$predicted[,1],  # Probability for positive class ("1")
  .pred_0 = calib_preds$predicted[,2],  # Probability for negative class ("0")
  class = calib_preds$class
) 
# probably::cal_estimate_* "The order of the identifiers will be considered the
# same as the order of the levels of the truth variable."

cat("Predictions generated on calibration partition.\n")
cat("Raw prediction range:", round(min(calib_preds$predicted[,1]), 4), 
    "to", round(max(calib_preds$predicted[,1]), 4), "\n")

raw_cal_plot_partitioncalib <- cal_plot_windowed(
  slice_sample(calib_results, n = 100e3), 
  truth = truth, 
  estimate = .pred_1,
  event_level = "first",
  include_rug = FALSE)
print(raw_cal_plot_partitioncalib)

ecdf_partitioncalib <- ggplot(calib_results, aes(x = .pred_1)) +
  stat_ecdf() +
  labs(y = "Cumulative Probability") +
  theme_bw()

cat("\nGenerating predictions on test set...\n")

# Prepare test data as data.frame
test_df <- test_data |>
  mutate(response = factor(response, levels = c("1", "0"))) |>
  as.data.frame()

# Ensure ar50 is properly handled as factor if it exists
if("ar50" %in% names(test_df)) {
  test_df$ar50 <- as.factor(as.character(test_df$ar50))
}

# Get test predictions
test_preds <- predict.rfsrc(final_model, test_df, importance = FALSE)

# Create test results data frame
test_results_raw <- tibble(
  truth = factor(test_data$response, levels = c("1", "0")),
  .pred_1 = test_preds$predicted[,1],
  .pred_0 = test_preds$predicted[,2],
  class = test_preds$class
)

cat("Predictions generated on test partition.\n")
cat("Raw prediction range:", round(min(test_preds$predicted[,1]), 4),
    "to", round(max(test_preds$predicted[,1]), 4), "\n")

raw_cal_plot_partitiontest <- cal_plot_windowed(
  slice_sample(test_results_raw, n = 100e3), 
  truth = truth, 
  estimate = .pred_1,
  event_level = "first",
  include_rug = FALSE)
print(raw_cal_plot_partitiontest)

ecdf_partitiontest <- ggplot(test_results_raw, aes(x = .pred_1)) +
  stat_ecdf() +
  labs(y = "Cumulative Probability") +
  theme_bw()

patchwork::wrap_plots(
  raw_cal_plot_partitioncalib, raw_cal_plot_partitiontest,
  ecdf_partitioncalib, ecdf_partitiontest,
  ncol = 2, nrow = 2) +
  patchwork::plot_annotation(
    title = "Calibration (left) vs Test (right) partitions",
    subtitle = "Raw predictions only"
  ) + patchwork::plot_layout(guides = "collect")

## Step 2: Cross-validation to select between beta calibration and no calibration ####

cat("\nUsing 5-fold CV within calibration set to compare beta calibration vs no calibration...\n")
cat("Evaluating with Brier score (selection criterion), ECE, and log-loss.\n")

# Create CV folds within calibration set only
set.seed(42)  # For reproducible results
cal_folds <- vfold_cv(calib_results, v = 5)

# Helper function to calculate Expected Calibration Error (ECE)
calculate_ece <- function(data, truth_col, prob_col, n_bins = 10) {
  # Convert truth to numeric (1 for positive class, 0 for negative)
  truth_numeric <- as.numeric(data[[truth_col]] == levels(data[[truth_col]])[1])
  
  # Create bins and calculate ECE using tidy approach
  ece_data <- tibble(
    prob = data[[prob_col]],
    truth = truth_numeric
  ) |>
    mutate(
      bin = cut(prob, breaks = seq(0, 1, length.out = n_bins + 1), 
                include.lowest = TRUE, right = TRUE)
    ) |>
    group_by(bin) |>
    summarise(
      bin_accuracy = mean(truth),
      bin_confidence = mean(prob),
      bin_count = n(),
      .groups = "drop"
    ) |>
    mutate(
      bin_weight = bin_count / sum(bin_count),
      bin_ece = bin_weight * abs(bin_accuracy - bin_confidence)
    )
  
  return(sum(ece_data$bin_ece, na.rm = TRUE))
}

# Helper function to find custom thresholds using ROCR
find_custom_thresholds <- function(data, truth_col, prob_col) {
  # Convert truth to numeric (1 for positive class, 0 for negative)
  truth_numeric <- as.numeric(data[[truth_col]] == levels(data[[truth_col]])[1])
  
  # Create ROCR prediction object
  pred <- prediction(data[[prob_col]], truth_numeric)
  
  # Calculate sensitivity and specificity at all thresholds
  sens <- performance(pred, "sens")
  spec <- performance(pred, "spec")
  
  # Get cutoffs and performance values
  cutoffs <- sens@x.values[[1]]  # These are the thresholds
  sens_values <- sens@y.values[[1]]  # Sensitivity values
  spec_values <- spec@y.values[[1]]  # Specificity values
  
  # Create a data frame for easier manipulation
  perf_data <- tibble(
    threshold = cutoffs,
    sensitivity = sens_values,
    specificity = spec_values
  ) |>
    # Remove any NA values and sort by threshold
    filter(!is.na(threshold), !is.na(sensitivity), !is.na(specificity)) |>
    arrange(threshold)
  
  # Find thresholds for 95% specificity and 95% sensitivity
  # For 95% specificity: find the lowest threshold where specificity >= 0.95
  spec_95_candidates <- perf_data |> filter(specificity >= 0.95)
  threshold_95spec <- if(nrow(spec_95_candidates) > 0) min(spec_95_candidates$threshold) else NA
  
  # For 95% sensitivity: find the highest threshold where sensitivity >= 0.95
  sens_95_candidates <- perf_data |> filter(sensitivity >= 0.95)
  threshold_95sens <- if(nrow(sens_95_candidates) > 0) max(sens_95_candidates$threshold) else NA
  
  return(list(
    threshold_95spec = threshold_95spec,
    threshold_95sens = threshold_95sens,
    n_thresholds = nrow(perf_data),
    max_spec = max(perf_data$specificity, na.rm = TRUE),
    max_sens = max(perf_data$sensitivity, na.rm = TRUE)
  ))
}

# Function to evaluate beta calibration vs no calibration
compare_methods <- function(split) {
  analysis_data <- analysis(split)
  assessment_data <- assessment(split)
  
  # Fit beta calibration method on analysis set
  tryCatch({
    beta_cal <- cal_estimate_beta(
      analysis_data,
      truth = truth, 
      estimate = dplyr::starts_with(".pred_"),
      event_level = "first"
    )
  }, error = function(e) beta_cal <<- NULL)
  
  # Initialize results
  results <- tibble()
  
  # Evaluate no calibration (raw predictions)
  raw_ece <- calculate_ece(assessment_data, "truth", ".pred_1")
  raw_brier <- brier_class(assessment_data, truth = truth, .pred_1, 
                           event_level = "first")$.estimate
  raw_logloss <- mn_log_loss(assessment_data, truth = truth, .pred_1, 
                             event_level = "first")$.estimate
  results <- bind_rows(results, 
                      tibble(method = "none", ece = raw_ece, brier = raw_brier, logloss = raw_logloss))
  
  # Evaluate beta calibration (if successful)
  if (!is.null(beta_cal)) {
    beta_preds <- cal_apply(assessment_data, beta_cal)
    beta_ece <- calculate_ece(beta_preds, "truth", ".pred_1")
    beta_brier <- brier_class(beta_preds, truth = truth, .pred_1,
                              event_level = "first")$.estimate
    beta_logloss <- mn_log_loss(beta_preds, truth = truth, .pred_1,
                                event_level = "first")$.estimate
    results <- bind_rows(results,
                        tibble(method = "beta", ece = beta_ece, brier = beta_brier, logloss = beta_logloss))
  }
  
  return(results)
}

# Run cross-validation comparison
cv_results <- map_dfr(cal_folds$splits, compare_methods, .id = "fold")

# Summarize results across folds
method_performance <- cv_results |>
  group_by(method) |>
  summarise(
    mean_brier = mean(brier, na.rm = TRUE), 
    sd_brier = sd(brier, na.rm = TRUE),
    mean_logloss = mean(logloss, na.rm = TRUE),
    sd_logloss = sd(logloss, na.rm = TRUE),
    mean_ece = mean(ece, na.rm = TRUE),
    sd_ece = sd(ece, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(mean_brier)  # Lower Brier score is better

cat("\n=== CALIBRATION METHOD COMPARISON (5-fold CV) ===\n")
print(method_performance, n = Inf)

# Select best method based on lowest Brier score
best_method <- method_performance$method[1]
cat(sprintf("\nSelected method: %s (lowest mean Brier score: %.4f, mean ECE: %.4f, mean log-loss: %.4f)\n", 
            best_method, method_performance$mean_brier[1], 
            method_performance$mean_ece[method_performance$method == best_method],
            method_performance$mean_logloss[method_performance$method == best_method]))

## Step 3: Fit selected method on full calibration set ####

if (best_method == "none") {
  cat("\nNo calibration selected - using raw predictions.\n")
  final_cal <- NULL
} else {
  cat(sprintf("\nFitting %s calibration on full calibration set...\n", best_method))
  
  # Fit the selected method on full calibration set
  if (best_method == "beta") {
    final_cal <- cal_estimate_beta(
      calib_results,
      truth = truth,
      estimate = dplyr::starts_with(".pred_"),
      event_level = "first"
    )
    cat(sprintf("%s calibration fitted successfully.\n", str_to_title(best_method)))
    print(final_cal)
  }
}

# Apply selected method to calibration partition for diagnostics
# Note: this is optimistic since we are applying to the same data used to fit

if (best_method == "none") {
  cat("Using raw predictions for calibration set diagnostics:\n")
  calib_results_calibrated <- calib_results
} else {
  cat("Applying selected calibration to calibration set for diagnostics:\n")
  cat("(this is optimistic since it's the training data for calibration)\n")
  calib_results_calibrated <- cal_apply(calib_results, final_cal)
}

cal_plot_windowed(
  slice_sample(calib_results_calibrated, n = 100e3), 
  truth = truth, 
  estimate = .pred_1,
  event_level = "first",
  include_rug = FALSE
)

## Step 4: Calculate custom thresholds from combined training and calibration partitions ####

cat("\nCalculating custom thresholds from combined training and calibration partitions...\n")

# First, apply the selected calibration method to training data if needed
if (best_method == "none") {
  train_results_calibrated <- train_results
} else {
  train_results_calibrated <- cal_apply(train_results, final_cal)
}

# Combine training and calibration predictions for threshold calculation
combined_results <- bind_rows(
  train_results_calibrated |> mutate(partition = "train"),
  calib_results_calibrated |> mutate(partition = "calib")
)

cat(sprintf("Combined dataset for threshold calculation: %d observations (%d train + %d calib)\n",
            nrow(combined_results), nrow(train_results_calibrated), nrow(calib_results_calibrated)))
cat(sprintf("Combined prevalence: %.2f%%\n", mean(combined_results$truth == "1") * 100))

# Find thresholds for 95% specificity and 95% sensitivity using calibrated probabilities from combined data
custom_thresholds <- find_custom_thresholds(combined_results, "truth", ".pred_1")

cat(sprintf("Number of thresholds evaluated: %d\n", custom_thresholds$n_thresholds))
cat(sprintf("Maximum achievable specificity: %.3f\n", custom_thresholds$max_spec))
cat(sprintf("Maximum achievable sensitivity: %.3f\n", custom_thresholds$max_sens))
cat(sprintf("Threshold for 95%% specificity: %s\n", 
            ifelse(is.na(custom_thresholds$threshold_95spec), "Not achievable", 
                   sprintf("%.4f", custom_thresholds$threshold_95spec))))
cat(sprintf("Threshold for 95%% sensitivity: %s\n", 
            ifelse(is.na(custom_thresholds$threshold_95sens), "Not achievable", 
                   sprintf("%.4f", custom_thresholds$threshold_95sens))))

## Step 5: Apply method to test set predictions ####

if (best_method == "none") {
  cat("\nUsing raw predictions on test set (no calibration)...\n")
  test_results_calibrated <- test_results_raw |> 
    mutate(
      estimate = class  # Use original model class predictions (prevalence-based threshold)
    )
} else {
  cat("\nCalibrating predictions on test set...\n")
  test_results_calibrated <- cal_apply(test_results_raw, final_cal) |> 
    mutate(
      estimate = if_else(.pred_1 >= 0.5, "1", "0"),  # Use 0.5 threshold for calibrated probabilities
      estimate = factor(estimate, levels = c("1", "0"))
    )
}

# Calculate positive predictions
uncalibrated_positives <- test_results_raw |> 
  filter(class == "1") |>
  nrow()
processed_positives <- test_results_calibrated |> 
  filter(estimate == "1") |>
  nrow()

cat("Test predictions processed.\n")
cat("Raw test prediction range:", round(min(test_preds$predicted[,1]), 4), 
    "to", round(max(test_preds$predicted[,1]), 4), "\n")
cat("Processed test prediction range:", 
    round(min(test_results_calibrated$.pred_1), 4), "to", 
    round(max(test_results_calibrated$.pred_1), 4), "\n")
if (best_method == "none") {
  cat("Raw vs processed positive predictions (model's prevalence-based threshold):\n", 
      uncalibrated_positives, "->", processed_positives, "\n")
} else {
  cat("Raw vs calibrated positive predictions (0.5 threshold on calibrated probabilities):\n", 
      uncalibrated_positives, "->", processed_positives, "\n")
}

## Step 6: Evaluate performance on test set ####

cat("\nEvaluating model performance on test set...\n")

# Calculate performance metrics for raw predictions
# Separate metrics by what they need: probabilities vs class predictions
prob_metrics <- metric_set(roc_auc, pr_auc)
class_metrics <- metric_set(accuracy, sens, spec, ppv, npv)

# Raw predictions: probability-based metrics
raw_prob_performance <- test_results_raw |>
  prob_metrics(truth = truth, .pred_1, event_level = "first") |>
  mutate(threshold = NA_real_, threshold_type = "prob_metrics")

# Raw predictions: class-based metrics (prevalence-based threshold)
raw_class_performance <- test_results_raw |>
  class_metrics(truth = truth, estimate = class, event_level = "first") |>
  mutate(threshold = mean(test_results_raw$truth == "1"), threshold_type = "prevalence")

# Calibrated predictions: probability-based metrics
calibrated_prob_performance <- test_results_calibrated |>
  prob_metrics(truth = truth, .pred_1, event_level = "first") |>
  mutate(threshold = NA_real_, threshold_type = "prob_metrics")

# Calibrated predictions: class-based metrics (0.5 threshold or prevalence for raw)
if (best_method == "none") {
  calibrated_class_performance <- test_results_calibrated |>
    class_metrics(truth = truth, estimate = estimate, event_level = "first") |>
    mutate(threshold = mean(test_results_calibrated$truth == "1"), threshold_type = "prevalence")
} else {
  calibrated_class_performance <- test_results_calibrated |>
    class_metrics(truth = truth, estimate = estimate, event_level = "first") |>
    mutate(threshold = 0.5, threshold_type = "fixed_0.5")
}

# Custom threshold evaluations using calibrated probabilities
custom_threshold_performance <- tibble()

if (!is.na(custom_thresholds$threshold_95spec)) {
  cat(sprintf("Evaluating at 95%% specificity threshold (%.4f)...\n", custom_thresholds$threshold_95spec))
  
  test_results_95spec <- test_results_calibrated |>
    mutate(
      estimate_95spec = if_else(.pred_1 >= custom_thresholds$threshold_95spec, "1", "0"),
      estimate_95spec = factor(estimate_95spec, levels = c("1", "0"))
    )
  
  custom_95spec_performance <- test_results_95spec |>
    class_metrics(truth = truth, estimate = estimate_95spec, event_level = "first") |>
    mutate(threshold = custom_thresholds$threshold_95spec, threshold_type = "95_specificity")
  
  custom_threshold_performance <- bind_rows(custom_threshold_performance, custom_95spec_performance)
}

if (!is.na(custom_thresholds$threshold_95sens)) {
  cat(sprintf("Evaluating at 95%% sensitivity threshold (%.4f)...\n", custom_thresholds$threshold_95sens))
  
  test_results_95sens <- test_results_calibrated |>
    mutate(
      estimate_95sens = if_else(.pred_1 >= custom_thresholds$threshold_95sens, "1", "0"),
      estimate_95sens = factor(estimate_95sens, levels = c("1", "0"))
    )
  
  custom_95sens_performance <- test_results_95sens |>
    class_metrics(truth = truth, estimate = estimate_95sens, event_level = "first") |>
    mutate(threshold = custom_thresholds$threshold_95sens, threshold_type = "95_sensitivity")
  
  custom_threshold_performance <- bind_rows(custom_threshold_performance, custom_95sens_performance)
}

# Combine all raw performance
raw_performance <- bind_rows(
  prob = raw_prob_performance, 
  class = raw_class_performance,
  .id = "metric_type")

# Combine all calibrated performance
calibrated_performance <- bind_rows(
  prob = calibrated_prob_performance, 
  class = calibrated_class_performance,
  class = custom_threshold_performance,
  .id = "metric_type")

cat("\n=== RAW PREDICTIONS PERFORMANCE ===\n")
print(raw_performance, n = Inf)

if (best_method == "none") {
  cat("\n=== PROCESSED PREDICTIONS PERFORMANCE (NO CALIBRATION) ===\n") 
} else {
  cat("\n=== CALIBRATED PREDICTIONS PERFORMANCE ===\n") 
}
print(calibrated_performance, n = Inf)

## Step 7: Calibration diagnostic plots ####

cat("\nGenerating calibration diagnostic plots...\n")

cal_plot <- cal_plot_windowed(
  slice_sample(test_results_calibrated, n = 100e3), 
  truth = truth, 
  estimate = .pred_1,
  event_level = "first", 
  include_rug = FALSE
)

cat("\nRaw predictions calibration plot:\n")
print(raw_cal_plot_partitiontest)

if (best_method == "none") {
  cat("\nProcessed predictions calibration plot (no calibration applied):\n")
  plot_title <- "Calibration Plots: Raw vs Processed Predictions"
  plot_subtitle <- "Method: No calibration (raw predictions)"
} else {
  cat("\nCalibrated predictions calibration plot:\n")
  plot_title <- "Calibration Plots: Raw vs Calibrated Predictions"  
  plot_subtitle <- sprintf("Calibration method: %s", str_to_title(best_method))
}
print(cal_plot)

patchwork::wrap_plots(raw_cal_plot_partitiontest, cal_plot, ncol = 2) +
  patchwork::plot_annotation(
    title = plot_title,
    subtitle = plot_subtitle
  )

## Step 8: Save outputs ####

cat("\nSaving results...\n")

# Save performance metrics in tidy format
performance_comparison <- bind_rows(
  raw_performance |> mutate(prediction_type = "raw"),
  calibrated_performance |> mutate(prediction_type = if_else(best_method == "none", "processed_no_calibration", "calibrated"))
) |>
  select(prediction_type, metric_type, .metric, .estimator, .estimate, threshold, threshold_type)

write_csv(performance_comparison, "output/test_performance_comparison.csv")
cat("Performance metrics saved to: output/test_performance_comparison.csv\n")

# Save calibration object if applicable
if (best_method != "none") {
  saveRDS(final_cal, "output/calibration60split.rds")
  cat("Calibration object saved to: output/calibration60split.rds\n")
} else {
  cat("No calibration object to save (raw predictions were optimal).\n")
}

# sessionInfo ####

sessioninfo::session_info()