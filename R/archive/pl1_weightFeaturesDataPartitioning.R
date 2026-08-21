# RF-based feature importance assessment at regional scale for data partitioning ####

library(readr)
library(dplyr)
library(terra)
library(tidyr)
library(ggplot2)

## Create training dataset ####
rf_global_current <- rast("output/rf_global_pred_regional_current.tif") 
names(rf_global_current) <- "rf_global"
preds_nor_250m_cur <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")
preds_train <- c(rf_global_current, preds_nor_250m_cur)

presence <- read_csv("output/presence_coords_regional.csv")
absence <- read_csv("output/absence_coords_regional.csv")

# Extract predictor values for presence points
presence_df <- terra::extract(preds_train, presence[c("x", "y")], xy=TRUE, ID=FALSE) |> 
  mutate(response = 1, .before = 1) |>
  drop_na() |> 
  tibble()

# Extract predictor values for absence points
absence_df <- terra::extract(preds_train, absence[c("x", "y")], xy=TRUE, ID=FALSE) |> 
  mutate(response = 0, .before = 1) |>
  drop_na() |> 
  tibble()

# Combine into training dataset (keep x,y coordinates for later use)
training_data_with_coords <- bind_rows(presence_df, absence_df) |>
  mutate(response = factor(response))

training_data <- training_data_with_coords |> 
  select(-x, -y) # Remove x and y coordinates for modeling

# Get final counts
n_final_presence <- sum(training_data$response == "1")
n_final_absence <- sum(training_data$response == "0")
n_final_prevalence <- n_final_presence / (n_final_presence + n_final_absence)

cat("Final training data - Presences:", n_final_presence, "Absence:", n_final_absence, "\n")
cat("Prevalence in training data:", round(n_final_prevalence * 100, 3), "%\n")

## Feature weights for feature space distance ####

training_data_MIAmaxent <- training_data |> 
  mutate(response = as.numeric(as.character(response))) |> # Convert factor to numeric
  as.data.frame()  # Convert tibble to data.frame for MIAmaxent
MIAmaxent::plotFOP(training_data_MIAmaxent, EV="rf_global")
MIAmaxent::plotFOP(training_data_MIAmaxent, EV="artype_60")
MIAmaxent::plotFOP(training_data_MIAmaxent, EV="artype_81")
MIAmaxent::plotFOP(training_data_MIAmaxent, EV="slope")

### With Balanced Random Forest (in ranger) ####

library(ranger)  # Faster for repeated runs

get_brf_weights <- function(outcome, predictors, n_reps = 1) {
  importance_matrix <- matrix(0, n_reps, ncol(predictors))

  full_train_data <- cbind(outcome, predictors) |> 
    mutate(outcome = as.factor(outcome))
  
  for(i in 1:n_reps) {
    small_train_data <- full_train_data |> 
      slice_sample(n = min(10e3, nrow(full_train_data)))
    
    # Calculate the fraction needed for balanced sampling
    n_pos <- sum(small_train_data$outcome == 1)
    n_neg <- sum(small_train_data$outcome == 0)
    
    # For balanced: sample ALL positives, and equal number of negatives
    frac_pos <- 1.0  # Use all positive cases
    frac_neg <- n_pos / n_neg  # Downsample negatives to match
    
    rf <- ranger(outcome ~ .,
                 data = small_train_data,
                 num.trees = 1000,
                 importance = "permutation",
                 scale.permutation.importance = TRUE,
                 sample.fraction = c(frac_pos, frac_neg),  # Order matters! Verify with keep.inbag = TRUE
                 replace = TRUE,
                 class.weights = NULL, # Don't use with sample.fraction
                 verbose = TRUE)  

    importance_matrix[i, ] <- rf$variable.importance
  }

  # Calculate median and sd for each feature
  weights_median <- apply(importance_matrix, 2, median)
  weights_sd <- apply(importance_matrix, 2, sd)
  names(weights_median) <- names(rf$variable.importance)
  names(weights_sd) <- names(rf$variable.importance)
  
  # Return all repetitions along with summary stats
  return(list(
    reps = importance_matrix,
    median = tibble::enframe(weights_median, name = "name", value = "median"),
    sd = tibble::enframe(weights_sd, name = "name", value = "sd")
  ))
}

weights_brf_results <- get_brf_weights(
  training_data$response,
  training_data[, -1],
  n_reps = 20)

weights_brf_results$median |> 
  arrange(desc(median)) |>
  print(n = 29)

### With Random Forest Q-classification (in rfsrc) ####

library(randomForestSRC)

get_rfq_weights <- function(outcome, predictors, n_reps = 1) {
  importance_matrix <- matrix(0, n_reps, ncol(predictors))
  full_train_data <- cbind(outcome, predictors) |> 
    mutate(outcome = as.factor(outcome))

  for(i in 1:n_reps) {
    small_train_data <- full_train_data |> 
      slice_sample(n = min(10e3, nrow(full_train_data))) |>
      as.data.frame()  # Convert tibble to data.frame for imbalanced()
    
    # Train final model with selected hyperparameters
    rf <- imbalanced(
      formula = outcome ~ .,
      data = small_train_data,
      # From docs: For computational speed, the default VIMP method has changed
      # from "permute" (Breiman-Cutler permutation) to "anti" (importance = "anti"
      # or importance = TRUE). While faster, this may be less accurate in settings
      # such as highly imbalanced classification. To revert to permutation VIMP,
      # use importance = "permute".
      importance = "permute",
      do.trace = 30
    )
    
    importance_matrix[i, ] <- rf$importance[,"all"]
  }

  # Calculate median and sd for each feature
  weights_median <- apply(importance_matrix, 2, median)
  weights_sd <- apply(importance_matrix, 2, sd)
  names(weights_median) <- rownames(rf$importance)
  names(weights_sd) <- rownames(rf$importance)
  
  # Return all repetitions along with summary stats
  return(list(
    reps = importance_matrix,
    median = tibble::enframe(weights_median, name = "name", value = "median"),
    sd = tibble::enframe(weights_sd, name = "name", value = "sd")
  ))
}

weights_rfq_results <- get_rfq_weights(
  training_data$response,
  training_data[, -1],
  n_reps = 20)

weights_rfq_results$median |> 
  arrange(desc(median)) |>
  print(n = 29)

## Plotting

# Combine median and sd data from both methods
weights_brf_summary <- left_join(
  weights_brf_results$median,
  weights_brf_results$sd,
  by = "name", suffix = c("_median", "_sd")
) |> 
  mutate(method = "Balanced Random Forest (ranger)")

weights_rfq_summary <- left_join(
  weights_rfq_results$median,
  weights_rfq_results$sd,
  by = "name", suffix = c("_median", "_sd")
) |> 
  mutate(method = "Random Forest Q-classification (rfsrc)")

# Combine both methods
weights_combined <- bind_rows(weights_brf_summary, weights_rfq_summary)

# Normalize within each method to make them comparable
weights_df <- weights_combined |> 
  group_by(method) |> 
  mutate(
    median_normalized = median / max(median),
    sd_normalized = sd / max(median)  # Scale SD by same factor as median
  ) |>
  ungroup() |>
  mutate(feature = forcats::fct_reorder(name, median_normalized, .fun = last))

# Create plot with error bars
ggplot(weights_df, aes(x = feature, y = median_normalized, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  geom_errorbar(aes(ymin = pmax(0, median_normalized - sd_normalized), 
                    ymax = median_normalized + sd_normalized),
                position = position_dodge(width = 0.9), width = 0.2) +
  coord_flip() +
  labs(title = "Feature Importance from Different Random Forest Methods",
       x = "Feature",
       y = "Normalized Importance",
       fill = "Method") +
  theme_minimal() +
  theme(legend.position = "bottom")

## Saving to file ####

# Save the combined data with median and SD
weights_df |> 
  select(feature = name, method, median, sd, median_normalized, sd_normalized) |>
  write_csv("output/pl1/weights_feature_data_partitioning.csv", append = FALSE)

# sessionInfo ####

sessioninfo::session_info()
