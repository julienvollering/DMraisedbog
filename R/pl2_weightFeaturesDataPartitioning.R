# RF-based feature importance assessment at regional scale for data partitioning ####

library(readr)
library(dplyr)
library(ggplot2)
library(ranger)
library(randomForestSRC)

## Read modeling frame ####
mf <- read_csv("output/pl2/modeling_frame_regional.csv")

mf_current <- mf |> 
  filter(scenario == "current") |>
  select(-scenario, -x, -y) # Remove x and y coordinates for modeling

n_presence <- sum(mf_current$response == "1")
n_absence <- sum(mf_current$response == "0")
prevalence <- n_presence / (n_presence + n_absence)

cat("Final training data - Presences:", n_presence, "Absence:", n_absence, "\n")
cat("Prevalence in training data:", round(prevalence * 100, 2), "%\n")

## Feature weights for feature space distance ####

### Explore ####
training_data_MIAmaxent <- mf_current |> 
  mutate(response = as.numeric(as.character(response))) |> # Convert factor to numeric
  as.data.frame()  # Convert tibble to data.frame for MIAmaxent
MIAmaxent::plotFOP(training_data_MIAmaxent, EV="rf_global")
MIAmaxent::plotFOP(training_data_MIAmaxent, EV="slope")
MIAmaxent::plotFOP(training_data_MIAmaxent, EV="elevation")

### Functions ####
get_brf_weights <- function(outcome, predictors, n_reps = 1) {
  importance_matrix <- matrix(0, n_reps, ncol(predictors))

  full_train_data <- cbind(outcome, predictors) |> 
    mutate(outcome = as.factor(outcome))
  
  for(i in 1:n_reps) {
    # Calculate the fraction needed for balanced sampling
    n_pos <- sum(full_train_data$outcome == 1)
    n_neg <- sum(full_train_data$outcome == 0)
    
    # For balanced: sample ALL positives, and equal number of negatives
    frac_pos <- 1.0  # Use all positive cases
    frac_neg <- n_pos / n_neg  # Downsample negatives to match
    
    rf <- ranger(outcome ~ .,
                 data = full_train_data,
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

get_rfq_weights <- function(outcome, predictors, n_train, n_reps = 1) {
  importance_matrix <- matrix(0, n_reps, ncol(predictors))
  full_train_data <- cbind(outcome, predictors) |> 
    mutate(outcome = as.factor(outcome))

  for(i in 1:n_reps) {
    # Randomly subsample training data for runtime efficiency
    small_train_data <- full_train_data |> 
      slice_sample(n = min(n_train, nrow(full_train_data))) |>
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
      do.trace = 60
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

### With Balanced Random Forest (in ranger) ####

weights_brf_results <- get_brf_weights(
  mf_current$response,
  mf_current[, -1],
  n_reps = 10)

weights_brf_results$median |> 
  arrange(desc(median)) |>
  print(n = 22)

### With Random Forest Q-classification (in rfsrc) ####

n_train <- 10e3
cat("Training RFQ models on", n_train, "observations...\n")
npresences95CI <- qbinom(
  p = c(0.025, 0.975), 
  size = n_train, 
  prob = prevalence)
cat("Data likely contain", npresences95CI[1], "to", npresences95CI[2], 
    "presences (95% CI, prevalence =", round(prevalence, 4), ")\n")

weights_rfq_results <- get_rfq_weights(
  mf_current$response,
  mf_current[, -1],
  n_train = n_train,
  n_reps = 10)

weights_rfq_results$median |> 
  arrange(desc(median)) |>
  print(n = 22)

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
  mutate(feature = forcats::fct_reorder(name, median_normalized, .fun = first))

# Create plot with error bars
ggplot(weights_df, aes(x = feature, y = median_normalized, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  geom_errorbar(aes(ymin = median_normalized - sd_normalized, 
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
  write_csv("output/pl2/weights_feature_data_partitioning.csv", append = FALSE)

# sessionInfo ####

sessioninfo::session_info()
