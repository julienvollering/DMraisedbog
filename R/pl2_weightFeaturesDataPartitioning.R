# Feature importance assessment at regional scale for data partitioning ####
# Compares three methods: RF Q-classification, Balanced RF, and Penalized Regression

library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(randomForestSRC)
library(randomForest)
library(glmnet)

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
  as.data.frame() # Convert tibble to data.frame for MIAmaxent
MIAmaxent::plotFOP(training_data_MIAmaxent, EV = "rf_global")
MIAmaxent::plotFOP(training_data_MIAmaxent, EV = "slope")
MIAmaxent::plotFOP(training_data_MIAmaxent, EV = "elevation")

### Functions ####
get_rfq_weights <- function(
  outcome,
  predictors,
  n_reps = 5,
  absence_ratio = 10
) {
  full_train_data <- cbind(outcome, predictors) |>
    mutate(outcome = as.factor(outcome))

  presences <- full_train_data |> filter(outcome == "1")
  absences <- full_train_data |> filter(outcome == "0")
  n_pres <- nrow(presences)
  n_abs_sample <- min(n_pres * absence_ratio, nrow(absences))

  map_dfr(
    1:n_reps,
    ~ {
      set.seed(1000 + .x)

      small_train_data <- bind_rows(
        presences, # All presences, every time
        absences |> slice_sample(n = n_abs_sample)
      ) |>
        as.data.frame()

      # Default for perf.type:
      # "Permutation-based VIMP is used by default in this setting,
      # in contrast to anti-VIMP which is the default for other families.
      # Empirical results suggest that permutation VIMP is more reliable
      # in highly imbalanced settings."
      rf <- imbalanced(
        formula = outcome ~ .,
        data = small_train_data,
        ntree = 500,
        importance = TRUE,
        do.trace = FALSE
      )

      importance_values <- tibble(
        rep = .x,
        feature = rownames(rf$importance),
        importance = rf$importance[, "all"]
      )

      rm(rf, small_train_data)
      gc()

      importance_values
    }
  ) |>
    group_by(feature) |>
    summarise(
      median = median(importance),
      sd = sd(importance),
      .groups = "drop"
    )
}

get_glmnet_weights <- function(
  outcome,
  predictors,
  n_reps = 5,
  absence_ratio = 10,
  alpha = 0 # Ridge by default to distribute weights among correlated features
) {
  # NOTE: We manually scale predictors so coefficients are on standardized scale
  # and comparable across features. glmnet's internal standardization returns
  # coefficients on original scale, making them incomparable (1-unit change in
  # bio03 vs. 1-unit change in gdd10 are not meaningful comparisons).

  full_train_data <- cbind(outcome, predictors) |>
    mutate(outcome = as.numeric(as.character(outcome)))

  presences <- full_train_data |> filter(outcome == 1)
  absences <- full_train_data |> filter(outcome == 0)
  n_pres <- nrow(presences)
  n_abs_sample <- min(n_pres * absence_ratio, nrow(absences))

  map_dfr(
    1:n_reps,
    ~ {
      set.seed(1000 + .x)

      small_train_data <- bind_rows(
        presences, # All presences, every time
        absences |> slice_sample(n = n_abs_sample)
      )

      # Prepare data for glmnet - EXPLICITLY SCALE for comparable coefficients
      X <- as.matrix(small_train_data[, -1])
      X_scaled <- scale(X)
      y <- small_train_data$outcome

      # Class weights to handle imbalance
      weights <- ifelse(y == 1, 1 / mean(y), 1 / mean(1 - y))

      # Fit with cross-validation on scaled data
      fit <- cv.glmnet(
        X_scaled,
        y,
        family = "binomial",
        weights = weights,
        alpha = alpha,
        standardize = FALSE # Already scaled manually
      )
      coefs <- coef(fit, s = "lambda.1se")[-1, 1] # Drop intercept, extract column

      importance_values <- tibble(
        rep = .x,
        feature = names(coefs),
        importance = abs(coefs)
      )

      rm(fit, X, X_scaled, y, small_train_data)
      gc()

      importance_values
    }
  ) |>
    group_by(feature) |>
    summarise(
      median = median(importance),
      sd = sd(importance),
      .groups = "drop"
    )
}

get_brf_weights <- function(
  outcome,
  predictors,
  n_reps = 5,
  absence_ratio = 10
) {
  full_train_data <- cbind(outcome, predictors) |>
    mutate(outcome = as.factor(outcome))

  presences <- full_train_data |> filter(outcome == "1")
  absences <- full_train_data |> filter(outcome == "0")
  n_pres <- nrow(presences)
  n_abs_sample <- min(n_pres * absence_ratio, nrow(absences))

  map_dfr(
    1:n_reps,
    ~ {
      set.seed(1000 + .x)

      small_train_data <- bind_rows(
        presences, # All presences, every time
        absences |> slice_sample(n = n_abs_sample)
      ) |>
        as.data.frame()

      # Balanced Random Forest: use sampsize to ensure balanced bootstrap samples
      n_per_class <- min(
        sum(small_train_data$outcome == "1"),
        sum(small_train_data$outcome == "0")
      )

      rf <- randomForest(
        formula = outcome ~ .,
        data = small_train_data,
        ntree = 500,
        sampsize = c("0" = n_per_class, "1" = n_per_class),
        importance = TRUE,
        replace = TRUE
      )

      # Use MeanDecreaseAccuracy as importance measure
      importance_values <- tibble(
        rep = .x,
        feature = rownames(rf$importance),
        importance = rf$importance[, "MeanDecreaseAccuracy"]
      )

      rm(rf, small_train_data)
      gc()

      importance_values
    }
  ) |>
    group_by(feature) |>
    summarise(
      median = median(importance),
      sd = sd(importance),
      .groups = "drop"
    )
}

### With Random Forest Q-classification (in rfsrc) ####

tictoc::tic()
weights_rfq_results <- get_rfq_weights(
  mf_current$response,
  mf_current[, -1],
  n_reps = 10,
  absence_ratio = 10
)
tictoc::toc()

weights_rfq_results |>
  arrange(desc(median)) |>
  print(n = Inf)

### With Penalized Regression (glmnet) ####

tictoc::tic()
weights_glmnet_results <- get_glmnet_weights(
  mf_current$response,
  mf_current[, -1],
  n_reps = 10,
  absence_ratio = 10,
  alpha = 0.5 # Elastic Net to balance between Ridge and Lasso
)
tictoc::toc()

weights_glmnet_results |>
  arrange(desc(median)) |>
  print(n = Inf)

### With Balanced Random Forest (randomForest) ####

tictoc::tic()
weights_brf_results <- get_brf_weights(
  mf_current$response,
  mf_current[, -1],
  n_reps = 10,
  absence_ratio = 10
)
tictoc::toc()

weights_brf_results |>
  arrange(desc(median)) |>
  print(n = Inf)

## Combine results ####

weights_df <- bind_rows(
  weights_rfq_results |>
    mutate(method = "Random Forest Q-classification"),
  weights_brf_results |>
    mutate(method = "Balanced Random Forest"),
  weights_glmnet_results |>
    mutate(method = "Penalized Regression (glmnet)")
)

## Plotting ####

# Faceted plot with independent y-axes (unscaled)
weights_df |>
  group_by(method) |>
  mutate(feature = forcats::fct_reorder(feature, median)) |>
  ungroup() |>
  ggplot(aes(x = feature, y = median)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(
    aes(ymin = pmax(0, median - sd), ymax = median + sd),
    width = 0.2
  ) +
  coord_flip() +
  facet_wrap(~method, scales = "free_x", ncol = 2) +
  labs(
    title = "Feature Importance Comparison (Original Scale)",
    x = "Feature",
    y = "Importance"
  ) +
  theme_minimal()

## Diagnostic: plotFOP for top and bottom features ####

cat("\n=== Diagnostic: Checking feature-response relationships ===\n")

# Get top 4 features for each method
diagnostic_features <- weights_df |>
  group_by(method) |>
  arrange(desc(median)) |>
  mutate(rank = row_number()) |>
  filter(rank <= 4) |>
  ungroup() |>
  select(method, feature, rank, median) |>
  arrange(method, rank, desc(median))

print(diagnostic_features)

# Plot FOP for top 4 features from each method

par(mfrow = c(2, 2))
for (method in unique(diagnostic_features$method)) {
  features_to_plot <- diagnostic_features |>
    filter(method == !!method) |>
    arrange(rank) |>
    pull(feature)

  cat("\n--- Method:", method, "---\n")
  pairs(training_data_MIAmaxent[, features_to_plot])
  for (feature in features_to_plot) {
    cat("Plotting feature:", feature, "\n")
    MIAmaxent::plotFOP(training_data_MIAmaxent, EV = feature)
  }
}
par(mfrow = c(1, 1)) # Reset plotting layout

## Saving to file ####

# Save both methods with median, SD, and scaled values
weights_df |>
  select(feature, method, median, sd) |>
  arrange(method, desc(median)) |>
  write_csv("output/pl2/weights_feature_data_partitioning.csv", append = FALSE)

# sessionInfo ####

sessioninfo::session_info()
