# RF-based feature importance assessment at regional scale for data partitioning ####

library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
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
get_rfq_weights <- function(outcome, predictors, n_train, n_reps = 1) {
  full_train_data <- cbind(outcome, predictors) |>
    mutate(outcome = as.factor(outcome))

  # Use purrr to run repetitions and collect results
  map_dfr(1:n_reps, ~ {
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
      do.trace = 30
    )

    # Return tibble with rep, feature, importance
    tibble(
      rep = .x,
      feature = rownames(rf$importance),
      importance = rf$importance[, "all"]
    )
  }) |>
    # Summarize across repetitions
    group_by(feature) |>
    summarise(
      median = median(importance),
      sd = sd(importance),
      .groups = "drop"
    )
}

### With Random Forest Q-classification (in rfsrc) ####

n_train_cap <- 20e3
n_train <- min(n_train_cap, nrow(mf_current))
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
  n_reps = 5)

weights_rfq_results |>
  arrange(desc(median)) |>
  print(n = 22)

## Plotting

weights_df <- weights_rfq_results |>
  mutate(
    method = "Random Forest Q-classification (rfsrc)"
  ) |>
  mutate(feature = forcats::fct_reorder(feature, median))

# Create plot with error bars
ggplot(weights_df, aes(x = feature, y = median)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = median - sd,
                    ymax = median + sd),
                width = 0.2) +
  coord_flip() +
  labs(title = "Feature Importance from Random Forest Q-classification",
       x = "Feature",
       y = "Importance") +
  theme_minimal()

## Saving to file ####

# Save the RFQ data with median and SD
weights_df |>
  select(feature, method, median, sd) |>
  write_csv("output/pl2/weights_feature_data_partitioning.csv", append = FALSE)

# sessionInfo ####

sessioninfo::session_info()
