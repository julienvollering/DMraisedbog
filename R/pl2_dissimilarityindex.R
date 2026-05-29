# Area of Applicability Assessment for Local-Scale Future Predictions
# Based on environmental distances between CV folds during training
# Using CAST::aoa() with variable importance weighting

# Load required packages
library(tidyverse)
library(terra)
library(CAST)
library(randomForestSRC)

source("R/functions.R")


## Read modeling frame ####
mf <- read_csv("output/pl2/modeling_frame_regional.csv")

raster_current <- rast("output/pl2/scenario_current.tif")

mf_current_with_coords <- mf |>
  filter(scenario == "current") |>
  select(-scenario)
mf_current <- mf_current_with_coords |>
  select(-x, -y)

toJoin <- mf |>
  filter(scenario == "future") |>
  select(-scenario)
mf_future_presence <- mf |>
  filter(response == 1) |>
  select(x, y) |>
  left_join(toJoin, by = c("x", "y")) |>
  select(-x, -y)

n_presence <- sum(mf_current$response == "1")
n_absence <- sum(mf_current$response == "0")
prevalence <- n_presence / (n_presence + n_absence)

cat("Training data - Presences:", n_presence, "Absences:", n_absence, "\n")
cat("Prevalence:", round(prevalence * 100, 2), "%\n")

## Load feature weights and identify most important feature ####
weights_features <- read_csv("output/pl2/weights_feature_data_partitioning.csv")
weighting_method <- "Balanced Random Forest"

weights <- weights_features |>
  filter(method == weighting_method) |>
  select(feature, median) |>
  arrange(desc(median)) |>
  tibble::deframe()

most_important_feature <- names(weights)[1]

cat("\nMost important feature:", most_important_feature, "\n")

## Compute Dissimilarity Index ####
?trainDI

### On (unspaced) 5-cut partitions based on most important feature ####

partitions_5cuts <- create_cut_based_partitions(
  data = mf_current,
  feature = most_important_feature,
  target_col = "response",
  n_cuts = 5,
  n_partitions = 5,
  seed = 42
)

if (!file.exists("output/pl2/dissimilarity_index_train_5cutsGDD10.rds")) {
  di_results <- trainDI(
    train = mf_current,
    variables = names(weights),
    weight = as.data.frame(t(weights)),
    CVtest = partitions$partitions,
    useWeight = TRUE,
    verbose = TRUE
  )
  saveRDS(di_results, "output/pl2/dissimilarity_index_train_5cutsGDD10.rds")
} else {
  di_results <- readRDS("output/pl2/dissimilarity_index_train_5cutsGDD10.rds")
}

aoa_results <- aoa(
  newdata = mf_future_presence,
  trainDI = di_results,
  variables = names(weights),
  weight = as.data.frame(t(weights)),
  useWeight = TRUE,
  verbose = TRUE
)

# Histograms of Dissimilarity Index faceted by scenario
di_train_df <- data.frame(
  DI = di_results$trainDI,
  scenario = "current"
)
di_future_df <- data.frame(
  DI = aoa_results$DI,
  scenario = "future"
)
di_combined_df <- rbind(di_train_df, di_future_df)
ggplot(di_combined_df, aes(x = DI, fill = scenario)) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.6,
    bins = 50
  ) +
  labs(
    title = "Dissimilarity Index Distribution by Scenario",
    subtitle = "Unspaced 5-cut partitions using most important feature",
    x = "Dissimilarity Index",
    y = "Density"
  ) +
  geom_vline(
    xintercept = di_results$threshold,
    linetype = "dashed",
    color = "black"
  ) +
  facet_grid(scenario ~ ., scales = "free_y") +
  theme_minimal() +
  scale_fill_manual(values = c("current" = "purple", "future" = "orange"))

### On spaced 5-cut partitions based on most important feature ####

partitions_5cutsspaced <- create_spacer_based_partitions(
  data = mf_current,
  feature = most_important_feature,
  target_col = "response",
  n_cuts = 1000,
  spacer_proportion = 0.5,
  n_partitions = 5,
  seed = 42
)

mf_current_spaced <- mf_current |>
  mutate(partition = partitions_5cutsspaced$partitions) |>
  filter(!is.na(partition))

if (!file.exists("output/pl2/dissimilarity_index_train_5cutsspacedGDD10.rds")) {
  di_results_spaced <- trainDI(
    train = mf_current_spaced,
    variables = names(weights),
    weight = as.data.frame(t(weights)),
    CVtest = mf_current_spaced$partition,
    useWeight = TRUE,
    verbose = TRUE
  )
  saveRDS(
    di_results_spaced,
    "output/pl2/dissimilarity_index_train_5cutsspacedGDD10.rds"
  )
} else {
  di_results_spaced <- readRDS(
    "output/pl2/dissimilarity_index_train_5cutsspacedGDD10.rds"
  )
}

aoa_results_spaced <- aoa(
  newdata = mf_future_presence,
  trainDI = di_results_spaced,
  variables = names(weights),
  weight = as.data.frame(t(weights)),
  useWeight = TRUE,
  verbose = TRUE
)

# Histograms of Dissimilarity Index faceted by scenario
di_train_spaced_df <- data.frame(
  DI = di_results_spaced$trainDI,
  scenario = "CV"
)
di_future_spaced_df <- data.frame(
  DI = aoa_results_spaced$DI,
  scenario = "current-future"
)
di_combined_spaced_df <- rbind(di_train_spaced_df, di_future_spaced_df)
p_spaced <- ggplot(di_combined_spaced_df, aes(x = DI, fill = scenario)) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.6,
    bins = 50
  ) +
  labs(
    title = "Five CV folds stratified on most important feature, with spacers",
    x = "Dissimilarity Index",
    y = "Density"
  ) +
  geom_vline(
    xintercept = di_results_spaced$threshold,
    linetype = "dashed",
    color = "black"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("CV" = "purple", "current-future" = "orange"))
p_spaced
ggsave(
  "meetings/20260325-Landskapskaffi/img/di_spaced.png",
  plot = p_spaced,
  width = 6,
  height = 4,
  dpi = 150
)

### On random partitions ####

set.seed(42)
partitions_random <- sample.int(5, nrow(mf_current), replace = TRUE)

if (!file.exists("output/pl2/dissimilarity_index_train_5random.rds")) {
  di_results_random <- trainDI(
    train = mf_current,
    variables = names(weights),
    weight = as.data.frame(t(weights)),
    CVtest = partitions_random,
    useWeight = TRUE,
    verbose = TRUE
  )
  saveRDS(di_results_random, "output/pl2/dissimilarity_index_train_5random.rds")
} else {
  di_results_random <- readRDS(
    "output/pl2/dissimilarity_index_train_5random.rds"
  )
}

aoa_results_random <- aoa(
  newdata = mf_future_presence,
  trainDI = di_results_random,
  variables = names(weights),
  weight = as.data.frame(t(weights)),
  useWeight = TRUE,
  verbose = TRUE
)

# Histograms of Dissimilarity Index faceted by scenario
di_train_df <- data.frame(
  DI = di_results_random$trainDI,
  scenario = "CV"
)
di_future_df <- data.frame(
  DI = aoa_results_random$DI,
  scenario = "current-future"
)
di_combined_df <- rbind(di_train_df, di_future_df)
p_random <- ggplot(di_combined_df, aes(x = DI, fill = scenario)) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.6,
    bins = 50
  ) +
  labs(
    title = "Five random CV folds",
    x = "Dissimilarity Index",
    y = "Density"
  ) +
  geom_vline(
    xintercept = di_results_random$threshold,
    linetype = "dashed",
    color = "black"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("CV" = "purple", "current-future" = "orange"))
p_random
ggsave(
  "meetings/20260325-Landskapskaffi/img/di_random.png",
  plot = p_random,
  width = 6,
  height = 4,
  dpi = 150
)


# sessionInfo

sessioninfo::session_info()
