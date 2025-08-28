# HDBSCAN Results Exploration ####

library(tidyverse)
library(terra)
library(cluster)
library(dbscan)
library(sf)

## Load and prepare data ####

rf_global_current <- rast("output/rf_global_pred_regional_current.tif") 
names(rf_global_current) <- "rf_global"
chelsa_nor_250m_cur <- rast("output/predictors_regional_250m_Norway_current_EPSG3035.tif")[[1:19]]
ar50_nor_250m <- rast("output/ar50_250m_EPSG3035.tif")
names(ar50_nor_250m) <- "ar50"
ar50_nor_250m <- as.factor(ar50_nor_250m)
preds_train <- c(rf_global_current, chelsa_nor_250m_cur, ar50_nor_250m)

presence <- read_csv("output/presence_coords_regional.csv")
absence <- read_csv("output/absence_coords_regional.csv")

presence_df <- terra::extract(preds_train, presence[c("x", "y")], xy=TRUE, ID=FALSE) |> 
  mutate(response = 1, .before = 1) |>
  drop_na() |> 
  tibble()

absence_df <- terra::extract(preds_train, absence[c("x", "y")], xy=TRUE, ID=FALSE) |> 
  mutate(response = 0, .before = 1) |>
  drop_na() |> 
  tibble()

training_data_with_coords <- bind_rows(presence_df, absence_df) |>
  mutate(response = factor(response), ar50 = factor(ar50))

training_data_with_coords |> 
  select(response, x, y) |> 
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
  sample_n(1e5) |> 
  arrange(response) |>
  ggplot() +
  geom_sf(aes(color = response), size = 0.5)

n_final_presence <- sum(training_data_with_coords$response == "1")
n_final_background <- sum(training_data_with_coords$response == "0")
n_final_prevalence <- n_final_presence / (n_final_presence + n_final_background)

cat("Final training data - Presences:", n_final_presence, "Background:", n_final_background, "\n")
cat("Prevalence in training data:", round(n_final_prevalence * 100, 3), "%\n")

weights_gini <- read_csv("output/feature_weights_gini.csv")
weights <- pull(weights_gini, weights_gini)
names(weights) <- pull(weights_gini, feature)

## Create clustering sample ####

absence_size <- 1e4 # 1e5 produces error in call to cluster::daisy()
set.seed(42)
training_sample_with_coords <- training_data_with_coords |> 
  filter(response == "0") |>
  sample_n(absence_size) |>
  bind_rows(training_data_with_coords |> filter(response == "1")) |>
  arrange(response)

sum(training_sample_with_coords$response == "1") / 
  (sum(training_sample_with_coords$response == "1") + 
     sum(training_sample_with_coords$response == "0"))

training_sample_with_coords |>
  select(response, x, y) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
  arrange(response) |>
  ggplot() +
  geom_sf(aes(color = response), size = 0.5)

training_sample <- training_sample_with_coords |>
  select(-x, -y)

## Run HDBSCAN clustering ####

feature_cols <- setdiff(names(training_sample), "response")

# gower_dist_weighted <- cluster::daisy(
#   x = training_sample[, feature_cols],
#   metric = "gower",
#   weights = weights # alt: 1
# )

gower_dist_unweighted <- cluster::daisy(
  x = training_sample[, feature_cols],
  metric = "gower",
  weights = 1 # alt: weights
)

hdb_result <- hdbscan(gower_dist_unweighted, minPts = 100)
hdb_result
plot(hdb_result, show_flat=TRUE)

training_sample_with_clusters <- training_sample_with_coords |>
  mutate(cluster = factor(hdb_result$cluster))

## Cluster statistics ####

training_sample_with_clusters |>
  group_by(cluster) |>
  summarise(
    n_points = n(),
    n_presence = sum(as.numeric(as.character(response))),
    prevalence = n_presence / n_points,
    .groups = 'drop'
  ) |>
  arrange(cluster) |> 
  print(n=50)

## Basic visualizations ####

# PCA for feature space
numeric_features <- training_sample |> 
  select(where(is.numeric)) |>
  na.omit()

pca_result <- prcomp(numeric_features, scale. = TRUE)
pca_data <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  cluster = factor(hdb_result$cluster[1:nrow(numeric_features)]),
  response = training_sample$response[1:nrow(numeric_features)]
)

# Feature space plot
plot(pca_data$PC1, pca_data$PC2, col = as.numeric(pca_data$cluster), 
     main = "HDBSCAN Clusters in PCA Space", xlab = "PC1", ylab = "PC2")

# Geographic plot using sf
training_sample_with_clusters |> 
  select(cluster, x, y) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(preds_train)) |>
  arrange(cluster) |>
  ggplot() +
  geom_sf(aes(color = cluster), size = 0.5) +
  labs(title = "HDBSCAN Clusters")

